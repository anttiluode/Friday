"""Mirror lock: five private bits stored in the plastic weights of Friday's graph.

Friday's 15 dark modes are not an accident of the spectrum. Each six-node ring
has a gauge reflection S_c (reflect j -> -j about the connecting node, then a
fixed diagonal phase) that commutes with the Laplacian, leaves every connecting
node untouched, and squares to the identity. The dark modes are exactly the
S = S_0 ... S_4 odd sector (15 = 5 x 3).

That symmetry applies to plastic edge MAGNITUDES as well as to states. Store
one bit per ring as the sign of the ring's mirror-antisymmetric weight pattern.
Then, for every input applied only at connecting nodes and every write rule
that is gauge-covariant:

  * public reads (drive and read at connecting nodes) are identical for the
    bit pattern and its mirror image, at every frequency, exactly;
  * public writes map a graph and its mirror image to mirror images, so they
    can never create a bit in a symmetric ring, and the whole public future is
    identical for all 32 stored words;
  * an interior "key" (drive or read at an interior node) can read and write.

This script measures how well the key works on the real 30-node graph, using
a frequency-domain version of SliderV3's Fourier write rule.

Run: python mirror_lock.py --output results/mirror_lock.json
"""
import argparse
import json
from pathlib import Path

import numpy as np

import audit

N = audit.N
PORTS = audit.PORTS
GAMMA = 0.2                            # damping of the driven linear system
SWEEP = np.linspace(0.05, 6.0, 48)     # probe frequencies; spectrum lies in [0, ~5.7]
RATE = 0.06                            # SliderV3 write rate (0.055) rounded
EXPOSURE_SCALE = 0.25                  # SliderV3 tanh scale
MAG_RANGE = (0.05, 3.0)
EPS = 0.08                             # stored asymmetry for static tests

# ring edges as (ring, j) -> nodes (6c+j, 6c+j+1); bridges between connecting nodes
RING_EDGES = [(6*c+j, 6*c+(j+1) % 6) for c in range(5) for j in range(6)]
BRIDGES = [(6*c, 6*((c+1) % 5)) for c in range(5)]
EDGES = RING_EDGES + BRIDGES
# mirror pairs inside a ring: (0,1)<->(5,0), (1,2)<->(4,5), (2,3)<->(3,4)
MIRROR_PAIRS = [[(6*c+j, 6*c+(j+1) % 6) for j in range(3)] for c in range(5)]
MIRROR_PARTNER = [[(6*c+(6-j-1) % 6, 6*c+(6-j) % 6) for j in range(3)] for c in range(5)]


def phases():
    """Unit-magnitude Hermitian hopping matrix of the golden configuration."""
    L = audit.laplacian(audit.bridge_phases("golden"))
    return np.diag(np.diag(L)) - L


PHASE = phases()


def mirror_operator():
    """Port-fixing gauge reflection S with S L S^-1 = L and S^2 = 1."""
    R = np.arange(N)
    for c in range(5):
        for j in range(6):
            R[6*c+j] = 6*c+(-j) % 6
    K = PHASE[np.ix_(R, R)]
    u = np.ones(N, complex)
    for c in range(5):
        for j in range(5):
            p, q = 6*c+j, 6*c+j+1
            u[q] = np.conj(PHASE[p, q]/(u[p]*K[p, q]))
    S = np.zeros((N, N), complex)
    S[np.arange(N), R] = u
    return S, R


S, REFLECT = mirror_operator()


def base_magnitudes():
    return {e: 1.0 for e in EDGES}


def laplacian(mag):
    H = np.zeros((N, N), complex)
    for (p, q), m in mag.items():
        H[p, q] = m*PHASE[p, q]
        H[q, p] = np.conj(H[p, q])
    return np.diag(np.abs(H).sum(1)) - H


def mirror(mag):
    """Magnitudes of the mirror-image graph (reflect every ring)."""
    out = {}
    for (p, q), m in mag.items():
        a, b = REFLECT[p], REFLECT[q]
        out[(a, b) if (a, b) in mag else (b, a)] = m
    return out


def antisym(mag):
    """Per-ring signed asymmetry: sum over mirror pairs of (m_e - m_mirror(e))."""
    return np.array([sum(mag[e]-mag[f] for e, f in zip(MIRROR_PAIRS[c], MIRROR_PARTNER[c]))
                     for c in range(5)])


def response(mag, drive, freqs=SWEEP):
    """Complex steady-state amplitudes A(w) of dx/dt = -iLx - gamma x + b e^{-iwt}.

    drive: dict node -> complex amplitude. Returns array (len(freqs), N).
    """
    L = laplacian(mag)
    b = np.zeros(N, complex)
    for node, amp in drive.items():
        b[node] = amp
    return np.array([np.linalg.solve(GAMMA*np.eye(N) + 1j*(L - w*np.eye(N)), b) for w in freqs])


def write(mag, drive, freq, mode="additive"):
    """SliderV3 Fourier rule, made gauge-covariant: exposure from |A_p - e^{i theta} A_q|.

    additive:       m <- clip(m + rate*E)     (SliderV3)
    multiplicative: m <- m*(1 + rate*E)       (ratios between twins survive uniform exposure)
    """
    A = response(mag, drive, [freq])[0]
    A = A/max(np.abs(A).max(), 1e-12)          # response normalized to unit peak
    out = {}
    for (p, q), m in mag.items():
        E = 0.48*np.tanh(abs(A[p] - PHASE[p, q]*A[q])**2/EXPOSURE_SCALE)
        if mode == "additive":
            out[(p, q)] = float(np.clip(m + RATE*E, *MAG_RANGE))
        else:
            out[(p, q)] = float(m*(1 + RATE*E))
    return out


def compete(mag, beta):
    """Soft winner-take-all between mirror twins, conserving each twin pair's sum.

    beta = 1 does nothing. beta > 1 makes the symmetric state unstable and the two
    mirror states attractors. It commutes with the mirror, so it decides nothing itself.
    """
    out = dict(mag)
    for c in range(5):
        for e, f in zip(MIRROR_PAIRS[c], MIRROR_PARTNER[c]):
            a, b = mag[e]**beta, mag[f]**beta
            s = (mag[e] + mag[f])/(a + b)
            out[e], out[f] = a*s, b*s
    return out


FLIPFLOP = {"beta": 1.05, "floor": 0.2}


def flipflop_step(mag, drive, freq, beta=FLIPFLOP["beta"], floor=FLIPFLOP["floor"]):
    """Multiplicative write, twin competition, ring renormalization, weight floor."""
    mag = renormalize(compete(write(mag, drive, freq, "multiplicative"), beta))
    return {e: max(m, floor) for e, m in mag.items()}


def renormalize(mag):
    """Multiplicative homeostasis: rescale each ring to mean magnitude 1 (mirror-equivariant)."""
    out = dict(mag)
    for c in range(5):
        ring = [e for e in RING_EDGES if e[0]//6 == c]
        s = np.mean([mag[e] for e in ring])
        for e in ring:
            out[e] = mag[e]/s
    return out


def pull(mag, rate=0.05):
    """Additive homeostasis toward the symmetric base (mirror-equivariant, but erodes bits)."""
    return {e: m + rate*(1.0-m) for e, m in mag.items()}


def stored(word, eps=EPS):
    mag = base_magnitudes()
    for c, bit in enumerate(word):
        s = eps if bit else -eps
        mag[MIRROR_PAIRS[c][0]] += s
        mag[MIRROR_PARTNER[c][0]] -= s
    return mag


def rel(x, y):
    return float(np.sqrt(np.mean(np.abs(x-y)**2)/((np.mean(np.abs(x)**2)+np.mean(np.abs(y)**2))/2)))


def public_reading(mag):
    """Drive each connecting node in turn, read all connecting nodes."""
    return np.concatenate([response(mag, {p: 1.0})[:, PORTS] for p in PORTS])


def key_contrast(mag, c, freqs=SWEEP):
    """Public drive at the ring's connecting node, private listening at its two interior twins.

    C(w) = (|A_1|^2 - |A_5|^2)/(|A_1|^2 + |A_5|^2). For a mirror-symmetric ring it is 0 at
    every w; mirroring the ring flips its sign at every w.
    """
    A = response(mag, {6*c: 1.0}, freqs)
    a1, a5 = np.abs(A[:, 6*c+1])**2, np.abs(A[:, 6*c+5])**2
    return (a1-a5)/(a1+a5)


def decode(mag, freqs=SWEEP):
    """Bit = sign of the band-averaged contrast. Margin = |average contrast|."""
    C = np.array([key_contrast(mag, c, freqs).mean() for c in range(5)])
    return (C > 0).astype(int), np.abs(C)


def random_public_history(rng, steps):
    hist = []
    for _ in range(steps):
        k = rng.integers(1, 4)
        ports = rng.choice(PORTS, size=k, replace=False)
        drive = {int(p): complex(rng.normal(), rng.normal()) for p in ports}
        hist.append((drive, float(rng.uniform(0.1, 5.8))))
    return hist


def run(mag, history, homeostasis=None):
    for drive, f in history:
        mag = write(mag, drive, f)
        if homeostasis == "renormalize":
            mag = renormalize(mag)
        elif homeostasis == "pull":
            mag = pull(mag)
    return mag


def key_write(mag, c, side, freqs=(0.6, 1.4, 2.5), homeostasis="renormalize"):
    """Write ring c's bit by driving one interior twin (side=+1 -> node 1, -1 -> node 5)."""
    node = 6*c + (1 if side > 0 else 5)
    for f in freqs:
        mag = write(mag, {node: 1.0}, f)
        if homeostasis == "renormalize":
            mag = renormalize(mag)
    return mag


def report(seed=0, public_steps=400):
    rng = np.random.default_rng(seed)
    out = {"protocol": {"graph": "Friday golden configuration, plastic edge magnitudes (all 35 edges)",
                        "dynamics": "dx/dt = -iLx - gamma x + drive, steady state", "gamma": GAMMA,
                        "sweep": [float(SWEEP[0]), float(SWEEP[-1]), len(SWEEP)],
                        "write_rule": "m += rate*0.48*tanh(|A_p - e^{i theta_pq} A_q|^2/0.25), clip [0.05,3]",
                        "rate": RATE, "stored_eps": EPS, "seed": seed}}

    # 0. structure: dark modes are the S-odd sector
    D, _ = audit.dark_basis()
    L0 = laplacian(base_magnitudes())
    out["symmetry"] = {
        "S_commutes_with_L": float(np.abs(S@L0-L0@S).max()),
        "S_squared_minus_identity": float(np.abs(S@S-np.eye(N)).max()),
        "S_fixes_connecting_nodes": bool(np.allclose(S[np.ix_(PORTS, PORTS)], np.eye(5))),
        "dark_modes_are_S_odd": float(np.abs(S@D+D).max()),
        "odd_sector_dimension": int(np.sum(np.linalg.eigvalsh((S+S.conj().T)/2) < 0)),
    }

    # 1. static: all 32 words look identical from the ports; the key decodes them
    words = [tuple((w >> c) & 1 for c in range(5)) for w in range(32)]
    ref = public_reading(stored(words[0]))
    pub = max(rel(public_reading(stored(w)), ref) for w in words[1:])
    correct, margins = 0, []
    for w in words:
        bits, m = decode(stored(w))
        correct += int(tuple(bits) == w)
        margins.append(m.min())
    out["static"] = {"max_public_difference_over_32_words": pub,
                     "key_decoded_words": f"{correct}/32", "min_key_margin": float(min(margins))}

    # 2. public traffic can not write a bit into symmetric rings
    hist = random_public_history(rng, public_steps)
    m_pub = run(base_magnitudes(), hist, "renormalize")
    out["public_cannot_write"] = {"steps": public_steps,
                                  "max_abs_ring_asymmetry": float(np.abs(antisym(m_pub)).max()),
                                  "mean_magnitude_after_raw_traffic":
                                      float(np.mean(list(run(base_magnitudes(), hist[:100]).values())))}

    # 3. key writes: write every word with the key into a blank graph, then bury it under traffic
    trials = []
    for w in words:
        m = base_magnitudes()
        for c, bit in enumerate(w):
            m = key_write(m, c, +1 if bit else -1)
        trials.append(m)
    # which sign does a node-1 key write produce? (calibration, same for all rings)
    sign = np.sign(antisym(key_write(base_magnitudes(), 0, +1))[0])
    asym_written = np.array([antisym(m) for m in trials])
    written_ok = sum(tuple(((sign*a) > 0).astype(int)) == w for a, w in zip(asym_written, words))
    crosstalk = float(np.abs(antisym(key_write(base_magnitudes(), 0, +1))[1:]).max())
    dec_ok = sum(tuple(decode(m)[0]) == w for m, w in zip(trials, words))
    out["key_write"] = {"node1_write_sign": int(sign),
                        "words_written_with_correct_signs": f"{written_ok}/32",
                        "words_key_decoded_after_write": f"{dec_ok}/32",
                        "min_written_asymmetry": float(np.abs(asym_written).min()),
                        "crosstalk_into_other_rings": crosstalk,
                        "public_difference_over_32_written_words":
                            max(rel(public_reading(m), public_reading(trials[0])) for m in trials[1:])}

    # 4. retention under public traffic, three homeostasis rules
    retention = {}
    sample = [words[i] for i in (5, 10, 19, 26)]
    for rule in ["none", "renormalize", "pull"]:
        rows = []
        for w in sample:
            m = base_magnitudes()
            for c, bit in enumerate(w):
                m = key_write(m, c, +1 if bit else -1)
            a0 = antisym(m)
            checkpoints = {}
            for k, (drive, f) in enumerate(hist, start=1):
                m = run(m, [(drive, f)], None if rule == "none" else rule)
                if k in (50, 200, public_steps):
                    a = antisym(m)
                    bits, margin = decode(m)
                    checkpoints[k] = {"asym_retained_fraction": float(np.min(np.abs(a)/np.abs(a0))),
                                      "sign_intact": bool(np.all(np.sign(a) == np.sign(a0))),
                                      "key_decoded": bool(tuple(bits) == w),
                                      "min_key_margin": float(margin.min()),
                                      "mean_magnitude": float(np.mean(list(m.values())))}
            rows.append(checkpoints)
        retention[rule] = {str(k): {key: (all(r[k][key] for r in rows) if isinstance(rows[0][k][key], bool)
                                          else float(min(r[k][key] for r in rows)) if key != "mean_magnitude"
                                          else float(np.mean([r[k][key] for r in rows])))
                                    for key in rows[0][k]} for k in rows[0]}
    out["retention"] = retention

    # 5. flip-flop rule: write, bury, rewrite with the lock open, bury again
    flip = {"rule": "multiplicative write, twin competition beta, ring renormalization, floor",
            **FLIPFLOP}
    retained = rewrites = intact = trials_n = 0
    needed = []
    for s_ in range(2):
        h = random_public_history(np.random.default_rng(seed + 10 + s_), 300)
        for w in words[::4]:
            trials_n += 1
            m = base_magnitudes()
            for c, bit in enumerate(w):
                m = key_write(m, c, +1 if bit else -1)
            for drive, f in h:
                m = flipflop_step(m, drive, f)
            retained += int(tuple(decode(m)[0]) == w)
            node, n = (1 if w[0] == 0 else 5), None
            for k in range(1, 301):
                m = flipflop_step(m, {node: 1.0}, (0.6, 1.4, 2.5)[k % 3], beta=1.0)
                if decode(m)[0][0] != w[0]:
                    n = k
                    break
            if n:
                rewrites += 1
                needed.append(n)
            for drive, f in h[:150]:
                m = flipflop_step(m, drive, f)
            intact += int(tuple(decode(m)[0]) == (1 - w[0],) + w[1:])
    sym = base_magnitudes()
    for drive, f in hist:
        sym = flipflop_step(sym, drive, f)
    m_lock = base_magnitudes()
    m_lock = key_write(m_lock, 0, +1)
    for drive, f in hist[:300]:
        m_lock = flipflop_step(m_lock, drive, f)
    flip.update({"words_retained_after_300_public_writes": f"{retained}/{trials_n}",
                 "ring0_rewritten_by_key_with_lock_open": f"{rewrites}/{trials_n}",
                 "median_key_writes_to_flip": float(np.median(needed)) if needed else None,
                 "new_word_intact_after_150_more_public_writes": f"{intact}/{trials_n}",
                 "symmetric_start_asymmetry_after_400_public_writes": float(np.abs(antisym(sym)).max()),
                 "ring0_twin_magnitudes_after_lock": [[round(m_lock[e], 3), round(m_lock[f], 3)]
                                                      for e, f in zip(MIRROR_PAIRS[0], MIRROR_PARTNER[0])]})
    out["flipflop"] = flip

    # 6. the whole public future is independent of the stored word
    m_a, m_b = stored(words[0]), stored(words[31])
    worst = 0.0
    for drive, f in hist[:100]:
        ra = response(m_a, drive, [f])[0][PORTS]
        rb = response(m_b, drive, [f])[0][PORTS]
        worst = max(worst, rel(ra, rb))
        m_a, m_b = write(m_a, drive, f), write(m_b, drive, f)
    out["public_future"] = {"steps": 100, "max_port_difference_word0_vs_word31": worst}
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--steps", type=int, default=400)
    args = parser.parse_args()
    text = json.dumps(report(args.seed, args.steps), indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
