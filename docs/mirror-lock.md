# The mirror lock: writing behind Friday's dark modes

*7 November 2026. A follow-up to [What survives, what can be read](what-survives.md). Written by Claude (Opus 5.5) from Antti Luode's Friday and [SliderV3](https://github.com/anttiluode/SliderV3) repositories. Code: [`mirror_lock.py`](../mirror_lock.py). Numbers: [`results/mirror_lock.json`](../results/mirror_lock.json).*

The research note ended with an open requirement: *a useful memory needs a controlled path between protection and access*. Opening a path was shown; writing, retaining and retrieving were not. SliderV3 had the opposite gap. It writes history into edge weights, but everything it writes is public and only weakly readable (about 2.9% probe difference). This note couples the two.

## 1. The dark modes are a mirror's odd sector

Number ring c's nodes 6c+j, j = 0…5, with j = 0 the connecting node. Let R reflect j → −j (mod 6) inside every ring, and let U be the diagonal phase matrix with, in every ring,

```
u_j = 1, e^{iπ/3}, e^{2iπ/3}, −1, e^{−2iπ/3}, e^{−iπ/3}      (j = 0…5)
```

Define (Sx)_i = u_i x_{R(i)}. Then:

- S L S⁻¹ = L for the original graph, with any outer phases (the bridges join fixed points of S);
- S² = 1, and S is the identity on the five connecting nodes;
- the −1 eigenspace of S has dimension 15, and every analytic dark mode lies in it.

The plain reflection R does not commute with L, because it reverses the orientation of each ring and conjugates its hopping phases. The ring flux π is reversal-invariant (−π ≡ π), so a gauge change U can repair it. Flux 0 and π are the only values for which this works. An unfrustrated (flux 0) ring would also have a mirror, but with u_3 = +1 its odd sector is two-dimensional. Flux π puts the far node 3 in the odd sector (u_3 = −1), which gives three dark modes per ring.

A connecting node is fixed by S. For a mirror-odd vector v, v_port = (Sv)_port = −v_port, so v_port = 0. Nothing more is needed to explain darkness at the interface. The original analysis's 15 = 3 × 5 is the dimension count of the odd sector.

## 2. The same symmetry acting on weights

Let every edge magnitude be plastic, with the hopping phases fixed. A configuration of magnitudes m defines L(m). The mirror acts on configurations by relabelling edges, m ↦ m^S, and L(m^S) = S L(m) S⁻¹ by the same gauge argument, because the gauge phases do not depend on magnitudes.

Store one bit per ring as the sign of its antisymmetric pattern, Σ over the three mirror pairs (0,1)|(5,0), (1,2)|(4,5), (2,3)|(3,4) of (m_e − m_mirror(e)). Flipping a ring's bit is applying that ring's S_c. Because S_c fixes every connecting node, three statements hold exactly.

**Reads.** For any drive b supported on connecting nodes, the steady-state amplitudes A = (γ + i(L − ω))⁻¹ b satisfy A(m^S) = S A(m). Restricted to connecting nodes, S is the identity. Every public reading is the same for a word and its mirror, at every frequency and for any number of simultaneous ports.

**Writes.** The write rule uses the gauge-covariant strain |A_p − e^{iθ_pq} A_q|. Under S, both endpoints pick up their gauge phases and the edge is relabelled, so the exposure on edge e in m^S equals the exposure on the mirror edge in m. Any rule built from these exposures by local operations, including clipping, renormalization over a mirror-invariant set, and twin competition, maps mirror pairs to mirror pairs. A symmetric ring (m^S = m) therefore stays symmetric under public traffic forever.

**The public future.** Combine the two. Graphs holding words w and w′ are related by a product of S_c's. Identical public traffic keeps them so related, and every port reading along the way is identical. Nothing an observer at the ports does, including writing, ever depends on the stored word.

**The key.** Anything at an interior node breaks the tie. The reader used here drives the connecting node (a public act) and listens privately at the two interior twins:

```
C(ω) = (|A_1|² − |A_5|²) / (|A_1|² + |A_5|²),     bit = sign(mean over the sweep of C)
```

For a symmetric ring, C(ω) = 0 at every ω; mirroring the ring negates C(ω) at every ω. The key writer drives one interior twin.

## 3. Measurements

Setup: Friday's golden configuration, all 35 magnitudes plastic and starting at 1. Steady-state damped dynamics with γ = 0.2, probe sweep of 48 frequencies in [0.05, 6]. Write rule: SliderV3's exposure 0.48·tanh(strain/0.25), with the response normalized to unit peak, rate 0.06. Public traffic: 1 to 3 random connecting nodes with random complex amplitudes, at a random frequency in [0.1, 5.8].

| Exact claims | Measured |
|---|---|
| ‖SL − LS‖, ‖S² − 1‖, ‖SD + D‖ | 1.7e-15, 1.7e-15, 8e-16 |
| Port difference across all 32 stored words | 7e-16 |
| Ring asymmetry after 400 public writes from symmetric | 2e-15 |
| Key-write crosstalk into the other rings | 2e-16 |
| Port difference, word 00000 vs 11111, over 100 shared public steps | 5e-15 |

| Key performance | Measured |
|---|---|
| Static words (ε = 0.08) decoded | 32/32, minimum margin 0.12 |
| Words written by key (3 drives per ring), correct sign / decoded | 32/32 / 32/32 |

### Retention: the public cannot read or write, but it can erase

With SliderV3's additive rule, public traffic erased the stored words. Minimum over four words:

| Homeostasis | After 50 public writes | After 200 |
|---|---|---|
| none | not decodable; magnitudes heading for the 3.0 clip | asymmetry 0 (all edges clipped) |
| multiplicative ring renormalization | decodable, margin 0.001 | 2% of asymmetry left, signs lost |
| additive pull to base | decodable, margin 0.0003 | 1e-5 left |

The mechanism is ordinary dilution. Additive increments that are roughly equal on two twins shrink their ratio, and a symmetric configuration is a mirror-symmetric target, so nothing forbids reaching it. Equivariance protects *which* of two mirror states the graph is in, but not the *distance* between them. When that distance goes to zero, the information is gone. This is the research note's point 4 ("isolation is not retention"), now for structure instead of state.

Some failures above are sign flips with nonzero asymmetry left. The antisymmetric pattern is a 3-vector per ring, and traffic can rotate it so that its sum changes sign. The pair (word, mirror word) is still distinct, but the fixed decoder no longer knows which is which.

### The mirror flip-flop

Three changes, all mirror-equivariant:

1. multiplicative writes, m ← m(1 + rate·E), so equal exposure leaves twin ratios unchanged;
2. twin competition, m ← m^β rescaled to keep each mirror pair's sum, with β = 1.05;
3. per-ring renormalization to mean 1, and a floor of 0.2.

| Test (16 trials: 8 words × 2 traffic seeds) | Result |
|---|---|
| Words retained after 300 public writes | 16/16 |
| Ring 0 flipped by driving the opposite twin, with competition off (β = 1) | 16/16, median 168 key writes |
| Flipped word intact after re-locking and 150 more public writes | 16/16 |

After locking, ring 0's twin pairs were (1.89, 0.20), (1.92, 0.20) and (1.66, 0.20): one side of the ring strong and the other at the floor. The stored bit is *which side of the ring was pruned*.

A wider sweep, done while developing and not in the committed script:

- β = 1.03 also saturated, but traffic rotated the pattern before it locked (0/33 decoded).
- β = 1.05 and 1.1 retained 33/33.
- Without the unlocking step (β kept above 1), no rewrite succeeded within 200 key writes. The locked flip-flop is write-once.
- With β = 1.1, the unlocked rewrite flipped ring 0 but the new word did not survive re-locking (0/16). I have not diagnosed this.

## 4. What this is and is not

**Prior art.** That symmetry about the input/output nodes makes a networked system uncontrollable and unobservable is known. Rahmani, Ji, Mesbahi and Egerstedt (*SIAM J. Control Optim.* 48, 2009) state it for leader-follower consensus networks, and the general version is standard invariant-subspace reasoning. The new part here is narrow:

- the payload is in plastic *weights*, so the write rule inherits the same blindness as the read;
- the public channel's entire future, including its own writes, is provably independent of the stored word;
- a bistable twin rule makes the hidden sector into storage that survives public traffic.

**Limits.**

- Exactness needs an exact mirror. Any fabrication asymmetry leaks the word to the ports at first order and gives public traffic a lever to write it.
- The flip-flop's symmetric state is unstable by design. Roundoff alone grew from 1e-15 to 7e-8 over 400 public writes from a perfectly symmetric start. Physical noise would write random bits into unwritten rings.
- The locked state is pruning, and rewriting is slow and needs a global "unlock" (β = 1) that is not itself a local mechanism.
- Steady-state linear responses replace time integration. SliderV3's real-valued time stepping showed the same exact invisibility on a 6-ring in a separate check, but that check is not in this repository.
- Five bits on a 30-node graph is a demonstration of a mechanism, not a memory with useful capacity. No noise, capacity or energy comparison against conventional storage is made.

## 5. Connections

- **Neurons.** A soma between two mirror sibling branches is exactly a fixed point of a mirror. It cannot know which sibling holds the stronger synapse. An input arriving off-axis, such as a localized tuft or a single branch, can. Twin competition with pruning is a familiar developmental pattern. This is an analogy, not a model of any measured cell.
- **SliderV3.** Its order memory (AB vs BA) is public. The mirror lock adds a private channel to the same kind of wave-written graph.
- **GhostWorld.** Its wave graph currently has no job. A public stream could write shared dynamics while a key holds identity bits the stream cannot overwrite. That is untested.
