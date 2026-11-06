"""Reproduce Friday's 30-node magnetic-graph results. Python + NumPy only.

This specifies an ideal complex-amplitude observer and, for finite-time
conditioning, unitary evolution. The original script specifies neither.
Run: python audit.py --output results/audit.json
"""
import argparse
import json
from pathlib import Path
import platform

import numpy as np

N = 30
PORTS = np.arange(0, N, 6)
EIGENVALUE_TOL = 1e-8
RANK_TOL = 1e-8


def bridge_phases(mode):
    c = np.arange(1, 6)
    if mode == "rational":
        return 2*np.pi*c/5
    if mode == "golden":
        return 2*np.pi*c/((1+np.sqrt(5))/2)**2
    if mode == "rational_11":
        return 2*np.pi*c/11
    raise ValueError(f"Unknown phase mode: {mode}")


def laplacian(phases, delta=0.0):
    """Original graph; delta changes edge (0,1) within every flower."""
    phases = np.asarray(phases, dtype=float)
    if phases.shape != (5,):
        raise ValueError("Expected five outer phases")
    H = np.zeros((N, N), complex)
    for c in range(5):
        for j in range(6):
            u, v = 6*c+j, 6*c+(j+1) % 6
            phase = 2*np.pi*j/6 + (delta if j == 0 else 0)
            H[u, v] = np.exp(1j*phase)
            H[v, u] = H[u, v].conjugate()
    for c, phase in enumerate(phases):
        u, v = 6*c, 6*((c+1) % 5)
        H[u, v] = np.exp(1j*phase)
        H[v, u] = H[u, v].conjugate()
    return np.diag(np.abs(H).sum(axis=1))-H


def dark_basis():
    """Closed-form basis; no eigensolver is used to construct these modes."""
    j = np.arange(6)
    D = np.zeros((N, 15), complex)
    values = []
    for r in range(3):
        alpha = (2*r+1)*np.pi/6
        local = np.exp(-1j*np.pi*j*(j-1)/6)*np.sin(alpha*j)/np.sqrt(3)
        for c in range(5):
            D[6*c:6*c+6, 5*r+c] = local
            values.append(2-2*np.cos(alpha))
    return D, np.array(values)


def observable_dimension(L, sensors):
    """PBH observable dimension for continuous-time generator -i L.

Sum the ranks of C restricted to each eigenspace. Grouping degenerate
eigenvalues is essential: counting eigenvectors with nonzero readings is wrong.
The numerical tolerances approximate algebraic rank; they are not noise levels.
"""
    values, V = np.linalg.eigh(L)
    groups = np.split(np.arange(N), np.where(np.diff(values) > EIGENVALUE_TOL)[0]+1)
    return int(sum(np.linalg.matrix_rank(V[sensors][:, g], tol=RANK_TOL) for g in groups))


def finite_singular_values(L, sensors, steps=80, dt=0.37):
    """Singular values of stacked C exp(-i k dt L) / sqrt(steps).

Computed in an orthonormal eigenbasis, preserving singular values. Node
readings are complex amplitudes, not intensities. Each sensor is read at
each time, including t=0. Normalization compares average squared read energy.
"""
    values, V = np.linalg.eigh(L)
    CV = V[sensors]
    O = np.concatenate([CV*np.exp(-1j*k*dt*values) for k in range(steps)])
    return np.linalg.svd(O/np.sqrt(steps), compute_uv=False)


def report():
    D, values = dark_basis()
    configurations = {}
    sensors = {"portals": PORTS, "interiors": PORTS+1,
               "both": np.sort(np.r_[PORTS, PORTS+1])}
    for mode in ["rational", "golden", "rational_11"]:
        L = laplacian(bridge_phases(mode))
        eig = np.linalg.eigvalsh(L)
        configurations[mode] = {
            "lambda_min": float(eig[0]),
            "spectral_spread": float(eig[-1]-eig[0]),
            "repeated_gaps_tol_1e_4": int(np.sum(np.diff(eig) < 1e-4)),
            "dark_eigenvector_max_residual": float(np.max(np.abs(L@D-D*values))),
            "observable_dimensions": {name: observable_dimension(L, s) for name, s in sensors.items()},
        }
    p = bridge_phases("golden")
    gauge_error = np.max(np.abs(np.linalg.eigvalsh(laplacian(p)) -
                               np.linalg.eigvalsh(laplacian([0, 0, 0, 0, sum(p)]))))
    opening = []
    for delta in [0, 0.0001, 0.01, 0.1, 0.5]:
        L = laplacian(p, delta)
        sv = finite_singular_values(L, PORTS)
        opening.append({"delta_radians": delta, "observable_dimension": observable_dimension(L, PORTS),
                        "normalized_sigma_min": float(sv[-1]), "normalized_sigma_max": float(sv[0])})
    return {
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "protocol": {"nodes": N, "eigenvalue_group_tolerance": EIGENVALUE_TOL,
                     "rank_tolerance": RANK_TOL, "finite_time_steps": 80, "dt": 0.37,
                     "evolution": "exp(-i t L)", "measurement": "complex node amplitude",
                     "normalization": "observability stack divided by sqrt(80)",
                     "note": "Ideal deterministic linear-state audit; no learned model or noise trial."},
        "dark_basis": {"dimension": 15, "eigenvalues": [float(x) for x in values[::5]],
                       "portal_max_residual": float(np.max(np.abs(D[PORTS]))),
                       "orthonormal_max_residual": float(np.max(np.abs(D.conj().T@D-np.eye(15))))},
        "configurations": configurations,
        "same_flux_max_spectral_error": float(gauge_error),
        "internal_opening": opening,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Also save the JSON report to this path")
    args = parser.parse_args()
    output = json.dumps(report(), indent=2, allow_nan=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()
