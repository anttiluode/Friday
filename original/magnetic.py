import numpy as np

def build_magnetic_laplacian(n_clusters=5, cluster_size=6, golden_angle_wiring=True):
    """
    Builds a Magnetic Graph Laplacian for 5 cyclic clusters (each of size cluster_size).
    Internal edges use a locked rational phase.
    Inter-cluster edges use either a rational phase or an irrational golden-angle phase.
    """
    N = n_clusters * cluster_size
    H = np.zeros((N, N), dtype=complex)
    
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    golden_phase = 2.0 * np.pi / (phi ** 2)
    
    # 1. Internal Cluster Edges (Resonant cyclic structure like FourierFlower petals)
    for c in range(n_clusters):
        base = c * cluster_size
        for i in range(cluster_size):
            u = base + i
            v = base + ((i + 1) % cluster_size)
            # Local cyclic weight with internal rational phase
            phase = 2.0 * np.pi * i / cluster_size
            H[u, v] = np.exp(1j * phase)
            H[v, u] = np.exp(-1j * phase) # Hermitian symmetry for Magnetic Laplacian
            
    # 2. Inter-Cluster Wiring
    for c in range(n_clusters):
        u = c * cluster_size
        v = ((c + 1) % n_clusters) * cluster_size # Connect ring of clusters
        
        if golden_angle_wiring:
            # Non-resonant irrational phase shift (Golden angle)
            link_phase = (c + 1) * golden_phase
        else:
            # Resonant rational phase shift
            link_phase = 2.0 * np.pi * (c + 1) / n_clusters
            
        H[u, v] = np.exp(1j * link_phase)
        H[v, u] = np.exp(-1j * link_phase)

    # Construct Magnetic Graph Laplacian: L_mag = D - H
    degrees = np.abs(H).sum(axis=1)
    D = np.diag(degrees)
    L_mag = D - H
    return L_mag

# Run comparison
for mode, is_golden in [("Rational (Resonant)", False), ("Golden-Angle (Non-Resonant)", True)]:
    L = build_magnetic_laplacian(n_clusters=5, cluster_size=6, golden_angle_wiring=is_golden)
    
    # Compute eigenvalues of the Hermitian Magnetic Laplacian
    eigenvals = np.linalg.eigvalsh(L)
    
    # Check spectral gap and smallest non-zero eigenvalue (proxy for observability / conditioning)
    non_zero_evals = eigenvals[eigenvals > 1e-5]
    min_eval = non_zero_evals[0] if len(non_zero_evals) > 0 else 0.0
    spectral_spread = eigenvals[-1] - eigenvals[0]
    
    print(f"--- {mode} Wiring ---")
    print(f"  Smallest non-zero eigenvalue (σ_min proxy): {min_eval:.5f}")
    print(f"  Spectral spread (λ_max - λ_min): {spectral_spread:.5f}")
    print(f"  Eigenvalue degeneracy count (tolerance 1e-4): {N_deg if (N_deg := np.sum(np.diff(np.sort(eigenvals)) < 1e-4)) > 0 else 0}")
    print()