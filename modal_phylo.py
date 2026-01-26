# -*- coding: utf-8 -*-
"""
Weighted Serial-Sample Neighbor-Joining (wsNJ) for ancient DNA phylogenetics.

See NEW_PHYLO.md for detailed algorithm explanation and theory.

Usage:
    from modal_phylo import run_wsnj_local, tree_to_linkage

    tree, mu, results, D_raw = run_wsnj_local(seq_matrix, tip_dates, weights)
    linkage = tree_to_linkage(tree)  # For scipy dendrogram plotting
"""

import numpy as np
from numba import njit, prange
from scipy.spatial.distance import squareform


# Floor value for corrected distances (keeps matrix positive-definite)
DISTANCE_FLOOR = 1e-9

# Fallback mutation rate from Soares et al. 2009 (mtDNA coding region)
FALLBACK_MU = 1.57e-8


@njit(parallel=True, cache=True)
def weighted_hamming_distance_numba(seq_matrix, weights, missing_val=4):
    """
    Compute weighted Hamming distance matrix.

    Parameters
    ----------
    seq_matrix : (n_samples, n_sites) int8 array
        Integer-encoded sequences (0=A, 1=C, 2=G, 3=T, 4=missing).
    weights : (n_sites,) float32 array
        Per-site weights (typically 1/contacts from protein structure).
    missing_val : int
        Code for missing data (skipped in computation).

    Returns
    -------
    D : (n_samples, n_samples) float64 array
        Symmetric pairwise distance matrix.
    """
    n_samples, n_sites = seq_matrix.shape
    D = np.zeros((n_samples, n_samples), dtype=np.float64)

    for i in prange(n_samples):
        for j in range(i + 1, n_samples):
            weighted_mismatches = 0.0
            weight_sum = 0.0

            for k in range(n_sites):
                si = seq_matrix[i, k]
                sj = seq_matrix[j, k]

                if si == missing_val or sj == missing_val:
                    continue

                w = weights[k]
                weight_sum += w

                if si != sj:
                    weighted_mismatches += w

            if weight_sum > 0:
                d = weighted_mismatches / weight_sum
            else:
                d = 1.0

            D[i, j] = d
            D[j, i] = d

    return D


@njit(parallel=True, cache=True)
def correct_distances_for_time_numba(D, tip_dates, mu):
    """
    Correct pairwise distances for sampling time differences.

    Computes d'_ij = d_ij - μ × |t_i - t_j|, floored at DISTANCE_FLOOR.

    Parameters
    ----------
    D : (n, n) float64 array
        Raw pairwise distance matrix.
    tip_dates : (n,) float64 array
        Sample ages in years BP.
    mu : float
        Mutation rate estimate.

    Returns
    -------
    D_corrected : (n, n) float64 array
        Time-corrected distance matrix.
    """
    n = len(tip_dates)
    D_corrected = np.zeros((n, n), dtype=np.float64)

    for i in prange(n):
        for j in range(i + 1, n):
            time_diff = abs(tip_dates[i] - tip_dates[j])
            d_corr = D[i, j] - mu * time_diff

            if d_corr < DISTANCE_FLOOR:
                d_corr = DISTANCE_FLOOR

            D_corrected[i, j] = d_corr
            D_corrected[j, i] = d_corr

    return D_corrected


def neighbor_joining(D, labels=None):
    """
    Neighbor-Joining tree construction (Saitou & Nei, 1987).

    Parameters
    ----------
    D : (n, n) float64 array
        Pairwise distance matrix.
    labels : list, optional
        Tip labels. Defaults to integer indices.

    Returns
    -------
    tree : dict
        Tree with keys: 'children', 'distances', 'labels', 'root', 'n_tips'.
    """
    n = D.shape[0]

    if labels is None:
        labels = list(range(n))

    # Pre-allocate for n tips + (n-1) internal nodes
    total_nodes = 2 * n - 1
    D_work = np.zeros((total_nodes, total_nodes), dtype=np.float64)
    D_work[:n, :n] = D

    active = list(range(n))
    children = {i: None for i in range(n)}
    distances = {i: None for i in range(n)}
    node_labels = {i: labels[i] for i in range(n)}
    next_node = n

    while len(active) > 2:
        m = len(active)

        # Compute row sums for Q-matrix
        row_sums = {}
        for i in active:
            row_sums[i] = sum(D_work[i, j] for j in active if j != i)

        # Find pair minimizing Q_ij = (m-2) × d_ij - r_i - r_j
        min_Q = float('inf')
        min_pair = None

        for idx_i, i in enumerate(active):
            for j in active[idx_i + 1:]:
                Q_ij = (m - 2) * D_work[i, j] - row_sums[i] - row_sums[j]
                if Q_ij < min_Q:
                    min_Q = Q_ij
                    min_pair = (i, j)

        i, j = min_pair

        # Compute branch lengths to new node
        if m > 2:
            d_iu = 0.5 * D_work[i, j] + (row_sums[i] - row_sums[j]) / (2 * (m - 2))
            d_ju = D_work[i, j] - d_iu
        else:
            d_iu = D_work[i, j] / 2
            d_ju = D_work[i, j] / 2

        d_iu = max(d_iu, 0)
        d_ju = max(d_ju, 0)

        # Create new internal node
        u = next_node
        next_node += 1

        children[u] = (i, j)
        distances[u] = (d_iu, d_ju)
        node_labels[u] = None

        # Update distances: d_uk = ½ × (d_ik + d_jk - d_ij)
        for k in active:
            if k != i and k != j:
                d_uk = 0.5 * (D_work[i, k] + D_work[j, k] - D_work[i, j])
                D_work[u, k] = d_uk
                D_work[k, u] = d_uk

        active.remove(i)
        active.remove(j)
        active.append(u)

    # Final join
    if len(active) == 2:
        i, j = active
        root = next_node
        children[root] = (i, j)
        distances[root] = (D_work[i, j] / 2, D_work[i, j] / 2)
        node_labels[root] = None
    else:
        root = active[0]

    return {
        'children': children,
        'distances': distances,
        'labels': node_labels,
        'root': root,
        'n_tips': n,
    }


def root_to_tip_distances(tree):
    """Compute root-to-tip distances for all tips."""
    children = tree['children']
    distances = tree['distances']
    root = tree['root']

    tip_dists = {}

    def traverse(node, dist_from_root):
        if children[node] is None:
            tip_dists[node] = dist_from_root
        else:
            left, right = children[node]
            d_left, d_right = distances[node]
            traverse(left, dist_from_root + d_left)
            traverse(right, dist_from_root + d_right)

    traverse(root, 0.0)
    return tip_dists


def estimate_rate_root_to_tip(tree, tip_dates, tip_indices=None):
    """
    Estimate mutation rate via root-to-tip regression.

    Returns
    -------
    mu : float
        Mutation rate (slope).
    intercept : float
        Y-intercept.
    r_squared : float
        Coefficient of determination.
    """
    tip_dists = root_to_tip_distances(tree)
    n_tips = tree['n_tips']

    if tip_indices is None:
        tip_indices = list(range(n_tips))

    dates = []
    dists = []

    for tip_idx in range(n_tips):
        if tip_idx in tip_dists:
            date_idx = tip_indices[tip_idx] if isinstance(tip_indices, dict) else tip_idx
            dates.append(tip_dates[date_idx])
            dists.append(tip_dists[tip_idx])

    dates = np.array(dates)
    dists = np.array(dists)

    # OLS regression
    X = np.column_stack([dates, np.ones(len(dates))])
    beta, _, _, _ = np.linalg.lstsq(X, dists, rcond=None)

    mu = beta[0]
    intercept = beta[1]

    # R-squared
    y_pred = X @ beta
    ss_res = ((dists - y_pred) ** 2).sum()
    ss_tot = ((dists - dists.mean()) ** 2).sum()
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    return mu, intercept, r_squared


def tree_to_newick(tree, include_distances=True):
    """Convert NJ tree to Newick format string."""
    children = tree['children']
    distances = tree['distances']
    labels = tree['labels']
    root = tree['root']

    def subtree_newick(node, branch_length=None):
        if children[node] is None:
            label = str(labels[node])
            if include_distances and branch_length is not None:
                return f"{label}:{branch_length:.6f}"
            return label
        else:
            left, right = children[node]
            d_left, d_right = distances[node]
            left_str = subtree_newick(left, d_left)
            right_str = subtree_newick(right, d_right)
            subtree = f"({left_str},{right_str})"
            if include_distances and branch_length is not None:
                return f"{subtree}:{branch_length:.6f}"
            return subtree

    return subtree_newick(root) + ";"


def tree_to_linkage(tree):
    """Convert NJ tree to scipy linkage matrix format for dendrogram plotting."""
    children = tree['children']
    distances = tree['distances']
    n_tips = tree['n_tips']

    linkage_rows = []
    node_heights = {i: 0.0 for i in range(n_tips)}
    node_counts = {i: 1 for i in range(n_tips)}

    internal_nodes = sorted([n for n in children if children[n] is not None])

    for node in internal_nodes:
        left, right = children[node]
        d_left, d_right = distances[node]

        height = max(node_heights[left] + d_left, node_heights[right] + d_right)
        node_heights[node] = height
        node_counts[node] = node_counts[left] + node_counts[right]

        linkage_rows.append([
            float(left),
            float(right),
            height,
            float(node_counts[node])
        ])

    return np.array(linkage_rows)


def estimate_rate_pairwise(D_raw, tip_dates):
    """
    Estimate mutation rate via pairwise regression (before tree construction).

    Returns mu, d0, r_squared.
    """
    n = len(tip_dates)

    distances = []
    time_diffs = []

    for i in range(n):
        for j in range(i + 1, n):
            if np.isnan(D_raw[i, j]) or np.isnan(tip_dates[i]) or np.isnan(tip_dates[j]):
                continue
            distances.append(D_raw[i, j])
            time_diffs.append(abs(tip_dates[i] - tip_dates[j]))

    distances = np.array(distances)
    time_diffs = np.array(time_diffs)

    X = np.column_stack([time_diffs, np.ones(len(time_diffs))])
    beta, _, _, _ = np.linalg.lstsq(X, distances, rcond=None)

    mu = beta[0]
    d0 = beta[1]

    y_pred = X @ beta
    ss_res = ((distances - y_pred) ** 2).sum()
    ss_tot = ((distances - distances.mean()) ** 2).sum()
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    return mu, d0, r_squared


def run_wsnj_local(seq_matrix, tip_dates, weights,
                   mu_init=None, max_iter=10, tol=1e-10,
                   verbose=True):
    """
    Weighted Serial-Sample Neighbor-Joining with iterative rate estimation.

    Parameters
    ----------
    seq_matrix : (n_samples, n_sites) int8 array
        Integer-encoded sequences.
    tip_dates : (n_samples,) float64 array
        Sample ages in years BP.
    weights : (n_sites,) float32 array
        Per-site weights.
    mu_init : float, optional
        Initial rate estimate.
    max_iter : int
        Maximum iterations.
    tol : float
        Convergence tolerance.
    verbose : bool
        Print progress.

    Returns
    -------
    tree : dict
        NJ tree structure.
    mu : float
        Final mutation rate estimate.
    results : dict
        Additional results (r_squared, n_negative, history, etc.).
    D_raw : array
        Raw distance matrix.
    """
    import time as time_module

    n_samples = len(seq_matrix)

    if verbose:
        print("=" * 60)
        print("WEIGHTED SERIAL-SAMPLE NEIGHBOR-JOINING")
        print("=" * 60)
        print(f"Samples: {n_samples}")
        print(f"Sites: {seq_matrix.shape[1]}")
        print()

    seq_matrix = np.ascontiguousarray(seq_matrix, dtype=np.int8)
    weights = np.ascontiguousarray(weights, dtype=np.float32)
    tip_dates = np.asarray(tip_dates, dtype=np.float64)

    # Step 1: Compute weighted distances
    if verbose:
        print("Step 1: Computing weighted pairwise distances...")
        t0 = time_module.time()

    D_raw = weighted_hamming_distance_numba(seq_matrix, weights, 4)

    if verbose:
        print(f"  Distance range: {D_raw.min():.6f} - {D_raw.max():.6f}")
        print(f"  Time: {time_module.time() - t0:.1f}s")
        print()

    # Step 2: Initial rate estimate
    if verbose:
        print("Step 2: Initial rate estimate (pairwise regression)...")

    mu_pairwise, d0_pairwise, r2_pairwise = estimate_rate_pairwise(D_raw, tip_dates)

    if verbose:
        print(f"  Pairwise: μ = {mu_pairwise:.4e}, d₀ = {d0_pairwise:.4f}, R² = {r2_pairwise:.4f}")

    if mu_init is not None:
        mu = mu_init
        if verbose:
            print(f"  Using provided initial μ = {mu:.4e}")
    elif mu_pairwise > 0:
        mu = mu_pairwise
    else:
        mu = FALLBACK_MU
        if verbose:
            print(f"  Warning: pairwise μ ≤ 0, using fallback {mu:.4e}")

    if verbose:
        print()

    # Step 3: Iterative refinement
    if verbose:
        print("Step 3: Iterative rate estimation (root-to-tip regression)...")

    history = {
        'iteration': [],
        'mu': [],
        'r_squared': [],
        'n_negative': [],
    }

    time_diffs_matrix = np.abs(np.subtract.outer(tip_dates, tip_dates))

    for iteration in range(max_iter):
        D_corrected = correct_distances_for_time_numba(D_raw, tip_dates, mu)
        n_negative = int(((D_raw - mu * time_diffs_matrix) < 0).sum() // 2)

        tree = neighbor_joining(D_corrected)
        mu_new, intercept, r_squared = estimate_rate_root_to_tip(tree, tip_dates)

        history['iteration'].append(iteration + 1)
        history['mu'].append(mu_new)
        history['r_squared'].append(r_squared)
        history['n_negative'].append(n_negative)

        if verbose:
            print(f"  Iter {iteration + 1}: μ = {mu_new:.4e}, R² = {r_squared:.4f}, neg = {n_negative}")

        if abs(mu_new - mu) < tol:
            if verbose:
                print(f"  Converged!")
            break

        if mu_new <= 0:
            if verbose:
                print(f"  Warning: negative rate, keeping previous value")
            break

        mu = mu_new

    if verbose:
        print()

    # Step 4: Final tree
    D_corrected = correct_distances_for_time_numba(D_raw, tip_dates, mu)
    tree = neighbor_joining(D_corrected)
    _, _, r_squared_final = estimate_rate_root_to_tip(tree, tip_dates)
    n_negative_final = int(((D_raw - mu * time_diffs_matrix) < 0).sum() // 2)

    if verbose:
        print("=" * 60)
        print("RESULTS")
        print("=" * 60)
        print(f"Final μ: {mu:.4e} subst/site/year")
        print(f"Root-to-tip R²: {r_squared_final:.4f}")
        print(f"Iterations: {len(history['iteration'])}")
        print(f"Negative corrections: {n_negative_final}")

    results = {
        'mu': mu,
        'd0': d0_pairwise,
        'r_squared': r_squared_final,
        'r_squared_pairwise': r2_pairwise,
        'n_negative': n_negative_final,
        'n_iterations': len(history['iteration']),
        'history': history,
    }

    return tree, mu, results, D_raw


def compute_distance_matrix_local(seq_matrix, weights, missing_val=4):
    """Compute weighted distance matrix locally."""
    seq_matrix = np.ascontiguousarray(seq_matrix, dtype=np.int8)
    weights = np.ascontiguousarray(weights, dtype=np.float32)
    return weighted_hamming_distance_numba(seq_matrix, weights, missing_val)


# =============================================================================
# MODAL CLOUD COMPUTE (optional)
# =============================================================================

try:
    import modal

    app = modal.App("phylo-nj")

    image = modal.Image.debian_slim(python_version="3.11").pip_install(
        "numpy",
        "numba",
        "scipy",
    )

    @app.function(image=image, cpu=8, memory=4096, timeout=600)
    def run_wsnj_modal(seq_matrix_bytes: bytes,
                       tip_dates_bytes: bytes,
                       weights_bytes: bytes,
                       seq_shape: tuple,
                       max_iter: int = 10) -> dict:
        """Run wsNJ on Modal cloud."""
        import numpy as np
        from numba import njit, prange

        seq_matrix = np.frombuffer(seq_matrix_bytes, dtype=np.int8).reshape(seq_shape)
        tip_dates = np.frombuffer(tip_dates_bytes, dtype=np.float64)
        weights = np.frombuffer(weights_bytes, dtype=np.float32)

        n_samples, n_sites = seq_matrix.shape

        @njit(parallel=True, cache=True)
        def _weighted_hamming(seq_matrix, weights, missing_val):
            n_samples, n_sites = seq_matrix.shape
            D = np.zeros((n_samples, n_samples), dtype=np.float64)
            for i in prange(n_samples):
                for j in range(i + 1, n_samples):
                    weighted_mismatches = 0.0
                    weight_sum = 0.0
                    for k in range(n_sites):
                        si = seq_matrix[i, k]
                        sj = seq_matrix[j, k]
                        if si == missing_val or sj == missing_val:
                            continue
                        w = weights[k]
                        weight_sum += w
                        if si != sj:
                            weighted_mismatches += w
                    if weight_sum > 0:
                        d = weighted_mismatches / weight_sum
                    else:
                        d = 1.0
                    D[i, j] = d
                    D[j, i] = d
            return D

        @njit(parallel=True, cache=True)
        def _correct_time(D, tip_dates, mu):
            n = len(tip_dates)
            D_corrected = np.zeros((n, n), dtype=np.float64)
            for i in prange(n):
                for j in range(i + 1, n):
                    time_diff = abs(tip_dates[i] - tip_dates[j])
                    d_corr = D[i, j] - mu * time_diff
                    if d_corr < 1e-9:
                        d_corr = 1e-9
                    D_corrected[i, j] = d_corr
                    D_corrected[j, i] = d_corr
            return D_corrected

        def _neighbor_joining(D):
            n = D.shape[0]
            total_nodes = 2 * n - 1
            D_work = np.zeros((total_nodes, total_nodes), dtype=np.float64)
            D_work[:n, :n] = D
            active = list(range(n))
            children = {i: None for i in range(n)}
            distances = {i: None for i in range(n)}
            next_node = n

            while len(active) > 2:
                m = len(active)
                row_sums = {i: sum(D_work[i, j] for j in active if j != i) for i in active}

                min_Q = float('inf')
                min_pair = None
                for idx_i, i in enumerate(active):
                    for j in active[idx_i + 1:]:
                        Q_ij = (m - 2) * D_work[i, j] - row_sums[i] - row_sums[j]
                        if Q_ij < min_Q:
                            min_Q = Q_ij
                            min_pair = (i, j)

                i, j = min_pair

                if m > 2:
                    d_iu = 0.5 * D_work[i, j] + (row_sums[i] - row_sums[j]) / (2 * (m - 2))
                    d_ju = D_work[i, j] - d_iu
                else:
                    d_iu = D_work[i, j] / 2
                    d_ju = D_work[i, j] / 2

                d_iu = max(d_iu, 0)
                d_ju = max(d_ju, 0)

                u = next_node
                next_node += 1
                children[u] = (i, j)
                distances[u] = (d_iu, d_ju)

                for k in active:
                    if k != i and k != j:
                        d_uk = 0.5 * (D_work[i, k] + D_work[j, k] - D_work[i, j])
                        D_work[u, k] = d_uk
                        D_work[k, u] = d_uk

                active.remove(i)
                active.remove(j)
                active.append(u)

            if len(active) == 2:
                i, j = active
                root = next_node
                children[root] = (i, j)
                distances[root] = (D_work[i, j] / 2, D_work[i, j] / 2)
            else:
                root = active[0]

            return {'children': children, 'distances': distances, 'root': root, 'n_tips': n}

        def _root_to_tip(tree, tip_dates):
            children = tree['children']
            distances = tree['distances']
            root = tree['root']
            n_tips = tree['n_tips']

            tip_dists = {}
            def traverse(node, dist):
                if children[node] is None:
                    tip_dists[node] = dist
                else:
                    left, right = children[node]
                    d_left, d_right = distances[node]
                    traverse(left, dist + d_left)
                    traverse(right, dist + d_right)
            traverse(root, 0.0)

            dates = np.array([tip_dates[i] for i in range(n_tips)])
            dists = np.array([tip_dists[i] for i in range(n_tips)])

            X = np.column_stack([dates, np.ones(len(dates))])
            beta, _, _, _ = np.linalg.lstsq(X, dists, rcond=None)
            mu = beta[0]

            y_pred = X @ beta
            ss_res = ((dists - y_pred) ** 2).sum()
            ss_tot = ((dists - dists.mean()) ** 2).sum()
            r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0

            return mu, r_squared

        # Run pipeline
        seq_matrix = np.ascontiguousarray(seq_matrix, dtype=np.int8)
        weights = np.ascontiguousarray(weights, dtype=np.float32)
        D_raw = _weighted_hamming(seq_matrix, weights, 4)

        # Initial rate (pairwise)
        distances = []
        time_diffs = []
        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                distances.append(D_raw[i, j])
                time_diffs.append(abs(tip_dates[i] - tip_dates[j]))
        distances = np.array(distances)
        time_diffs = np.array(time_diffs)
        X = np.column_stack([time_diffs, np.ones(len(time_diffs))])
        beta, _, _, _ = np.linalg.lstsq(X, distances, rcond=None)
        mu = float(beta[0]) if beta[0] > 0 else 1.57e-8
        d0 = float(beta[1])

        # Iterate
        for iteration in range(max_iter):
            D_corrected = _correct_time(D_raw, tip_dates, mu)
            tree = _neighbor_joining(D_corrected)
            mu_new, r_squared = _root_to_tip(tree, tip_dates)

            if abs(mu_new - mu) < 1e-10 or mu_new <= 0:
                break
            mu = mu_new

        # Final
        D_corrected = _correct_time(D_raw, tip_dates, mu)
        tree = _neighbor_joining(D_corrected)
        _, r_squared_final = _root_to_tip(tree, tip_dates)

        time_diffs_matrix = np.abs(np.subtract.outer(tip_dates, tip_dates))
        n_negative = int(((D_raw - mu * time_diffs_matrix) < 0).sum() // 2)

        # Convert tree to linkage format
        linkage_rows = []
        node_heights = {i: 0.0 for i in range(tree['n_tips'])}
        node_counts = {i: 1 for i in range(tree['n_tips'])}

        internal_nodes = sorted([n for n in tree['children'] if tree['children'][n] is not None])
        for node in internal_nodes:
            left, right = tree['children'][node]
            d_left, d_right = tree['distances'][node]
            height = max(node_heights[left] + d_left, node_heights[right] + d_right)
            node_heights[node] = height
            node_counts[node] = node_counts[left] + node_counts[right]
            linkage_rows.append([float(left), float(right), height, float(node_counts[node])])

        linkage_matrix = np.array(linkage_rows)

        return {
            "tree_bytes": linkage_matrix.tobytes(),
            "tree_shape": linkage_matrix.shape,
            "D_raw_bytes": D_raw.tobytes(),
            "D_raw_shape": D_raw.shape,
            "mu": mu,
            "d0": d0,
            "r_squared": r_squared_final,
            "n_negative": n_negative,
        }

    def run_wsnj_on_modal(seq_matrix, tip_dates, weights, verbose=True):
        """Client function to run wsNJ on Modal cloud."""
        if verbose:
            print("Submitting job to Modal...")

        seq_matrix = np.ascontiguousarray(seq_matrix, dtype=np.int8)
        weights = np.ascontiguousarray(weights, dtype=np.float32)
        tip_dates = np.asarray(tip_dates, dtype=np.float64)

        with app.run():
            result = run_wsnj_modal.remote(
                seq_matrix.tobytes(),
                tip_dates.tobytes(),
                weights.tobytes(),
                seq_matrix.shape
            )

        tree_linkage = np.frombuffer(result["tree_bytes"], dtype=np.float64).reshape(result["tree_shape"])
        D_raw = np.frombuffer(result["D_raw_bytes"], dtype=np.float64).reshape(result["D_raw_shape"])

        results = {
            "mu": result["mu"],
            "d0": result["d0"],
            "r_squared": result["r_squared"],
            "n_negative": result["n_negative"],
        }

        if verbose:
            print("=" * 60)
            print("MODAL RESULTS")
            print("=" * 60)
            print(f"Estimated μ: {result['mu']:.4e} subst/site/year")
            print(f"Root-to-tip R²: {result['r_squared']:.4f}")
            print(f"Negative corrections: {result['n_negative']}")

        return tree_linkage, result["mu"], results, D_raw

    MODAL_AVAILABLE = True

except ImportError:
    MODAL_AVAILABLE = False

    def run_wsnj_on_modal(*args, **kwargs):
        raise ImportError("Modal not installed. Run: pip install modal")


# Legacy aliases
run_wsupgma_local = run_wsnj_local
run_wsupgma_on_modal = run_wsnj_on_modal if MODAL_AVAILABLE else lambda *a, **k: None


if __name__ == "__main__":
    import argparse
    import time

    parser = argparse.ArgumentParser(description="Weighted Serial-Sample NJ")
    parser.add_argument("--local", action="store_true", help="Run local benchmark")
    parser.add_argument("--test", action="store_true", help="Run with small test data")
    args = parser.parse_args()

    if args.local or args.test:
        n_samples = 100 if args.test else 500
        n_sites = 1000 if args.test else 5000

        print(f"Generating test data: {n_samples} samples × {n_sites} sites")
        np.random.seed(42)
        seq_matrix = np.random.randint(0, 4, (n_samples, n_sites), dtype=np.int8)
        missing_mask = np.random.random((n_samples, n_sites)) < 0.01
        seq_matrix[missing_mask] = 4

        weights = np.random.uniform(0.05, 0.35, n_sites).astype(np.float32)
        tip_dates = np.random.uniform(200, 10000, n_samples).astype(np.float64)

        print(f"\nRunning wsNJ...")
        t0 = time.time()

        tree, mu, results, D_raw = run_wsnj_local(
            seq_matrix, tip_dates, weights, verbose=True
        )

        elapsed = time.time() - t0
        print(f"\nTotal time: {elapsed:.1f}s")
