#!/usr/bin/env python3
"""Six-system C-alpha shared-basis PCA, hierarchical clustering and FEL.
Requires: mdtraj, numpy, scipy, scikit-learn, matplotlib, pandas.
Reads CA-only XTC/GRO generated from genuine GROMACS trajectories.
"""
from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault('MPLBACKEND', 'Agg')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.ndimage import gaussian_filter
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

SYSTEMS = ('WT_BIO', 'R170C_BIO', 'WT_BLD', 'R170C_BLD', 'WT_BLI', 'R170C_BLI')
SHORT = {'BIO': 'Biocytin', 'BLD': 'Bilirubin diglucuronide',
         'BLI': 'Unconjugated bilirubin'}


def savefig(fig, folder, stem):
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(folder / f'{stem}.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)


def load_data(prepared, minimum_ns=50, maximum_ns=200):
    trjs = {}
    for s in SYSTEMS:
        f = prepared / s / 'ca_50_200ns.xtc'
        top = prepared / s / 'ca_reference_50ns.gro'
        if not f.is_file() or not top.is_file():
            raise FileNotFoundError(f'Missing {s} CA-only XTC/GRO: {f} or {top}')
        traj = md.load_xtc(str(f), top=str(top))
        if traj.n_frames < 15:
            raise RuntimeError(f'Insufficient frames for {s}: {traj.n_frames}')
        if traj.n_atoms < 100 or any(a.name != 'CA' for a in traj.top.atoms):
            raise ValueError(f'{s}: topology must contain only protein C-alpha atoms')
        if np.any(np.diff(traj.time) <= 0):
            raise ValueError(f'{s}: non-increasing trajectory timestamps')
        ns = traj.time / 1000.0
        keep = np.flatnonzero((ns >= minimum_ns - 1e-4) & (ns <= maximum_ns + 1e-4))
        if len(keep) < 15 or abs(ns[keep[0]] - minimum_ns) > 0.5 or abs(ns[keep[-1]]- maximum_ns) > 0.5:
            raise RuntimeError(f'{s}: trajectory does not span {minimum_ns}–{maximum_ns} ns; range={ns.min():.3f}–{ns.max():.3f}')
        traj = traj.slice(keep)
        trjs[s] = traj
        print(f'LOADED {s}: {traj.n_frames} frames, {traj.n_atoms} CA atoms, {traj.time[0]/1000:.2f}–{traj.time[-1]/1000:.2f} ns', flush=True)
    na = trjs[SYSTEMS[0]].n_atoms
    for s in SYSTEMS:
        if trjs[s].n_atoms != na:
            raise ValueError('CA atom counts differ across systems; cannot compute shared PCA')
    return trjs


def analyze(project):
    work = project / '08_analysis' / 'mdtraj_50_200ns'
    prepared = work / 'prepared'
    out = work / 'results' / 'PCA_clustering_FEL_English'
    out.mkdir(parents=True, exist_ok=True)
    trjs = load_data(prepared)
    na = trjs[SYSTEMS[0]].n_atoms
    # Globally superpose every system to the *same* WT_BIO frame.
    # Exclude 10 residues from each terminus for the least-squares fit.
    ref = trjs['WT_BIO'][0]
    fit = np.arange(10, na - 10, dtype=int)
    if len(fit) < 100:
        raise ValueError('Unexpectedly short protein chain')
    for s in SYSTEMS:
        trjs[s].superpose(ref, frame=0, atom_indices=fit, ref_atom_indices=fit)

    count = {s: trjs[s].n_frames for s in SYSTEMS}
    x = np.concatenate([trjs[s].xyz.reshape(count[s], -1) for s in SYSTEMS], axis=0)
    # Subtract the *shared* mean; per-system centering would erase inter-system shifts.
    pca = PCA(n_components=10, svd_solver='randomized', random_state=314159)
    score = pca.fit_transform(x)
    explained = pca.explained_variance_ratio_
    offsets = {}
    start = 0
    for s in SYSTEMS:
        offsets[s] = score[start:start + count[s]]
        start += count[s]
    with open(out / 'global_PCA_variance.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['PC', 'variance_explained_fraction', 'cumulative_fraction'])
        for i, v in enumerate(explained):
            w.writerow([f'PC{i+1}', float(v), float(np.sum(explained[:i+1]))])

    labels = {'WT': 'WT', 'R170C': 'R170C'}
    fig, ax = plt.subplots(figsize=(9.0, 7.0))
    markers = {'WT': 'o', 'R170C': '^'}
    for s in SYSTEMS:
        p = offsets[s]
        lig = s.split('_')[1]
        strain = s.split('_')[0]
        ax.scatter(p[::5, 0], p[::5, 1], s=9, alpha=0.5,
                   marker=markers[strain], label=f'{strain} / {lig}')
    ax.set_xlabel(f'PC1 ({100*explained[0]:.1f}% variance)')
    ax.set_ylabel(f'PC2 ({100*explained[1]:.1f}% variance)')
    ax.set_title('Shared-basis Cα PCA of six biotinidase systems (50–200 ns)')
    ax.legend(fontsize=8, frameon=False, ncol=2)
    savefig(fig, out, 'Figure_PCA_shared_6systems_English')

    fig, axes = plt.subplots(3, 1, figsize=(9.5, 11.5), sharex=True, sharey=True, layout='constrained')
    for ax, lig in zip(axes, ('BIO', 'BLD', 'BLI')):
        for genotype in ('WT', 'R170C'):
            s = f'{genotype}_{lig}'
            p = offsets[s]
            ax.scatter(p[::4, 0], p[::4, 1], s=10, alpha=0.45,
                       marker=markers[genotype], label=genotype)
        ax.set_title(SHORT[lig]); ax.set_ylabel('PC2 (nm)')
        ax.legend(frameon=False)
    axes[-1].set_xlabel('PC1 (nm)')
    fig.suptitle('Cα conformational sampling: WT versus R170C (shared PCs)')
    savefig(fig, out, 'Figure_PCA_WT_vs_R170C_by_ligand_English')

    # FEL = -ln(P/P_max) in kBT units, from smoothed 2D occupancy.
    # A descriptive single-trajectory projection; not a converged thermodynamic PMF.
    all_pc = score[:, :2]
    lower = all_pc.min(axis=0)
    upper = all_pc.max(axis=0)
    pad = np.maximum((upper - lower) * 0.04, 0.01)
    bounds = [[lower[0]-pad[0], upper[0]+pad[0]], [lower[1]-pad[1], upper[1]+pad[1]]]
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.5), sharex=True, sharey=True, layout='constrained')
    norm = matplotlib.colors.Normalize(vmin=0, vmax=6)
    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad('#f4f4f4')
    for ax, s in zip(axes.flat, SYSTEMS):
        pts = offsets[s]
        h, xe, ye = np.histogram2d(pts[:, 0], pts[:, 1], bins=40, range=bounds)
        sm = gaussian_filter(h, sigma=1.1)
        rel = sm / sm.max()
        eng = -np.log(np.clip(rel, 1e-10, None))
        eng = np.ma.masked_where(rel < np.exp(-6.0), eng)
        ax.pcolormesh(xe, ye, eng.T, shading='auto', cmap=cmap, norm=norm, rasterized=True)
        ax.set_title(s)
        ax.set_xlabel('PC1 (nm)'); ax.set_ylabel('PC2 (nm)')
    smap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(smap, ax=list(axes.flat), label='Relative free energy (kBT)', shrink=0.82)
    fig.suptitle('Projected conformational landscapes (50–200 ns; shared PCA axes)')
    savefig(fig, out, 'Figure_FEL_sixsystems_English')

    cluster_rows, representative_rows = [], []
    for s in SYSTEMS:
        times_ns = trjs[s].time / 1000
        dt_ns = float(np.median(np.diff(times_ns)))
        step = max(1, int(round(1.0 / dt_ns)))  # ~1ns sampling to limit near-duplicate frames
        ids = np.arange(0, count[s], step, dtype=int)
        sample = offsets[s][ids, :]
        cum = np.cumsum(explained)
        ndim = int(np.clip(np.searchsorted(cum, 0.8)+1, 3, 10))
        coords = sample[:, :ndim]
        tree = linkage(coords, method='ward', metric='euclidean')
        possible = range(2, min(7, len(ids)-1))
        options = []
        for k in possible:
            lab = fcluster(tree, k, criterion='maxclust')
            if len(set(lab)) != k:
                continue
            options.append((silhouette_score(coords, lab), k, lab))
        if not options:
            raise RuntimeError(f'Clustering did not yield valid labels for {s}')
        best_sil, kbest, labels_best = max(options, key=lambda v: v[0])
        # Representative sampled frame: medoid in PCA subspace among frames of this cluster.
        for cluster_id in sorted(np.unique(labels_best)):
            member_idxs = np.flatnonzero(labels_best == cluster_id)
            sub = coords[member_idxs]
            distances = np.linalg.norm(sub[:, None, :] - sub[None, :, :], axis=2)
            medoid_local = member_idxs[int(np.argmin(distances.sum(axis=1)))]
            frame_global = int(ids[medoid_local])
            population = len(member_idxs) / len(ids)
            cluster_rows.append(dict(system=s, cluster=int(cluster_id), n_sampled=int(len(member_idxs)),
                                     sample_fraction=float(population), n_clusters=int(kbest),
                                     silhouette=float(best_sil), sampling_interval_ns=float(step * dt_ns)))
            representative_rows.append(dict(system=s, cluster=int(cluster_id),
                                             population_sampled_pct=float(100*population),
                                             frame_index=frame_global,
                                             time_ns=float(times_ns[frame_global]),
                                             time_ps=float(trjs[s].time[frame_global]),
                                             n_pcs=int(ndim), silhouette=float(best_sil)))
        fig, ax = plt.subplots(figsize=(9.0, 4.7))
        dendrogram(tree, ax=ax, no_labels=True, truncate_mode='lastp', p=22,
                   color_threshold=0, above_threshold_color='gray')
        ax.set_title(f'{s}: hierarchical clustering (Ward, ~1 ns sampling)')
        ax.set_xlabel('Truncated dendrogram leaves')
        ax.set_ylabel('Euclidean distance in shared PCA space (nm)')
        savefig(fig, out, f'Figure_Dendrogram_{s}_English')
        print(f'CLUSTER {s}: k={kbest}, silhouette={best_sil:.3f}, PCA dimensions={ndim}, sampled={len(ids)}', flush=True)

    pd.DataFrame(cluster_rows).to_csv(out / 'hierarchical_cluster_populations.csv', index=False)
    pd.DataFrame(representative_rows).to_csv(out / 'representative_frames_for_GROMACS.csv', index=False)
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), layout='constrained')
    for ax, s in zip(axes.flat, SYSTEMS):
        rows = [r for r in cluster_rows if r['system'] == s]
        ax.bar([str(r['cluster']) for r in rows], [100*r['sample_fraction'] for r in rows])
        ax.set_title(s); ax.set_xlabel('Cluster')
        ax.set_ylabel('Sampled frames (%)')
        ax.set_ylim(0, 100)
    fig.suptitle('Hierarchical cluster populations (1 ns subsampling)')
    savefig(fig, out, 'Figure_Cluster_populations_sixsystems_English')

    meta = {
        'systems': list(SYSTEMS), 'n_atoms_CA': na,
        'n_frames_per_system': count,
        'time_window_ns': [50, 200], 'global_fit_reference': 'WT_BIO at 50 ns',
        'superpose_core': f'CA indices 11 through {na-10} (terminal 10 excluded per side)',
        'PCA': 'shared PCA fitted to all aligned CA coordinates of all six systems',
        'explained_variance_ratio_first_10': explained.tolist(),
        'FEL': 'dimensionless -ln(relative 2D histogram density), 40 bins, Gaussian sigma 1.1, not a converged binding free energy',
        'clustering': 'Ward in first 3–10 shared PCs; frame thinning ~1ns; silhouette selected k=2..6',
        'limitations': [
          'one simulation per system, correlated MD frames and no between-replicate uncertainty',
          'PCA/FEL are descriptive; cluster populations not independent observations',
          'CA-only snapshots not suitable for MM/PBSA/MMGBSA or ligand interactions',
          'shared reference may be influenced by structural alignment; verify core fit and PBC visually',
        ]
    }
    (out/'methodology_and_qc.json').write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding='utf-8')
    with open(out / 'RESULTS_OVERVIEW.txt', 'w') as file:
        print('SHARED PCA 50–200 ns', file=file)
        print('PC1 variance explained: %.2f%%' % (explained[0]*100), file=file)
        print('PC2 variance explained: %.2f%%' % (explained[1]*100), file=file)
        print('PC1+PC2: %.2f%%' % (explained[:2].sum()*100), file=file)
        print('\nCLUSTERS (1ns-subsampled frames):', file=file)
        for s in SYSTEMS:
            rows = [r for r in cluster_rows if r['system']==s]
            print(f"{s}: {len(rows)} clusters; silhouette={rows[0]['silhouette']:.3f}; largest population={max(r['sample_fraction'] for r in rows):.2%}", file=file)
        print('\nNo mechanistic/inhibitory conclusions from PCA/FEL alone.', file=file)
    print('DONE: results:', out, flush=True)
    print((out / 'RESULTS_OVERVIEW.txt').read_text(), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('project', type=Path)
    args = parser.parse_args()
    analyze(args.project)
