#!/usr/bin/env python3
"""RMSF core sensitivity check using already calculated genuine GROMACS .xvg files."""
import argparse, os
from pathlib import Path
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument('project', type=Path)
args = p.parse_args()
root = args.project / '08_analysis' / 'mdtraj_50_200ns' / 'results'
output = root / 'RMSF_core_QC_English'
output.mkdir(parents=True, exist_ok=True)
rows=[]
fig, axes = plt.subplots(3, 1, figsize=(11, 8.5), sharex=True, layout='constrained')
for ax, lig in zip(axes, ['BIO','BLD','BLI']):
    data = {}
    for g in ['WT','R170C']:
        fn = root / f'{g}_{lig}_RMSF_CA_50_200ns.xvg'
        a = np.loadtxt(fn, comments=('#', '@'))
        if a.ndim != 2 or a.shape[1] < 2:
            raise RuntimeError(f'Unexpected XVG content: {fn}')
        data[g] = a
    if data['WT'].shape != data['R170C'].shape or not np.array_equal(data['WT'][:,0],data['R170C'][:,0]):
        raise RuntimeError(f'Residue index mismatch for {lig}. Do not compare without residue mapping.')
    residues = data['WT'][:,0].astype(int)
    n = len(residues)
    core = np.zeros(n, dtype=bool)
    core[10:n-10] = True
    diff = data['R170C'][:,1] - data['WT'][:,1]
    ax.plot(residues, diff, linewidth=1.1, label='R170C − WT')
    ax.axhline(0, color='gray', linewidth=.8)
    ax.axvline(170, color='gray', linestyle='--', linewidth=.8)
    ax.axvspan(float(residues.min()), float(residues[9]), color='gray', alpha=.10)
    ax.axvspan(float(residues[-10]), float(residues.max()), color='gray', alpha=.10)
    ax.set_title(lig)
    ax.set_ylabel('Δ Cα RMSF (nm)')
    for g in ('WT', 'R170C'):
        rmsf = data[g][:,1]
        topidx = np.argsort(rmsf)[::-1][:8]
        rows.append(dict(system=f'{g}_{lig}', n_residues=n,
                         mean_all_nm=float(rmsf.mean()),
                         mean_core_nm=float(rmsf[core].mean()),
                         median_core_nm=float(np.median(rmsf[core])),
                         max_core_nm=float(rmsf[core].max()),
                         max_core_resid=int(residues[core][rmsf[core].argmax()]),
                         top_8_residues=','.join(map(str,residues[topidx]))))
    pd.DataFrame({'residue_index':residues, 'WT_RMSF_nm':data['WT'][:,1],
                  'R170C_RMSF_nm':data['R170C'][:,1],
                  'R170C_minus_WT_RMSF_nm':diff,
                  'core_excluding_terminal_10':core}).to_csv(output/f'RMSF_difference_{lig}.csv',index=False)
axes[-1].set_xlabel('Mature biotinidase residue index')
fig.suptitle('Residue-level RMSF differences, mutant minus WT (50–200 ns)')
for ext in ('png','svg','pdf'):
    fig.savefig(output/f'Figure_RMSF_mutant_minus_WT_terminal_QC.{ext}',dpi=300,bbox_inches='tight')
plt.close(fig)
df=pd.DataFrame(rows)
df.to_csv(output/'RMSF_core_versus_full_summary.csv', index=False)
print('RMSF core QC completed:')
print(df[['system','mean_all_nm','mean_core_nm','max_core_resid']].to_string(index=False))
