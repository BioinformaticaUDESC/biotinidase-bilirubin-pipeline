# Biotinidase–bilirubin pipeline

Research repository associated with the manuscript **Conformational dynamics of human biotinidase and the Arg211Cys variant in the presence of biocytin and bilirubins**.

## Study scope

The study combines multimethod molecular docking with molecular dynamics to examine human biotinidase, biocytin, unconjugated bilirubin, and bilirubin diglucuronide. The p.Arg211Cys variant is designated R170C in the 502-residue structural model used for the analyses.

Six protein–ligand conditions were simulated for 200 ns each, with one trajectory per condition and an aggregate simulation time of 1.2 microseconds.

| System | Protein | Ligand |
| --- | --- | --- |
| WT_BIO | Wild type | Biocytin |
| WT_BLI | Wild type | Unconjugated bilirubin |
| WT_BLD | Wild type | Bilirubin diglucuronide |
| R170C_BIO | Arg211Cys variant | Biocytin |
| R170C_BLI | Arg211Cys variant | Unconjugated bilirubin |
| R170C_BLD | Arg211Cys variant | Bilirubin diglucuronide |

## Computational workflow

1. Structural preparation, residue mapping, and variant generation.
2. Molecular docking with AutoDock Vina, GNINA, and AutoDock-GPU.
3. Pose comparison and geometric contact analysis.
4. Ligand parameterization, solvation, equilibration, and molecular dynamics with GROMACS.
5. Analysis of RMSD, RMSF, radius of gyration, solvent-accessible surface area, ligand contacts, PCA, and conformational clustering.

Software versions documented in the manuscript include AutoDock Vina 1.2.5, GNINA 1.3.2, and GROMACS 2024.2.

## Repository status

The original post-processing package BTD_PCA_Clustering_FEL_Z8 is available as a ZIP archive and as individual source files. It covers shared-basis C-alpha PCA, Ward hierarchical clustering, descriptive occupancy-based FEL projections, and RMSF sensitivity analysis excluding terminal residues.

Docking preparation scripts, MD production configurations, the upstream trajectory-preparation scripts, and simulation datasets have not yet been deposited. This is a partial source deposit, rather than a complete end-to-end reproduction package.

Reported observations are conditional on one trajectory per condition and do not establish substrate competition or enzymatic inhibition.

## Available files

| File | Purpose |
| --- | --- |
| `BTD_PCA_Clustering_FEL_Z8.zip` | Archived analysis package |
| `analyze_pca_fel.py` | Shared PCA, hierarchical clustering, representative frames, and exploratory FEL projections |
| `rmsf_core_qc.py` | RMSF comparison with terminal-residue sensitivity analysis |
| `run_on_Z8.sh` | GROMACS C-alpha extraction and analysis launcher |
| `README.txt` | Original package instructions |

The original package includes exploratory FEL plots; these are descriptive occupancy projections and should not be interpreted as converged binding free energies.

## Required inputs and execution

The launcher expects prepared trajectories for the six systems listed above, with the following project layout:

```text
PROJECT_ROOT/07_gromacs_md/systems/SYSTEM/md_200ns_rep01.tpr
PROJECT_ROOT/08_analysis/mdtraj_50_200ns/prepared/SYSTEM/aligned_50_200ns.xtc
PROJECT_ROOT/08_analysis/mdtraj_50_200ns/prepared/SYSTEM/groups.ndx
PROJECT_ROOT/08_analysis/mdtraj_50_200ns/results/SYSTEM_RMSF_CA_50_200ns.xvg
```

The original launcher assumes existing Conda environments named `gromacs_cpu` and `btd_interactions`. Python dependencies are MDTraj, NumPy, pandas, SciPy, scikit-learn, and Matplotlib. GROMACS must be available as `gmx` before launching.

From the directory containing the three source scripts:

```bash
conda activate gromacs_cpu
bash run_on_Z8.sh /path/to/PROJECT_ROOT
```

Outputs are written under `08_analysis/mdtraj_50_200ns/results/PCA_clustering_FEL_English/` and `08_analysis/mdtraj_50_200ns/results/RMSF_core_QC_English/`.

Syntax checks passed for the archived Python and shell scripts. No trajectory-based execution was performed as part of this repository deposit because the required simulation inputs are not included.

## Authors

Alessandra da Silva Santos; Juan Philippe Teixeira; Daniel Ferreira de Lima Neto; Chiara Ascencio Brancalion; Gustavo Felippe da Silva; Maria de Lourdes Borba Magalhaes; Rosiléia Marinho de Quadros; Miklos Maximiliano Bajay.

## Funding

Juan Philippe Teixeira received a scholarship from the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior–Brasil (CAPES), Finance Code 001. Alessandra da Silva Santos received a scholarship from PROMOP at the Universidade do Estado de Santa Catarina (UDESC).

## Contact

Juan Philippe Teixeira — juan.bioinformatics@gmail.com
