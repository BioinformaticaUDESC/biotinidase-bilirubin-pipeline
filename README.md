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
6. Complementary MM/PBSA calculations with residue decomposition, cluster-window sensitivity analyses, and higher-resolution shared-basis PCA/FEL analysis.

Software versions documented in the manuscript include AutoDock Vina 1.2.5, GNINA 1.3.2, and GROMACS 2024.2.

## Repository status

The original post-processing package `BTD_PCA_Clustering_FEL_Z8.zip` remains available as an archive and as individual source files. It covers shared-basis C-alpha PCA, Ward hierarchical clustering, descriptive occupancy-based FEL projections, and RMSF sensitivity analysis excluding terminal residues.

On 16 November 2026, the original source export `BTD_V3_20261116_162456_jc5hltty_DADOS_PARA_PAPER.zip` was added with the complementary numerical results used in manuscript V4. Its V3 filename is preserved for provenance; it is not the separately assembled `SupplementaryData1_BTD_V4.zip` manuscript supplement. The source export includes derived data, selected analysis code, and the original plotting outputs. Source-figure layouts were subsequently revised for manuscript V4 without changing these numerical data.

Original molecular dynamics trajectories, full production `.tpr` files, complete simulation topologies, docking preparation scripts, and the upstream trajectory-preparation scripts are not included. This deposit supports verification of the exported post-processing results and is not a complete end-to-end simulation reproduction package.

Reported observations are conditional on one trajectory per condition and do not establish substrate competition or enzymatic inhibition.

## Available files

| File | Purpose |
| --- | --- |
| [BTD_PCA_Clustering_FEL_Z8.zip](BTD_PCA_Clustering_FEL_Z8.zip) | Original archived analysis package |
| [BTD_V3_20261116_162456_jc5hltty_DADOS_PARA_PAPER.zip](BTD_V3_20261116_162456_jc5hltty_DADOS_PARA_PAPER.zip) | Complementary derived source data used in manuscript V4, including the nested PCA/FEL V2 results package |
| `analyze_pca_fel.py` | Original shared PCA, hierarchical clustering, representative frames, and exploratory FEL projections |
| `rmsf_core_qc.py` | RMSF comparison with terminal-residue sensitivity analysis |
| `run_on_Z8.sh` | Original GROMACS C-alpha extraction and analysis launcher |
| `README.txt` | Original package instructions |

## Complementary data used in manuscript V4

The [source archive at deposit commit 8daa414](https://github.com/BioinformaticaUDESC/biotinidase-bilirubin-pipeline/blob/8daa414ea004002a261f73848eb6d81e00fd7053/BTD_V3_20261116_162456_jc5hltty_DADOS_PARA_PAPER.zip) contains:

| Archive location | Contents |
| --- | --- |
| `configuracoes_e_resultados/` | MM/PBSA input files and final summary/per-frame outputs for all six systems under both `serial_` and `cluster_` sampling schemes |
| `metricas_50_200ns/` | Rg, SASA, receptor–ligand COM distance, and minimum atom-distance XVG series for all six systems |
| `tabelas_primarias/` | MM/PBSA component and residue-decomposition summary tables, including model Arg497 (canonical Arg538) |
| `resultados_primarios/` | Copies of the primary MM/PBSA summary outputs |
| `Table_Arg497_BTD_V3.csv` | Exported Arg497 component and total values |
| `proveniencia_BTD_V3.json` | Source-file provenance records |
| `PCA_FEL_V2/BTD_PCA_FEL_V2_20261116_124531_8dk2_xgn_RESULTADOS.zip` | Shared PCA basis, projections, FEL grids, methods, code, metadata, and plotting outputs |

Extract the outer archive first, then extract the nested PCA/FEL ZIP. That nested package includes `BTD_PCA_base_comum_V2.npz`, `BTD_FEL_grades_comuns_V2.npz`, six PC1–PC2 projection tables, variance and quality-control tables, provenance metadata, GROMACS logs, and `btd_pca_fel_comum_v2.py`.

### Shared PCA and occupancy-based FEL

The complementary analysis uses 1,501 frames per system from 50–200 ns at 100 ps spacing (9,006 frames in total) and 502 corresponding C-alpha atoms. Coordinates were fitted by unweighted Kabsch alignment to a common reference, the first WT_BIO frame at 50 ns. A single Cartesian mean and covariance matrix define the common PCA basis. Coordinates were not mass weighted or z-score standardized.

PC1 and PC2 explain 48.61% and 10.51% of the pooled variance, respectively (59.12% combined). Every system is projected onto the same eigenvectors. The original, more sparsely sampled PCA/clustering analysis used 151 frames per system and remains a separate analysis; its variance percentages should not be substituted for those of the complementary PCA.

FEL values use `Delta F = -R T ln[P(PC1, PC2)/Pmax]` at 310 K on a common 38 × 38 histogram grid. Each system has its own zero at its most populated bin. Empty bins are masked; no pseudocount or kernel-density smoothing is applied. The plotted common color range is 0–4 kcal/mol. These are descriptive conformational occupancy landscapes, not converged binding free energies. Separate per-system zeros do not establish absolute free-energy differences between systems.

### MM/PBSA sampling and uncertainty

The primary `serial_` calculations sample the common 50–200 ns interval at 1 ns spacing (151 frames per system). The `cluster_` calculations use the system-specific continuous intervals below at the same spacing. Energies exclude conformational entropy. The MM/PBSA input temperature is 298.15 K; the MD reference temperature and occupancy-based FEL temperature are 310 K. These settings are documented separately and were not retrospectively changed.

| System | Primary 50–200 ns mean ± block SEM (kcal/mol) | Cluster interval (ns) | Cluster frames | Cluster mean ± block SEM (kcal/mol) |
| --- | ---: | ---: | ---: | ---: |
| WT_BIO | -34.87 ± 1.50 | 102–200 | 99 | -37.15 ± 0.80 |
| WT_BLI | -16.29 ± 0.24 | 50–145 | 96 | -15.94 ± 0.42 |
| WT_BLD | -16.04 ± 2.26 | 80–200 | 121 | -18.03 ± 2.14 |
| R170C_BIO | -37.60 ± 0.36 | 81–200 | 120 | -37.82 ± 0.47 |
| R170C_BLI | -14.89 ± 0.62 | 50–107 | 58 | -13.96 ± 1.21 |
| R170C_BLD | -14.03 ± 1.65 | 84–200 | 117 | -12.11 ± 1.66 |

Block SEM uses eight equal consecutive blocks of length `floor(n/8)`. Trailing remainder frames are excluded only from this block-SEM calculation; reported means retain all sampled frames. For 151 frames, the SEM uses eight blocks of 18 frames. This within-trajectory uncertainty does not represent independent simulation replicates. Because cluster intervals differ across systems, their differences are sampling-sensitivity comparisons, not matched-window mutation free-energy estimates.

### Proximity and time-series interpretation

COM distance is measured between the entire receptor COM and the ligand COM, not the center of a predefined pocket. Minimum distances refer to the selected receptor and ligand atoms and include hydrogens. Proximity plots use identical COM axes (0–3 nm) and identical minimum-distance axes across systems. Structural-descriptor moving means use five frames sampled every 100 ps; MM/PBSA moving means use five energy samples sampled every 1 ns. Boundary handling avoids zero-padding artifacts.

A smaller COM distance alone does not demonstrate stronger binding or persistent occupation of a particular pocket. Per-frame energies and proximity are complementary descriptive measures; neither establishes a dissociation constant, inhibitory activity, or kinetic residence time.

## Required inputs and execution of the original package

The original `run_on_Z8.sh` launcher expects prepared trajectories for the six systems listed above, with the following project layout:

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

Syntax checks passed for the original archived Python and shell scripts. The original repository deposit did not execute them against trajectories because those required simulation inputs are not included. The complementary archive now supplies exported numerical results and available code, but still does not include the full coordinates needed to independently repeat trajectory extraction and alignment.

## Data and code availability

Post-processing scripts and complementary derived numerical data for this study are publicly available at https://github.com/BioinformaticaUDESC/biotinidase-bilirubin-pipeline. The complementary source export deposited at commit `8daa414ea004002a261f73848eb6d81e00fd7053` contains the shared PCA mean, fitting reference and eigenvectors, per-frame PC1–PC2 projections, FEL histogram counts and relative energies, Rg/SASA and proximity series, MM/PBSA inputs and per-frame results for both sampling schemes, residue-decomposition summaries, and source-provenance records. Original MD trajectories, full production TPR files, and complete simulation topologies are not included. The deposited data support verification of the exported post-processing results but do not permit an independent rerun of the original simulations.

## Authors

Alessandra da Silva Santos; Juan Philippe Teixeira; Bárbara Miranda-Vilela; Daniel Ferreira de Lima Neto; Chiara Ascencio Brancalion; Gustavo Felippe da Silva; Maria de Lourdes Borba Magalhaes; Rosiléia Marinho de Quadros; Miklos Maximiliano Bajay.

## Funding

Juan Philippe Teixeira received a scholarship from the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior–Brasil (CAPES), Finance Code 001. Alessandra da Silva Santos received a scholarship from PROMOP at the Universidade do Estado de Santa Catarina (UDESC).

## Contact

Juan Philippe Teixeira — juan.teixeira0863@edu.udesc.br (institutional); juan.bioinformatics@gmail.com (additional contact).
