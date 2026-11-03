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

The study overview is available in this initial repository. Original execution scripts, configuration files, and reproducibility instructions are pending upload and verification. No executable pipeline release or simulation dataset has been deposited here yet.

Reported observations are conditional on one trajectory per condition and do not establish substrate competition or enzymatic inhibition.

## Authors

Alessandra da Silva Santos; Juan Philippe Teixeira; Daniel Ferreira de Lima Neto; Chiara Ascencio Brancalion; Gustavo Felippe da Silva; Maria de Lourdes Borba Magalhaes; Rosiléia Marinho de Quadros; Miklos Maximiliano Bajay.

## Funding

Juan Philippe Teixeira received a scholarship from the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior–Brasil (CAPES), Finance Code 001. Alessandra da Silva Santos received a scholarship from PROMOP at the Universidade do Estado de Santa Catarina (UDESC).

## Contact

Juan Philippe Teixeira — juan.bioinformatics@gmail.com
