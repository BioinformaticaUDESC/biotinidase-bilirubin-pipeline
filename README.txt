BTD: Genuine trajectory analysis (Z8) — PCA, hierarchical clustering, FEL

Inputs confirmed to exist on the user's Z8:
- /dados/projects/biotinidase_bilirrubina_importado/07_gromacs_md/systems/{SYSTEM}/md_200ns_rep01.tpr
- /dados/projects/biotinidase_bilirrubina_importado/08_analysis/mdtraj_50_200ns/prepared/{SYSTEM}/aligned_50_200ns.xtc
- corresponding groups.ndx

Run in terminal on Z8:
  cd /home/juan
  unzip -o BTD_PCA_Clustering_FEL_Z8.zip
  conda activate gromacs_cpu
  bash /home/juan/BTD_PCA_Clustering_FEL_Z8/run_on_Z8.sh /dados/projects/biotinidase_bilirrubina_importado

If Python reports a missing package:
  conda activate btd_interactions
  conda install -c conda-forge mdtraj scipy scikit-learn matplotlib pandas numpy -y
  conda activate gromacs_cpu
  bash /home/juan/BTD_PCA_Clustering_FEL_Z8/run_on_Z8.sh /dados/projects/biotinidase_bilirrubina_importado

Outputs (in English):
  08_analysis/mdtraj_50_200ns/results/PCA_clustering_FEL_English/
  - Figure_PCA_shared_6systems_English.{png,svg,pdf}
  - Figure_PCA_WT_vs_R170C_by_ligand_English.{png,svg,pdf}
  - Figure_FEL_sixsystems_English.{png,svg,pdf}
  - Figure_Dendrogram_{SYSTEM}_English.{png,svg,pdf}
  - Figure_Cluster_populations_sixsystems_English.{png,svg,pdf}
  - global_PCA_variance.csv
  - hierarchical_cluster_populations.csv
  - representative_frames_for_GROMACS.csv
  - RESULTS_OVERVIEW.txt
  - methodology_and_qc.json

Scientific caveats:
- PCA is global shared basis on aligned C-alpha atoms, starting at 50ns.
- Relative FEL (-ln occupancy, kBT units) is descriptive, not evidence of converged free energy or competitive inhibition.
- Cluster population estimates are based on one trajectory per system; they are not independent biological replicates.
- Representative frames must be extracted from original aligned full-system trajectories for molecular depiction or MM/GBSA.
- RMSF terminal hotspots (residues 1/3) require interpretation with termini excluded sensitivity analysis.
