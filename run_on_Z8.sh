#!/usr/bin/env bash
# Entrar com conda activate gromacs_cpu para GROMACS, script ativa btd_interactions depois.
set -Eeuo pipefail
PROJ="${1:-/dados/projects/biotinidase_bilirrubina_importado}"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
OUT="$PROJ/08_analysis/mdtraj_50_200ns"
mkdir -p "$OUT/prepared" "$OUT/results" "$OUT/logs"

if ! command -v gmx >/dev/null 2>&1; then
    echo 'ERRO: gmx ausente. Execute: conda activate gromacs_cpu' >&2
    exit 2
fi
for s in WT_BIO R170C_BIO WT_BLD R170C_BLD WT_BLI R170C_BLI; do
    pre="$OUT/prepared/$s"
    mkdir -p "$pre"
    ndx="$pre/groups.ndx"
    input="$pre/aligned_50_200ns.xtc"
    tpr="$PROJ/07_gromacs_md/systems/$s/md_200ns_rep01.tpr"
    # These inputs were created by the RMSF pipeline. Re-use them, do not reprocess full XTC.
    [[ -s "$input" && -s "$ndx" && -s "$tpr" ]] || {
        echo "ERRO: faltam aligned_50_200ns.xtc, groups.ndx ou TPR em $s" >&2; exit 3;
    }
    ca_xtc="$pre/ca_50_200ns.xtc"
    ca_gro="$pre/ca_reference_50ns.gro"
    if [[ ! -s "$ca_xtc" ]]; then
        printf 'C-alpha\n' | gmx trjconv -s "$tpr" -f "$input" -n "$ndx" -o "$ca_xtc" \
            > "$OUT/logs/${s}_CA_only.log" 2>&1
    fi
    if [[ ! -s "$ca_gro" ]]; then
        printf 'C-alpha\n' | gmx trjconv -s "$tpr" -f "$input" -n "$ndx" -dump 50000 -o "$ca_gro" \
            > "$OUT/logs/${s}_CA_reference.log" 2>&1
    fi
    echo "CA PREPARED: $s"
done

# Use conda bash initialization even if invoking this from a zsh session.
if ! command -v conda >/dev/null 2>&1; then
   echo 'ERRO: conda indisponível' >&2; exit 4
fi
BASE="$(conda info --base)"
# shellcheck disable=SC1090
source "$BASE/etc/profile.d/conda.sh"
conda activate btd_interactions
python - <<'PY'
try:
    import mdtraj, sklearn, scipy, pandas, numpy, matplotlib
except ImportError as exc:
    raise SystemExit(f'Falta dependência Python: {exc}; instale no ambiente btd_interactions')
print('Python analysis dependencies OK')
PY
python "$SCRIPT_DIR/analyze_pca_fel.py" "$PROJ" \
   2>&1 | tee "$OUT/logs/PCA_FEL_clustering_execution.log"
python "$SCRIPT_DIR/rmsf_core_qc.py" "$PROJ" \
   2>&1 | tee "$OUT/logs/RMSF_core_QC_execution.log"
