#!/usr/bin/env bash
# =============================================================================
#  emu_classification.sh
# -----------------------------------------------------------------------------
#  Clasificacion taxonomica con EMU de los FASTQ ya libres de quimeras, y
#  combinacion de los perfiles por muestra en tablas por nivel taxonomico.
#
#  Uso: emu_classification.sh <NOMBRE_CORRIDA> [THREADS]
#
#  Estructura esperada (rutas genericas, editar BASE segun corresponda):
#    $BASE/<NOMBRE_CORRIDA>/final/<barcode>.sin_quimeras.fastq.gz   (entrada)
#    $BASE/<NOMBRE_CORRIDA>/emu/                                    (salida)
#
#  Decisiones de diseno:
#    1) Activa el entorno de conda por su cuenta: 'emu' no esta en 'base'. Al
#       lanzar con nohup o en shell no interactivo, la activacion no se hereda.
#    2) No usa '| head' en ningun sitio: con 'set -o pipefail', un head que corta
#       la tuberia manda SIGPIPE al productor (salida 141) -> fallo intermitente.
#    3) Es reanudable: salta cualquier muestra que ya tenga su _rel-abundance.tsv.
# =============================================================================
set -euo pipefail

CORRIDA="${1:-}"; THREADS="${2:-8}"
BASE="${BASE:-$HOME/project/runs}"      # ruta generica; editar segun corresponda

if [ -z "$CORRIDA" ]; then
    echo "Uso: emu_classification.sh <NOMBRE_CORRIDA> [THREADS]"
    echo "Corridas disponibles:"
    ls -1 "$BASE" 2>/dev/null || echo "  (ninguna)"
    exit 1
fi

IN="$BASE/$CORRIDA/final"
OUT="$BASE/$CORRIDA/emu"
[ -d "$IN" ] || { echo "ERROR: no existe $IN"; exit 1; }

CONDA_ENV="${CONDA_ENV:-emu}"
if ! command -v emu >/dev/null 2>&1; then
    CB="$(conda info --base 2>/dev/null || echo "$HOME/miniforge3")"
    set +u
    source "$CB/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
    set -u
fi
command -v emu >/dev/null 2>&1 \
    || { echo "ERROR: no encuentro 'emu' ni en PATH ni en el entorno '$CONDA_ENV'"; exit 1; }

[ -n "${EMU_DATABASE_DIR:-}" ] || EMU_DATABASE_DIR="$HOME/emu_db"
export EMU_DATABASE_DIR
[ -s "$EMU_DATABASE_DIR/species_taxid.fasta" ] \
    || { echo "ERROR: no encuentro species_taxid.fasta en $EMU_DATABASE_DIR"; exit 1; }

mkdir -p "$OUT"

EMU_VER="$(emu --version 2>&1 || true)"; EMU_VER="${EMU_VER%%$'\n'*}"

{   echo "Fecha   : $(date '+%Y-%m-%d %H:%M:%S')"
    echo "emu     : $EMU_VER"
    echo "entorno : $CONDA_ENV"
    echo "DB      : $EMU_DATABASE_DIR"
    ls -la "$EMU_DATABASE_DIR"
} > "$OUT/emu_entorno.txt"

tot=$(find "$IN" -maxdepth 1 -name '*.fastq.gz' | wc -l)
[ "$tot" -gt 0 ] || { echo "ERROR: no hay FASTQ en $IN"; exit 1; }

echo "=== $CORRIDA : $tot muestras, $THREADS hilos, DB=$EMU_DATABASE_DIR ==="
n=0
for f in "$IN"/*.fastq.gz; do
    b=$(basename "$f"); b="${b%.sin_quimeras.fastq.gz}"; n=$((n+1))
    if [ -s "$OUT/${b}_rel-abundance.tsv" ]; then
        echo "[$n/$tot] $b ya hecho, se salta"; continue
    fi
    echo "[$n/$tot] $b ..."
    emu abundance "$f" --type map-ont --db "$EMU_DATABASE_DIR" \
        --output-dir "$OUT" --output-basename "$b" --threads "$THREADS" --keep-counts
done

echo "=== combinando tablas por nivel taxonomico ==="
for rank in phylum class order family genus species; do
    emu combine-outputs "$OUT" "$rank" --counts
    emu combine-outputs "$OUT" "$rank"
done

echo ">>> Listo. Resultados en: $OUT"
ls -1 "$OUT"
