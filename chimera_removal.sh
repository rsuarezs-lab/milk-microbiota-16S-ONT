#!/usr/bin/env bash
# =============================================================================
#  chimera_removal.sh
# -----------------------------------------------------------------------------
#  Preparacion de FASTQ para EMU: subsampling reproducible + deteccion y
#  eliminacion de quimeras con VSEARCH/UCHIME, para 16S full-length de Nanopore.
#
#  FLUJO:
#    FASTQ original
#      -> seqkit (subsampling reproducible, <= MAX_READS)
#      -> FASTQ subsampleado                        [SE CONSERVA]
#      -> FASTA temporal (solo para VSEARCH)
#      -> [opcional] vsearch --orient (solo el FASTA)
#      -> vsearch --derep_fulllength   (contabilidad de duplicados exactos)
#      -> vsearch --cluster_size       (genera estructura de ABUNDANCIA real)
#      -> vsearch --uchime_denovo      (sobre los centroides)
#      -> centroides quimericos -> TODOS sus reads miembros (via .uc)
#      -> seqkit grep -v sobre el FASTQ subsampleado
#      -> FASTQ final sin quimeras                  [ENTRADA DE EMU]
#
#  GARANTIA: el FASTQ final es EXACTAMENTE el subsampleado menos las lecturas
#  quimericas. No se modifican secuencias, IDs ni quality scores. Se valida.
#
#  NOTA METODOLOGICA IMPORTANTE:
#    NO se usa derep_fulllength como entrada de uchime_denovo. En Nanopore las
#    lecturas casi nunca son identicas byte a byte, por lo que la dereplicacion
#    produce ~1 secuencia unica por read (todas con size=1). Como uchime_denovo
#    exige que los padres sean >= --abskew (2x) mas abundantes que el candidato,
#    con todo a size=1 NUNCA se detectaria ninguna quimera (0% falso).
#    Por eso se usa clustering al CLUSTER_ID de identidad, que si genera
#    abundancias reales, y despues se propaga el resultado a todos los miembros.
# =============================================================================

set -Eeuo pipefail
export LC_ALL=C          # imprescindible: sort/join/comm deben usar el mismo orden

# =============================================================================
#  1. CONFIGURACION  (todo se puede sobreescribir con variables de entorno)
# =============================================================================

# Rutas genericas: reemplazar por las rutas reales o pasarlas como variables de entorno.
INPUT_DIR="${INPUT_DIR:-/path/to/data/barcodes}"
OUTPUT_DIR="${OUTPUT_DIR:-/path/to/results/chimera}"

# --- Subsampling -------------------------------------------------------------
MAX_READS="${MAX_READS:-50000}"        # techo de reads por muestra
SEED="${SEED:-100}"                    # semilla fija -> reproducible
OVERSAMPLE="${OVERSAMPLE:-1.20}"       # margen del muestreo Bernoulli antes de recortar

# --- VSEARCH: clustering -----------------------------------------------------
# 0.97 es razonable para R10.4.1 con basecalling sup (~99% de precision).
# Para quimica antigua/ruidosa (R9.4.1, modelo fast) baja a 0.95 o 0.93.
CLUSTER_ID="${CLUSTER_ID:-0.97}"

# --- VSEARCH: UCHIME ---------------------------------------------------------
ABSKEW="${ABSKEW:-2.0}"                # asimetria minima de abundancia padre/quimera
MINH="${MINH:-0.28}"                   # score minimo (valor por defecto de vsearch)
CHIMERA_MODE="${CHIMERA_MODE:-denovo}" # denovo | ref | both
REF_DB="${REF_DB:-}"                   # FASTA de referencia (solo para ref/both)
ORIENT_DB="${ORIENT_DB:-}"             # FASTA para orientar reads (opcional, recomendado)
REMOVE_BORDERLINE="${REMOVE_BORDERLINE:-0}"  # 1 = eliminar tambien las "borderline" (?)

# --- Comportamiento ----------------------------------------------------------
THREADS="${THREADS:-4}"
FASTQ_QMAX="${FASTQ_QMAX:-93}"         # Nanopore supera el default 41 de vsearch
GZIP_OUTPUT="${GZIP_OUTPUT:-1}"        # 1 = salidas .fastq.gz (EMU/minimap2 las leen)
OVERWRITE="${OVERWRITE:-0}"            # 0 = no sobreescribir muestras ya procesadas
STOP_ON_ERROR="${STOP_ON_ERROR:-1}"    # 1 = abortar todo si una muestra falla
KEEP_TMP="${KEEP_TMP:-1}"              # 1 = conservar intermedios (auditoria)

# --- Estructura de salida ----------------------------------------------------
DIR_SUB="$OUTPUT_DIR/subsampled"
DIR_VS="$OUTPUT_DIR/vsearch"
DIR_FIN="$OUTPUT_DIR/final"
DIR_RES="$OUTPUT_DIR/resumen"
DIR_LOG="$OUTPUT_DIR/logs"
DIR_TMP="$OUTPUT_DIR/tmp"

CSV="$DIR_RES/resumen_quimeras.csv"
MANIFEST="$DIR_RES/manifiesto_entrada.tsv"
VERSIONS="$DIR_RES/versiones.txt"
PARAMS="$DIR_RES/parametros.txt"

EXT=".fastq"
[ "$GZIP_OUTPUT" = "1" ] && EXT=".fastq.gz"

# =============================================================================
#  2. UTILIDADES
# =============================================================================

C_OK=$'\033[32m'; C_WARN=$'\033[33m'; C_ERR=$'\033[31m'; C_INF=$'\033[36m'; C_0=$'\033[0m'
if [ ! -t 1 ]; then C_OK=""; C_WARN=""; C_ERR=""; C_INF=""; C_0=""; fi

log()  { printf '%s[%s]%s %s\n' "$C_INF" "$(date '+%H:%M:%S')" "$C_0" "$*"; }
ok()   { printf '%s[%s] OK%s %s\n' "$C_OK" "$(date '+%H:%M:%S')" "$C_0" "$*"; }
warn() { printf '%s[%s] AVISO%s %s\n' "$C_WARN" "$(date '+%H:%M:%S')" "$C_0" "$*" >&2; }
die()  { printf '%s[%s] ERROR%s %s\n' "$C_ERR" "$(date '+%H:%M:%S')" "$C_0" "$*" >&2; exit 1; }

trap 'die "fallo inesperado en la linea $LINENO (comando: $BASH_COMMAND)"' ERR

nlines() { local n=0; if [ -s "$1" ]; then n=$(wc -l < "$1"); fi; echo "$n"; }

count_reads() {
    local f="$1" n
    n=$(seqkit stats -T -j "$THREADS" "$f" 2>/dev/null | awk 'NR==2 {print $4}')
    if [ -z "$n" ]; then n=0; fi
    echo "$n"
}

vsearch_supports() {
    local opt="$1" out
    out="$(vsearch --help 2>&1 || true)"
    case "$out" in
        *"$opt"*) return 0 ;;
        *)        return 1 ;;
    esac
}

primera_linea() {
    local out
    out="$("$@" 2>&1 || true)"
    printf '%s' "${out%%$'\n'*}"
}

fasta_ids() {
    if [ -s "$1" ]; then
        grep '^>' "$1" | sed -e 's/^>//' -e 's/[[:space:]].*$//' -e 's/;size=[0-9]*;*//g' | sort -u
    fi
}

# =============================================================================
#  3. COMPROBACION DE DEPENDENCIAS Y VERSIONES
# =============================================================================

log "Comprobando dependencias..."
for prog in seqkit vsearch awk sort join comm grep sed find gzip; do
    command -v "$prog" >/dev/null 2>&1 || die "'$prog' no esta instalado o no esta en el PATH."
done

SEQKIT_V="$(primera_linea seqkit version)"
VSEARCH_V="$(primera_linea vsearch --version)"

vsearch_supports '--fastx_filter' || die "Tu vsearch no soporta --fastx_filter. Actualiza a >= 2.8.0."

HAS_ORIENT=0
if vsearch_supports '--orient'; then HAS_ORIENT=1; fi

[ -d "$INPUT_DIR" ] || die "INPUT_DIR no existe: $INPUT_DIR"

case "$CHIMERA_MODE" in
    denovo) ;;
    ref|both)
        [ -n "$REF_DB" ] || die "CHIMERA_MODE=$CHIMERA_MODE requiere REF_DB=<fasta de referencia>."
        [ -s "$REF_DB" ] || die "REF_DB no existe o esta vacio: $REF_DB" ;;
    *) die "CHIMERA_MODE debe ser: denovo | ref | both (recibido: '$CHIMERA_MODE')" ;;
esac

STRAND="both"
if [ -n "$ORIENT_DB" ]; then
    [ -s "$ORIENT_DB" ] || die "ORIENT_DB no existe o esta vacio: $ORIENT_DB"
    if [ "$HAS_ORIENT" = "1" ]; then
        STRAND="plus"
    else
        warn "Tu vsearch no soporta --orient (requiere >= 2.22.0). Se ignora ORIENT_DB."
        ORIENT_DB=""
    fi
fi

mkdir -p "$DIR_SUB" "$DIR_VS" "$DIR_FIN" "$DIR_RES" "$DIR_LOG" "$DIR_TMP"

{
    echo "Fecha de ejecucion : $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "Host               : $(uname -a)"
    echo "seqkit             : $SEQKIT_V"
    echo "vsearch            : $VSEARCH_V"
    if [ "$HAS_ORIENT" = "1" ]; then
        echo "vsearch --orient   : disponible"
    else
        echo "vsearch --orient   : NO disponible"
    fi
} > "$VERSIONS"

{
    echo "INPUT_DIR         = $INPUT_DIR"
    echo "OUTPUT_DIR        = $OUTPUT_DIR"
    echo "MAX_READS         = $MAX_READS"
    echo "SEED              = $SEED"
    echo "OVERSAMPLE        = $OVERSAMPLE"
    echo "CLUSTER_ID        = $CLUSTER_ID"
    echo "ABSKEW            = $ABSKEW"
    echo "MINH              = $MINH"
    echo "CHIMERA_MODE      = $CHIMERA_MODE"
    echo "REF_DB            = ${REF_DB:-<ninguna>}"
    echo "ORIENT_DB         = ${ORIENT_DB:-<ninguna>}"
    echo "STRAND (cluster)  = $STRAND"
    echo "REMOVE_BORDERLINE = $REMOVE_BORDERLINE"
    echo "FASTQ_QMAX        = $FASTQ_QMAX"
    echo "THREADS           = $THREADS"
} > "$PARAMS"

echo
echo "==============================================================="
echo " Pipeline de quimeras - Nanopore 16S full-length"
echo "==============================================================="
cat "$VERSIONS"
echo "---------------------------------------------------------------"
cat "$PARAMS"
echo "==============================================================="
echo

case "$OUTPUT_DIR" in
    *OneDrive*)
        warn "OUTPUT_DIR esta dentro de OneDrive: la sincronizacion puede ralentizar mucho la E/S. Considera pausar la sincronizacion durante la ejecucion." ;;
esac

# =============================================================================
#  4. DESCUBRIMIENTO DE MUESTRAS
#     a) INPUT_DIR/barcodeXX/*.fastq.gz   (un subdirectorio por barcode)
#     b) INPUT_DIR/barcodeXX.fastq.gz     (un fichero por barcode)
# =============================================================================

log "Descubriendo muestras en: $INPUT_DIR"
: > "$MANIFEST"

while IFS= read -r d; do
    [ -n "$d" ] || continue
    files="$(find "$d" -type f \( -name '*.fastq' -o -name '*.fq' -o -name '*.fastq.gz' -o -name '*.fq.gz' \) | sort | tr '\n' '|')"
    [ -n "$files" ] || continue
    printf '%s\t%s\n' "$(basename "$d")" "${files%|}" >> "$MANIFEST"
done < <(find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

while IFS= read -r f; do
    [ -n "$f" ] || continue
    n="$(basename "$f")"; n="${n%.gz}"; n="${n%.fastq}"; n="${n%.fq}"
    printf '%s\t%s\n' "$n" "$f" >> "$MANIFEST"
done < <(find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type f \( -name '*.fastq' -o -name '*.fq' -o -name '*.fastq.gz' -o -name '*.fq.gz' \) | sort)

N_SAMPLES=$(nlines "$MANIFEST")
[ "$N_SAMPLES" -gt 0 ] || die "No se encontro ningun FASTQ (.fastq/.fq/.fastq.gz/.fq.gz) en $INPUT_DIR"

DUPS="$(cut -f1 "$MANIFEST" | sort | uniq -d)"
[ -z "$DUPS" ] || die "Nombres de muestra duplicados (se mezclarian muestras): $DUPS"

ok "$N_SAMPLES muestra(s) detectada(s):"
awk -F'\t' '{n=split($2,a,"|"); printf "      - %-24s (%d fichero(s))\n", $1, n}' "$MANIFEST"
echo

# =============================================================================
#  5. CABECERA DEL CSV RESUMEN
# =============================================================================

if [ ! -s "$CSV" ]; then
    echo 'sample,reads_originales,reads_subsampleados,reads_evaluados,secuencias_unicas,clusters_total,clusters_singleton,clusters_quimericos,clusters_borderline,reads_quimericos,reads_borderline,porcentaje_quimeras,reads_finales,estado,advertencias' > "$CSV"
fi

# =============================================================================
#  6. PROCESADO DE UNA MUESTRA
# =============================================================================

procesar_muestra() {
    local sample="$1" filelist="$2"
    local vsd="$DIR_VS/$sample"
    local tmp="$DIR_TMP/$sample"
    local sub="$DIR_SUB/${sample}.subsampled${EXT}"
    local fin="$DIR_FIN/${sample}.sin_quimeras${EXT}"
    local advert=""

    mkdir -p "$vsd" "$tmp"

    local src="$tmp/entrada.fastq.gz"
    local -a files=()
    IFS='|' read -r -a files <<< "$filelist"

    if [ "${#files[@]}" -eq 1 ]; then
        src="${files[0]}"
    else
        log "  [$sample] uniendo ${#files[@]} ficheros de entrada..."
        : > "$src"
        local f
        for f in "${files[@]}"; do
            case "$f" in
                *.gz) cat "$f" >> "$src" ;;
                *)    gzip -c "$f" >> "$src" ;;
            esac
        done
    fi

    local n_orig
    n_orig=$(count_reads "$src")
    [ "$n_orig" -gt 0 ] || die "[$sample] el FASTQ de entrada no contiene lecturas."
    log "  [$sample] reads originales: $n_orig"

    if [ "$n_orig" -le "$MAX_READS" ]; then
        log "  [$sample] <= $MAX_READS reads: se conservan todas las lecturas (sin submuestreo)."
        seqkit seq "$src" -o "$sub"
    else
        local prop
        prop=$(awk -v m="$MAX_READS" -v t="$n_orig" -v f="$OVERSAMPLE" \
               'BEGIN{p=m*f/t; if(p>1)p=1; if(p<1e-6)p=1e-6; printf "%.8f", p}')
        log "  [$sample] submuestreando a $MAX_READS reads (semilla=$SEED, p=$prop)..."
        seqkit sample --rand-seed "$SEED" --proportion "$prop" "$src" -o "$tmp/sobremuestra.fastq.gz"
        seqkit shuffle --rand-seed "$SEED" "$tmp/sobremuestra.fastq.gz" -o "$tmp/barajado.fastq.gz"
        seqkit head -n "$MAX_READS" "$tmp/barajado.fastq.gz" -o "$sub"
        rm -f "$tmp/sobremuestra.fastq.gz" "$tmp/barajado.fastq.gz"
    fi

    local n_sub
    n_sub=$(count_reads "$sub")
    log "  [$sample] reads subsampleados: $n_sub"

    if [ "$n_orig" -gt "$MAX_READS" ] && [ "$n_sub" -lt "$MAX_READS" ]; then
        advert="${advert}submuestreo por debajo del objetivo (sube OVERSAMPLE); "
        warn "  [$sample] se obtuvieron $n_sub < $MAX_READS reads; sube OVERSAMPLE."
    fi

    local fa="$vsd/${sample}.reads.fasta"
    vsearch --fastx_filter "$sub" \
            --fastq_qmax "$FASTQ_QMAX" \
            --fasta_width 0 \
            --fastaout "$fa" --quiet

    local fa_use="$fa"
    if [ -n "$ORIENT_DB" ]; then
        log "  [$sample] orientando lecturas (solo en el FASTA temporal)..."
        vsearch --orient "$fa" --db "$ORIENT_DB" \
                --fastaout "$vsd/${sample}.oriented.fasta" \
                --notmatched "$vsd/${sample}.orient_notmatched.fasta" \
                --tabbedout "$vsd/${sample}.orient.tsv" \
                --fasta_width 0 --threads "$THREADS" --quiet
        cat "$vsd/${sample}.oriented.fasta" "$vsd/${sample}.orient_notmatched.fasta" \
            > "$vsd/${sample}.forvsearch.fasta"
        fa_use="$vsd/${sample}.forvsearch.fasta"
    fi

    log "  [$sample] dereplicando (contabilidad de duplicados exactos)..."
    vsearch --derep_fulllength "$fa_use" \
            --output "$vsd/${sample}.derep.fasta" \
            --uc "$vsd/${sample}.derep.uc" \
            --sizeout --minseqlength 1 --fasta_width 0 --strand plus --quiet

    local n_uniq
    n_uniq=$(grep -c '^>' "$vsd/${sample}.derep.fasta") || n_uniq=0

    log "  [$sample] clustering al $CLUSTER_ID de identidad (strand=$STRAND)..."
    vsearch --cluster_size "$vsd/${sample}.derep.fasta" \
            --id "$CLUSTER_ID" --strand "$STRAND" \
            --sizein --sizeout \
            --centroids "$vsd/${sample}.centroids.fasta" \
            --uc "$vsd/${sample}.clusters.uc" \
            --minseqlength 1 --fasta_width 0 --threads "$THREADS" --quiet

    vsearch --sortbysize "$vsd/${sample}.centroids.fasta" \
            --output "$vsd/${sample}.centroids.sorted.fasta" \
            --minseqlength 1 --fasta_width 0 --quiet

    awk -F'\t' 'BEGIN{OFS="\t"}
        $1=="S"{print $9,$9}
        $1=="H"{print $10,$9}' "$vsd/${sample}.derep.uc" \
      | sed 's/;size=[0-9]*;*//g' | sort -k1,1 > "$tmp/rep2read.tsv"

    awk -F'\t' 'BEGIN{OFS="\t"}
        $1=="S"{print $9,$9}
        $1=="H"{print $9,$10}' "$vsd/${sample}.clusters.uc" \
      | sed 's/;size=[0-9]*;*//g' | sort -k1,1 > "$tmp/rep2centroid.tsv"

    join -t "$(printf '\t')" -1 1 -2 1 -o 1.2,2.2 \
         "$tmp/rep2centroid.tsv" "$tmp/rep2read.tsv" \
      | sort -k1,1 > "$vsd/${sample}.centroide2read.tsv"

    local n_eval
    n_eval=$(nlines "$vsd/${sample}.centroide2read.tsv")

    cut -f1 "$vsd/${sample}.centroide2read.tsv" | uniq -c \
        | awk '{print $2"\t"$1}' > "$vsd/${sample}.cluster_sizes.tsv"
    local n_clust n_single
    n_clust=$(nlines "$vsd/${sample}.cluster_sizes.tsv")
    n_single=$(awk -F'\t' '$2==1' "$vsd/${sample}.cluster_sizes.tsv" | wc -l)

    log "  [$sample] reads evaluados=$n_eval | secuencias unicas=$n_uniq | clusters=$n_clust (singletons=$n_single)"

    if [ "$n_eval" -gt 0 ]; then
        local pct_reads_single
        pct_reads_single=$(awk -v s="$n_single" -v e="$n_eval" 'BEGIN{printf "%.1f", 100*s/e}')
        log "  [$sample] lecturas en clusters singleton: ${pct_reads_single}%"
        if awk -v p="$pct_reads_single" 'BEGIN{exit !(p>50)}'; then
            advert="${advert}${pct_reads_single}% de las LECTURAS estan en clusters de 1 read: uchime_denovo pierde poder; "
            warn "  [$sample] ${pct_reads_single}% de las lecturas caen en clusters de un solo read. Considera bajar CLUSTER_ID (0.95 / 0.93) o usar CHIMERA_MODE=ref."
        fi
    fi

    local chim_reads="$vsd/${sample}.reads_quimericos.txt"
    local bord_reads="$vsd/${sample}.reads_borderline.txt"
    : > "$chim_reads"; : > "$bord_reads"
    local n_cl_chim=0 n_cl_bord=0

    if [ "$CHIMERA_MODE" = "denovo" ] || [ "$CHIMERA_MODE" = "both" ]; then
        log "  [$sample] ejecutando vsearch --uchime_denovo (abskew=$ABSKEW, minh=$MINH)..."
        vsearch --uchime_denovo "$vsd/${sample}.centroids.sorted.fasta" \
                --abskew "$ABSKEW" --minh "$MINH" \
                --chimeras    "$vsd/${sample}.denovo_chimeras.fasta" \
                --nonchimeras "$vsd/${sample}.denovo_nonchimeras.fasta" \
                --borderline  "$vsd/${sample}.denovo_borderline.fasta" \
                --uchimeout   "$vsd/${sample}.uchime_denovo.tsv" \
                --fasta_width 0 --quiet

        fasta_ids "$vsd/${sample}.denovo_chimeras.fasta" > "$tmp/chim_centroids.txt"
        n_cl_chim=$(nlines "$tmp/chim_centroids.txt")
        if [ "$n_cl_chim" -gt 0 ]; then
            join -t "$(printf '\t')" -1 1 -2 1 -o 2.2 \
                 "$tmp/chim_centroids.txt" "$vsd/${sample}.centroide2read.tsv" \
              >> "$chim_reads"
        fi

        fasta_ids "$vsd/${sample}.denovo_borderline.fasta" > "$tmp/bord_centroids.txt"
        n_cl_bord=$(nlines "$tmp/bord_centroids.txt")
        if [ "$n_cl_bord" -gt 0 ]; then
            join -t "$(printf '\t')" -1 1 -2 1 -o 2.2 \
                 "$tmp/bord_centroids.txt" "$vsd/${sample}.centroide2read.tsv" \
              >> "$bord_reads"
        fi
    fi

    if [ "$CHIMERA_MODE" = "ref" ] || [ "$CHIMERA_MODE" = "both" ]; then
        log "  [$sample] ejecutando vsearch --uchime_ref contra $REF_DB (puede tardar)..."
        vsearch --uchime_ref "$fa_use" --db "$REF_DB" \
                --chimeras    "$vsd/${sample}.ref_chimeras.fasta" \
                --nonchimeras "$vsd/${sample}.ref_nonchimeras.fasta" \
                --borderline  "$vsd/${sample}.ref_borderline.fasta" \
                --uchimeout   "$vsd/${sample}.uchime_ref.tsv" \
                --fasta_width 0 --threads "$THREADS" --quiet

        fasta_ids "$vsd/${sample}.ref_chimeras.fasta"   >> "$chim_reads"
        fasta_ids "$vsd/${sample}.ref_borderline.fasta" >> "$bord_reads"
    fi

    sort -u "$chim_reads" -o "$chim_reads"
    sort -u "$bord_reads" -o "$bord_reads"
    comm -23 "$bord_reads" "$chim_reads" > "$tmp/bord_only.txt"
    mv "$tmp/bord_only.txt" "$bord_reads"

    local n_bord
    n_bord=$(nlines "$bord_reads")

    if [ "$REMOVE_BORDERLINE" = "1" ] && [ "$n_bord" -gt 0 ]; then
        cat "$bord_reads" >> "$chim_reads"
        sort -u "$chim_reads" -o "$chim_reads"
        advert="${advert}borderline eliminadas; "
    fi

    local n_chim
    n_chim=$(nlines "$chim_reads")

    if [ "$n_chim" -gt 0 ]; then
        log "  [$sample] eliminando $n_chim lecturas quimericas del FASTQ..."
        seqkit grep -v -f "$chim_reads" "$sub" -o "$fin" 2>/dev/null
    else
        log "  [$sample] no se detectaron quimeras: el FASTQ final es identico al subsampleado."
        cp "$sub" "$fin"
    fi

    local n_fin
    n_fin=$(count_reads "$fin")

    log "  [$sample] validando..."

    [ "$n_orig" -ge "$n_sub" ] || die "[$sample] reads_originales ($n_orig) < reads_subsampleados ($n_sub)"
    [ "$n_sub"  -ge "$n_fin" ] || die "[$sample] reads_subsampleados ($n_sub) < reads_finales ($n_fin)"

    local esperado=$(( n_sub - n_fin ))
    [ "$esperado" -eq "$n_chim" ] \
        || die "[$sample] descuadre: reads_sub - reads_fin = $esperado pero reads_quimericos = $n_chim"

    local fmt avgq
    seqkit stats -T -a -j "$THREADS" "$fin" > "$tmp/stats_final.tsv"
    fmt=$(awk -F'\t' 'NR==2{print $2}' "$tmp/stats_final.tsv")
    [ "$fmt" = "FASTQ" ] || die "[$sample] el fichero final no es FASTQ (formato detectado: '$fmt')"
    if [ "$n_fin" -gt 0 ]; then
        avgq=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="AvgQual") c=i}
                           NR==2{ if(c) print $c; else print "NA" }' "$tmp/stats_final.tsv")
        if [ "$avgq" = "NA" ]; then
            warn "  [$sample] esta version de seqkit no reporta AvgQual; se omite esa comprobacion."
        else
            awk -v q="$avgq" 'BEGIN{exit !(q+0>0)}' \
                || die "[$sample] calidad media = $avgq: se han perdido los quality scores."
        fi
    fi

    seqkit fx2tab -i "$sub" 2>/dev/null | cut -f1-3 | sort > "$tmp/sub.tsv"
    seqkit fx2tab -i "$fin" 2>/dev/null | cut -f1-3 | sort > "$tmp/fin.tsv"
    comm -23 "$tmp/fin.tsv" "$tmp/sub.tsv" > "$tmp/no_encontrados.tsv"
    [ ! -s "$tmp/no_encontrados.tsv" ] \
        || die "[$sample] $(nlines "$tmp/no_encontrados.tsv") registros del FASTQ final NO coinciden con el subsampleado."

    if [ "$n_chim" -gt 0 ]; then
        cut -f1 "$tmp/fin.tsv" | sort -u > "$tmp/fin_ids.txt"
        comm -12 "$tmp/fin_ids.txt" "$chim_reads" > "$tmp/supervivientes.txt"
        [ ! -s "$tmp/supervivientes.txt" ] \
            || die "[$sample] $(nlines "$tmp/supervivientes.txt") lecturas quimericas siguen presentes en el FASTQ final."
    fi

    ok "  [$sample] todas las validaciones superadas."

    local pct estado
    pct=$(awk -v c="$n_chim" -v s="$n_sub" 'BEGIN{ if(s>0) printf "%.4f", 100*c/s; else printf "NA" }')
    if [ "$n_chim" -gt 0 ]; then estado="OK"; else estado="OK_SIN_QUIMERAS"; fi
    if [ -z "$advert" ]; then advert="-"; fi

    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,"%s"\n' \
        "$sample" "$n_orig" "$n_sub" "$n_eval" "$n_uniq" "$n_clust" "$n_single" \
        "$n_cl_chim" "$n_cl_bord" "$n_chim" "$n_bord" "$pct" "$n_fin" "$estado" "$advert" \
        >> "$CSV"

    ok "  [$sample] $n_orig -> $n_sub -> quimeras=$n_chim (${pct}%) -> final=$n_fin"

    if [ "$KEEP_TMP" != "1" ]; then rm -rf "$tmp"; fi
    return 0
}

# =============================================================================
#  7. BUCLE PRINCIPAL
# =============================================================================

N_OK=0; N_SKIP=0; N_FAIL=0; IDX=0

while IFS=$'\t' read -r sample filelist; do
    [ -n "$sample" ] || continue
    IDX=$(( IDX + 1 ))
    echo "---------------------------------------------------------------"
    log "[$IDX/$N_SAMPLES] MUESTRA: $sample"

    fin_check="$DIR_FIN/${sample}.sin_quimeras${EXT}"
    if [ -s "$fin_check" ] && [ "$OVERWRITE" != "1" ]; then
        warn "  [$sample] ya existe $(basename "$fin_check"): se omite (usa OVERWRITE=1 para rehacerla)."
        N_SKIP=$(( N_SKIP + 1 ))
        continue
    fi

    trap - ERR
    set +e
    ( set -Eeuo pipefail
      trap 'die "fallo inesperado en la linea $LINENO (comando: $BASH_COMMAND)"' ERR
      procesar_muestra "$sample" "$filelist" ) 2>&1 | tee "$DIR_LOG/${sample}.log"
    rc=${PIPESTATUS[0]}
    set -e
    trap 'die "fallo inesperado en la linea $LINENO (comando: $BASH_COMMAND)"' ERR

    if [ "$rc" -eq 0 ]; then
        N_OK=$(( N_OK + 1 ))
    else
        N_FAIL=$(( N_FAIL + 1 ))
        printf '%s,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,ERROR,"ver logs/%s.log"\n' "$sample" "$sample" >> "$CSV"
        warn "  [$sample] FALLO (detalles en $DIR_LOG/${sample}.log)"
        if [ "$STOP_ON_ERROR" = "1" ]; then
            die "Se detiene el pipeline por error en la muestra '$sample' (usa STOP_ON_ERROR=0 para continuar con las demas)."
        fi
    fi
done < "$MANIFEST"

# =============================================================================
#  8. CIERRE
# =============================================================================

echo
echo "==============================================================="
ok "Procesadas: $N_OK | Omitidas: $N_SKIP | Con error: $N_FAIL"
echo
echo "Resumen : $CSV"
echo "Finales : $DIR_FIN/*.sin_quimeras${EXT}   <-- ENTRADA PARA EMU"
echo "==============================================================="
echo
column -s, -t < "$CSV" 2>/dev/null || cat "$CSV"
echo

[ "$N_FAIL" -eq 0 ] || exit 1
exit 0
