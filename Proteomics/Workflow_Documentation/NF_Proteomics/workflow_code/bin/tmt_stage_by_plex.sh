#!/usr/bin/env bash
# Reorganize TMT mzML files into plex folders and create annotation.txt per plex.
# FragPipe expects: plex_folder/mzML_files + plex_folder/annotation.txt
#
# FragPipe annotation notes (https://fragpipe.nesvilab.org/docs/tutorial_fragpipe.html):
#   - One annotation.txt per folder; unused channels = sample "NA"
#
# Usage: cd DIR_WITH_MZML_FILES && tmt_stage_by_plex.sh [MANIFEST] [EXP_ANNOT] [TMT_FORMAT]
#   Run from dir containing the mzML files. Creates TMTa_1/, TMTb_1/, etc. in cwd.
#   MANIFEST   - manifest.tsv (file, experiment, bioreplicate, data_type; no header). Default: manifest.tsv
#   EXP_ANNOT  - experiment_annotation.tsv. Default: experiment_annotation.tsv
#   TMT_FORMAT - TMT6 | TMT10 | TMT11 | TMT16 | TMT18. Required.
#
# Output: ./{plex}_{biorep}/*.mzML, ./{plex}_{biorep}/annotation.txt

set -euo pipefail

BASE="$(pwd)"
MANIFEST="${1:-$BASE/manifest.tsv}"
EXP_ANNOT="${2:-$BASE/experiment_annotation.tsv}"
TMT_FORMAT="${3:-}"

[[ -f "$MANIFEST" ]] || { echo "ERROR: Manifest not found: $MANIFEST" >&2; exit 1; }
[[ -f "$EXP_ANNOT" ]] || { echo "ERROR: Experiment annotation not found: $EXP_ANNOT" >&2; exit 1; }
[[ -n "$TMT_FORMAT" ]] || { echo "ERROR: TMT_FORMAT required (TMT6, TMT10, TMT11, TMT16, TMT18)" >&2; exit 1; }

echo "BASE=$BASE"
echo "MANIFEST=$MANIFEST"
echo "EXP_ANNOT=$EXP_ANNOT"
echo "TMT_FORMAT=$TMT_FORMAT"

# Create plex dirs and move files
# Manifest format: file, experiment (plex), bioreplicate, data_type
while IFS=$'\t' read -r fname experiment bioreplicate data_type; do
    [[ -z "$fname" ]] && continue
    plex="$experiment"
    # Include bioreplicate in folder name: TMTa_1, TMTb_1, etc.
    if [[ -n "$bioreplicate" ]]; then
        plex_dir="${plex}_${bioreplicate}"
    else
        plex_dir="$plex"
    fi
    mkdir -p "$BASE/$plex_dir"
    # Use basename for lookup: manifest Path should match staged files (sample_name.mzML) in work dir
    fname_base="$(basename "$fname")"
    src="$BASE/$fname_base"
    if [[ -f "$src" ]]; then
        mv "$src" "$BASE/$plex_dir/"
        echo "  mv $fname_base -> $plex_dir/"
    else
        echo "WARN: $src not found, skipping" >&2
    fi
done < "$MANIFEST"

# Generate annotation.txt per plex
# TMT channel order (must match FragPipe's expected order from QuantLabel.java)
case "$TMT_FORMAT" in
    TMT6)  CHANNEL_ORDER="126 127N 128C 129N 130C 131N" ;;
    TMT10) CHANNEL_ORDER="126 127N 127C 128N 128C 129N 129C 130N 130C 131N" ;;
    TMT11) CHANNEL_ORDER="126 127N 127C 128N 128C 129N 129C 130N 130C 131N 131C" ;;
    TMT16) CHANNEL_ORDER="126 127N 127C 128N 128C 129N 129C 130N 130C 131N 131C 132N 132C 133N 133C 134N" ;;
    TMT18) CHANNEL_ORDER="126 127N 127C 128N 128C 129N 129C 130N 130C 131N 131C 132N 132C 133N 133C 134N 134C 135N" ;;
    *) echo "ERROR: Unknown TMT_FORMAT: $TMT_FORMAT (use TMT6, TMT10, TMT11, TMT16, TMT18)" >&2; exit 1 ;;
esac

# Get unique plex_dir names (Experiment_Bioreplicate) from manifest
# Folders: TMTa_1, TMTa_1_1 (with TechRepMixture), TMTb_1, etc.
# Experiment annotation plex column must exactly match folder name
plex_dirs=$(awk -F'\t' '{if ($2 && $3) print $2 "_" $3; else if ($2) print $2}' "$MANIFEST" | sort -u)
for plex_dir in $plex_dirs; do
    annot="$BASE/$plex_dir/annotation.txt"
    
    # Extract channels for this folder: exact match on plex column (plex_dir = Experiment_Bioreplicate)
    tail -n +2 "$EXP_ANNOT" | awk -v plex_dir="$plex_dir" -F'\t' '$1 == plex_dir { print $2 "\t" $3 }' > "$annot.tmp"
    
    # Build channel map
    declare -A channel_map
    while IFS=$'\t' read -r channel sample; do
        [[ -n "$channel" ]] && channel_map["$channel"]="$sample"
    done < "$annot.tmp"
    
    # Write annotation.txt with ALL channels in order (FragPipe requires rows == label type channel count).
    # Unused channels get sample name "NA" per FragPipe docs: "To ignore certain channel... set sample name to NA"
    > "$annot"
    for channel in $CHANNEL_ORDER; do
        sample="${channel_map[$channel]:-NA}"
        echo -e "$channel\t$sample" >> "$annot"
    done
    
    rm -f "$annot.tmp"
    [[ -s "$annot" ]] || echo "WARN: No annotation for plex $plex_dir" >&2
done

# Update manifest with plex_dir paths so FragPipe finds files
# Manifest format: file, experiment (plex), bioreplicate, data_type
while IFS=$'\t' read -r fname experiment bioreplicate data_type; do
    [[ -z "$fname" ]] && continue
    fname_base="$(basename "$fname")"
    # Build plex_dir name (plex_biorep)
    if [[ -n "$bioreplicate" ]]; then
        plex_dir="${experiment}_${bioreplicate}"
    else
        plex_dir="$experiment"
    fi
    echo -e "${plex_dir}/${fname_base}\t${experiment}\t${bioreplicate}\t${data_type}"
done < "$MANIFEST" > "$MANIFEST.tmp" && mv "$MANIFEST.tmp" "$MANIFEST"

echo "Done. Plex dirs in $BASE"
