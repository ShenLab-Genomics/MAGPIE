#!/bin/bash
MODE="pred"
FILE_STATE="annotated"
VISUALIZATION=0
ANNOVAR_DIR="annovar/"
TEMP_DIR="data/temp/"
BPCA_DIR="BPCA/"
ANNOVAR_DATA_DIR="data/output/annovar/"
SPLICEAI_DATA_DIR="data/output/spliceai/"
TEST_FILE="data/datasets/test.csv"
TRAIN_FILE="data/datasets/denovo.csv"
MODEL_FILE="data/result/MAGPIE.model"
FEATURE_FILE="data/result/openFE.features"
SELECTION_FILE="data/result/selection.csv"

if [ $# -eq 0 ]; then
  echo "No arguments provided. Please use 'bash magpie.sh --help' to get more info."
  exit 1
else
  SHORT_OPTS="p:"
  LONG_OPTS="mode:,input_file:,train_file:,test_file:,model_file:,feature:,selection:,file_state:,help,visualization" ## 選項後有冒號代表包含參數

  opt=$(getopt -o $SHORT_OPTS --long $LONG_OPTS --name "$(basename "$0")" -- "$@")
  eval set --"$opt"

  while true; do
    case "$1" in
      --mode) MODE="$2"; shift 2;;
      --file_state) FILE_STATE="$2"; shift 2;;
      --input_file) TRAIN_FILE=$(readlink -f "$2"); shift 2;;
      --train_file) TRAIN_FILE=$(readlink -f "$2"); shift 2;;
      --test_file) TEST_FILE=$(readlink -f "$2"); shift 2;;
      --model_file) MODEL_FILE=$(readlink -f "$2"); shift 2;;
      --feature) FEATURE_FILE=$(readlink -f "$2"); shift 2;;
      --selection) SELECTION_FILE=$(readlink -f "$2"); shift 2;;
      --visualization) VISUALIZATION=1; shift;;
      --help) HELP=1; shift;;
      --) shift; break;;
      *) echo "Invalid option: $1"; exit 1;;
    esac
  done

  if [ -n "$HELP" ]; then
    cat << EOF
Usage: $(basename "$0") [--mode {running mode of MAGPIE}] | default: pred | pred & train supported
                 [--file_state {file state of input file.}] | default: annotated | required when mode is pred | annotated & unannotated supported
                 [--input_file {path of train file}] | required when mode is train
                 [--train_file {path of train file of trained model}] | required when mode is pred
                 [--test_file {path of test file}] | required when mode is pred
                 [--model_file {path of model file}] | required when mode is pred
                 [--feature_file {path of feature list file}] | required when mode is pred
                 [--selection_file {path of selection list file}] | required when mode is pred
                 [--visualization] | Visualize MAGPIE prediction results or not.
                 [--help] | display usage of MAGPIE script
EOF
    exit
  fi
fi

TRAIN_FILE_NAME=$(basename "$TRAIN_FILE")
TRAIN_FILE_NAME=${TRAIN_FILE_NAME%.*}
TEST_FILE_NAME=$(basename "$TEST_FILE")
TEST_FILE_NAME=${TEST_FILE_NAME%.*}

annotate_file() {
  local FULL_FILE_NAME="$1"
  local FILE_NAME="$2"
  python python/magpie.py --mode prepare --input_file "${FULL_FILE_NAME}"
  perl ${ANNOVAR_DIR}table_annovar.pl ${ANNOVAR_DATA_DIR}"${FILE_NAME}".avinput ${ANNOVAR_DATA_DIR}humandb/ -buildver hg38 -out ${ANNOVAR_DATA_DIR}"${FILE_NAME}" -remove -protocol refGene,phastConsElements100way,gnomad30_genome,dbnsfp33a,dbnsfp42a -operation g,r,f,f,f -csvout
  conda activate spliceai
  spliceai -I ${SPLICEAI_DATA_DIR}"${FILE_NAME}".vcf -O ${SPLICEAI_DATA_DIR}"${FILE_NAME}"_out.vcf -R ${SPLICEAI_DATA_DIR}hg38.fa -A grch38
  conda deactivate
  python python/magpie.py --mode merge --input_file "${ANNOVAR_DATA_DIR}${FILE_NAME}.hg38_multianno.csv" --spliceai_out "${SPLICEAI_DATA_DIR}/${FILE_NAME}_out.vcf"
  matlab -nodesktop -nosplash -r "filename=${TEMP_DIR}${FILE_NAME}.csv" ${BPCA_DIR}MAGPIE_fill.m
}

if [ "$MODE" = "pred" ]; then
  if [ "$FILE_STATE" = "annotated" ]; then
    if [ "$VISUALIZATION" = 1 ]; then
      python python/magpie.py --mode pred --test_file "$TEST_FILE" --model_file "$MODEL_FILE" --feature "$FEATURE_FILE" --selection "$SELECTION_FILE" --file_state "$FILE_STATE" --visualization
    else
      python python/magpie.py --mode pred --test_file "$TEST_FILE" --model_file "$MODEL_FILE" --feature "$FEATURE_FILE" --selection "$SELECTION_FILE" --file_state "$FILE_STATE"
    fi
  elif [ "$FILE_STATE" = "unannotated" ]; then
    annotate_file "$TEST_FILE" "$TEST_FILE_NAME"
    if [ "$VISUALIZATION" = 1 ]; then
      python python/magpie.py --mode pred --test_file "${TEMP_DIR}${TEST_FILE_NAME}.csv" --model_file "$MODEL_FILE" --feature "$FEATURE_FILE" --selection "$SELECTION_FILE" --file_state "$FILE_STATE" --visualization
    else
      python python/magpie.py --mode pred --test_file "${TEMP_DIR}${TEST_FILE_NAME}.csv" --model_file "$MODEL_FILE" --feature "$FEATURE_FILE" --selection "$SELECTION_FILE" --file_state "$FILE_STATE"
    fi
    rm "${TEMP_DIR}${TEST_FILE_NAME}.csv" "${TEMP_DIR}${TEST_FILE_NAME}_bpca.csv" 
  fi
elif [ "$MODE" = 'train' ]; then
  annotate_file "$TRAIN_FILE" "$TRAIN_FILE_NAME"
  python python/magpie.py --mode train --input_file "${TEMP_DIR}${TRAIN_FILE_NAME}.csv"
  rm "${TEMP_DIR}${TRAIN_FILE_NAME}.csv"
fi
