#!/bin/bash

# 원본 pQTL 데이터가 있는 폴더
PQTL_DIR="/BiO/hae/000005_MRSPI/12_pqtl_sig_clump_each_GWAS"

# Clumping 결과를 저장할 폴더
OUTPUT_BASE_DIR="/BiO/hae/000005_MRSPI/13_pqtl_clump"

# PLINK 실행 파일 경로
PLINK_PATH="plink"

# 1000G Reference 데이터
REF_DATA="/BiO/hae/000006_ref_1000G/ref"

# GWAS 폴더 리스트 가져오기
for GWAS_FOLDER in "$PQTL_DIR"/*; do
    # GWAS ID 추출 (예: GCST002245)
    GWAS_ID=$(basename "$GWAS_FOLDER")

    # Clumping 결과 저장할 폴더 생성
    OUTPUT_DIR="$OUTPUT_BASE_DIR/$GWAS_ID"
    mkdir -p "$OUTPUT_DIR"

    echo "🔹 Processing GWAS: $GWAS_ID"

    # GWAS 폴더 내 pQTL 파일들 찾기
    for PQTL_FILE in "$GWAS_FOLDER"/*_v2.txt; do
        # 파일명에서 확장자 제거
        FILE_NAME=$(basename "$PQTL_FILE" .txt)

        # Clumping을 위한 input 파일 생성
        CLUMP_INPUT="$OUTPUT_DIR/${FILE_NAME}_clump_input.txt"
        awk 'NR==1 || !seen[$2]++' "$PQTL_FILE" > "$CLUMP_INPUT"

        # PLINK Clumping 실행
        PLINK_OUTPUT="$OUTPUT_DIR/${FILE_NAME}_clumped"
        "$PLINK_PATH" --bfile "$REF_DATA" \
                      --clump "$CLUMP_INPUT" \
                      --clump-kb 1000 \
                      --clump-r2 0.01 \
                      --clump-p1 1 \
                      --clump-p2 1 \
                      --out "$PLINK_OUTPUT"

        # Clumping 결과 확인 후 저장
        if [ -f "${PLINK_OUTPUT}.clumped" ]; then
            mv "${PLINK_OUTPUT}.clumped" "${PLINK_OUTPUT}_clumped.tsv"
            echo "✅ Clumping 완료: $PLINK_OUTPUT"
        else
            echo "⚠️ Clumping 결과 없음: $PQTL_FILE"
        fi
    done
done
