#!/bin/bash

# Diretórios (substitua pelos caminhos corretos)
HAP_DIR="/caminho/para/haplotypebased"
OUTPUT_DIR="/caminho/para/DMR_hap"
REF="/caminho/para/genoma/genoma_referencia.fa"

# Número de threads
THREADS=12

# Criar diretório de saída se não existir
mkdir -p "$OUTPUT_DIR"

# Percorrer todos os arquivos *_haplotyped_1.bed.gz no diretório
for HAP1 in "$HAP_DIR"/*_haplotyped_1.bed.gz; do
    # Extrair o nome da amostra (removendo diretório e sufixo)
    SAMPLE_NAME=$(basename "$HAP1" _haplotyped_1.bed.gz)

    # Definir o arquivo do segundo haplótipo
    HAP2="${HAP_DIR}/${SAMPLE_NAME}_haplotyped_2.bed.gz"

    # Definir o caminho dos arquivos de saída
    OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.bed"
    LOG_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_dmr.log"

    # Verificar se ambos os arquivos de haplótipos existem
    if [[ -f "$HAP1" && -f "$HAP2" ]]; then
        echo "Processando DMR para $SAMPLE_NAME..."
        modkit dmr pair -a "$HAP1" -b "$HAP2" -o "$OUTPUT_FILE" --ref "$REF" --base C --threads "$THREADS" --log-filepath "$LOG_FILE"
    else
        echo "Aviso: Arquivo faltando para $SAMPLE_NAME, pulando..."
    fi
done

echo "Processamento concluído!"
