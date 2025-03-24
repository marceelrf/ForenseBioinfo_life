#!/bin/bash

# Diretórios (substitua pelos caminhos corretos)
BAM_DIR="/caminho/para/bams"
OUTPUT_DIR="/caminho/para/haplotypebased"
REF="/caminho/para/genoma/genoma_referencia.fa"

# Número de threads
THREADS=12

# Criar diretório de saída se não existir
mkdir -p "$OUTPUT_DIR"

# Percorrer todos os arquivos BAM no diretório
for BAM_FILE in "$BAM_DIR"/*.bam; do
    # Extrair o nome da amostra (removendo diretório e extensão)
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)

    # Definir prefixo para saída
    PREFIX="${SAMPLE_NAME}_haplotyped"

    # Executar modkit pileup
    echo "Processando $SAMPLE_NAME..."
    modkit pileup "$BAM_FILE" "$OUTPUT_DIR" --cpg --ref "$REF" --partition-tag HP --prefix "$PREFIX" -t "$THREADS"
done

echo "Processamento concluído!"
