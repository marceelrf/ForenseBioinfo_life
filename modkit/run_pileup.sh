#!/bin/bash

# Diretórios (substitua pelos caminhos corretos)
BAM_DIR="/caminho/para/bams"
BED_DIR="/caminho/para/methylbeds"
REF="/caminho/para/genoma/genoma_referencia.fa"

# Número de threads
THREADS=12

# Criar diretório de saída se não existir
mkdir -p "$BED_DIR"

# Percorrer todos os arquivos BAM no diretório
for BAM_FILE in "$BAM_DIR"/*.bam; do
    # Extrair o nome da amostra (removendo diretório e extensão)
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)

    # Definir o caminho do arquivo BED de saída
    BED_FILE="$BED_DIR/${SAMPLE_NAME}.bed"

    # Executar modkit pileup
    echo "Processando $SAMPLE_NAME..."
    modkit pileup "$BAM_FILE" "$BED_FILE" --cpg --ref "$REF" --combine-strands -t "$THREADS"
done

echo "Processamento concluído!"
