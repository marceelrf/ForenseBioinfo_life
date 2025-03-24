#!/bin/bash

# Diretórios (substitua pelos caminhos corretos)
BAM_DIR="/caminho/para/bams"
OUTPUT_DIR="/caminho/para/MethylationEntropy"
REGIONS="/caminho/para/BEDs/Genes.bed"
REF="/caminho/para/genoma/genoma_referencia.fa"
LOG_FILE="${OUTPUT_DIR}/modkit_entropy.log"

# Número de threads
THREADS=12

# Criar diretório de saída se não existir
mkdir -p "$OUTPUT_DIR"

# Percorrer todos os arquivos BAM no diretório
for BAM_FILE in "$BAM_DIR"/*.bam; do
    # Extrair o nome da amostra (removendo diretório e extensão)
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)

    # Executar modkit entropy
    echo "Calculando entropia de metilação para $SAMPLE_NAME..."
    modkit entropy --in-bam "$BAM_FILE" -o "$OUTPUT_DIR" --prefix "$SAMPLE_NAME" --regions "$REGIONS" --cpg --ref "$REF" --threads "$THREADS" --log-filepath "$LOG_FILE"
done

echo "Processamento concluído!"
