#!/bin/bash

# Diretório onde os arquivos SAM estão localizados
SAM_DIR="/media/lab/Data1/longReads/R10_1KG/skin/assembly/retro_BAMs/"  # Altere para o diretório correto onde os arquivos SAM estão

# Diretório onde os arquivos BAM serão salvos
OUTPUT_DIR="/media/lab/Data1/longReads/R10_1KG/skin/assembly/retro_BAMs/"  # Altere para o diretório onde você quer salvar os arquivos BAM

# Verifica se o diretório de saída existe, caso contrário, cria
mkdir -p "$OUTPUT_DIR"

# Loop através de todos os arquivos .sam no diretório especificado
for sam_file in "$SAM_DIR"/*.sam; do
    # Extrai o nome base do arquivo (sem o caminho e extensão)
    base_name=$(basename "$sam_file" .sam)
    
    # Define o nome do arquivo BAM de saída
    bam_file="$OUTPUT_DIR/$base_name.bam"
    
    # Aplica o comando samtools sort para ordenar e converter para BAM
    echo "Ordenando e convertendo $sam_file para $bam_file..."
    samtools sort -o "$bam_file" "$sam_file"
    
    echo "Arquivo $bam_file criado com sucesso."
done

echo "Processo completo."
