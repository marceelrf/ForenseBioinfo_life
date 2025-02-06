#!/bin/bash

# Caminho para a pasta gVCF
gvcf_dir="path/to/clair3/gVCF"

# Percorrer todos os arquivos .gvcf.gz na pasta gVCF
for gvcf_file in "$gvcf_dir"/*.gvcf.gz; do
    # Pega o nome da amostra a partir do nome do arquivo
    # Assume que o nome da amostra está presente antes do sufixo ".gvcf.gz"
    amostra=$(basename "$gvcf_file" .gvcf.gz)

    # Verifica se o nome da amostra é diferente de "SAMPLE"
    if [[ "$amostra" != "SAMPLE" ]]; then
        # Substitui o nome da amostra no arquivo VCF usando bcftools
        bcftools reheader -s <(echo "$amostra") -o "${gvcf_file%.gvcf.gz}_renomeado.gvcf.gz" "$gvcf_file"

        # Opcional: Remove o arquivo original, se desejar manter apenas os renomeados
        # rm "$gvcf_file"

        echo "Arquivo renomeado: $gvcf_file para ${gvcf_file%.gvcf.gz}_renomeado.gvcf.gz"
    fi
done
