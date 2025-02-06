#!/bin/bash

# Caminho para a pasta gVCF
gvcf_dir="gVCF"

# Percorrer todas as subpastas de amostras (presumindo que são pastas de amostras dentro de clair3)
for amostra_dir in */; do
    # Verifica se é uma pasta de amostra (ignora a pasta gVCF)
    if [[ "$amostra_dir" != "$gvcf_dir/" ]]; then
        # Pega o nome da amostra a partir do nome da pasta
        amostra=$(basename "$amostra_dir")
        
        # Caminho completo dos arquivos merge_output.gvcf.gz e merge_output.gvcf.gz.tbi
        gvcf_file="$amostra_dir/merge_output.gvcf.gz"
        tbi_file="$amostra_dir/merge_output.gvcf.gz.tbi"
        
        # Verifica se os arquivos existem antes de movê-los
        if [[ -f "$gvcf_file" && -f "$tbi_file" ]]; then
            # Move e renomeia os arquivos para a pasta gVCF com o nome da amostra
            mv "$gvcf_file" "$gvcf_dir/$amostra.gvcf.gz"
            mv "$tbi_file" "$gvcf_dir/$amostra.gvcf.gz.tbi"
            echo "Arquivos de amostra $amostra movidos e renomeados com sucesso."
        else
            echo "Arquivos de amostra $amostra não encontrados."
        fi
    fi
done
