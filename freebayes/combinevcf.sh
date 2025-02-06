for file in /media/lab/Data1/longReads/R10_1KG/skin/Freebayes_genes/combined/Freebayes_genes/*.vcf; do
    bgzip -c "$file" > "$file.gz"
    tabix -p vcf "$file.gz"
done


bcftools concat --threads 14 --allow-overlaps --remove-duplicates -o ../merged.vcf.gz -Oz $(ls *.vcf.gz)
