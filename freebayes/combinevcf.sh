for file in path/to/data/*.vcf; do
    bgzip -c "$file" > "$file.gz"
    tabix -p vcf "$file.gz"
done


bcftools concat --threads 14 --allow-overlaps --remove-duplicates -o ../merged.vcf.gz -Oz $(ls *.vcf.gz)
