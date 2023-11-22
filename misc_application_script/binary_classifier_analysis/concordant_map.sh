while read -r line;do
	conda run -n r4 Rscript concordant_map.R $line
done < 'special_species'
