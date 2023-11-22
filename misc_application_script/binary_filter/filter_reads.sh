while read -r line;do
	conda run -n r4 Rscript filter_reads.R $line
done < 'filter_species'
