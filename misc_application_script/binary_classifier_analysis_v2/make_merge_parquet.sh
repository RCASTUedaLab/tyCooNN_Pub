while read -r line;do
	python make_merge_parquet.py $line
done < 'special_species_size'
