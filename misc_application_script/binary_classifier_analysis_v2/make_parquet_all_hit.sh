while read -r line;do
	python make_parquet_all_hit.py $line
done < 'species'
