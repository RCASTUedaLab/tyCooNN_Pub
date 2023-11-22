while read -r line;do
	python make_parquet_v0.py $line
done < 'species'
