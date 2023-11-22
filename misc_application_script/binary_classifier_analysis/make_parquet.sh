while read -r line;do
	python make_parquet.py $line
done < 'species'
