while read -r line;do python make_parquet_with_filter.py $line;done < 'species'
