while read -r line;do
	python make_label_from_fast5.py $line
done < 'species'
