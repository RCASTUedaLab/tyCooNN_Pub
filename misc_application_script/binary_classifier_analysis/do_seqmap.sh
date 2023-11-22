while read -r line;do
	python seqmap.py $line
done < 'species_confusion'
