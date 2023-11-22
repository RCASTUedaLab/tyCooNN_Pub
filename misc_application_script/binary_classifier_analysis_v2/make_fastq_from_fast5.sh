while read -r line;do
	python make_fastq_from_fast5.py $line
done < 'species'
