while read -r line;do
	bash distribute_reads.sh $line
done < 'species'
