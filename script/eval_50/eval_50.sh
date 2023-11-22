
while read -r species;do
	sed "s/input_species/${species}/g" gen.yaml > test.yaml
	echo python ../Infer_on_pq.py test.yaml
	python ../Infer_on_pq.py test.yaml
	rm test.yaml
done < 'rcc_species'
