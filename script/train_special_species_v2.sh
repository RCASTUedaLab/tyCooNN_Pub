while read -r line;do
	comment=`echo $line | grep "^#" | wc -l`
	if [[ $comment == "1" ]];then
		echo "Not Doing $line"
	else
		echo "Doing $line"
		python train_binary_v2.py $line
		python train_binary_with_aug_v2.py $line
	fi
done < 'special_species'
