inp=$1

cd $inp

grep -v 'None' ${inp}.label > ${inp}_hit.label

cut -f5 -d" " ${inp}_hit.label | sort | uniq | while read -r name;do
	grep "$name" ${inp}_hit.label | cut -f1,5 -d" " > ${name}.tmp
done

cat *.tmp > add.tmp

nl=`wc -l add.tmp | cut -f1 -d" "`

yes default | head -n $nl > default.tmp

paste -d" " add.tmp default.tmp > ${inp}_distr.label

rm ${inp}_hit.label *.tmp
