ls infer_filter/*.yaml | while read -r line;do python Infer_on_fast5.py $line;done
