# This is for only 46 tRNA inference. We already run 44 tRNA inference while making input yaml file, basecalling.
# That means files are already prepared, and now just run.
while read -r line;do
  
  echo $line

  python Infer_on_fast5_with_io.py $line

done < 'dbname/lb'
