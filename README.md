# tyCooNN        <sub><sup>[tīˈko͞on] </sup></sub>

 About : tRNA expression Analysis using 1D-CNN


## Commands

  Usage: tyCooNN.py [OPTIONS] COMMAND [ARGS]...

   Options:
     --help  Show this message and exit.
   
   Commands:
     analysis
     evaluatetest
     makeparquetall
     makeparqueteach
     train

### 1.　makeparquet
  - Prepare training dataset for each tRNA,Read fast5, format(trimmed) and convert it to parquet file
     
    python tyCooNN.py makeparqueteach     
     
        -l, --tRNAlabel
        -i, --indir
        -o, --outpq
        -c, --takeCount,default=12000
        -p, --paramPath,default='settings.yaml'

     or

    python tyCooNN.py makeparquetall

        -ls, --listOfIOPath
        -p, --paramPath,default='settings.yaml'

![image](https://user-images.githubusercontent.com/70622849/140273121-1d7312ee-d1e9-4891-aa3d-dc2a4a86853d.png)

### 2. train 
  - train CNN model using isolated tRNA data sets and save weight
     
    python tyCooNN.py train
        
        -i, --input
        -o, --outdir
        -e, --epoch,default=50
        -a, --data_argument,default=50        

![image](https://user-images.githubusercontent.com/70622849/140274886-1758e556-b769-4088-b1be-a6ff77659b8f.png)

### 3. evaluatetest
  - test accuracy of model using isolated tRNA data sets
     
    python tyCooNN.py evaluatetest
    
       -i, --input
       -o, --outdir
       -c, --csvout


![image](https://user-images.githubusercontent.com/70622849/140274997-45208886-d4f7-4a21-846f-7a7ab42dd6d2.png)

### 4. analysis

   - Inference using data of total tRNA and classify
     
    
    python tyCooNN.py analysis   
   
       -p, --paramPath,default='settings.yaml'
       -i, --indir
       -c, --configdir
       -o, --outpath
       -f, --fasta


![image](https://user-images.githubusercontent.com/70622849/140275201-811d7f05-112e-4609-acfd-bb800f088a83.png)
    
 
## Parameter file
-- example of parameter file

   max_core: 5  #maxcore　<br>
   #for filterling　<br>
   qval_min: 4  #q-value filter  threshold　<br> 
   delta_min: 40 <br> 
   delta_max: 65 <br>
   readlen_min: 50　<br>
   readlen_max: 200　<br>
   signallen_max: 20000　<br>
   duratio_rate_max: 1.1　<br>
   #for trimming by HMM　<br>
   meantoSet: 92　<br>
   adap1thery: 75　<br>
   adap2thery: 125　<br>
   #for trriming by Aligment　<br>
   firstAdaptor: GUAUCCUCCCAGAGAGAGAGAGAGAGAGAGAGA　<br>
   #trimlen　<br>
   trimlen: 8192　<br>

## fast5 output format

  Analysis/Basecall_1D_099/BaseCalled_template/Fast5 ： tRNA reference sequence will be add for each read
  
  Analysis/Basecall_1D_099　:following data will be add as a attribute to the group <br>
  
            tRNA : tRNA label　
            tRNAIndex :tRNA Index　
            value: inference value (0-1)　
            filterpass: filter pass　
            filterflg: filter flg 
            trimSuccess: trimSuccess　
  
  
## Additional information

  [Readme more](./README_ADDITIONAL)

