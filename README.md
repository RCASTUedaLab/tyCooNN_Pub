# tyCooNN        <sub><sup>[tīˈko͞on] </sup></sub>

 About : tRNA expression analysis using 1D-CNN

## Commands

  Usage: tyCooNN.py [OPTIONS] COMMAND [ARGS]...

   Options:
     --help  Show this message and exit.
   
   Commands:
     analysis
     infer
     makeparquetAll
     makeparquetEach
     train
     
   Additional commands: 
     evaluate 
     handyevaluate 
     train2

### 1.　makeparquet
  - Prepare training dataset for each tRNA,Read fast5, format (trimming of custom Adapter) and convert it to parquet file
     

    python tyCooNN.py makeparquetEach     
     
        -l, --tRNAlabel
        -i, --indir
        -o, --outpq
        -c, --takeCount,default=12000
        -p, --paramPath,default='settings.yaml'

  or

    python tyCooNN.py makeparquetAll

        -ls, --listOfIOPath
        -p, --paramPath,default='settings.yaml'

![image](https://user-images.githubusercontent.com/70622849/140273121-1d7312ee-d1e9-4891-aa3d-dc2a4a86853d.png)

### 2. train 
  - train CNN model using isolated tRNA data sets and save weight
     

    python tyCooNN.py train
    
        -i, --input
        -o, --outdir
        -e, --epoch,default=50
        -a, --data_augment,default=3        

![image](https://user-images.githubusercontent.com/70622849/140274886-1758e556-b769-4088-b1be-a6ff77659b8f.png)

### 3. evaluate
  - test accuracy of model using isolated tRNA data sets
    
    optionally use a post-filter threshold on soft-max probabilities
    

    python tyCooNN.py evaluate
    
       -i, --input
       -o, --outdir
       -c, --csvout
       -t, --threshold [default: 0.75]


![image](https://user-images.githubusercontent.com/70622849/140274997-45208886-d4f7-4a21-846f-7a7ab42dd6d2.png)

### 4. infer

   - Inference using data of total tRNA and classify

     optionally use a post-filter threshold on soft-max probabilities

    python tyCooNN.py infer   
   
       -p, --paramPath,default='settings.yaml'
       -i, --input     # input directory for fast5 files
       -m, --modeldir  # location of trained model
       -o, --outpath   # location of output data
       -r, --ref       # reference fasta file for each tRNA species
       -f, --fast5fmt  # 'S' or 'M' for single-fast5 or multi-fast5 output
       -t, --threshold # post-filter threshold


![image](https://user-images.githubusercontent.com/70622849/140275201-811d7f05-112e-4609-acfd-bb800f088a83.png)
    
### 5. handyevaluate
  - Inference on any parquet file (just for a handy feature)
    
    note this will not generate any fast5 file however it can be applied
    to a pre-made parquet file
    

    python tyCooNN.py handyevaluate   

       -i, --input     # input directory for fast5 files
       -m, --modeldir  # location of trained model
       -o, --outcsv    # file-name for abundances of species
       -t, --threshold # post-filter threshold

### 5. train2

- misc. other training when we pre-separate data into test and train


    python tyCooNN.py train2
    
        -t, --train   # location of parquet files holding train data
        -v, --test    # location of parquet files holding test data
        -l, --label   # text file including names of tRNAs for parquet files
        -o, --outdir  # trained model directory
        -e, --epoch,default=100
        -a, --data_augment,default=0 

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
            filterflg: filter flag 
            trimSuccess: trimSuccess　

