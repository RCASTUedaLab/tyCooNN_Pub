# tyCooNN        <sub><sup>[tīˈko͞on] </sup></sub>

 About : tRNA expression analysis using 1D-CNN

## Commands

  Usage: tyCooNN.py [OPTIONS] COMMAND [ARGS]...

   Options:
     --help  Show this message and exit.
   
   Commands:
     infer
     makeparquetAll
     makeparquetEach
     train
     evaluate
     
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

![image](https://github.com/RCASTUedaLab/tyCooNN_Pub/blob/013f4717bc4de9e511323606930f9933e75e1370/.doc/doc_image_01.jpg)

### 2. train 
  - train CNN model using isolated tRNA data sets and save weight
     

    python tyCooNN.py train
    
        -i, --input
        -o, --outdir
        -e, --epoch,default=50
        -a, --data_augment,default=3        

![image](https://github.com/RCASTUedaLab/tyCooNN_Pub/blob/e41c3ae29571abbd252a89b720a2880f01ce2759/.doc/doc_image_02.jpg)

### 3. Evaluate
  - test accuracy of model using isolated tRNA data sets


    

    python tyCooNN.py evaluate
    
       -i, --input
       -o, --outdir
       -c, --csvout
       -t, --threshold [default: 0.75]


![image](https://github.com/RCASTUedaLab/tyCooNN_Pub/blob/32e80c118a4288618e2917fef73b27a2a97d0c18/.doc/doc_image_03.jpg)


### 4. Inference using data of total tRNA and classify

    optionally use a post-filter threshold on soft-max probabilities

    python tyCooNN.py infer   
   
    -p, --paramPath,default='settings.yaml'
       -i, --input     # input directory for fast5 files
       -m, --modeldir  # location of trained model
       -o, --outpath   # location of output data

       -f, --fast5fmt  # 'S' or 'M' for single-fast5 or multi-fast5 output
       -t, --threshold # post-filter threshold


![image](https://github.com/RCASTUedaLab/tyCooNN_Pub/blob/13a4d929f20dd05f503c5c1f7aa806c524863f48/.doc/doc_image_04.jpg)
    

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

