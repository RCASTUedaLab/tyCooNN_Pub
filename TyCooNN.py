import click
import os.path
import sys
import preprocess.MakeTrainingPq as mkpq
import inference.Inference as inference
import inference.Evaluate as evaluate
import training.Training as training

@click.group()
def cmd():
    pass

def main():
    cmd()

# Making parquet files for training for each tRNA
@cmd.command(name='makeParquetEach')
@click.option('-l', '--tRNAlabel')
@click.option('-i', '--indir')
@click.option('-o', '--outpq')
@click.option('-c', '--takeCount',default=12000)
@click.option('-s','--outstat')
@click.option('-p', '--paramPath',default='settings.yaml')
def makeParquetEach(paramPath, tRNALabel, indir, outpq, outstat, takeCount):

    mkpq.genaratePqForTraining(paramPath, tRNALabel, indir, outpq, outstat, takeCount)

# Making parquet files for training for all 44+2 tRNA (rcc dataset)
@cmd.command(name='makeParquetAll')
@click.option('-i','--listOfIOPath',"listOfIOPath")
@click.option('-c', '--takeCount',"takeCount",default=12000)
@click.option('-p', '--paramPath',"paramPath",default='settings.yaml')
def makeparquetall(listOfIOPath, takeCount, paramPath):

    mkpq.generatePqForTrainingAll(paramPath, listOfIOPath, takeCount)

# Train model based on rcc data (full, 12000 reads/tRNA). Internally it splits
# data into train and test (10000 reads/tRNA + 2000 reads/tRNA)
@cmd.command(name='train')
@click.option('-i', '--input')
@click.option('-o', '--outdir')
@click.option('-e', '--epoch', default=100)
@click.option('-a', '--data_augment',default=0)
def train(input, outdir,  epoch, data_augment):

    training.train(input, outdir, epoch, data_augment)

# infer on a fast5 dataset
@cmd.command(name='infer')
@click.option('-i', '--input')
@click.option('-m', '--modeldir')
@click.option('-o', '--outpath')
@click.option('-r', '--ref')
@click.option('-p', '--parampath', default='settings.yaml')
@click.option('-f', '--fast5fmt', required=False, default="M")
@click.option('-t', '--threshold', default=0.75)
def infer(input, modeldir, outpath, ref, parampath, fast5fmt, threshold):

    inference.infer(input, modeldir, outpath, ref, fast5fmt, threshold, parampath)

# infer on validation split on full data (as used in train)
# see train above.
# note this is for validation and this will not generate any fast5 file
@cmd.command(name='evaluate')
@click.option('-i', '--input') # input parquet directory
@click.option('-m', '--modeldir')
@click.option('-o', '--outcsv')
@click.option('-t', '--threshold', default=0.75)
def evaluate(input, modeldir, outcsv, threshold):
    filename,extension = os.path.splitext(outcsv)
    filename2 = filename + "_thresh"
    outcsv2 = filename2 + extension
    evaluate.evaluate(input, modeldir, outcsv, outcsv2, threshold)

if __name__ == '__main__':
    main()
