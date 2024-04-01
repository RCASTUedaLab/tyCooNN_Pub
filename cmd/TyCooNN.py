import click

import sys
sys.path.append("../")
import preprocess.MakeTrainingPq as mkpq

import inference.Inference as inference

@click.group()
def cmd():
    pass

def main():
    cmd()

@cmd.command(name='makeParquetEach')
@click.option('-l', '--tRNAlabel')
@click.option('-i', '--indir')
@click.option('-o', '--outpq')
@click.option('-c', '--takeCount',default=12000)
@click.option('-s','--outstat')
@click.option('-p', '--paramPath',default='settings.yaml')
def makeParquetEach(paramPath,tRNALabel,indir,outpq,outstat,takeCount):

    mkpq.genaratePqForTraining(paramPath,tRNALabel,indir,outpq,outstat,takeCount)


@cmd.command(name='makeParquetAll')
@click.option('-i','--listOfIOPath',"listOfIOPath")
@click.option('-c', '--takeCount',"takeCount",default=12000)
@click.option('-p', '--paramPath',"paramPath",default='settings.yaml')
def makeparquetall(listOfIOPath,takeCount,paramPath):

    mkpq.generatePqForTrainingAll(paramPath,listOfIOPath,takeCount)


@cmd.command(name='train')
@click.option('-i', '--input')
@click.option('-o', '--outdir')
@click.option('-e', '--epoch',default=50)
@click.option('-a', '--data_argument',default=50)
def train(input, outdir,  epoch,data_argument):

    traning.train(input, outdir, epoch,data_argument)


@cmd.command(name='evaluate')
@click.option('-i', '--input')
@click.option('-o', '--outdir')
@click.option('-c', '--csvout')
def evaluateTest(input, outdir, csvout):

    evaluate.evaluate(input, outdir, csvout)


@cmd.command(name='analysis')
@click.option('-p', '--paramPath',default='settings.yaml')
@click.option('-i', '--indir')
@click.option('-c', '--configdir')
@click.option('-o', '--outpath')
@click.option('-f', '--fasta')
@click.option('-f5', '--fasta5out',required=False,default="m")
def analysis(paramPath,indirs,configdir,outpath,fasta,fasta5out):

    inference.evaluate(paramPath,indirs,configdir,outpath,fasta,fasta5out)


if __name__ == '__main__':
    main()
