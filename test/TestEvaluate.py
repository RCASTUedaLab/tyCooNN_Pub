import glob
import os
import sys
sys.path.append("../")
import yaml
import inference.Evaluate as ev

def testEvaluate():

    input     = "/path/to/input/parquet/files/directory"
    modeldir  = "/path/to/trained/model/"
    outcsv    = "/output/path/to/save/tRNA/profile/without/post-filter/threshold"
    outcsv2   = "/output/path/to/save/tRNA/profile/with/post-filter/threshold"
    threshold = 0.70
    ev.evaluate(input, modeldir, outcsv, outcsv2, threshold):

testEvaluate()

