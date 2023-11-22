import glob
import os
import sys
sys.path.append("../")
import yaml
import training.Evaluate as evaluate
import training.GenaratePqForTrainning as pqg

def testEvaluate(opts):

    evaluate.evaluate(opts)

input_options = sys.argv[1]
with open(input_options) as f:
    opts = yaml.safe_load(f)
    testEvaluate(opts)

