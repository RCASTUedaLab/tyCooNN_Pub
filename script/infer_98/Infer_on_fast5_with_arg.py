import sys
sys.path.append("../")
import yaml
import inference.Inference as inference

import os

def testEvaluate(opts):
    inference.evaluate(opts)

input_options = sys.argv[1]
inp_loc = sys.argv[2]
out_loc = sys.argv[3]

with open(input_options) as f:
    opts = yaml.safe_load(f)
    opts['inp_loc'] = inp_loc
    opts['out_loc'] = out_loc
    testEvaluate(opts) 
