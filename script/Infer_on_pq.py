import sys
sys.path.append("/home/bhaskar/work/tyCooNN_2206/")
import yaml
import inference.Inference_v2 as inference

import os

def testEvaluate(opts):

    inference.evaluate_on_pq(opts)

input_options = sys.argv[1]
with open(input_options) as f:
    opts = yaml.safe_load(f)
    testEvaluate(opts) 
