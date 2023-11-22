import sys,os,yaml
sys.path.append("../")
import inference.Inference as inference

input_options = sys.argv[1]
with open(input_options) as f:
    opts = yaml.safe_load(f)
    inference.evaluate(opts)
