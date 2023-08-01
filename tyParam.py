import multiprocessing


class tyParam():

    def __init__(self, param):

        self.meantoSet = int(param['meantoSet'])
        self.adap1thery = int(param['adap1thery'])
        self.adap2thery = int(param['adap2thery'])
        self.qval_min = int(param['qval_min'])
        self.delta_min= int(param['delta_min'])
        self.delta_max = int(param['delta_max'])
        self.readlen_min = int(param['readlen_min'])
        self.readlen_max = int(param['readlen_max'])

        self.signallen_max = int(param['signallen_max'])
        self.duratio_rate_max = float(param['duratio_rate_max'])

        self.ncore = get_number_of_core(int(param['max_core']))

        # for trriming by Aligment
        self.firstAdaptor = param['firstAdaptor']
        self.trimlen = param['trimlen']


def get_number_of_core(MAX_CORE:int):
    ncore = multiprocessing.cpu_count()
    if ncore > MAX_CORE:
        ncore = MAX_CORE
    return ncore
