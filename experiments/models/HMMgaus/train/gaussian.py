import os
import sys
is_exp_script = True
script_dir = os.path.split(os.path.realpath(__file__))[0]
import util.myconfig as myconfig
import util.utility as util
if is_exp_script:
    root_path = os.path.join(script_dir, "../../../")
    configs = myconfig.config(root_path, log_dir = script_dir)
else:
    root_path = os.path.join(script_dir, "../../")
    configs = myconfig.config(root_path)

import numpy as np
from pomegranate import *
resolution = 200
assay_list = configs.chromHMM_assay
E = len(assay_list)
pilot_index_region = {}


def read_pilot_region():
    '''
    open the pilot file
    record the pilot location
    :return:
    '''
    global pilot_index_region
    with open(configs.pilot_index_region(resolution), "r") as pilot_index_region_f: 
        while True:
            line = pilot_index_region_f.readline()
            if not line:
                break
            inf = line.strip('\n').split()
            chromosome = inf[0]
            start = int(inf[1])
            end = int(inf[2])
            if inf[0] not in pilot_index_region:
                pilot_index_region[chromosome] = [(start, end)]
            else:
                pilot_index_region[chromosome].append((start, end))
    # # for test
    # pilot_index_region = {"chr1":pilot_index_region["chr1"][:2], "chr12":pilot_index_region["chr12"][:2]}

def fetch_data(chromosome, start_index, end_index):
    resolution_fs = []
    for i in range(len(assay_list)):
        resolution_file = os.path.join(configs.processed_pilot_data_path, assay_list[i] + "/resolution-" + str(resolution) + "bp/", chromosome + ".begGraph")
        resolution_fs.append(open(resolution_file, 'r'))
    data_array = [[] for f in resolution_fs] # shape (E, n)
    done = False
    while not done:
        for i in range(len(resolution_fs)):
            resolution_f = resolution_fs[i]
            line = resolution_f.readline().strip('\n')
            if not line:
                print(chromosome, start_index, end_index, "-", index)
                raise Exception("invalid line")
            inf = line.split()
            index = int(inf[0])
            signal = float(inf[1])
            if index < start_index:
                continue
            data_array[i].append(np.arcsinh(signal))
            if index == end_index:
                done = True
    return np.asarray(data_array).T

def main(argv):
    K = int(argv[0])
    # start a log
    configs.openLog()
    configs.toLog("E {}".format(E))
    configs.toLog("K {}".format(K))
    configs.toLog("assay {}".format(" ".join(assay_list)))
    configs.toLog("resolution {}".format(resolution))
    configs.closeLog()
    # fill up the data
    sequences = [] # shape(m, (n,E))
    read_pilot_region()
    for chromosome in pilot_index_region:
        for region in pilot_index_region[chromosome]:
            sequence = fetch_data(chromosome, region[0], region[1])
            print(sequence.shape)
            sequences.append(sequence)
    mymodel=HiddenMarkovModel.from_samples(MultivariateGaussianDistribution, n_components=K,X=sequences.copy(), init='first-k', max_iterations=15, batch_size=1000, batches_per_epoch=10, n_jobs=4, verbose=True)
    with open(os.path.join(script_dir, "param.json"), "w") as param_f:
        print(mymodel.to_json(), file=param_f)
    return 0
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
