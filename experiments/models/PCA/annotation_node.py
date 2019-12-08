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
import subprocess
from sklearn.decomposition import PCA
import json
def main(argvs):
    process_id = int(argvs[1])
    state_n = int(argvs[2])
    chromosome = argvs[3]
    start_index = int(argvs[4])
    end_index = int(argvs[5])
    resolution_size = int(argvs[6])
    annotation_result_dir = argvs[7]
    print(process_id, state_n, chromosome, start_index, end_index, resolution_size, annotation_result_dir)

    E = len(configs.chromHMM_assay)
    assay_list = configs.chromHMM_assay
    K = state_n
    window_size = 100
    # load pca model
    pca = PCA(n_components=K)
    with open(os.path.join(script_dir, "../train/param_{}.json".format(K))) as param_f:
        param = json.load(param_f)
    pca.components_ = np.array(param['components_'])
    pca.mean_ = np.array(param['mean_'])
    
    if not os.path.exists(annotation_result_dir):
        try:
            os.makedirs(annotation_result_dir)
        except(FileExistsError):
            pass


    annotation_result_file = os.path.join(annotation_result_dir, chromosome + "_" + str(process_id) + ".bed")
    with open(annotation_result_file, 'w') as annotation_f:
        # start annotation
        resol_fs = []
        for assay in assay_list:
            resol_fs.append(open(os.path.join(
                                        configs.blacklist_rm_data_path,
                                        assay + "/resolution-" + str(resolution_size) + "bp/",
                                        chromosome + ".bed"), "r")
                            )
        # read in genome data 
        done = False
        while not done:
            data_array = [[] for i in range(E)]
            start_annotation_index = None
            for i in range(E):
                resol_f = resol_fs[i]
                g = 0
                while g < window_size:
                    line = resol_f.readline()
                    if not line:
                        done = True
                        # raise(Exception("Reach end of file, wroing End index."))
                        break
                    inf = line.strip('\n').split()
                    index = int(inf[0])
                    signal = float(inf[1])
                    if index < start_index:
                        continue     
                    if start_annotation_index is None:
                        start_annotation_index = index
                    signal = np.arcsinh(signal)
                    data_array[i].append(signal)
                    print(index)
                    print("append" + assay_list[i] + ":" + str(index))
                    g += 1
                    if index >= end_index:
                        done = True
                        break
            # currently skip the data when it is too small
            # a better implementation is to fetch the data from former position
            # test if data array lines up
            is_lined_up = True
            for myarray in data_array:
                if len(myarray) != len(data_array[0]):
                    is_lined_up = False
                    for array in data_array:
                        print(len(array))
                    break    
            if g >= window_size/2 and is_lined_up:
                x_m = np.asarray(data_array, dtype=np.int8)
                result = pca.transform(x_m.T) # shape (G,E)
                index = start_annotation_index
                print(result.shape)
                for i in range(g):
                    print("{} {} ".format(index*resolution_size, (index+1)*resolution_size), end="", file=annotation_f)
                    for k in range(K):
                        print("{} ".format(result[i, k]),end="",file=annotation_f)
                    print(file=annotation_f)
                    index += 1

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
