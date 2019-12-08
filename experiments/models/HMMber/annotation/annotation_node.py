import os
import sys
is_exp_script = True
script_dir = os.path.split(os.path.realpath(__file__))[0]
#sys.path.append("../../../src")
import util.myconfig as myconfig
import util.utility as util
data_type = "fc"
if is_exp_script:
    root_path = os.path.join(script_dir, "../../../")
    configs = myconfig.config(root_path, log_dir = script_dir, data_type=data_type)
else:
    root_path = os.path.join(script_dir, "../../")
    configs = myconfig.config(root_path, data_type=data_type)

import models.state_space_model.ssm as ssm
import numpy as np
import subprocess
from pomegranate import *

def main(argvs):
    process_id = int(argvs[1])
    train_dir = argvs[2]
    chromosome = argvs[3]
    start_index = int(argvs[4])
    end_index = int(argvs[5])
    print(process_id, train_dir, chromosome, start_index, end_index)
    # set up parameter file
    with open(os.path.join(train_dir, "log.txt"), 'r') as param_f:
        E = int(param_f.readline().strip('\n').split()[-1])
        K = int(param_f.readline().strip('\n').split()[-1])
        assay_list = param_f.readline().strip('\n').split()[1:]
        # only need is resolution size here
        resolution_size = int(param_f.readline().strip('\n').split()[-1])
        if data_type == "pval":
            threshold = int(param_f.readline().strip('\n').split()[-1])
        elif data_type == "fc":
            next(param_f) # skip a line
            threshold = np.squeeze(np.asarray(util.read_matrix(param_f, 1, E))) # reshape to 1D vector
    with open(os.path.join(train_dir, "param.json"), 'r') as json_f:
        param_json = json_f.read()
    mymodel = HiddenMarkovModel.from_json(param_json)
    window_size = 100
    # copy the log file as reference 
    result = subprocess.call(" ".join(["cp", os.path.join(train_dir, "log.txt"), os.path.join(script_dir, "log.txt")]), shell=True)
    dis_annotation_result_dir = os.path.join(script_dir, "dis-annotation-result/")
    con_annotation_result_dir = os.path.join(script_dir, "con-annotation-result/")

    if not os.path.exists(dis_annotation_result_dir):
        try:
            os.makedirs(dis_annotation_result_dir)
        except(FileExistsError):
            pass
    
    if not os.path.exists(con_annotation_result_dir):
        try:
            os.makedirs(con_annotation_result_dir)
        except(FileExistsError):
            pass


    dis_annotation_file = os.path.join(dis_annotation_result_dir, chromosome + "_" + str(process_id) + ".bed")
    con_annotation_file = os.path.join(con_annotation_result_dir, chromosome + "_" + str(process_id) + ".bed")
    with open(dis_annotation_file, 'w') as dis_annotation_f, open(con_annotation_file, 'w') as con_annotation_f:
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
                if type(threshold) == int: # in the case of pval data
                    thr = threshold
                else:
                    thr = threshold[i] # in the case of fc data
                while g < window_size:
                    line = resol_f.readline()
                    if not line:
                        done = True
                        # raise(Exception("Reach end of file, wroing End index."))
                        break
                    inf = line.strip('\n').split()
                    index = int(inf[0])
                    print(assay_list[i] + ":" + str(index))
                    signal = float(inf[1])
                    if index < start_index:
                        continue     
                    if start_annotation_index is None:
                        start_annotation_index = index
                    signal = int(signal >= thr)
                    data_array[i].append(signal)
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
                x_m = np.asarray(data_array, dtype=np.int8).T
                dis_result = mymodel.predict(x_m.copy()) # 
                con_result = mymodel.predict_proba(x_m.copy())
                # discrete result
                index = start_annotation_index
                for i in range(g):
                    print("{} {} {}".format(index*resolution_size, (index+1)*resolution_size, dis_result[i]), file=dis_annotation_f)
                    index += 1
                # continuous result
                index = start_annotation_index
                for i in range(g):
                    print("{} {} ".format(index*resolution_size, (index+1)*resolution_size), end="", file=con_annotation_f)
                    for k in range(K):
                        print("{} ".format(con_result[i, k]),end="",file=con_annotation_f)
                    print(file=con_annotation_f)
                    index += 1

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
