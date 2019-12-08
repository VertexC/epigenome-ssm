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

import models.state_space_model.ssmKalman as ssm
import numpy as np
import subprocess

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
        window_size = int(param_f.readline().strip('\n').split()[-1])
        K = int(param_f.readline().strip('\n').split()[-1])
        positive_state = int(param_f.readline().strip('\n').split()[-1])
        positive_em = int(param_f.readline().strip('\n').split()[-1])
        lambda_1 = float(param_f.readline().strip('\n').split()[-1])
        lambda_2 = float(param_f.readline().strip('\n').split()[-1])
        lambda_3 = float(param_f.readline().strip('\n').split()[-1])
        assay_list = param_f.readline().strip('\n').split()[1:]
        resolution_size = int(param_f.readline().strip('\n').split()[-1])
    print(assay_list)
    theta_m = util.read_matrix(os.path.join(train_dir, "theta.txt"))
    lambda_m = util.read_matrix(os.path.join(train_dir, "lambda.txt"))
    # copy the log file as reference 
    result = subprocess.call(" ".join(["cp", os.path.join(train_dir, "log.txt"), os.path.join(script_dir, "log.txt")]), shell=True)

    annotation_result_dir = os.path.join(script_dir, "annotation-result/")
    if not os.path.exists(annotation_result_dir):
        try:
            os.makedirs(annotation_result_dir)
        except(FileExistsError):
            pass
    with open(os.path.join(annotation_result_dir, chromosome + "_" + str(process_id) + ".bed"), 'w') as annotation_f:
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
                    print(assay_list[i] + ":" + str(index))
                    signal = float(inf[1])
                    if index < start_index:
                        continue     
                    if start_annotation_index is None:
                        start_annotation_index = index
                    signal = np.arcsinh(signal)
                    data_array[i].append(signal)
                    print("append" + assay_list[i] + ":" + str(index))
                    g += 1
                    if index >= end_index:
                        done = True
                        break
            # currently skip the data when it is too small
            # a better implementation is to fetch the data from former position
            iteration = 1
            # test if data array lines up
            is_lined_up = True
            for myarray in data_array:
                if len(myarray) != len(data_array[0]):
                    is_lined_up = False
                    for array in data_array:
                        print(len(array))
                    break    
            if g >= window_size/2 and is_lined_up:
                y_m = np.asmatrix(np.asanyarray(data_array))
                mymodel = ssm.ssmKalman(P=E, T=g, M=K, positive_state=positive_state, positive_em=positive_em, lambda_1=lambda_1, lambda_2=lambda_2, lambda_3=lambda_3)
                mymodel.set_y(y_m)
                mymodel.set_z(theta_m)
                mymodel.set_t(lambda_m)
                mymodel.optimization(iteration=iteration)
                index = start_annotation_index
                for g in range(mymodel.T):
                    print("{} {} ".format(index*resolution_size, (index+1)*resolution_size), end="", file=annotation_f)
                    for k in range(mymodel.M):
                        print("{} ".format(mymodel.y_m[k, g]), end="", file=annotation_f)
                    print("\n", end="", file=annotation_f)
                    index += 1
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
