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
    # set argvs
    chromosome = argvs[1]
    start_index = int(argvs[2])
    end_index = int(argvs[3])
    skip = int(argvs[4])

    resolution = int(argvs[5])
    assay_list = list((argvs[6]).split(','))

    E = int(argvs[7])
    G = int(argvs[8])
    K= int(argvs[9])
    positive_state = int(argvs[10])
    positive_em = int(argvs[11])
    lambda_1 = float(argvs[12])
    lambda_2 = float(argvs[13])
    lambda_3 = float(argvs[14])
    process_id = int(argvs[15])
    print(chromosome, start_index, end_index, skip, resolution, ','.join(assay_list), E, G, K, positive_state, positive_em, lambda_1, lambda_2, lambda_3, process_id)

    # apply ssm
    # read in data to SS_model
    theta_lhs = np.asmatrix(np.zeros((K, K)))
    theta_rhs = np.asmatrix(np.zeros((K, E)))
    lambda_lhs = np.asmatrix(np.zeros((K, K)))
    lambda_rhs = np.asmatrix(np.zeros((K, K)))
    window_size = 100
    done = False
    print("Train ssm on" + chromosome + "(" +str(start_index)+ "," +str(end_index) + ")")

    resolution_files = []
    skipped_files = []
    for i in range(len(assay_list)):
        resolution_file = os.path.join(configs.processed_pilot_data_path, assay_list[i] + "/resolution-" + str(resolution) + "bp/", chromosome + ".begGraph")
        resolution_files.append(resolution_file)
        skipped_file = os.path.join(script_dir, assay_list[i] + chromosome + str(skip) + ".txt")
        skipped_files.append(skipped_file)
        result = subprocess.call(" ".join(["tail", "-n", "+" + str(skip+1), resolution_file, ">", skipped_file]), shell=True)
        print(skipped_file, result)
    # get resolution chromosome signal from raw data
    skipped_fs = []
    for skipped_file in skipped_files:
        skipped_fs.append(open(skipped_file, 'r'))
    # read data
    is_first_window = True
    is_first_index = True # for test
    while not done:
        g = 0
        data_array = [[] for i in range(len(assay_list))]
        while True:
            for i in range(len(skipped_fs)):
                skipped_f = skipped_fs[i]
                line = skipped_f.readline().strip('\n')
                if not line:
                    raise Exception("invalid line")
                inf = line.split()
                index = int(inf[0])
                signal = float(inf[1])
                if is_first_index:
                    if index != start_index:
                        raise Exception("start index not begin with 0")
                    is_first_index = False
                data_array[i].append(np.arcsinh(signal))
            g += 1
            if index == end_index:
                done = True
                break
            if g == window_size:
                break
        # apply ssm
        if g >= window_size/2:
            iteration = 30
            print("=========Another region=============")
            model = ssm.ssmKalman(P = E, T = g, M=K, positive_state=positive_state, positive_em=positive_em,  lambda_1=lambda_1, lambda_2=lambda_2,lambda_3=lambda_3,verbose=False)
            # if is_first_window:
            #    is_first_window = False
            # else:
            #    model.update_z(theta_lhs, theta_rhs)
            #    model.set_t(np.asmatrix(np.linalg.lstsq(lambda_lhs, lambda_rhs, rcond=1)[0]).T)
            # import pdb; pdb.set_trace();
            print('data mean\n', np.mean(data_array))
            model.set_y(np.asmatrix(data_array))
            print('model y mean\n', np.mean(model.y_m))
            model.optimization(iteration=iteration)
            print('state average\n', np.mean(model.alpha_m, axis=1))
            print('error\n', model.error_m)
            temp_theta_lhs, temp_theta_rhs = model.get_z_msg()
            temp_lambda_lhs, temp_lambda_rhs = model.get_t_msg()
            theta_lhs += temp_theta_lhs
            theta_rhs += temp_theta_rhs
            lambda_lhs += temp_lambda_lhs
            lambda_rhs += temp_lambda_rhs
    model = ssm.ssmKalman(P = E, T = g, M=K, positive_state=positive_state, positive_em=positive_em,  lambda_1=lambda_1, lambda_2=lambda_2,lambda_3=lambda_3)
    model.update_z(theta_lhs, theta_rhs)
    theta_m = model.z_m
    lambda_m = np.asmatrix(np.linalg.lstsq(lambda_lhs, lambda_rhs, rcond=1)[0]).T
    # print to lambda and theta record
    util.output_matrix(theta_m, os.path.join(script_dir, "theta_" + str(process_id) + ".txt"))
    util.output_matrix(lambda_m, os.path.join(script_dir, "lambda_" + str(process_id) + ".txt"))
    # clean up the skipped_files
    for skipped_f in skipped_fs:
        skipped_f.close()
    for skipped_file in skipped_files:
        os.remove(skipped_file)
    open(os.path.join(script_dir, str(process_id) + ".done"), 'w')
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
