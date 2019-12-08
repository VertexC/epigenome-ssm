"""
This script is acting as the header of parrell trainning.
The main functionality is to:
1. assign chunked region to trainning
2. average the trainned theta and lambda

not completed
"""
import os
import sys
is_exp_script = True
script_dir = os.path.split(os.path.realpath(__file__))[0]
import util.myconfig as myconfig
import util.utility as util
if is_exp_script:
    root_path = os.path.join(script_dir, "../../../")
else:
    root_path = os.path.join(script_dir, "../../../")
configs = myconfig.config(root_path, log_dir=script_dir)
import numpy as np
import subprocess
import time
pilot_index_region = {}
# parameter of ssm
resolution = 200
assay_list = configs.chromHMM_assay

E = len(assay_list)
G = 100
K = 5
positive_state = 1
positive_em = 1
lambda_1 = 0.1 
lambda_2 = 0.1
lambda_3 = 0.1
configs.openLog()
configs.toLog("E {}".format(E))
configs.toLog("G {}".format(G))
configs.toLog("K {}".format(K))
configs.toLog("positive_state {}".format(positive_state))
configs.toLog("positive_em {}".format(positive_em))
configs.toLog("lambda_1 {}".format(lambda_1))
configs.toLog("lambda_2 {}".format(lambda_2))
configs.toLog("lambda_3 {}".format(lambda_3))
configs.toLog("assay {}".format(" ".join(assay_list)))
configs.toLog("resolution {}".format(resolution))
configs.closeLog()

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

def master():
    child_id = 0
    for chromosome in pilot_index_region:
        skip = 0
        for region in pilot_index_region[chromosome]:
            start = region[0]
            end = region[1]

            if child_id == 4:
                child = subprocess.Popen(["python", "train_node.py", chromosome, str(start), str(end), str(skip), str(resolution), ','.join(assay_list), str(E),str(G),str(K),str(positive_state),str(positive_em),str(lambda_1),str(lambda_2),str(lambda_3), str(child_id)])
            
            child_id += 1
            skip += end - start + 1
    # test if all sub process done
    while True:
        all_done = True
        for i in range(child_id):
            if not os.path.exists(os.path.join(script_dir, str(i) + ".done")):
                all_done = False
                break
        print("done!")
        if all_done:
            break
        else:
            time.sleep(60)
    time.sleep(30)
    # clean temp file
    for i in range(child_id):
        os.remove(os.path.join(script_dir, str(i) + ".done"))
    files = []
    for (_, _, filenames) in os.walk(script_dir):
         files.extend(filenames)
         break
    theta_m = np.asmatrix(np.zeros((E,K))) 
    lambda_m = np.asmatrix(np.zeros((K,K)))
    # assemble theta
    for myfile in files:
        if myfile[:5] == "theta":
            theta_m += util.read_matrix(os.path.join(script_dir, myfile))
            os.remove(os.path.join(script_dir,myfile))
    # assemble lambda
    for myfile in files:
        if myfile[:6] == "lambda":
            lambda_m += util.read_matrix(os.path.join(script_dir, myfile))
            os.remove(os.path.join(script_dir,myfile))
    theta_m = theta_m / child_id
    lambda_m = lambda_m / child_id
    with open("theta.txt", 'w') as theta_f:
        print(util.matrix_to_str(theta_m), file=theta_f)
    with open("lambda.txt", 'w') as lambda_f:
        print(util.matrix_to_str(lambda_m), file=lambda_f)
    
def main():
    read_pilot_region()
    master()
    return 0

if __name__ == "__main__":
    sys.exit(main())

