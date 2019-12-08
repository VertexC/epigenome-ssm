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
import matplotlib as matplt
matplt.use('TkAgg') 
import sklearn

import matplotlib.pyplot as plt
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score

region = {} #{chromosome:[[start,end,response]}
test_set = {}
train_set = {}
state_n = None
SEARCH_REGION = 10 * 1000 # 10kb
exp_id = None
def readin_region():
    with open(configs.gexp_file, 'r') as gexp_f:
        while True:
            line = gexp_f.readline().strip('\n')
            if not line:
                break
            inf = line.split()
            chromosome = inf[0]
            start = int(inf[1])
            end = int(inf[2])
            value = np.arcsinh(float(inf[3]))
            strand = int(inf[4])
            if chromosome not in region:
                region[chromosome] = [(start, end, value, strand)]
            else:
                region[chromosome].append((start, end, value, strand))

def bin_vector(state):
    global state_n
    result = np.zeros((1, state_n))
    result[0,state] = 1
    return result

def split():
    """
    90% to train
    10% to test
    """
    flag = 0
    for chromosome in region:
        for inf in region[chromosome]:
            if flag == 0:
                if chromosome not in test_set:
                    test_set[chromosome] = [inf]
                else:
                    test_set[chromosome].append(inf)
            else:
                if chromosome not in train_set:
                    train_set[chromosome] = [inf]
                else:
                    train_set[chromosome].append(inf)

            flag += 1
            flag %= 10

def fectch_predictor_avg(chromosome, base, final):
    global exp_id
    index = 0
    result = np.zeros((1,state_n))
    count = 0
    # configs.toLog("base:{} final:{}".format(base,final))
    with open(os.path.join(configs.annotation_result_dir,exp_id, chromosome + ".bed"), 'r') as annotation_result_f:
        while True:
            line = annotation_result_f.readline()
            if not line:
                break
            inf = line.strip('\n').split()
            start = int(inf[0])
            end = int(inf[1])
            state = int(inf[2])
            # within the frame
            if end <= base:
                continue
            if start >= final:
                break  
            if start < base and end > base:
                start = base
            if start < final and end > final:
                end = final
            # for chromHMM's result
            if exp_id == "exp26-1":
                state -= 1
            result += (end - start) * bin_vector(state)
            count += (end - start)
    if count == 0:
        return None
    else:
        return result/count

def validation():
    """ 
    do train and test on both 
    """
    global SEARCH_REGION
    predictor = None
    response = []
    count = 0
    for chromosome in train_set:
        if chromosome not in configs.chromosome_list:
            continue
        for inf in train_set[chromosome]:
            strand = inf[3]
            if strand == 1:
                base = inf[0] - SEARCH_REGION
                final = inf[0] 
            else:
                base = inf[1]
                final = inf[1] + SEARCH_REGION
            value = inf[2]
            if base < 0:
                continue
            result = fectch_predictor_avg(chromosome, base, final)
            if result is None:
                continue
            response.append(value)
            if predictor is None:
                predictor = result
            else:
                predictor = np.vstack((predictor, result))
            count += 1
    print("in train:", predictor.shape)
    response = np.asarray(response).T
    regr = linear_model.LinearRegression()
    regr.fit(predictor, response)
    
    pre_response = regr.predict(predictor)
    adj_r2 = util.adj_r2_score(response, pre_response,count,state_n)
    r2 = sklearn.metrics.r2_score(response, pre_response)
    configs.toLog("train r2:{}".format(r2))
    configs.toLog("train adjr2:{}".format(adj_r2))

    predictor = None
    response = [] 
    count = 0
    for chromosome in test_set:
        if chromosome not in configs.chromosome_list:
            continue
        for inf in test_set[chromosome]:
            strand = inf[3]
            if strand == 1:
                base = inf[0] - SEARCH_REGION
                final = inf[0] 
            else:
                base = inf[1]
                final = inf[1] + SEARCH_REGION
            value = inf[2]
            if base < 0:
                continue
            result = fectch_predictor_avg(chromosome, base, final)
            if result is None:
                continue
            response.append(value)
            if predictor is None:
                predictor = result
            else:
                predictor = np.vstack((predictor, result))
            count += 1
    print("in test:", predictor.shape)
    pre_response = regr.predict(predictor)
    adj_r2 = util.adj_r2_score(response, pre_response, count, state_n)
    r2 = sklearn.metrics.r2_score(response, pre_response)
    configs.toLog("test r2:{}".format(r2))
    configs.toLog("test adjr2:{}".format(adj_r2))

def main(args):
    global exp_id
    global state_n
    exp_id = args[0]
    state_n = int(args[1])
    configs.openLog(log_name=exp_id)
    readin_region()
    split()
    validation()
    configs.closeLog()
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
