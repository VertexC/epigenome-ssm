"""
This experiment is to generate the r2 score of linear regression result of 
annotation result to predict the level of gene expression
"""
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

import math
import numpy as np
import matplotlib as matplt
matplt.use('TkAgg') 
import sklearn

import matplotlib.pyplot as plt
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score

gexp_region = {} #{chromosome:[[TTS,TSS,response]}
test_gexp = {}
train_gexp = {}
SEARCH_REGION = 110 * 1000 # 110kb
resolution = 1000 # 1kb
state_n = None
inside_interval_n = 10
exp_id = None
def readin_gexp_region():
    with open(configs.gexp_file, 'r') as gexp_f:
        while True:
            line = gexp_f.readline().strip('\n')
            if not line:
                break
            inf = line.split()
            chromosome = inf[0]
            start = int(inf[1])
            end = int(inf[2])
            rpkm = np.arcsinh(float(inf[3]))
            strand = int(inf[4])
            if chromosome not in gexp_region:
                gexp_region[chromosome] = [(start, end, rpkm, strand)]
            else:
                gexp_region[chromosome].append((start, end, rpkm, strand))

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
    for chromosome in gexp_region:
        for inf in gexp_region[chromosome]:
            if flag == 0:
                if chromosome not in test_gexp:
                    test_gexp[chromosome] = [inf]
                else:
                    test_gexp[chromosome].append(inf)
            else:
                if chromosome not in train_gexp:
                    train_gexp[chromosome] = [inf]
                else:
                    train_gexp[chromosome].append(inf)

            flag += 1
            flag %= 10

def fectch_predictor_per_index(chromosome, base, final, resolution):
    """
    TODO: add the weight according to the position
    """
    index = 0
    remain = np.zeros((1,state_n))
    remain_count = 0
    result = None
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
            if exp_id == "exp26-1":
                state -= 1
            # configs.toLog("\nstart:{} end:{}".format(start, end))
            # within the frame
            if end <= base:
                continue
            if start >= final:
                break  # finish the resolution on this frame
            if start < base and end > base:
                start = base
            if start < final and end > final:
                end = final
            if  (index + 1) * resolution + base  >= end:
                # configs.toLog("1. Larger than end")
                remain += (end - start) * bin_vector(state)
                remain_count += (end - start)
            else:
                # configs.toLog("2. Smaller than end")
                while (index + 1) * resolution + base < end: 
                    if (index+1) * resolution + base > start:
                        valid_bp = (index + 1) * resolution + base - max(index * resolution + base, start)
                        temp_result = (remain + valid_bp * bin_vector(state)) / (valid_bp + remain_count) 
                        if result is None:
                            result = temp_result
                        else:
                            result = np.vstack((result, temp_result))
                        remain = np.zeros((1,state_n))
                        remain_count = 0
                        index += 1
                    else:
                        return None
                    # configs.toLog("index {} Increase {},{}".format(index, index * resolution + base, (index+1)*resolution + base))
                remain = (end - max(start, index * resolution + base)) * bin_vector(state)
                remain_count = (end - max(start, index * resolution + base))
        # from (index * resolution + base, (index+1)*resolution + base
        if remain_count != 0:
            if result is None:
                result = remain / remain_count
            else:
                result = np.vstack((result, remain / remain_count))
    return result


def gexp_validation(start_index, end_index, start_offset, end_offset, index_n):
    """
    do train and test on both 
    """
    predictors = []
    response = []
    regrs = []
    for chromosome in train_gexp:
        if chromosome not in configs.chromosome_list:
            continue
        print(chromosome)
        for inf in train_gexp[chromosome]:
            base = inf[start_index] + start_offset
            final = inf[end_index] + end_offset
            rpkm = inf[2]
            if base < 0:
                continue
            resolution = math.floor((final + 1 - base) / index_n)
            result = fectch_predictor_per_index(chromosome, base, final, resolution)
            if result is None:
                continue
            m,n = result.shape
            if m >= index_n:
                # print(m)
                result = result[:index_n,:]
            else:
                continue
            strand = inf[3]
            if strand == -1:
                result = result[::-1, :]
            response.append(rpkm)
            predictors.append(result)

    response = np.asarray(response).T
    for index in range(index_n):
        regr_per = linear_model.LinearRegression()
        predictor_per = None
        for predictor in predictors:
            if predictor_per is None:
                predictor_per = predictor[index, :]
            else:
                predictor_per = np.vstack((predictor_per, predictor[index, :]))
        regr_per.fit(predictor_per, response)
        r2 = sklearn.metrics.r2_score(response, regr_per.predict(predictor_per))
        # configs.tolog(r2) 
        regrs.append(regr_per)
    
    predictors = []
    response = []
    for chromosome in test_gexp:
        if chromosome not in configs.chromosome_list:
            continue
        print(chromosome)
        for inf in test_gexp[chromosome]:
            base = inf[start_index] + start_offset
            final = inf[end_index] + end_offset
            rpkm = inf[2]
            if base < 0:
                continue
            resolution = math.floor((final + 1 - base) / index_n)
            result = fectch_predictor_per_index(chromosome, base, final, resolution)
            if result is None:
                continue
            m,n = result.shape
            if m >= index_n:
                result = result[:index_n, :]
            else:
                continue
            if strand == -1:
                result = result[::-1, :]
            response.append(rpkm)
            predictors.append(result)
    response = np.asarray(response).T
    for index in range(index_n):
        regr_per = linear_model.LinearRegression()
        predictor_per = None
        for predictor in predictors:
            if predictor_per is None:
                predictor_per = predictor[index, :]
            else:
                predictor_per = np.vstack((predictor_per, predictor[index, :]))
        r2 = sklearn.metrics.r2_score(response, regrs[index].predict(predictor_per))
        configs.toLog(r2)

def main(args):
    global inside_interval_n
    global exp_id
    global state_n
    if len(args) == 0:
        raise(Exception("No given exp_id"))
    exp_id = args[0]
    state_n = int(args[1])
    readin_gexp_region()
    split()
    configs.openLog(exp_id + "_summary")
    configs.toLog("upstream")
    gexp_validation(0, 0, -SEARCH_REGION, 0, int(SEARCH_REGION/resolution))
    configs.toLog("between")
    gexp_validation(0, 1, 0, 0, inside_interval_n) # start site
    configs.toLog("downstream")
    gexp_validation(1, 1, 0, SEARCH_REGION, int(SEARCH_REGION/resolution))
    configs.closeLog()
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))