'''
    intersection of two tuple
'''
import math
import numpy as np
import sklearn
from sklearn.metrics import mean_squared_error, r2_score

def intersect(a: tuple, b: tuple):
    if a[0] < b[0]:
        if a[1] <= b[0]:
            return None
        elif a[1] > b[1]:
            return b
        elif a[1] <= b[1]:
            return b[0], a[1]
    else:
        if a[0] >= b[1]:
            return None
        elif a[1] > b[1]:
            return a[0], b[1]
        elif a[1] <= b[1]:
            return a

def test_interset():
    print(intersect((1, 4), (0, 2)),
          intersect((1, 4), (0, 0.5)),
          intersect((1, 4), (2, 3)),
          intersect((1, 4), (3, 5)),
          intersect((1, 4), (5, 6)),
          intersect((1, 4), (0, 6)))

def rm_intersect(a:tuple, b:tuple, exclude=False):
    '''
    :param a: region A
    :param b: region B
    :return: an array of region from A - A intersect B
    '''
    result = []
    inter = intersect(a, b)
    if inter is not None:
        if inter[0] != a[0]:
            if exclude:
                if inter[0]-1 >= a[0]:
                    result.append((a[0], inter[0]-1))
            else:
                result.append((a[0], inter[0]))
        if inter[1] != a[1]:
            if exclude:
                if inter[1]+1 <= a[1]:
                    result.append((inter[1]+1, a[1]))
            else:
                result.append((inter[1], a[1]))
    else:
        result.append(a)
    return result

def test_rm_interset():
    print(rm_intersect((1, 4), (0, 2)),
          rm_intersect((1, 4), (0, 0.5)),
          rm_intersect((1, 4), (1, 4)),
          rm_intersect((1, 4), (2, 3)),
          rm_intersect((1, 4), (3, 5)),
          rm_intersect((1, 4), (5, 6)),
          rm_intersect((1, 4), (0, 6)))
    print(rm_intersect((1, 4), (0, 2), exclude=True),
          rm_intersect((1, 4), (0, 0.5), exclude=True),
          rm_intersect((1, 4), (1, 4), exclude=True),
          rm_intersect((1, 4), (2, 3), exclude=True),
          rm_intersect((1, 4), (3, 5), exclude=True),
          rm_intersect((1, 4), (5, 6), exclude=True),
          rm_intersect((1, 4), (0, 6), exclude=True))
'''
    probability  
'''


def log_add(p1, p2, base=math.e):
    return p1 + math.log(1 + math.exp(p2 - p1), base)


def log_multi(p1, p2, base=math.e):
    return math.log(p1, base) + math.log(p2, base)


def test_on_log_pro():
    p1 = math.exp(-30)
    p2 = 2 * math.exp(-30)
    pass


def create_resolution(in_file, out_file, base, final, resolution=100):
    '''
    make the chrs into resolution(window), default as 100bp
    base,final for start and end of frame
    :param in_file: input file
    :param out_file: output resolution file
    :param base: start position
    :param final: end position
    :param resolution: size of resolution
    :return:
    '''
    in_f = open(in_file, 'r')
    out_f = open(out_file, 'a')
    index = 0
    remain = 0
    cut = resolution
    while True:
        line = in_f.readline()
        if not line:
            break
        inf = line.strip('\n').split()
        start = int(inf[0])
        end = int(inf[1])

        # within the frame
        if end <= base:
            continue
        elif start >= final:
            break  # finish the resolution on this frame
        elif start < base and end > base:
            start = base
        elif start < final and end > final:
            end = final

        signal = float(inf[2])
        if cut >= end - start:
            remain += (end - start) * signal
            cut = (index + 1) * resolution + base - end
        else:
            while (index + 1) * resolution + base < end:  # careful here
                avg_signal = (remain + cut * signal) / resolution
                # print(index * resolution + base, base + (index + 1) * resolution, avg_signal, file=out_f)
                print(index, avg_signal, file=out_f)
                remain = 0
                cut = resolution
                index += 1
            remain = (end - index * resolution - base) * signal
            cut = (index + 1) * resolution + base - end
    # print(index * resolution + base, (index + 1) * resolution +base, remain / resolution, file=out_f)
    print(index, remain / resolution, file=out_f)
    out_f.close()
    in_f.close()

def create_resolution_plus(in_file, out_file, base, final, resolution=100, out_mod='w'):
    '''
    another version to create resolution
    allow to create resolution direct from raw data by specify chromosome
    :param in_file: input file
    :param out_file: output resolution file
    :param chr_resol: chromosome name
    :param base: start position
    :param final: end position
    :param resolution: size of resolution
    :return:
    '''
    in_f = open(in_file, 'r')
    out_f = open(out_file, out_mod)
    index = 0
    remain = 0
    cut = resolution
    while True:
        line = in_f.readline()
        if not line:
            break
        inf = line.strip('\n').split()
        start = int(inf[0])
        end = int(inf[1])
        signal = float(inf[2])
 
        same_chr = True
        # within the frame
        if end <= base:
            continue
        elif start >= final:
            break  # finish the resolution on this frame
        elif start < base and end > base:
            start = base
        elif start < final and end > final:
            end = final

        if cut >= end - start:
            remain += (end - start) * signal
            cut = (index + 1) * resolution + base - end
        else:
            while (index + 1) * resolution + base < end:  # careful here
                avg_signal = (remain + cut * signal) / resolution
                # print(index * resolution + base, base + (index + 1) * resolution, avg_signal, file=out_f)
                print(index, avg_signal, file=out_f)
                remain = 0
                cut = resolution
                index += 1
            remain = (end - index * resolution - base) * signal
            cut = (index + 1) * resolution + base - end
    # print(index * resolution + base, (index + 1) * resolution +base, remain / resolution, file=out_f)
    print(index, remain / resolution, file=out_f)
    out_f.close()
    in_f.close()


def output_matrix(matrix, filename):
    '''
    output a matrix data to a file
    :param filename:
    :return:
    '''
    with open(filename, 'w') as f:
        print(matrix_to_str(matrix), file=f)

def read_matrix(filename):
    '''
    read in a format
    :param filename:
    :return: numpy matrix
    '''
    in_f = open(filename, "r")
    inf = in_f.readline().strip('\n').split()
    m = int(inf[0])
    n = int(inf[1])
    matrix = np.asmatrix(np.zeros((m,n)))
    for i in range(m):
        inf = in_f.readline().strip('\n').split()
        for j in range(n):
            matrix[i,j] = float(inf[j])
    in_f.close()
    return matrix

def matrix_to_str(m):
    """
    :param m: numpy matrix
    :return s: nicely formatted matrix data in string
    """
    s = ""
    row,col = m.shape
    s += ("{} {}\n".format(row, col))
    for i in range(row):
        for j in range(col):
            s += "{} ".format(m[i,j])
        s += "\n"
    return s

def test_matrix_readinout():
    a = np.asmatrix([[1, 2, 3], [3, 3, 4]])
    output_matrix(a, "out.txt")
    print(readin_matrix("out.txt"))

def read_trained_ssm(in_file):
    '''
    read in theta param for ssm
    :param in_file: input theta file
    :return positive_flag, penalize_flag, factor_count, lamda1, assay_count, theta_m:parameter of ssm model
    update 2018-July-01: currently abandoned as the parameter changes in trainning
    '''
    in_f = open(in_file)
    positive_flag = bool(in_f.readline().strip('\n').split()[1])
    penalize_flag = bool(in_f.readline().strip('\n').split()[1])
    lamda1 = float(in_f.readline().strip('\n').split()[1])
    factor_count, assay_count = [int(element) for element in in_f.readline().strip('\n').split()]
    theta_m = np.asmatrix(np.zeros((factor_count, assay_count)))
    for k in range(factor_count):
        inf = in_f.readline().strip('\n').split()
        for e in range(assay_count):
            theta_m[k, e] = inf[e]
    return positive_flag, penalize_flag, factor_count, lamda1, assay_count, theta_m

def resolution_index_locate(pos, resolution_size, start=0):
    """
    :param pos: bp position 
    :return index: index after resolution applied
    """
    return math.floor((pos - start)*1.0 / resolution_size)

def adj_r2_score(real_data, predict_data, n, p):
    """
    :param n: number of samples
    :param p: number of states
    """
    r2 = sklearn.metrics.r2_score(real_data, predict_data)
    adj_r2 = 1 - (1-r2)*(n-1)/(n-p-1)
    return adj_r2


def ncc_score(a, b, l):
    a = np.asarray(a)[0]
    b = np.asarray(b)[0]
    # return np.correlate(a, b)
    return np.dot(a-np.mean(a), (b-np.mean(b))) / np.sqrt(np.var(a) * np.var(b)) / l


def merge_em(a_m, b_m, count=1):
    """
    merge emission matrix by take two column vector most correlated (ncc score)
    a_m, b_m: state_n by assay_n, both are emission matrix
    """
    state_n, assay_n = a_m.shape
    state_dic = {}
    merge_m = np.asmatrix(np.zeros((state_n, assay_n)))
    for i in range(state_n):
        max_score = 0
        max_j = 0
        a = a_m[i,:]
        for j in range(state_n):
            b = b_m[j,:]
            score = ncc_score(a / count, b, assay_n)
            if score > max_score:
                max_j = j
                max_score = score
        state_dic[i] = max_j
        merge_m[i, :] = a_m[i, :] + b_m[max_j, :]
    print("[Merge]a_m:\n", a_m)
    print("[Nerge]b_m:\n", b_m)
    print("[Merge]merged_m:\n", merge_m)
    print("[Merge]dict:\n", state_dic)
    return merge_m, state_dic

def merge_tr(a_m, b_m, state_dic):
    n_state, _ = a_m.shape
    merge_m =  np.asmatrix(np.zeros((n_state, n_state)))
    for i in range(n_state):
        for j in range(n_state):
            merge_m[i, j] = (a_m[i, j] + b_m[state_dic[i], state_dic[j]])
    return merge_m

def test_merge_em():
    a_m = np.asmatrix([[0.0,0.1,0.3], [0.1,0.0,0.2], [0.1, 0.2, 0.0]])
    b_m = np.asmatrix([[0.0,0.12,0.35], [0.1, 0.2, 0.0], [0.1,0.0,0.25]])
    new_m = merge_em(a_m, b_m)
    print(new_m)

if __name__ == '__main__':
    # test_rm_interset()
    test_merge_em()
