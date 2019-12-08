"""
this script is to generate blacklist free resoluted data for annotation
"""
import os
import sys
is_exp_script = True
script_dir = os.path.split(os.path.realpath(__file__))[0]
if is_exp_script:
    sys.path.insert(0, os.path.join(script_dir, "../../../src/"))
else:
    sys.path.insert(0, os.path.join(script_dir, ""))

import util.utility as util
import util.myconfig as myconfig
import copy
if is_exp_script:
    root_path = os.path.join(script_dir, "../../../")
else:
    root_path = os.path.join(script_dir, "../")
configs = myconfig.config(root_path, log_dir=script_dir)

assay_list = configs.chromHMM_assay

chunk_region = {}
blacklist_index_region = {}
resolution_size = 200

def locate_blacklist_index_region():
    """
    locatet the index of blacklist region
    """
    global blacklist_index_region
    with open(configs.blacklist_region_file, 'r') as blacklist_region_f:
        while True:
            line = blacklist_region_f.readline()
            if not line:
                break
            inf = line.strip('\n').split()
            chromosome = inf[0]
            start = int(inf[1])
            end = int(inf[2])
            start_index = util.resolution_index_locate(start, resolution_size)
            end_index = util.resolution_index_locate(end, resolution_size)
            if chromosome not in blacklist_index_region:
                blacklist_index_region[chromosome] = [(start_index, end_index)]
            else:
                blacklist_index_region[chromosome].append((start_index, end_index))
        # sort the blacklist_inde_region
        for chromosome in blacklist_index_region:
            blacklist_index_region[chromosome] = list(sorted(blacklist_index_region[chromosome], key=lambda x:x[0]))

        # test
        for chromosome in blacklist_index_region:
            for region in blacklist_index_region[chromosome]:
                configs.toLog("{} {} {}".format(chromosome, region[0], region[1]))
def assay_chunk(assay):
    '''
    :param assay: prefix name of assay, like for E126-H3K36me3
    :return:
    '''
    chunked_dir = os.path.join(configs.blacklist_rm_data_path, assay, "chunk")
    if not os.path.exists(chunked_dir):
        os.makedirs(chunked_dir)
        print("Chunk assay " + assay + "...")
        chunked_fs = {}
        for chromosome in configs.chromosome_list:
            chunked_fs[chromosome] = open(os.path.join(chunked_dir, chromosome + ".bed"), 'w')
        with open(configs.assay_data_path + assay + ".begGraph", "r") as assay_f:
            last_line = None
            while True:
                line = assay_f.readline()
                if not line:
                    break
                inf = line.strip('\n').split()
                chromosome = inf[0]
                start = int(inf[1])
                end = int(inf[2])
                signal = float(inf[3])
                if chromosome not in chunked_fs:
                    continue
                print(start, end, signal, file=chunked_fs[chromosome])
    else:
        print("Chunk assay" + assay + "already exists")

    
def assay_resolution(assay):
    global resolution_size
    print("Building resolution on assay " + assay + "...")
    # build resolution
    resolution_dir = configs.blacklist_rm_data_path + assay + "/resolution-" + str(resolution_size) + "bp/"
    if not os.path.exists(resolution_dir):
        os.makedirs(resolution_dir)
        configs.toLog(assay)
        temp_blacklist_index_region = copy.deepcopy(blacklist_index_region)
        for chromosome in configs.chromosome_list:
            with open(os.path.join(resolution_dir, chromosome + ".bed"), 'w') as resol_f, open(os.path.join(configs.blacklist_rm_data_path, assay, "chunk", chromosome + ".bed"), 'r') as chunked_f:
                index = 0
                remain = 0
                cut = resolution_size
                while True:
                    line = chunked_f.readline()
                    if not line:
                        break
                    inf = line.strip('\n').split()
                    start = int(inf[0])
                    end = int(inf[1])
                    signal = float(inf[2])
            
                    same_chr = True
                    # within the frame
                    if cut >= end - start:
                        remain += (end - start) * signal
                        cut = (index + 1) * resolution_size - end
                    else:
                        while (index + 1) * resolution_size < end:  
                            avg_signal = (remain + cut * signal) / resolution_size
                            if chromosome in temp_blacklist_index_region:
                                while len(temp_blacklist_index_region[chromosome]) > 0 and index > temp_blacklist_index_region[chromosome][0][1]:
                                    temp_blacklist_index_region[chromosome] = temp_blacklist_index_region[chromosome][1:]
                                if len(temp_blacklist_index_region[chromosome]) == 0 or not (index >=  temp_blacklist_index_region[chromosome][0][0] and index <= temp_blacklist_index_region[chromosome][0][1]):
                                    print(index, avg_signal, file=resol_f)
                                else:
                                    configs.toLog("chromosome{} index{} b_region{}".format(chromosome, index, temp_blacklist_index_region[chromosome][0]))
                            else:
                                print(index, avg_signal, file=resol_f)
                            remain = 0
                            cut = resolution_size
                            index += 1
                        remain = (end - index * resolution_size) * signal
                        cut = (index + 1) * resolution_size - end
                if chromosome in temp_blacklist_index_region and len(temp_blacklist_index_region[chromosome]) > 0:
                    if not (index >=  temp_blacklist_index_region[chromosome][0][0] and index <= temp_blacklist_index_region[chromosome][0][1]):
                        print(index, remain / resolution_size, file=resol_f)
                    else:
                        configs.toLog("chromosome{} index{} b_region{}".format(chromosome, index, temp_blacklist_index_region[chromosome][0]))
                else:
                    print(index, remain / resolution_size, file=resol_f)
    else:
        print("resolution-" + str(resolution_size) + "bp already exists")

def process_data():
    # chunk assay 
    for assay in assay_list:
        assay_chunk(assay)
    # apply resolution
    for assay in assay_list:
        assay_resolution(assay)

def main():
    locate_blacklist_index_region()
    process_data()
    return 0
    
if __name__ == "__main__":
    configs.openLog()
    sys.exit(main())
    configs.closeLog()
