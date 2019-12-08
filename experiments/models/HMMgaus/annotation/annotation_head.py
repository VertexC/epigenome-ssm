"""
This the script to serve as head to distribute annotation task to node script.
"""
import os
import sys
is_exp_script = True
script_dir = os.path.split(os.path.realpath(__file__))[0]

import util.utility as util
import util.myconfig as myconfig
if is_exp_script:
    root_path = os.path.join(script_dir, "../../../")
else:
    root_path = os.path.join(script_dir, "../")
configs = myconfig.config(root_path, log_dir = script_dir)
configs.openLog()

import subprocess
import re
import time
def main():
    train_dir = os.path.join(configs.experiment_path, "57-1020-18-pilotTrainHMMgaus/exp1")
    with open(os.path.join(train_dir, "log.txt"), 'r') as param_f:
        E = int(param_f.readline().strip('\n').split()[-1])
        K = int(param_f.readline().strip('\n').split()[-1])
        assay_list = param_f.readline().strip('\n').split()[1:]
        # only need is resolution size here
        resolution_size = int(param_f.readline().strip('\n').split()[-1])

    index_region = {}
    with open(configs.assay_region_file, 'r') as assay_region_f:
        while True:
            line = assay_region_f.readline()
            if not line:
                break
            inf = line.strip('\n').split()
            chromosome = inf[0]
            start = int(inf[1])
            end = int(inf[2])
            start_index = util.resolution_index_locate(start, resolution_size)
            end_index = util.resolution_index_locate(end, resolution_size)
            if chromosome in index_region:
                raise(Exception("Multi chromosome region"))
            index_region[chromosome] = [(start_index, end_index)]

        # # for test
        # for chromosome in index_region:
        #     print(chromosome, index_region[chromosome])

    # remove blacklist_index_region
    blacklist_index_region = {}
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

    # apply exclude intersection remove, assay_index_region - blacklist_index_region
    for chromosome in blacklist_index_region:
        for region_blacklist in blacklist_index_region[chromosome]:
            # configs.toLog("\nblack{}".format(region_blacklist))
            # configs.toLog(index_region[chromosome])
            # like a queue
            for temp in index_region[chromosome]:
                region_assay = index_region[chromosome][0]
                index_region[chromosome].remove(region_assay)
                # configs.toLog("after remove{}".format(index_region[chromosome]))
                index_region[chromosome] += util.rm_intersect(region_assay, region_blacklist, exclude=True)
                # configs.toLog("after pending{}".format(index_region[chromosome]))
                # configs.toLog("assay{}".format(region_assay))
                # configs.toLog("result{}".format(util.rm_intersect(region_assay, region_blacklist, exclude=True)))
    
    for chromosome in index_region:
        index_region[chromosome] = list(sorted(index_region[chromosome], key=lambda x:x[0]))
        # tested
        for region in index_region[chromosome]:
            configs.toLog("{} {}".format(chromosome, region))
    
   # assign jobs to subprocess
    child_list = []
    child_id = 0
    for chromosome in index_region:
        for region in index_region[chromosome]:
            start_index = region[0]
            end_index = region[1]
            success = False
            while not success:
                try:
                    child_process = subprocess.Popen(["sbatch", "node.sh", str(child_id), train_dir, chromosome, str(start_index), str(end_index)], stdout=subprocess.PIPE)
                    # # for test
                    # end_index = int(end_index/1000)
                    # subprocess.call(["./node.sh", str(child_id), train_dir, chromosome, str(start_index), str(end_index)])
                    # exit(0)
                except Exception:
                    time.sleep(60)
                    continue # rebatch after some sleep
                child_list.append(child_process)
                success = True
            # sleep after each job submission, to avoid scoket time out error
            time.sleep(40)
            child_id += 1
            
    job_ids = []
    for child in child_list:
        inf = child.stdout.readline().strip(b'\n').split()
        job_id = inf[-1]
        job_ids.append(job_id)

    while True:
        all_job_done = True
        # get job states 
        time.sleep(20)
        try:
            result = subprocess.check_output(["squeue", "-u", "vertexc"]).split(b'\n')  
        except Exception:
            time.sleep(60)
            continue
        header = result[0].split()
        job_index = 0
        st_index = 0
        for i in range(len(header)):
            if header[i] == b'JOBID':
                job_index = i
            if header[i] == b'ST':
                st_index = i
        for inf in result[1:]:
            if not inf:
                break
            if inf.split()[job_index] in job_ids:
                all_job_done = False
                break
        if all_job_done:
            break
    # assemble annotation file
    dis_annotation_result_dir = os.path.join(script_dir, "dis-annotation-result")
    con_annotation_result_dir = os.path.join(script_dir, "con-annotation-result")
    
    def collect(annotation_result_dir):
        sub_annotation = {}
        files = []
        for (_, _, filenames) in os.walk(annotation_result_dir):
            files.extend(filenames)
            break
        for myfile in files:
            file_id = int(re.search(r'_\d+', myfile).group()[1:])
            chromosome = re.search(r'.*_', myfile).group()[:-1]
            if chromosome not in sub_annotation:
                sub_annotation[chromosome] = [(myfile, file_id)]
            else:
                sub_annotation[chromosome].append((myfile, file_id))
                
        for chromosome in sub_annotation:
            sub_annotation[chromosome] = list(sorted(sub_annotation[chromosome], key=lambda x:x[1]))
            final_annotation_file = os.path.join(annotation_result_dir, chromosome + ".bed")
            open(final_annotation_file, 'w')
            for sub_annotation_file,_ in sub_annotation[chromosome]:
                result = subprocess.call(" ".join(["cat", os.path.join(annotation_result_dir, sub_annotation_file), ">>", final_annotation_file]), shell=True)
                os.remove(os.path.join(annotation_result_dir, sub_annotation_file))
    collect(dis_annotation_result_dir)
    collect(con_annotation_result_dir)
    return 0

if __name__ == "__main__":
    sys.exit(main())
