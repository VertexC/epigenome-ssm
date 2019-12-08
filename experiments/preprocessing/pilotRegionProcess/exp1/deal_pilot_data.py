import os
import sys
is_exp_script = True
script_dir = os.path.split(os.path.realpath(__file__))[0]
import util.myconfig as myconfig
import util.utility as util

if is_exp_script:
    root_path = os.path.join(script_dir, "../../../")
else:
    root_path = os.path.join(script_dir, "../")
configs = myconfig.config(root_path, log_dir=script_dir, data_type="p")

import numpy 
assay_list = configs.chromHMM_assay
resolution = 200

pilot_inf = {}  # {#chromosome: (start_location, end_location)}
pilot_index_region = {}  # {#chromosome: (start_location, end_location)}
pilot_index_region_file = configs.pilot_index_region(resolution)
first_resolution = True


def read_in_pilot():
    '''
    open the pilot file
    record the pilot location
    :return: 
    '''
    global pilot_inf
    print("read in pilot inf...")
    pilot_file = open(configs.pilot_region_file, "r")
    while True:
        line = pilot_file.readline()
        if not line:
            break
        inf = line.strip('\n').split()
        chr_name = inf[0]
        start = int(inf[1])
        end = int(inf[2])
        if inf[0] not in pilot_inf:
            pilot_inf[chr_name] = [(start, end)]
        else:
            pilot_inf[chr_name].append((start, end))

    pilot_file.close()

def assay_chunk(assay):
    print("Chunk assay " + assay + "...")
    assay_f = open(configs.assay_data_path + assay + ".begGraph", "r")
    last_line = assay_f.readline()
    new_data = False
    while True:
        if not last_line:
            break
        inf = last_line.strip('\n').split()
        # look up in pilot_inf
        chr_name = inf[0]
        if chr_name not in pilot_inf:
            last_line = assay_f.readline()
            continue
        else:
            # print("begin")
            current_chr_name = inf[0]
            assay_name = assay
            if not os.path.exists(configs.processed_pilot_data_path + assay_name + "/chunk/"):
                os.makedirs(configs.processed_pilot_data_path + assay_name + "/chunk/")
                new_data = True
            if not new_data:
                print("chunk of " + assay_name + " already exists")
                break
            pilot_f = open(configs.processed_pilot_data_path + assay_name + "/chunk/" + current_chr_name + ".begGraph", "w")
           
            frame_index = 0
            while True:
                inf = last_line.strip('\n').split()
                chr_name = inf[0]
                start = int(inf[1])
                end = int(inf[2])
                signal = float(inf[3])
                if chr_name != current_chr_name:
                    break
                if frame_index > len(pilot_inf[current_chr_name]) - 1:
                    # read until next chromosome
                    while True:
                        last_line = assay_f.readline()
                        if not last_line:
                            break
                        inf = last_line.strip('\n').split()
                        chr_name = inf[0]
                        if chr_name != current_chr_name:
                            break
                    break
                if not last_line:
                    break

                base = pilot_inf[current_chr_name][frame_index][0]
                final = pilot_inf[current_chr_name][frame_index][1]

                if end <= base:
                    last_line = assay_f.readline()
                    continue
                elif start >= final:
                    last_line = assay_f.readline()
                    frame_index += 1
                    continue
                elif start <= base < final <= end:
                    # no need to read next line
                    start = base
                    end = final
                    frame_index += 1
                    print(str(start), str(end), str(signal), file=pilot_f)
                    continue
                elif start <= base < end <= final:
                    start = base
                    print(str(start), str(end), str(signal), file=pilot_f)
                    last_line = assay_f.readline()
                elif base <= start < final <= end:
                    end = final
                    print(str(start), str(end), str(signal), file=pilot_f)
                    frame_index += 1
                elif base <= start < end <= final:
                    print(str(start), str(end), str(signal), file=pilot_f)
                    last_line = assay_f.readline()
                else:
                    # for test
                    print(start, end)
                    print(base, final)
                    print("Warining! Chunk unexpected case!")
                    exit(0)
            pilot_f.close()
    assay_f.close()

def assay_resolution(assay):
    global first_resolution
    print("Building resolution on assay " + assay + "...")
    # build resolution
    if not os.path.exists(configs.processed_pilot_data_path + assay + "/resolution-" + str(resolution) + "bp/"):
        os.makedirs(configs.processed_pilot_data_path + assay + "/resolution-" + str(resolution) + "bp/")
        for chromosome in pilot_inf:
            for region in pilot_inf[chromosome]:
                base = region[0]
                final = region[1]
                # extract from util to store the index
                in_f = open(configs.processed_pilot_data_path + assay + "/chunk/" + chromosome + ".begGraph", 'r')
                out_f = open(configs.processed_pilot_data_path + assay + "/resolution-" + str(resolution) + "bp/" + chromosome + ".begGraph", 'a')
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
                            print(index, avg_signal, file=out_f)
                            remain = 0
                            cut = resolution
                            index += 1
                        remain = (end - index * resolution - base) * signal
                        cut = (index + 1) * resolution + base - end
                print(index, remain / resolution, file=out_f)
                out_f.close()
                in_f.close()
                # to record region index
                if first_resolution:
                    if chromosome not in pilot_index_region:
                        pilot_index_region[chromosome] = [(0, index)]
                    else:
                        pilot_index_region[chromosome].append((0, index))
        if first_resolution:
            with open(pilot_index_region_file, 'w') as pilot_index_region_f:
                for chromosome in pilot_index_region:
                    for region in pilot_index_region[chromosome]:
                        print("{} {} {}".format(chromosome, region[0], region[1]), file=pilot_index_region_f)
            first_resolution = False
    else:
        print("resolution-" + str(resolution) + "bp already exists")

def process_data():
    # process pilot data, chunk and resolution
    for assay in assay_list:
        assay_chunk(assay)

    # build resolution
    for assay in assay_list:
        assay_resolution(assay)

def main():
    # read in pilot information
    read_in_pilot()
    process_data()
    return 0

if __name__ == "__main__":
    sys.exit(main())

    
