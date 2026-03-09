#################################################
#  File Name:get_bin_segments.py
#  Author: Pengwei.Xing
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se,xpw1992@gmail.com
#  Created Time: Fri 23 Feb 2024 01:26:33 PM CET
#################################################


from pybedtools import BedTool
import numpy as np
import pandas as pd
import sys
import copy
import re

windowfile = sys.argv[1]

f = open('list','r')
f1 = open(windowfile,'r')##C2_peaks.sort.bed
sample_file = [line.strip() for line in f.readlines()]
sample_name = [line.strip().split('MG')[0] for line in sample_file]
all_bins = [line.strip() for line in f1.readlines()]
all_bins_name = []
for each in all_bins:
    each2 = each.split('\t')
    each3 = '_'.join(each2)
    all_bins_name.append(each3)
allpeaks = BedTool(windowfile)## merged bed file
matrix = pd.DataFrame(index=all_bins_name, columns=sample_name)
matrix = matrix.astype(str)

# Create an empty dictionary to store matrix data
matrix_dict = {}
for each_sample_name, file in zip(sample_name, sample_file):
    print(each_sample_name)
    U3068 = BedTool(file)
    intersections = allpeaks.intersect(U3068, wa=True, wb=True)
    print('finished %s' % each_sample_name)
    
    for a in intersections:
        CHR,ss,ee = a[0:3]
        ss2,ee2,states = a[4:7]
        LEN = int(ee2) - int(ss2)
        INDEX = CHR + '_' + ss + '_' + ee

        if INDEX not in matrix_dict:
            matrix_dict[INDEX] = {each_sample_name: {states: LEN}}
        else:
            if each_sample_name not in matrix_dict[INDEX]:
                matrix_dict[INDEX][each_sample_name] = {states: LEN}
            else:
                if states not in matrix_dict[INDEX][each_sample_name]:
                    matrix_dict[INDEX][each_sample_name][states] =LEN
                else:
                    matrix_dict[INDEX][each_sample_name][states] += LEN
    print('finished %s' % each_sample_name)

#keys_to_delete = {'E6', 'E8','E12','E13'}

# Iterate over the data dictionary
#for key, inner_dict in matrix_dict.items():
    # Iterate over the inner dictionaries
#    for inner_key, inner_inner_dict in inner_dict.items():
#        # Delete keys_to_delete from inner_inner_dict
#        for k in keys_to_delete:
#            inner_inner_dict.pop(k, None)
            
for key, inner_dict in matrix_dict.items():
    # Get keys with empty values
    empty_keys = [inner_key for inner_key, inner_inner_value in inner_dict.items() if not inner_inner_value or (isinstance(inner_inner_value, set) and not inner_inner_value)]
    # Remove keys with empty values
    for empty_key in empty_keys:
        del inner_dict[empty_key]


result = {}
for key, inner_dict in matrix_dict.items():
    result[key] = {}
    for inner_key, inner_inner_dict in inner_dict.items():
        max_key = max(inner_inner_dict, key=inner_inner_dict.get)
        result[key][inner_key] = max_key

mapping = {
    r"\bE5\b": "ActiveTSS",
    r"\bE1\b": "TssFlank",
    r"\bE2\b": "TssFlankD",
    r"\bE3\b": "TxFlank",
    r"\bE4\b": "EnhA1",
    r"\bE10\b": "EnhA2",
    r"\bE7\b": "EnhPoised1",
    r"\bE11\b": "EnhPrimed",
    r"\bE9\b": "EnhPoised2",
    r"\bE15\b": "EnhBivalent",
    r"\bE14\b": "ReprPC",
    r"\bE6\b": "Low1",
    r"\bE12\b": "Low2",
    r"\bE13\b": "Low3",
    r"\bE8\b": "Low4"
}

for index in result:
    for each_sample_name in result[index]:
        states = result[index][each_sample_name]
        for old_value, new_value in mapping.items():
            states = re.sub(old_value, new_value, states)
        result[index][each_sample_name] = states
# Convert the dictionary to DataFrame
result = pd.DataFrame(result).T.fillna('')

# Remove rows with all empty values
result2 = result.dropna(how='all')

result2.to_csv(sys.argv[2],sep="\t")
