#############################################################################################
# Data Storage for Error Correction Benchmarking
#   project zar-lab ucla
#   4/5/18

#  Functions Contained: store_base_data, store_read_data, baseline, data, compression
#############################################################################################


import os
import logging
import sys
import csv


def check_existence(base_dir_join, storage_type, file_name):
    if os.path.exists(base_dir_join + file_name + storage_type + '.txt'):
        return 'a'  # append if already exists
    else:
        return 'w'  # make a new file if not

def my_log(base_dir_join, file_name, message):
    storage_ext = '_log'
    type = check_existence(base_dir_join, storage_ext, file_name)
    with open(base_dir_join + file_name + storage_ext + '.txt', type) as logout:
        print (message)
        logout.write(str(message) + "\n")

def store_base_data(base_dir_join, file_name, read_name, length, base_counts):
    # stats_dict = {'TN': 0, 'TP': 0, 'FN': 0, 'FN WRONG': 0, 'FP': 0, 'FP INDEL': 0, 'FP TRIM': 0, 'TP TRIM': 0}
    storage_ext = '_base_data'
    type = check_existence(base_dir_join, storage_ext, file_name)
    with open(base_dir_join + file_name + storage_ext + '.txt', type) as baseout:
        baseout = csv.writer(baseout, delimiter=',')
        baseout.writerow([read_name, length, base_counts['TP'], base_counts['FN'], base_counts['FN WRONG'], base_counts['FP'], base_counts['FP INDEL'], base_counts['FP TRIM'], base_counts['TP TRIM']])


def store_read_data(base_dir_join, file_name, read_name, read_classification):
    storage_ext = '_read_data'
    type = check_existence(base_dir_join, storage_ext, file_name)
    with open(base_dir_join + file_name + storage_ext + '.txt', type) as readout:
        readout = csv.writer(readout, delimiter=',')
        if read_classification != "TN":
            readout.writerow([read_name, read_classification])


def baseline(base_dir_join, file_name, read_name, length, base_counts):
    # stats_dict = {'TN': 0, 'TP': 0, 'FN': 0, 'FN WRONG': 0, 'FP': 0, 'FP INDEL': 0, 'FP TRIM': 0, 'TP TRIM': 0}
    storage_ext = '_baseline_data'
    type = check_existence(base_dir_join, storage_ext, file_name)
    with open(base_dir_join + file_name + storage_ext + '.txt', type) as baseline_out:
        baseline_out = csv.writer(baseline_out, delimiter=',')
        baseline_out.writerow([read_name, length, base_counts['TP'], base_counts['FP'], base_counts['FP INDEL']])


def data_compression(base_dir_join, file_name, read_name, length, position_calls):
    storage_ext = '_data_compression'
    type = check_existence(base_dir_join, storage_ext, file_name)
    with open(base_dir_join + file_name + storage_ext + '.txt', type) as data_comp_out:
        data_comp = csv.writer(data_comp_out, delimiter=',')
        for key, value in sorted(position_calls.items()):
            position = key
            call, true_bp, raw_bp, ec_bp = value[0], value[1], value[2], value[3]
            data_comp.writerow([read_name, length, position, call, true_bp, raw_bp, ec_bp])
