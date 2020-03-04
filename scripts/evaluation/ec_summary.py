#####################################################################################
# Compression Summary for Error Correction Benchmarking
#   zar-lab ucla
#   5/17/18
#   supervisor: Serghei Mangul
#   author: Keith Mitchell
"""
  Functions Contained
    check_existence: checks if the filename supplied has already been written to (should append: 'a') else we use 'w'.
    append_summary: appends the data from the read level and base level compression to a master or summary file.
"""
#####################################################################################

from __future__ import division
import argparse
import csv
import os

def check_existence(out_dir, file_name):
    if os.path.exists(out_dir + file_name):
        return 'a'  # append if already exists
    else:
        return 'w'  # make a new file if not


def append_summary(wrapper, kmer, read_data, base_data, ec_name, outdir, filename):
    read_sum = {"TN": 0, "TP": 0, "FN": 0, 'FN WRONG': 0, 'FP': 0, 'FP INDEL': 0, "Read Total": 0, "Not TN": 0, "Trimmed Total": 0}
    base_sum = {"Total Bases": 0, 'TP': 0, 'FN': 0, 'FN WRONG': 0, 'FP': 0, 'FP INDEL': 0, 'FP TRIM': 0, 'TP TRIM': 0, 'TN': 0}
    with open(read_data, 'r') as reads_csv:
        with open(base_data, 'r') as base_csv:
            reads_reader = csv.reader(reads_csv, delimiter=',')
            # readout.writerow([read_name, read_classification])
            for row in reads_reader:
                read_sum[str(row[1])] += 1
                read_sum["Not TN"] += 1

            base_reader = csv.reader(base_csv, delimiter=',')
            for row in base_reader:
                # NOTE: we take the "Read Total" here since we dont record TN reads but we record everything for
                # base evaluation so we know this is the number of reads that were analyzed
                read_sum["Read Total"] += 1
                if int(row[7]) != 0 or int(row[8]) != 0:
                    read_sum["Trimmed Total"] += 1

                # Aggregrate the base data per file
                base_sum["Total Bases"] += int(row[1])
                base_sum["TP"] += int(row[2])
                base_sum["FN"] += int(row[3])
                base_sum["FN WRONG"] += int(row[4])
                base_sum["FP"] += int(row[5])
                base_sum["FP INDEL"] += int(row[6])
                base_sum["FP TRIM"] += int(row[7])
                base_sum["TP TRIM"] += int(row[8])
                base_sum["TN"] += int(row[1]) - (int(row[2]) + int(row[3]) + int(row[4])
                                                 + int(row[5]) + int(row[6]) + int(row[7]) + int(row[8]))
            read_sum["TN"] = read_sum["Read Total"] - read_sum["Not TN"]


    type = check_existence(outdir, filename)
    with open(outdir + filename, type) as summaryout:
        summaryout = csv.writer(summaryout, delimiter=',')
        if type == 'w':
            write_header=["EC Filename", "Wrapper Name", "Kmer Size",
                          "Read - TP", "Read - TN", "Read - FN", "Read - FN WRONG", "Read - FP", "Read - FP INDEL",
                          "Total Trimmed Reads", "Total Reads",
                          "Base - TP", "Base - TN", "Base - FN", "Base - FN WRONG", "Base - FP", "Base - FP INDEL",
                          "Base - FP TRIM", "Base - TP TRIM", "Total Bases"]
            summaryout.writerow(write_header)

        ec_cleaned = ec_name.split("/")
        ec_cleaned = ec_cleaned[-1]
        wrapper_cleaned = wrapper.split("/")
        wrapper_cleaned = wrapper_cleaned[-1]

        write_data = [ec_cleaned, wrapper_cleaned, kmer,
                      read_sum["TP"], read_sum["TN"], read_sum["FN"], read_sum["FN WRONG"], read_sum["FP"],
                      read_sum["FP INDEL"], read_sum["Trimmed Total"], read_sum["Read Total"],
                      base_sum["TP"], base_sum["TN"], base_sum["FN"], base_sum["FN WRONG"], base_sum["FP"],
                      base_sum["FP INDEL"], base_sum["FP TRIM"], base_sum["TP TRIM"], base_sum["Total Bases"]]

        summaryout.writerow(write_data)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Produces EC Evaluation for either paired-end or single-end sequences.\
                                                 takes either 6, 5, or 4 arguments.( will adjust if needed)')
    #TODO: should we specify the filename
    parser.add_argument('-wrapper', help='wrapper that was ran', required=True)
    parser.add_argument('-kmer', help='kmer size used to run the EC tool', required=True)
    parser.add_argument('-read_data', help='read data file from ec_evaluation', required=True)
    parser.add_argument('-base_data', help='base data file from ec_evaluation', required=True)
    parser.add_argument('-ec_name', help='name of the error corrected file produced from the wrapper', required=True)
    parser.add_argument('-outdir', help='outdir for the filename to be placed', required=True)
    # parser.add_argument('-filename', help='output file name', required=True)

    args = vars(parser.parse_args())
    wrapper = args['wrapper']
    kmer = args['kmer']
    read_data = args['read_data']
    base_data = args['base_data']
    ec_name = args['ec_name']
    out_dir = args['outdir']
    # filename = args['outdir']

    filename='master_summary.txt'
    append_summary(wrapper, kmer, read_data, base_data, ec_name, out_dir, filename)
