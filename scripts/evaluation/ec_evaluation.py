#####################################################################################
# Evaluation Code for Error Correction Benchmarking
#   zar-lab ucla
#   4/20/18
#   supervisor: Serghei Mangul
#   author: Keith Mitchell
"""
Functions Contained:
    msa - performs a multiple sequence alignment between the true, raw, and ec read
          returns: a dict of each individual sequence from the alignment ("TRUE", "RAW", "EC")
    analyze_bases - iterates through
    analyze_reads -
    handle_sequences -
"""
#####################################################################################

from Bio import pairwise2, SeqIO
from ec_data_compression import store_base_data, store_read_data, baseline, data_compression, my_log
import sys
import argparse
import os
import csv
import uuid
import subprocess
import time


def msa(true_sequence, raw_sequence, ec_sequence, description):
    #TODO
    # create a temporary fasta for the msa to be performed on
    id_temp = uuid.uuid4()
    # create a temporary fasta with id for the msa output
    id_msa = uuid.uuid4()

    # alter directory based on where temporary msas can go (this can use a lot of memory)
    with open('./%s.fasta' %id_temp, 'w') as fasta:
        fasta.write(">TRUE\n")
        fasta.write(str(true_sequence) + '\n')
        fasta.write(">RAW\n")
        fasta.write(str(raw_sequence) + '\n')
        fasta.write(">EC\n")
        fasta.write(str(ec_sequence) + '\n')

    # run the msa
    subprocess.call(["./muscle3.8.31_i86darwin64","-quiet",
                     "-in", "./%s.fasta" %id_temp,
                     "-out", "./%s.fasta" %id_msa])

    # remove the temporary file from prior to the msa
    os.remove("./%s.fasta" % id_temp)

    # parse the msa file to get each sequences we care about
    aligned = {"RAW": '',"TRUE": '', "EC": ''}
    fasta_msa = SeqIO.parse("./%s.fasta" % id_msa, 'fasta')
    for rec in fasta_msa:
        aligned[rec.description] = rec.seq

    # remove the temporary file from after the msa
    os.remove("./%s.fasta" % id_msa)


    #print (description)
    #print ("True:", aligned["TRUE"])
    #print ("Raw: ", aligned["RAW"])
    #print ("EC:  ", aligned["EC"])

    return aligned["TRUE"], aligned["RAW"], aligned["EC"]

def resolve_trim(true, raw, ec):
    fp_trim = 0
    tp_trim = 0
    position_calls_trim = {}
    position = 1

    # start from the left
    for true_bp, raw_bp, ec_bp in zip(true, raw, ec):
        if ec_bp == '-':
            if raw_bp == true_bp:
                position_calls_trim[position] = ['FP TRIM', true_bp, raw_bp, ec_bp]
                fp_trim += 1
            else:
                position_calls_trim[position] = ['TP TRIM', true_bp, raw_bp, ec_bp]
                tp_trim += 1
            position += 1
        else:
            break

    # start from the right
    position = len(true)
    for true_bp, raw_bp, ec_bp in zip(reversed(true), reversed(raw), reversed(ec)):
        if ec_bp == '-':
            if raw_bp == true_bp:
                position_calls_trim[position] = ['FP TRIM', true_bp, raw_bp, ec_bp]
                fp_trim += 1
            else:
                position_calls_trim[position] = ['TP TRIM', true_bp, raw_bp, ec_bp]
                tp_trim += 1
            position -= 1
        else:
            break

    return fp_trim, tp_trim, position_calls_trim


def analyze_bases(true_read, raw_read, ec_read):

    #TODO this function here
    # base level statistics now that there is an "MSA" that we can compare.
    stats_dict = {'TN': 0, 'TP': 0, 'FN': 0, 'FN WRONG': 0, 'FP': 0, 'FP INDEL': 0, 'FP TRIM': 0, 'TP TRIM': 0}
    position_calls = {}

    trim_results = resolve_trim(true_read, raw_read, ec_read)
    trim_positions = trim_results[2]
    tp_trim_count = trim_results[1]
    fp_trim_count = trim_results[0]
    stats_dict['FP TRIM'] += fp_trim_count
    stats_dict['TP TRIM'] += tp_trim_count

    pos_marker = 1
    ec_bp, true_bp, raw_bp = "", "", ""
    for true_bp, raw_bp, ec_bp in zip(true_read, raw_read, ec_read):
        if pos_marker in trim_positions.keys():
            pos_marker += 1
            continue

        elif true_bp == raw_bp == ec_bp:
            stats_dict['TN'] += 1

        # if the tool changed the BP when it did not need to be fixed
        elif true_bp == raw_bp != ec_bp:
            if raw_bp == '-' or true_bp == '-' or ec_bp == '-':
                stats_dict['FP INDEL'] += 1
                position_calls[pos_marker] = ['FP INDEL', true_bp, raw_bp, ec_bp]
            else:
                stats_dict['FP'] += 1
                position_calls[pos_marker] = ['FP', true_bp, raw_bp, ec_bp]

        # if the EC was correct
        elif ec_bp == true_bp != raw_bp:
            stats_dict['TP'] += 1
            position_calls[pos_marker] = ['TP', true_bp, raw_bp, ec_bp]

        # if none of the bases are equal or the tool did not make a fix when one should have been performed
        elif (true_bp != raw_bp == ec_bp) or (true_bp != raw_bp != ec_bp and true_bp != ec_bp != raw_bp):
            if true_bp != raw_bp == ec_bp:
                stats_dict['FN'] += 1
                position_calls[pos_marker] = ['FN', true_bp, raw_bp, ec_bp]
            else:
                stats_dict['FN WRONG'] += 1
                position_calls[pos_marker] = ['FN WRONG', true_bp, raw_bp, ec_bp]

        else:
            message = "ERROR Base Evaluation Case Was Not Found."
            #print (message)
            my_log(base_dir_join, cleaned_filename, message)
        pos_marker += 1

    length = len(true_read)
    if length != sum(stats_dict.values()):
        message = "ERROR in Base Level Evaluation Summary"
        #print (message)
        my_log(base_dir_join, cleaned_filename, message)
    position_calls.update(trim_positions)
    return stats_dict, length, position_calls


def analyze_read(stats_dict, length):
    # stats_dict = {'TN': 0, 'TP': 0, 'FN': 0, 'FN WRONG': 0, 'FP': 0, 'FP INDEL': 0, 'FP TRIM': 0, 'TP TRIM': 0}
    non_trim_summary = 0
    for key in stats_dict.keys():
        if key != "FP TRIM" and key != "TP TRIM":
            non_trim_summary += stats_dict[key]

    # If all bases are the TN then the read is a TN.
    if stats_dict['TN'] == non_trim_summary:
        #print ("TN")
        return "TN"

    # If there are normal FP but no other FP(INDEL/trim) then the read is a normal FP.
    elif stats_dict['FP'] != 0 and stats_dict['FP INDEL'] == 0:
        #print ("FP")
        return "FP"

    # Next if there are FP from INDEL then the read is as well.
    elif stats_dict['FP INDEL'] != 0:
        #print ("FP INDEL")
        return "FP INDEL"

    # If there are any TP bases and no FN bases then the read is a TP
    elif stats_dict['TP'] != 0 and stats_dict['FN'] == 0 and stats_dict['FN WRONG'] == 0:
        #print ("TP")
        return "TP"

    # Otherwise the read is a FN
    elif stats_dict['FN'] != 0 and stats_dict['FN WRONG'] == 0:
        #print ("FN")
        return "FN"

    else:
        #print ("FN WRONG")
        return "FN WRONG"


def find_match(fastq, true):
    try:
        sequence = fastq[true]
        return sequence
    except:
        return None


def handle_sequences(true_check, true_rec, two_raw, fastq_raw1_parser, fastq_raw2_parser, fastq_ec1_parser):

    if two_raw is True:
        raw_rec = find_match(fastq_raw1_parser, true_check[0])
        if raw_rec is None:
            raw_rec = find_match(fastq_raw2_parser, true_check[0])
    else:
        raw_rec = find_match(fastq_raw1_parser, true_check[0])

    if raw_rec is not None:
        ec_rec = find_match(fastq_ec1_parser, true_check[0])

        if ec_rec is not None:
            alignment = msa(true_rec, raw_rec, ec_rec, true_check[0])
            base_counts = analyze_bases(alignment[0], alignment[1], alignment[2])

            if base_counts is None:
                message = "FAILURE: Base count == 'None'(improper MSA) %s" % true_check[0]
                my_log(base_dir_join, cleaned_filename, message)

            else:
                position_calls = base_counts[2]
                length = base_counts[1]
                base_stats = base_counts[0]
                #print (base_stats)
                read_class = analyze_read(base_stats, length)


                store_base_data(base_dir_join, cleaned_filename, true_check[0], length, base_stats)
                store_read_data(base_dir_join, cleaned_filename, true_check[0], read_class)
                baseline(base_dir_join, cleaned_filename, true_check[0], length, base_stats)
                data_compression(base_dir_join, cleaned_filename, true_check[0], length, position_calls)
    else:
        message = 'It seems that no match was found for the sequence %s (Paired End) for the %s lookup' % (true_check[0], "RAW")
        my_log(base_dir_join, cleaned_filename, message)


def make_dict(parsed_fastq):
    dict = {}
    for rec in parsed_fastq:
        rec_id = rec.description.split()
        #print (rec_id[0])
        dict[rec_id[0]] = rec.seq
    return dict


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Produces EC Evaluation for either paired-end or single-end sequences.\
                                                 takes either 6, 5, or 4 arguments.( will adjust if needed)')

    parser.add_argument('-base_dir', help='This is the directory to produce the raw sequences to', required=True)
    parser.add_argument('-true_1', help='This is the true file to be evaluated.', required=True)
    parser.add_argument('-true_2', help='This is the second true file to be evaluated.', required=False)
    parser.add_argument('-raw_1', help='This is the raw file to be evaluated.', required=True)
    parser.add_argument('-raw_2', help='This is the second raw file, if paired-end, to be evaluated.', required=False)
    parser.add_argument('-ec_1', help='This is the ec file to be evaluated.', required=True)


    args = vars(parser.parse_args())
    base_dir = args['base_dir']
    true1_filename = args['true_1']
    true2_filename = args['true_2']
    raw1_filename = args['raw_1']
    raw2_filename = args['raw_2']
    ec1_filename = args['ec_1']

    end_file_directory = ec1_filename.split('/')
    global base_dir_join
    base_dir_join = os.path.join(str(base_dir))

    global cleaned_filename
    cleaned_filename = end_file_directory[-1].replace(".fastq", "")

    fastq_true1_parser = SeqIO.parse(os.path.join(str(true1_filename)), 'fastq')


    two_raw = False
    fastq_raw2_parser = ''
    if raw2_filename is not None and raw2_filename != '':
        two_raw = True

    two_true = False
    if true2_filename is not None and true2_filename != '':
        fastq_true2_parser = SeqIO.parse(os.path.join(str(true2_filename)), 'fastq')
        two_true = True

    if two_true == True:
        fastq_ec1_parser = SeqIO.parse(os.path.join(str(ec1_filename)), 'fastq')
        fastq_raw1_parser = SeqIO.parse(os.path.join(str(raw1_filename)), 'fastq')
        fastq_raw2_parser = SeqIO.parse(os.path.join(str(raw2_filename)), 'fastq')

        #print ("EC")
        fastq_ec1 = make_dict(fastq_ec1_parser)
        #print ("RAW1")
        fastq_raw1 = make_dict(fastq_raw1_parser)
        #print ("RAW2")
        fastq_raw2 = make_dict(fastq_raw2_parser)

        for true1_rec, true2_rec in zip(fastq_true1_parser, fastq_true2_parser):

            true_check1 = true1_rec.description.split()
            true_check2 = true2_rec.description.split()

            if true1_rec:
                handle_sequences(true_check1, true1_rec.seq, two_raw, fastq_raw1, fastq_raw2, fastq_ec1)

            if true2_rec:
                handle_sequences(true_check2, true2_rec.seq, two_raw, fastq_raw1, fastq_raw2, fastq_ec1)
    else:

            fastq_ec1_parser = SeqIO.parse(os.path.join(str(ec1_filename)), 'fastq')
            fastq_raw1_parser = SeqIO.parse(os.path.join(str(raw1_filename)), 'fastq')
            fastq_raw2_parser = SeqIO.parse(os.path.join(str(raw2_filename)), 'fastq')

            #Lines added to transform generator pobject into dictionaries
            fastq_ec1 = make_dict(fastq_ec1_parser)
            fastq_raw1 = make_dict(fastq_raw1_parser)

            for true_rec in fastq_true1_parser:
                true_check = true_rec.description.split(' ')


                #Changed the names of the variables fastq_raw1_parser to fastq_raw1 and fastq_ec1_parser to fastq_ec1
                handle_sequences(true_check, true_rec.seq, two_raw, fastq_raw1, fastq_raw2_parser, fastq_ec1)


    message = "DONE: %s" % (cleaned_filename)
    my_log(base_dir_join, cleaned_filename, message)
