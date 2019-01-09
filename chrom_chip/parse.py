#!/usr/bin/env python3

import os
import pickle
from datetime import datetime

import yaml
import numpy as np

from chrom_chip import Experiment
from chrom_chip.common import val_folder, output, read_pd, glob_check, make_folder, version


def parse_config(config_file):
    '''
    Parse experimental info from yaml file
    '''

    with open(config_file, 'r') as file:
        yml = yaml.safe_load(file)

    # Make a new experimental object
    exp = Experiment()

    # Project
    exp.project = yml['LSF_Project']

    # Setting Scratch folder
    exp.scratch = f'{os.getcwd()}/{yml["Name"]}_tmp/' if yml["Scratch_folder"] is None else f'{val_folder(yml["Scratch_folder"])}{yml["Name"]}/'
    os.makedirs(exp.scratch, exist_ok=True)

    # check whether experiment has been attempted
    exp.name = yml['Name']
    exp.out_dir = make_folder(yml['Output_directory'])
    filename = f'{exp.scratch}{exp.name}_incomplete.pkl'

    if os.path.isfile(filename):
        if yml['Restart'] is False:
            with open(filename, 'rb') as experiment:
                exp = pickle.load(experiment)
            os.remove(filename)

            # set new date
            exp.date = f'{datetime.now():%Y-%m-%d}'

            # For output of R logs into job_log_folder
            os.chdir(exp.job_folder)

            output(f'\n#############\nRestarting pipeline on {datetime.now():%Y-%m-%d %H:%M:%S}, from last completed step.', exp.log_file)

            return exp
        else:
            os.remove(filename)

    # Passing paramters to new object
    exp.date = f'{datetime.now():%Y-%m-%d}'

    # Log file
    exp.log_file = f'{exp.out_dir}{exp.name}-{exp.date}.log'

    output(f'Pipeline version {version()} run on {exp.date} \n', exp.log_file)
    output(f'Beginning ChIPseq Analysis: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)
    output('Reading experimental file...\n', exp.log_file)
    output(f"Pipeline output folder: {exp.out_dir}\n", exp.log_file)

    # Setting Job Folder
    exp.job_folder = f'{val_folder(exp.scratch)}logs/'
    os.makedirs(exp.job_folder, exist_ok=True)

    # Load sample info
    exp.sample_df = read_pd(yml['Sample_file'])

    # Make Sample Name
    exp.sample_df.replace([np.nan], 'none', inplace=True)
    exp.sample_df['Sample_Name'] = exp.sample_df.Condition + '_' + exp.sample_df.Replicate
    output(f'Processing samples:\n{exp.sample_df}')

    # Paired
    exp.sample_df['paired'] = [x != 'none' for x in exp.sample_df.File2.tolist()]

    exp.IPs = exp.sample_df[exp.sample_df['Background Sample'] != 'none']
    sample_dict = exp.sample_df.Sample_Name.to_dict()
    exp.IPs['Background_Name'] = exp.IPs['Background Sample'].map(sample_dict)
    exp.samples = exp.IPs.Sample_Name.tolist()

    # Make out directory if it doesn't exist
    exp.out_dir = f'{val_folder(yml["Output_directory"])}{exp.name}/'
    os.makedirs(exp.out_dir, exist_ok=True)

    # Lab specific files
    exp.genome_indicies['spike_index'] = yml['Spike_index']

    # Locating genome indicies
    tsvs = yml['Genome_tsv'].split(',')
    genomes = ['hg38', 'hg19', 'mm10']
    for tsv in tsvs:
        glob_check(tsv)
        exp.genome_indicies['encode_tsv'] = {**exp.genome_indicies['encode_tsv'],
                                             **{genome: tsv for genome in genomes if genome in tsv}
                                             }

    exp.encode3_folder = val_folder(yml['ENCODE3_folder'])

    # Initialized Process Complete List
    exp._parsed = True

    output(f'Experiment file parsed: {datetime.now():%Y-%m-%d %H:%M:%S}\n', exp.log_file)

    return exp
