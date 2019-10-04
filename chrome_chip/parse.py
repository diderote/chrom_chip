#!/usr/bin/env python3

import os
import pickle
from datetime import datetime

import yaml
import numpy as np

from chrome_chip import Experiment
from chrome_chip.common import val_folder, output, read_pd, glob_check, make_folder, version


def parse_config(config_file, run_main=False):
    '''
    Parse experimental info from yaml file
    '''

    with open(config_file, 'r') as file:
        yml = yaml.safe_load(file)

    # Make a new experimental object
    exp = Experiment()

    # Project
    exp.project = yml['LSF_Project']

    # Check if running as pipeline
    exp.run_main = run_main

    # Setting Scratch folder
    exp.scratch = f'{os.getcwd()}/{yml["Name"]}_tmp/' if yml["Scratch_folder"] is None else f'{val_folder(yml["Scratch_folder"])}{yml["Name"]}/'
    os.makedirs(exp.scratch, exist_ok=True)

    # check whether experiment has been attempted
    exp.name = yml['Name']
    exp.out_dir = make_folder(f"{val_folder(yml['Output_directory'])}{exp.name}/")
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

            output(f'\n#############\nRestarting pipeline on {datetime.now():%Y-%m-%d %H:%M:%S}, from last completed step.', log_file=exp.log_file, run_main=exp.run_main)

            return exp
        else:
            os.remove(filename)

    # Passing paramters to new object
    exp.date = f'{datetime.now():%Y-%m-%d}'

    # Log file
    exp.log_file = f'{exp.out_dir}{exp.name}-{exp.date}.log'

    output(f'Pipeline version {version()} run on {exp.date} \n', log_file=exp.log_file, run_main=run_main)
    output(f'Beginning ChIPseq Analysis: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log_file=exp.log_file, run_main=run_main)
    output('Reading experimental file...\n', log_file=exp.log_file, run_main=run_main)
    output(f"Pipeline output folder: {exp.out_dir}\n", log_file=exp.log_file, run_main=run_main)

    # Setting Job Folder
    exp.job_folder = f'{val_folder(exp.scratch)}logs/'
    os.makedirs(exp.job_folder, exist_ok=True)

    # Load sample info
    exp.sample_df = read_pd(yml['Sample_file'])

    # Make Sample Name
    exp.sample_df.replace([np.nan], 'none', inplace=True)
    exp.sample_df['Sample_Name'] = exp.sample_df.Condition + '_' + exp.sample_df.Replicate
    output(f'Processing samples:\n{exp.sample_df}', log_file=exp.log_file, run_main=run_main)

    # Paired
    exp.sample_df['paired'] = [x != 'none' for x in exp.sample_df.File2.tolist()]

    exp.IPs = exp.sample_df[exp.sample_df['Background Sample'] != 'none'].copy()
    sample_dict = exp.sample_df.Sample_Name.to_dict()
    exp.IPs['Background_Name'] = exp.IPs['Background Sample'].map(sample_dict)
    exp.samples = exp.IPs.Sample_Name.tolist()

    # Convert Comparisons to a column of lists, then make unique comparisons
    exp.IPs['Comparisons'] = exp.IPs.Comparisons.apply(lambda x: [x.replace(' ', '') for x in x.split(',')])
    exp.IPs['Comparison_names'] = exp.IPs[['Condition', 'Comparisons']].apply(lambda x: ['_v_'.join(sorted([x[0], y])) for y in x[1] if x[1][0] != 'none'], axis=1)

    comparisons = []
    for comparison in exp.IPs.Comparison_names.tolist():
        comparisons += comparison
    exp.overlaps = {o_name: o_name.split('_v_') for o_name in set(comparisons)}

    exp.IPs['UMI'] = [x.lower() for x in exp.IPs.UMI.tolist()]

    # Spike-in comparisons
    exp.IPs['Spike-in Comparisons'] = exp.IPs['Spike-in Comparisons'].apply(lambda x: [x.replace(' ', '') for x in x.split(',')])
    exp.IPs['Spike_names'] = exp.IPs[['Condition', 'Spike-in Comparisons']].apply(lambda x: ['_v_'.join(sorted([x[0], y])) for y in x[1] if x[1][0] != 'none'], axis=1)

    sp_comparisons = [comparison for subls in exp.IPs.Spike_names.tolist() for comparison in subls]
    exp.spike_comparisons = {s_name: s_name.split('_v_') for s_name in set(sp_comparisons)}

    spike_samples = [condition for subls in exp.spike_comparisons.values() for condition in subls]
    exp.spike_samples = exp.IPs[exp.IPs.Condition.isin(spike_samples)].Sample_Name.tolist()

    # Make out directory if it doesn't exist
    exp.out_dir = make_folder(f'{val_folder(yml["Output_directory"])}{exp.name}/')

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

    output(f'Experiment file parsed: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log_file=exp.log_file, run_main=run_main)

    return exp
