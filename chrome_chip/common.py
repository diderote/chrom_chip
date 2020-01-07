#!/usr/bin/env python3

'''
University of Miami - Pegasus LSF Cluster ChIPseq Pipeline

Reads an experimental design yaml file (Version 0.7).
Requires a conda environment 'ChIPseq' made from environment.yml
www.github.com/diderote/LSF-ChIPseq/

To do:
    - remove unwanted bams
    - add diff binding
    - combine overlap two and three
    - make an install script\
    - change away from os.system\
    - add spike comparisons
    - change to cromwell
    - add custom options for cromwell in parser

'''
import os
import re
import glob
import pickle
import random
import time
import shutil
from datetime import datetime

__author__ = 'Daniel L. Karl'
__license__ = 'MIT'
__version__ = '0.1'


def version():
    return __version__


def html_header(version=__version__, author=__author__, license=__license__):
    return ''.join(['<h1>ChIPseq Analysis Notebook</h1>',
                    f'<body><b>Experiment Date: {datetime.now():%Y-%m-%d}<br>',
                    f'Pipeline version: {version}</b><br>',
                    '<a href="http://www.github.com/diderote/chrome_chip">Pipeline Code</a><br>',
                    f'License: {license} <br> Author: {author}'
                    ])


def val_folder(folder):
    folder = folder if folder.endswith('/') else f'{folder}/'
    return f'{os.getcwd()}/' if folder == '/' else folder


def rout_write(rout):
    '''
    function for setting r_out to print to file
    '''
    print(rout, file=open(f'{os.getcwd()}/R_out_{datetime.now():%Y-%m-%d}.txt', 'a'))


def is_fastq(file):
    ends = ['.fastq.gz', '.fasta.gz', 'fq.gz', 'fa.gz', 'fastq', 'fasta', 'fa', 'fq']
    for end in ends:
        if end in file:
            return True
    return False


def read_pd(file):
    import pandas as pd

    glob_check(file)
    if (file.endswith('txt')) or (file.endswith('tab')):
        return pd.read_table(file, header=0, index_col=0)
    elif (file.endswith('xls')) or (file.endswith('xlsx')):
        return pd.read_excel(file, index_col=0)
    elif file.endswith('eak.gz'):
        return pd.read_table(file, compression='gzip', header=None, index_col=None)
    elif file.endswith('bed'):
        return pd.read_table(file, header=None, index_col=None).iloc[:, :3]
    else:
        raise IOError("Cannot parse file.  Make sure it is .txt, .xls, .xlsx, .bed, or *Peak.gz")


def output(text, log_file=None, run_main=False):
    if run_main:
        print(text, file=open(log_file, 'a'))
    else:
        print(text)


def glob_remove(path):
    files = glob.glob(path)
    for file in files:
        if os.path.isfile(file):
            os.remove(file)


def make_folder(folder):
    folder = val_folder(folder)
    os.makedirs(folder, exist_ok=True)
    return folder


def image_display(file):
    from IPython.display import display, Image

    display(Image(file))


def load_bedtool(file):
    import pandas as pd
    from pybedtools import BedTool

    if file.endswith('.gz'):
        return BedTool.from_dataframe(pd.read_csv(file, compression='gzip', index_col=None, header=None, sep="\t").iloc[:, :3])
    else:
        return BedTool(file)


def glob_check(path):
    file = glob.glob(path)
    return 'none' if len(file) == 0 else file[0]


def bed2df(bed):
    if len(bed) == 0:
        return None
    else:
        return bed.to_dataframe()


def txt_replace(string):
    return string.replace('.txt.bz2', '.fastq.gz')


def out_result(image, text, run_main):
    from IPython.display import HTML, display

    if not run_main:
        if os.path.isfile(image):
            display(HTML(f'<h2>{text}</h2>'))
            image_display(image)
        else:
            display(HTML(f'<body>No result for {text} found.</body>'))


def close_out(task, exp):
        output(f'Error in {task}.', log_file=exp.log_file, run_main=exp.run_main)
        filename = f'{exp.scratch}{exp.name}_incomplete.pkl'
        with open(filename, 'wb') as experiment:
            pickle.dump(exp, experiment)
        raise RuntimeError(f'Error in {task}. Fix problem then resubmit with same command to continue from last completed step.')


def send_job(command_list, job_name, job_log_folder, q, mem, log_file, project, cores=1, submit=False, run_main=False):
    '''
    Sends job to LSF pegasus.ccs.miami.edu
    '''

    os.makedirs(job_log_folder, exist_ok=True)

    rand_id = str(random.randint(0, 100000))

    cmd = '''\n'''.join(['#!/bin/bash',
                         f'#BSUB -J ID_{rand_id}_JOB_{job_name}',
                         f'#BSUB -R "rusage[mem={mem}]"',
                         f'#BSUB -R "span[ptile={cores}]"',
                         f'#BSUB -o {job_log_folder}{job_name}_logs_{rand_id}.stdout.%J',
                         f'#BSUB -e {job_log_folder}{job_name}_logs_{rand_id}.stderr.%J',
                         '#BSUB -W 120:00',
                         f'#BSUB -n {cores}',
                         f'#BSUB -q {q}',
                         f'#BSUB -P {project}',
                         '',
                         '\n'.join(command_list)
                         ])

    job_path_name = f'{job_log_folder}{job_name}.sh'
    with open(job_path_name, 'w') as file:
        file.write(cmd)
    os.system(f'bsub < {job_path_name}')
    output(f'sending {job_name} as ID_{rand_id}...', log_file=log_file, run_main=run_main)
    time.sleep(2)  # too many conda activations at once sometimes leads to inability to activate during a job.

    return rand_id


def job_wait(id_list, log_file, run_main=False):
    '''
    Waits for jobs sent by send job to finish.
    '''
    waiting = True
    while waiting:
        with os.popen('bhist -w') as stream:
            job_list = stream.read()
        current = []
        for rand_id in id_list:
            if len([j for j in re.findall(r'ID_(\d+)', job_list) if j == rand_id]) != 0:
                current.append(rand_id)
        if len(current) == 0:
            waiting = False
        else:
            output(f'Waiting for jobs to finish... {datetime.now():%Y-%m-%d %H:%M:%S}', log_file=log_file, run_main=run_main)
            time.sleep(60)


def move_file(file, dest, dest_type, log_file, run_main):
    if dest_type.lower() == 'folder':
        os.makedirs(dest, exist_ok=True)
    if (os.path.isfile(file)) and (os.path.isfile(dest) is False):
        try:
            shutil.move(file, dest)
        except FileNotFoundError:
            print(f'{file} not found.')
    else:
        output(f'{file} cannot be moved to {dest}', log_file=log_file, run_main=run_main)


def rename(file, name):
    if os.path.isfile(file):
        os.rename(file, name)


def job_pending(job, log_file, run_main=False):
    '''
    Waits for jobs sent by send job to start running.
    '''
    waiting = True
    while waiting:
        with os.popen('bjobs -p') as stream:
            job_list = stream.read()
        if len([j for j in re.findall(r'ID_(\d+)', job_list) if j == job]) != 0:
            output(f'Waiting for jobs to start running... {datetime.now():%Y-%m-%d %H:%M:%S}', log_file=log_file, run_main=run_main)
        else:
            waiting = False
        time.sleep(60)


def validated_run(task, func, exp):
    try:
        if task in exp.tasks_complete:
            output(f'Skipping {task}...', log_file=exp.log_file, run_main=exp.run_main)
            return exp
        else:
            return func(exp)
    except:
        close_out(task, exp)


def clean_encode_folder(exp):
    encode_dir = f'{exp.scratch}ENCODE3/'

    mk = ['bams', 'sample_peaks', 'nodup_bam', 'overlap_peaks', 'idr_peaks', 'reports', 'bws', 'jsons']

    for folder in mk:
        make_folder(f'{exp.scratch}{folder}')

    for sample in exp.IPs.Condition.unique().tolist():
        crom_dir = f'{encode_dir}{sample}/cromwell-executions/chip/*/'
        rep_number = len(exp.IPs[exp.IPs.Condition == sample])

        rename(glob_check(f'{crom_dir}call-reproducibility_idr/execution/optimal_peak.narrowPeak.gz'),
               f'{exp.scratch}idr_peaks/{sample}.idr.optimal_peak.narrowPeak.gz')
        rename(glob_check(f'{crom_dir}call-reproducibility_idr/execution/conservative_peak.narrowPeak.gz'),
               f'{exp.scratch}idr_peaks/{sample}.idr.conservative_peak.narrowPeak.gz')
        rename(glob_check(f'{crom_dir}call-reproducibility_overlap/execution/optimal_peak.narrowPeak.gz'),
               f'{exp.scratch}overlap_peaks/{sample}.overlap.optimal_peak.narrowPeak.gz')
        rename(glob_check(f'{crom_dir}call-reproducibility_overlap/execution/conservative_peak.narrowPeak.gz'),
               f'{exp.scratch}overlap_peaks/{sample}.overlap.conservative_peak.narrowPeak.gz')
        rename(glob_check(f'{crom_dir}call-qc_report/execution/qc.html'), f'{exp.scratch}reports/{sample}_qc_report.html')
        rename(glob_check(f'{encode_dir}{sample}/*.json'), f'{exp.scratch}jsons/{sample}_ENCODE3.json')

        for rep in range(rep_number):
            n = rep + 1
            rename(glob_check(f'{crom_dir}call-macs2/shard-{rep}/execution/*fc.signal.bigwig'),
                   f'{exp.scratch}bws/{sample}_Rep{n}.fc.signal.bigwig')
            rename(glob_check(f'{crom_dir}call-macs2/shard-{rep}/execution/*bfilt.narrowPeak.gz'),
                   f'{exp.scratch}sample_peaks/{sample}_Rep{n}.bfilt.narrowPeak.gz')
            rename(glob_check(f'{crom_dir}call-macs2/shard{rep}/execution/*broadPeak'),
                   f'{exp.scratch}sample_peaks/{sample}_Rep{n}.broadPeak')
            rename(glob_check(f'{crom_dir}call-filter/shard-{rep}/execution/*.nodup.bam'),
                   f'{exp.scratch}nodup_bam/{sample}_Rep{n}.nodup.bam')
            rename(glob_check(f'{crom_dir}call-filter/shard-{rep}/execution/*.nodup.bam.bai'),
                   f'{exp.scratch}nodup_bam/{sample}_Rep{n}.nodup.bam.bai')
            rename(glob_check(f'{crom_dir}call-filter_ctl/shard-{rep}/execution/*nodup.bam'),
                   f'{exp.scratch}nodup_bam/{sample}_Rep{n}_background.nodup.bam')
            rename(glob_check(f'{crom_dir}call-filter_ctl/shard-{rep}/execution/*nodup.bam.bai'),
                   f'{exp.scratch}nodup_bam/{sample}_Rep{n}_background.nodup.bam.bai')
            rename(glob_check(f'{crom_dir}call-bwa/shard-{rep}/execution/*.bam'), f'{exp.scratch}bams/{sample}_Rep{n}.bam')
            rename(glob_check(f'{crom_dir}call-bwa/shard-{rep}/execution/*.bam.bai'), f'{exp.scratch}bams/{sample}_Rep{n}.bam.bai')
            rename(glob_check(f'{crom_dir}call-bwa_ctl/shard-{rep}/execution/*.bam'), f'{exp.scratch}bams/{sample}_Rep{n}_background.bam')
            rename(glob_check(f'{crom_dir}call-bwa_ctl/shard-{rep}/execution/*.bam.bai'), f'{exp.scratch}bams/{sample}_Rep{n}_background.bam.bai')

    folders = ['raw_data']
    for folder in folders:
        glob_folders = glob.glob(f'{exp.scratch}{folder}')
        for glob_folder in glob_folders:
            if os.path.isdir(glob_folder):
                shutil.rmtree(glob_folder)
            else:
                output(f'Could not find {glob_folder}', exp.log_file, exp.run_main)


def extract_ENCODE_report_data(exp):
    '''
    Inputs
    -----
    base_folder:  AQUAS results folder.  Will use subfolders for sample name and look for report in those subfolders.
    report_type: 'AQUAS' or 'cromwell'
    replicate: Whether the ChIPseq was performed as a repliate or not.
    Returns
    -----
    DataFrame of results
    '''
    import pandas as pd
    from chrome_chip.plot import plot_col

    base_folder = val_folder(f'{exp.scratch}reports')

    reports = glob.glob(f'{base_folder}*report.html')

    results_df = pd.DataFrame(index=['Percent_mapped',
                                     'Filtered_Uniquely_Mapped_Reads',
                                     'Fraction_Duplicated',
                                     'S_JS_Distance',
                                     'PBC1',
                                     'RSC',
                                     'Overlap_Optimal_Peak_Number',
                                     'Overlap_Conservative_Peak_Number',
                                     'FrIP_IDR',
                                     'IDR_Peak_Number'])

    for file in reports:
        name = file.split('/')[-1].split('_qc_')[0]
        report = pd.read_html(file)
        if len(report) < 6:
            output(f'Skipping {name} from qc plots.', exp.log_file, exp.run_main)
            continue
        series = pd.Series()
        series['Percent_mapped'] = report[0].iloc[7, 1]
        series['Filtered_Uniquely_Mapped_Reads'] = report[3].iloc[5, 1]
        series['Fraction_Duplicated'] = report[1].iloc[7, 1]
        series['PBC1'] = report[2].iloc[6, 1]
        series['RSC'] = report[5].iloc[9, 1]
        series['Overlap_Optimal_Peak_Number'] = report[4].iloc[4, 1]
        series['Overlap_Conservative_Peak_Number'] = report[4].iloc[5, 1]

        chip_type = exp.IPs[exp.IPs.Condition == name]['ChIP Type'].unique().tolist()
        if 'tf' in [x.lower() for x in chip_type]:
            series['FrIP_IDR'] = report[7].iloc[1, 1]
            series['IDR_Peak_Number'] = report[4].iloc[4, 2]
            series['S_JS_Distance'] = report[8].iloc[8, 1]
        else:
            series['S_JS_Distance'] = report[7].iloc[8, 1]
        results_df[name] = series

    for index in results_df.index.tolist():
        plot_col(results_df.loc[index], out=base_folder, title=f'{index}', ylabel=index.replace('_', ' '), plot_type=['violin', 'swarm'], log_file=exp.log_file, run_main=exp.run_main)

    return results_df


def submission_prepend(submission=None, source='chrome_chip', conda=None, module_list=['python']):
    '''
    Prepends a string for a submission script
    Edit these defualts to optimize for another HPC or other environment
    '''
    prepend = f'module rm {" ".join(module_list)}\n' if len(module_list) > 0 else ''

    if source:
        prepend += f'source activate {source}\n'

    elif conda:
        prepend += f'conda activate {conda}\n'

    return prepend if submission is None else prepend + submission
