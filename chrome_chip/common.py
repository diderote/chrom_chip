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
from datetime import datetime

from IPython.display import HTML, display, Image
import pandas as pd

from . import __name__

__author__ = 'Daniel L. Karl'
__license__ = 'MIT'
__version__ = '0.1'


def run_main():
    return True if __name__ == '__main__' else False


def version():
    return __version__


def html_header(version=__version__, author=__author__, license=__license__):
    return ''.join(['<h1>ChIPseq Analysis Notebook</h1>',
                    f'<body><b>Experiment Date: {datetime.now():%Y-%m-%d}<br>',
                    f'Pipeline version: {version}</b><br>',
                    '<a href="http://www.github.com/diderote/LSF-aquaFor">Pipeline Code</a><br>',
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
    glob_check(file)
    if (file.split('.')[-1] == 'txt') or (file.split('.')[-1] == 'tab'):
        return pd.read_table(file, header=0, index_col=0)
    elif (file.split('.')[-1] == 'xls') or (file.split('.')[-1] == 'xlsx'):
        return pd.read_excel(file, index_col=0)
    else:
        raise IOError("Cannot parse file.  Make sure it is .txt, .xls, or .xlsx")


def output(text, log_file=None, run_main=run_main()):
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
    display(Image(file))


def glob_check(path):
    file = glob.glob(path)
    return 'none' if len(file) == 0 else file[0]


def gunzip_return(path):
    if path == 'none':
        return path
    elif path.endswith('.gz'):
        os.system(f'gunzip {path}')
        return path[:-3]
    else:
        return path


def txt_replace(string):
    return string.replace('.txt.bz2', '.fastq.gz')


def out_result(image, text, run_main):
    if not run_main:
        if os.path.isfile(image):
            display(HTML(f'<h2>{text}</h2>'))
            image_display(image)
        else:
            display(HTML(f'<body>No result for {text} found.</body>'))


def close_out(task, exp):
        output(f'Error in {task}.', exp.log_file)
        filename = f'{exp.scratch}{exp.name}_incomplete.pkl'
        with open(filename, 'wb') as experiment:
            pickle.dump(exp, experiment)
        raise RuntimeError(f'Error in {task}. Fix problem then resubmit with same command to continue from last completed step.')


def send_job(command_list, job_name, job_log_folder, q, mem, log_file, project, cores=1, submit=False):
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
    output(f'sending {job_name} as ID_{rand_id}...', log_file)
    time.sleep(2)  # too many conda activations at once sometimes leads to inability to activate during a job.

    return rand_id


def job_wait(id_list, log_file):
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
            output(f'Waiting for jobs to finish... {datetime.now():%Y-%m-%d %H:%M:%S}', log_file)
            time.sleep(60)


def job_pending(job, log_file):
    '''
    Waits for jobs sent by send job to start running.
    '''
    waiting = True
    while waiting:
        with os.popen('bjobs -p') as stream:
            job_list = stream.read()
        if len([j for j in re.findall(r'ID_(\d+)', job_list) if j == job]) != 0:
            output(f'Waiting for jobs to start running... {datetime.now():%Y-%m-%d %H:%M:%S}', log_file)
        else:
            waiting = False
        time.sleep(60)


def validated_run(task, func, exp):
    try:
        if task in exp.tasks_complete:
            output(f'Skipping {task}...', exp.log_file)
            return exp
        else:
            return func(exp)
    except:
        close_out(task, exp)
