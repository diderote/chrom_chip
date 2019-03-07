#!/usr/bin/env python3

import glob
import os
import time
from datetime import datetime
from shutil import copy2

from chrome_chip.common import output, val_folder, send_job, job_wait, make_folder, is_fastq, txt_replace, submission_prepend


def stage(exp):
    '''
    Stages files in scratch folder
    '''
    output(f'Staging in {exp.scratch}\n', log_file=exp.log_file, run_main=exp.run_main)
    exp.data_folder = make_folder(f'{exp.scratch}raw_data/')

    Scratch_File1 = []
    Scratch_File2 = []

    # join multiple files
    for sample in exp.sample_df.Sample_Name.tolist():
        index = exp.sample_df['Sample_Name'] == sample

        paired = exp.sample_df.loc[index, 'paired'].values[0]
        R1_list = ','.join(exp.sample_df.loc[index, 'File1']).split(',')
        R2_list = ','.join(exp.sample_df.loc[index, 'File2']).split(',')

        main_file = R1_list[0]

        # convert form bz2 to fastq
        if main_file.endswith('.tzt.bz2'):
            for file in R1_list + R2_list:
                newfile = file.replace(".txt.bz2", ".fastq.gz")
                os.system(f'bunzip2 -c < {file} | gzip -c > {newfile}')
            R1_list = [txt_replace(x) for x in R1_list]
            R2_list = [txt_replace(x) for x in R2_list]
            main_file = txt_replace(main_file)

        if main_file.endswith('.fastq.gz'):
            fileend = '_R1.fastq.gz' if paired else '.fastq.gz'
            filename = f'{exp.data_folder}{sample}{fileend}'
            os.system(f'cat {" ".join(R1_list)} > {filename}')
            Scratch_File1.append(filename)

            if paired:
                fileend = '_R2.fastq.gz'
                filename = f'{exp.data_folder}{sample}{fileend}'
                os.system(f'cat {" ".join(R2_list)} > {filename}')
                Scratch_File2.append(filename)
            else:
                Scratch_File2.append('none')

        elif main_file.endswith('.bam'):
            copy2(main_file, f'{exp.data_folder}{sample}.bam')

        else:
            raise IOError('Filetype not recognized.')

    exp.sample_df['Scratch_File1'] = Scratch_File1
    exp.sample_df['Scratch_File2'] = Scratch_File2
    exp.sample_df.replace([f'{exp.data_folder}none'], 'none', inplace=True)

    exp.tasks_complete.append('Stage')
    output(f'Staging complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log_file=exp.log_file, run_main=exp.run_main)

    return exp


def fastqc(exp):
    '''
    Performs fastq spec analysis with FastQC
    '''
    output('Assessing fastq quality. \n', log_file=exp.log_file, run_main=exp.run_main)

    # Make QC folder
    exp.qc_folder = make_folder(f'{exp.scratch}QC/')

    all_samples = exp.sample_df.Scratch_File1.tolist() + exp.sample_df.Scratch_File2.tolist()
    samples = [file for file in all_samples if is_fastq(file)]

    for sample in samples:
        command_list = [submission_prepend(f'fastqc {sample}')]

        exp.job_id.append(send_job(command_list=command_list,
                                   job_name=f'{sample.split("/")[-1]}_fastqc',
                                   job_log_folder=exp.job_folder,
                                   q='general',
                                   mem=5000,
                                   log_file=exp.log_file,
                                   project=exp.project,
                                   run_main=exp.run_main
                                   ))

    # Wait for jobs to finish
    job_wait(exp.job_id, exp.log_file)

    # move to qc folder
    fastqc_files = glob.glob(f'{exp.data_folder}*.zip')
    fastqc_files = fastqc_files + glob.glob(f'{exp.data_folder}*.html')
    for f in fastqc_files:
        copy2(f, exp.qc_folder)
        os.remove(f)

    exp.tasks_complete.append('FastQC')
    output(f'FastQC complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log_file=exp.log_file, run_main=exp.run_main)

    return exp


def fastq_screen(exp):
    '''
    Checks fastq files for contamination with alternative genomes using Bowtie2
    '''

    output(f'Screening for contamination during sequencing: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log_file=exp.log_file, run_main=exp.run_main)

    # Make QC folder
    exp.qc_folder = make_folder(f'{exp.scratch}QC/')

    cwd = val_folder(os.getcwd())
    os.chdir(exp.data_folder)

    samples = [file for file in exp.sample_df.Scratch_File1.tolist() if is_fastq(file)]

    # Submit fastqc and fastq_screen jobs for each sample
    for sample in samples:
        command_list = [submission_prepend(f'fastq_screen --threads 4 --aligner bowtie2 {sample}')]

        exp.job_id.append(send_job(command_list=command_list,
                                   job_name=f'{sample.split("/")[-1]}_fastq_screen',
                                   job_log_folder=exp.job_folder,
                                   q='general',
                                   mem=3000,
                                   log_file=exp.log_file,
                                   project=exp.project,
                                   cores=2,
                                   run_main=exp.run_main
                                   ))
        time.sleep(1)

    # Wait for jobs to finish
    job_wait(exp.job_id, exp.log_file)

    # move to qc folder
    fastqs_files = glob.glob(f'{exp.data_folder}*screen*')
    for f in fastqs_files:
        copy2(f, exp.qc_folder)
        os.remove(f)

    # change to experimental directory in scratch
    os.chdir(cwd)

    exp.tasks_complete.append('Fastq_screen')
    output(f'Screening complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log_file=exp.log_file, run_main=exp.run_main)

    return exp


def trim(exp):
    '''
    Trimming based on standard UM SCCC Core Nextseq 500 technical errors.
    Cudadapt can hard clip both ends, but may ignore 3' in future.
    '''

    output(f'Beginning fastq trimming: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log_file=exp.log_file, run_main=exp.run_main)

    for sample_dict in exp.sample_df[['Scratch_File1', 'Scratch_File2', 'Sequencer', 'Sample_Name']].to_dict(orient='records'):

        quality = '--nextseq-trim=20' if sample_dict['Sequencer'].lower() == 'nextseq' else '-q 20'
        seq_type = 'single' if sample_dict['Scratch_File2'] == 'none' else 'paired'

        sample = sample_dict["Sample_Name"]

        paired = f'{exp.data_folder}{sample}_trim_R2.fastq.gz'
        single = f'{exp.data_folder}{sample}_trim.fastq.gz'
        data_files = glob.glob(f'{exp.data_folder}*.gz')

        if (single in data_files) or (paired in data_files):
            continue
        else:
            output(f'Trimming {sample}: ', log_file=exp.log_file, run_main=exp.run_main)

            if seq_type == 'paired':
                cutadapt = f'cutadapt -j 4 -a AGATCGGAAGAGC -A AGATCGGAAGAGC --cores=10 {quality} -m 18 '
                cutadapt += f'-o {exp.data_folder}{sample}_trim_R1.fastq.gz -p {exp.data_folder}{sample}_trim_R2.fastq.gz '
                cutadapt += f'{sample_dict["Scratch_File1"]} {sample_dict["Scratch_File2"]}'
            elif seq_type == 'single':
                cutadapt = f'cutadapt -j 4 -a AGATCGGAAGAGC --cores=10 {quality} -m 18 '
                cutadapt += f'-o {exp.data_folder}{sample}_trim_R1.fastq.gz {sample_dict["Scratch_File1"]}'

            command_list = [submission_prepend(cutadapt)]

            exp.job_id.append(send_job(command_list=command_list,
                                       job_name=f"{sample}_trim",
                                       job_log_folder=exp.job_folder,
                                       q='general',
                                       mem=1000,
                                       log_file=exp.log_file,
                                       project=exp.project,
                                       cores=2,
                                       run_main=exp.run_main
                                       ))

    # Wait for jobs to finish
    job_wait(exp.job_id, exp.log_file, exp.run_main)

    # move logs to qc folder
    output('\nTrimming logs are found in stdout files from bsub.  Cutadapt does not handle log files in multi-core mode.', log_file=exp.log_file, run_main=exp.run_main)

    exp.tasks_complete.append('Trim')
    output(f'Trimming complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log_file=exp.log_file, run_main=exp.run_main)

    return exp
