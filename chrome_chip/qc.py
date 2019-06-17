#!/usr/bin/env python3

import os
from datetime import datetime
from shutil import copy2, copytree, rmtree
import pickle

import yaml
from IPython.display import HTML, display

from chrome_chip.common import output, send_job, read_pd, close_out, make_folder, submission_prepend, clean_encode_folder, extract_ENCODE_report_data
from chrome_chip.plot import plot_col


def preseq(exp):

    output('\nRunning QC plots: library complexity extrapolation, signal correlation and pca plots.', log_file=exp.log_file, run_main=exp.run_main)

    for sample in exp.samples:

        out_dir = make_folder(f'{exp.scratch}QC/preseq/{sample}/')

        command_list = [submission_prepend(f'preseq lc_extrap -bam -output {out_dir}{sample}_preseq.txt {exp.sample_files[sample]["bam"]}')]

        exp.job_id.append(send_job(command_list=command_list,
                                   job_name=f"{sample}_preseq",
                                   job_log_folder=exp.job_folder,
                                   q='general',
                                   mem=5000,
                                   log_file=exp.log_file,
                                   project=exp.project,
                                   cores=1,
                                   run_main=exp.run_main
                                   ))

    exp.tasks_complete.append('preseq')

    return exp


def principal_component_analysis(exp):

    out_dir = make_folder(f'{exp.scratch}PCA/')

    bigwigs = {sample: exp.sample_files[sample]['bw'] for sample in exp.samples if len(exp.sample_files[sample]['bw']) != 0}
    multibw_command = f"multiBigwigSummary bins -b {' '.join(list(bigwigs.values()))} -l {' '.join(list(bigwigs.keys()))} -p 4 --chromosomesToSkip chrM,chrX,chrY -o {out_dir}{exp.name}_bwsummary.npz"

    correlation_command = f'plotCorrelation --corData {out_dir}{exp.name}_bwsummary.npz --corMethod pearson --whatToPlot heatmap --skipZeros --plotTitle "{exp.name} Binned Pearson Correlation Heatmap" --plotFileFormat png --outFileCorMatrix {out_dir}{exp.name}_CorMatrix.tab --colorMap Purples -o {out_dir}{exp.name}_CorHeatmap.png'

    pca_command = f'plotPCA --corData {out_dir}{exp.name}_bwsummary.npz --plotTitle "{exp.name} PCA Plot" --plotFileFormat png --outFileNameData {out_dir}{exp.name}_PCA_data.tab --log2 -o {out_dir}{exp.name}_PCA_Plot.png'

    command_list = [submission_prepend(),
                    multibw_command,
                    correlation_command,
                    pca_command
                    ]

    exp.job_id.append(send_job(command_list=command_list,
                               job_name=f"{exp.name}_Cor_PCA",
                               job_log_folder=exp.job_folder,
                               q='general',
                               mem=4000,
                               log_file=exp.log_file,
                               project=exp.project,
                               cores=5,
                               run_main=exp.run_main
                               ))

    exp.tasks_complete.append('PCA')

    return exp


def final_qc(exp):

    ''' add preseq '''
    try:
        output(f'Beginning final qc: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log_file=exp.log_file, run_main=exp.run_main)

        if os.path.isdir(f'{exp.scratch}QC/multiqc_data') is False:
            os.system(f'multiqc {exp.scratch}* -o {exp.scratch}QC/')

        # Summary plots for FastQC data
        fastqc_file = f'{exp.scratch}/QC/multiqc_data/multiqc_fastqc.txt'
        if os.path.isfile(fastqc_file):
            gen_stats = read_pd(f'{exp.scratch}/QC/multiqc_data/multiqc_general_stats.txt')
            samples = exp.sample_df.Sample_Name.tolist()
            plot_col(df=gen_stats.loc[samples, 'FastQC_mqc-generalstats-fastqc-total_sequences'] / 1e6,
                     title='Total Sequencer Reads per Sample',
                     ylabel='Reads (Millions)',
                     log_file=exp.log_file,
                     run_main=exp.run_main
                     )

            plot_col(df=gen_stats.loc[samples, 'FastQC_mqc-generalstats-fastqc-percent_gc'],
                     title='Percent GC Content per Sample',
                     ylabel='Percentage of Reads with GC Content',
                     log_file=exp.log_file,
                     run_main=exp.run_main
                     )

            if os.path.isdir('plots/'):
                copytree('plots/', f'{exp.scratch}QC/plots/')
                rmtree('plots')

        display(HTML('<h1>Final QC Summary</h1>'))
        display(HTML(f'{exp.scratch}/multiqc_report.html'))

        exp.tasks_complete.append('MultiQC')

        return exp

    except:
        close_out('final qc', exp)


def finish(exp):
    try:
        output('Cleaning and reorganizing ENCODE3 files...\n')
        clean_encode_folder(exp)
        extract_ENCODE_report_data(exp)

        output(f'\nConda environment file: {exp.job_folder}{exp.name}_environmnet.yml\nPackage versions: ', log_file=exp.log_file, run_main=exp.run_main)
        os.system(f'conda env export > {exp.job_folder}{exp.name}_environmnet.yml')
        with open(f'{exp.job_folder}{exp.name}_environmnet.yml', 'r') as fp:
            versions = yaml.load(fp)
        for package in versions['dependencies']:
            output(package, log_file=exp.log_file, run_main=exp.run_main)

        output(f'\n{exp.name} analysis complete! \n', log_file=exp.log_file, run_main=exp.run_main)
        output(f'Copying all results into {exp.out_dir}: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log_file=exp.log_file, run_main=exp.run_main)

        scratch_log = f'{exp.scratch}{exp.log_file.split("/")[-1]}'
        if exp.run_main:
            copy2(exp.log_file, scratch_log)

        '''for sample in exp.sample_files.values():
            for file in exp.sample_files[sample].values():
                if exp.scratch in file:
                    file.replace(exp.scratch, exp.out_dir)'''

        rmtree(exp.out_dir)
        copytree(exp.scratch, exp.out_dir)
        rmtree(exp.scratch)

        exp.tasks_complete.append('Finished')

        filename = f'{exp.out_dir}{exp.name}_{exp.date}.pkl'
        with open(filename, 'wb') as experiment:
            pickle.dump(exp, experiment)

        output(f'Python Experiment: \n{exp}', log_file=exp.log_file, run_main=exp.run_main)
        output(f'Moved all files into {exp.out_dir}: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log_file=exp.log_file, run_main=exp.run_main)

        return exp

    except:
        close_out('finishing pipeline', exp)
