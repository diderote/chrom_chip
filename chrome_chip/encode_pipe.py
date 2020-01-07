#!/usr/bin/env python3

import json
import shutil
import re
import os
from datetime import datetime

from chrome_chip.common import output, send_job, job_wait, job_pending, glob_check, make_folder, glob_remove, submission_prepend, out_result
from chrome_chip.preprocess import stage
from chrome_chip.plot import spike_in_plot


def encode3(exp):

    if 'Stage' not in exp.tasks_complete:
        output('Files not staged.\n', log_file=exp.log_file)
        exp = stage(exp)

    output('Running alignment and peak calling using ENCODE3 standards.', log_file=exp.log_file, run_main=exp.run_main)
    output('ENCODE3 cromwell pipeline.', log_file=exp.log_file, run_main=exp.run_main)

    out_dir = make_folder(f'{exp.scratch}ENCODE3/')

    IPs = exp.IPs

    end_types = {'q.gz': 'fastq',
                 '.bam': 'bam'
                 }

    for experiment in IPs.Condition.unique().tolist():

        exp_dir = make_folder(f'{out_dir}{experiment}/')

        IP_sample_indicies = [(rep, index) for rep, index in enumerate(IPs[IPs.Condition == experiment].index.tolist(), start=1)]

        if len(IP_sample_indicies) > 6:
            raise IOError('Pipeline cannot handle more than 6 replicates.')

        seq_type = False if 'none' in IPs[IPs.Condition == experiment]['File2'].tolist() else True
        final_stage = 'align' if 'align' in IPs[IPs.Condition == experiment]['Final Stage'].tolist() else 'all'

        UMI_list = [x.lower() for x in IPs[IPs.Condition == experiment]['UMI'].unique().tolist()]
        if len(set(UMI_list)) > 1:
            raise IOError('All samples must be UMI processed or not for each condition.')
        UMI = True if UMI_list[0] == 'yes' else False

        try:
            file_type = end_types[exp.sample_df[exp.sample_df.Condition == experiment]['Scratch_File1'].tolist()[0][-4:]]
        except KeyError:
            output(f"{exp.sample_df[exp.sample_df.Condition == experiment]['Scratch_File1'].tolist()[0]} not a valid file type for this pipeline.", log_file=exp.log_file, run_main=exp.run_main)

        file_type = 'bam' if (UMI is True) & ('UMI' in exp.tasks_complete) else file_type

        genome = IPs[IPs.Condition == experiment]['Genome'].unique().tolist()
        if len(genome) > 1:
            raise IOError('Cannot align to more than one genome per condition.')

        chip_type = IPs[IPs.Condition == experiment]['ChIP Type'].unique().tolist()
        if len(chip_type) > 1:
            raise IOError('Cannot have more than one chip type (histone or TF) for a condition.')
        chip_type = 'histone' if chip_type[0].lower() == 'histone' else 'tf'

        json_file = {'chip.pipeline_type': chip_type,
                     'chip.paired_end': seq_type,
                     'chip.genome_tsv': exp.genome_indicies['encode_tsv'][genome[0]],
                     'chip.bwa.mem_mb': 30000,
                     'chip.macs2_mem_mb': 30000,
                     'chip.peak_caller': 'macs2',
                     "chip.true_rep_only": False,
                     "chip.dup_marker": "picard",
                     "chip.mapq_thresh": 30,
                     "chip.regex_filter_reads": "chrM",
                     "chip.subsample_reads": 0,
                     "chip.ctl_subsample_reads": 0,
                     "chip.xcor_subsample_reads": 15000000,
                     "chip.keep_irregular_chr_in_bfilt_peak": False,
                     "chip.always_use_pooled_ctl": False,
                     "chip.ctl_depth_ratio": 1.2,
                     "chip.macs2_cap_num_peak": 500000,
                     "chip.pval_thresh": 0.01,
                     "chip.idr_thresh": 0.05,
                     "chip.bwa_cpu": 4,
                     "chip.bwa_mem_mb": 20000,
                     "chip.bwa_time_hr": 48,
                     "chip.filter_cpu": 2,
                     "chip.filter_mem_mb": 20000,
                     "chip.filter_time_hr": 24,
                     "chip.bam2ta_cpu": 2,
                     "chip.bam2ta_mem_mb": 10000,
                     "chip.bam2ta_time_hr": 6,
                     "chip.fingerprint_cpu": 2,
                     "chip.fingerprint_mem_mb": 12000,
                     "chip.fingerprint_time_hr": 6,
                     "chip.xcor_cpu": 2,
                     "chip.xcor_mem_mb": 16000,
                     "chip.xcor_time_hr": 24,
                     "chip.macs2_time_hr": 24,
                     "chip.spr_mem_mb": 16000
                     }
        bams = []
        ctl_bams = []

        for rep, index in IP_sample_indicies:
            sample = exp.sample_df.loc[index, 'Sample_Name']
            input_sample = IPs.loc[index, 'Background_Name']

            if file_type == 'fastq':
                json_file[f'chip.fastqs_rep{rep}_R1'] = [f'{exp.data_folder}{sample}_trim_R1.fastq.gz']
                json_file[f'chip.ctl_fastqs_rep{rep}_R1'] = [f'{exp.data_folder}{input_sample}_trim_R1.fastq.gz']
                if seq_type:
                    json_file[f'chip.fastqs_rep{rep}_R2'] = [f'{exp.data_folder}{sample}_trim_R2.fastq.gz']
                    json_file[f'chip.ctl_fastqs_rep{rep}_R2'] = [f'{exp.data_folder}{input_sample}_trim_R2.fastq.gz']
            else:
                bams.append(f'{exp.data_folder}{sample}.bam')
                ctl_bams.append(f'{exp.data_folder}{input_sample}.bam')

        if file_type == 'bam':
            json_file[f'chip.bams'] = bams
            json_file[f'chip.ctl_bams'] = ctl_bams

        json_file['chip.align_only'] = True if UMI & (file_type == 'fastq') else False
        json_file['chip.align_only'] = True if final_stage == 'align' else json_file['chip.align_only']

        json_file['chip.no_dup_removal'] = True if UMI else False
        json_file['chip.title'] = f'{experiment}_postUMI_dedup' if UMI & (file_type == 'bam') else experiment
        json_file["chip.description"] = f"Cromwell ENCODE3 {experiment}: {'paired-end' if seq_type else 'single-end'} {chip_type}."

        encode_file = f'{exp_dir}{experiment}_ENCODE3.json'
        with open(encode_file, 'w') as file:
            json.dump(json_file, file, indent=4, sort_keys=True)

        pythonpath = shutil.which('python')
        miniconda = [x for x in pythonpath.split('/') if 'miniconda' in x]
        cromwell_jar = re.sub(r'{}/.*'.format(miniconda), '{}/envs/chrome_chip/share/cromwell/cromwell.jar'.format(miniconda), pythonpath)
        jar = cromwell_jar if os.path.isfile(cromwell_jar) else '~/miniconda3/envs/chrome_chip/share/cromwell/cromwell.jar'

        command_list = [submission_prepend(source='encode-chip-seq-pipeline'),
                        f'cd {exp_dir}',
                        f'java -jar -Dconfig.file={exp.encode3_folder}backends/backend.conf -Dbackend.default=Local {jar} run {exp.encode3_folder}chip.wdl -i {encode_file}'
                        ]

        sent_job = send_job(command_list=command_list,
                            job_name=f"{experiment}_ENCODE3",
                            job_log_folder=exp.job_folder,
                            q='bigmem',
                            mem=35000,
                            log_file=exp.log_file,
                            project=exp.project,
                            cores=1,
                            run_main=exp.run_main
                            )

        exp.job_id.append(sent_job)
        job_pending(sent_job, exp.log_file)

    # Wait for jobs to finish
    job_wait(exp.job_id, exp.log_file)

    # Check fraglength and resubmit with set 200 fraglen for macs2 if xcor error
    for experiment in exp.IPs.Condition.unique().tolist():
        rep_number = len(exp.IPs[exp.IPs.Condition == experiment])
        
        frag_list = []

        for rep in range(rep_number):
            file = glob_check(f'{exp.scratch}ENCODE3/{experiment}/cromwell-executions/chip/*/call-xcor/shard-{rep}/execution/*fraglen.txt')
            with open(file, 'r') as f:
                frag_list.append(f.read().split()[0])

        if '-' in [x[0] for x in frag_list]:
            output(f'Xcor failed for {experiment}.  Resubmitting with fragment length set to 200 for failed sample/s', 
                   log_file=exp.log_file, run_main=exp.run_main)

            frag_list = [x if x[0] != '-' else '200' for x in frag_list]
            exp_dir = f'{exp.scratch}ENCODE3/{experiment}/'
            encode_file = f'{exp_dir}{experiment}_ENCODE3.json'

            with open(encode_file, 'r') as file:
                json_file = json.load(file)

            json_file["chip.fraglen"] = frag_list

            resubmit_file = f'{exp_dir}/{experiment}_ENCODE3_setfraglenth.json'
            with open(resubmit_file, 'w') as file:
                json.dump(json_file, file, indent=4, sort_keys=True)

            pythonpath = shutil.which('python')
            miniconda = [x for x in pythonpath.split('/') if 'miniconda' in x]
            cromwell_jar = re.sub(r'{}/.*'.format(miniconda), '{}/envs/chrome_chip/share/cromwell/cromwell.jar'.format(miniconda), pythonpath)
            jar = cromwell_jar if os.path.isfile(cromwell_jar) else '~/miniconda3/envs/chrome_chip/share/cromwell/cromwell.jar'

            command_list = [submission_prepend(source='encode-chip-seq-pipeline'),
                            f'cd {exp_dir}',
                            f'java -jar -Dconfig.file={exp.encode3_folder}backends/backend.conf -Dbackend.default=Local {jar} run {exp.encode3_folder}chip.wdl -i {resubmit_file}'
                            ]

            sent_job = send_job(command_list=command_list,
                                job_name=f"{experiment}_ENCODE3_resubmission",
                                job_log_folder=exp.job_folder,
                                q='bigmem',
                                mem=35000,
                                log_file=exp.log_file,
                                project=exp.project,
                                cores=1,
                                run_main=exp.run_main
                                )

            exp.job_id.append(sent_job)
            job_pending(sent_job, exp.log_file)

    # Wait for jobs to finish
    job_wait(exp.job_id, exp.log_file)

    exp = encode_results(exp)

    exp.tasks_complete.append('ENCODE3')

    return exp


def UMI(exp):

    # exp.data_type = 'bam'

    IPs = exp.IPs

    for experiment in IPs.Condition.unique().tolist():
        UMI = True if 'yes' in IPs[IPs.Condition == experiment]['UMI'].tolist() else False

        if not UMI:
            return exp

        else:
            out_dir = make_folder(f'{exp.scratch}UMI/')
            output('Deduplicating bam files using UMIs with UMI-tools.', log_file=exp.log_file, run_main=exp.run_main)

            for index in IPs[IPs.Condition == experiment].index.tolist():
                sample = IPs.loc[index, 'Sample_Name']
                input_sample = IPs.loc[index, 'Background_Name']

                bam = exp.sample_files[sample]['bam']
                input_bam = exp.sample_files[input_sample]['bam']
                nodup_bam = f'{out_dir}{sample}.UMI.dedup.bam'
                nodup_input = f'{out_dir}{input_sample}.UMI.dedup.bam'

                umi_string = 'umi_tools dedup --umi-separator=":" --output-stats={out_dir}{sample}deduplicated.qc -I {inbam} -S {outbam} -L {out_dir}{sample}.UMI.log'

                seq_type = False if 'none' in IPs[IPs.Condition == experiment]['Scratch_File2'].tolist() else True
                if seq_type == 'paired':
                    umi_string += ' --paired'

                command_list = [submission_prepend(),
                                f'samtools index {bam}',
                                f'samtools index {input_bam}',
                                umi_string.format(inbam=bam, outbam=nodup_bam, sample=sample, out_dir=out_dir),
                                umi_string.format(inbam=input_bam, outbam=nodup_input, sample=input_sample, out_dir=out_dir)
                                ]

                exp.job_id.append(send_job(command_list=command_list,
                                           job_name=f"{sample}_UMI_dedup",
                                           job_log_folder=exp.job_folder,
                                           q='bigmem',
                                           mem=40000,
                                           log_file=exp.log_file,
                                           project=exp.project,
                                           cores=1,
                                           run_main=exp.run_main
                                           ))

                exp.sample_files[sample]['nodup_bam'] = nodup_bam
                exp.sample_files[input_sample]['nodup_bam'] = nodup_input

    job_wait(exp.job_id, exp.log_file)

    output('Dedplication complete.  Submitting deduplicated files for the remainder of processing.', log_file=exp.log_file, run_main=exp.run_main)
    exp.tasks_complete.append('UMI')

    return encode3(exp)


def encode_results(exp):

    # exract file locations

    IPs = exp.IPs

    for experiment in IPs.Condition.unique().tolist():

        exp_dir = f'{exp.scratch}ENCODE3/{experiment}/'
        cromwell_folder = '{}cromwell-executions/chip/*/call-{}/shard-*/execution/'
        non_shard_folder = '{}cromwell-executions/chip/*/call-{}/execution/'
        output(f'Extracting results from {experiment}', log_file=exp.log_file, run_main=exp.run_main)

        samples = exp.IPs[exp.IPs.Condition == experiment].Sample_Name.tolist()

        for sample in samples:

            input_index = int(IPs[IPs.Sample_Name == sample]['Background Sample'])
            input_sample = exp.sample_df.loc[input_index, 'Sample_Name']

            exp.sample_files[sample] = {}
            exp.sample_files[input_sample] = {}

            exp.sample_files[sample]['bam'] = glob_check(f"{cromwell_folder.format(exp_dir, 'bwa')}*{sample}*.bam")
            exp.sample_files[sample]['all_flagstat'] = glob_check(f"{cromwell_folder.format(exp_dir,'bwa')}*{sample}*.flagstat.qc")

            exp.sample_files[sample]['nodup_bam'] = glob_check(f"{cromwell_folder.format(exp_dir,'filter')}*{sample}*.nodup.bam")
            exp.sample_files[sample]['nodup_flagstat'] = glob_check(f"{cromwell_folder.format(exp_dir,'filter')}*{sample}*.nodup.flagstat.qc")
            exp.sample_files[sample]['dup_qc'] = glob_check(f"{cromwell_folder.format(exp_dir,'filter')}*{sample}*.dup.qc")
            exp.sample_files[sample]['pbs_qc'] = glob_check(f"{cromwell_folder.format(exp_dir,'filter')}*{sample}*.pbc.qc")

            exp.sample_files[sample]['jsd_qc'] = glob_check(f"{non_shard_folder.format(exp_dir,'fingerprint')}*{sample}*.jsd.qc")
            exp.sample_files[sample]['fingerprint_plot'] = glob_check(f"{non_shard_folder.format(exp_dir,'fingerprint')}*{sample}*.png")

            exp.sample_files[sample]['xcor_qc'] = glob_check(f"{cromwell_folder.format(exp_dir,'xcor')}*{sample}*.cc.qc")
            exp.sample_files[sample]['xcor_png'] = glob_check(f"{cromwell_folder.format(exp_dir,'xcor')}*{sample}*.png")
            exp.sample_files[sample]['xcor_fraglen'] = glob_check(f"{cromwell_folder.format(exp_dir,'xcor')}*{sample}*.fraglen.txt")

            exp.sample_files[sample]['narrowPeak'] = glob_check(f"{cromwell_folder.format(exp_dir,'macs2')}*{sample}*.bfilt.narrowPeak.gz")
            exp.sample_files[sample]['bw'] = glob_check(f"{cromwell_folder.format(exp_dir,'macs2')}*{sample}*.fc.signal.bigwig")
            exp.sample_files[sample]['frip'] = glob_check(f"{cromwell_folder.format(exp_dir,'macs2')}*{sample}*.frip.qc")

            exp.sample_files[input_sample]['bam'] = glob_check(f"{cromwell_folder.format(exp_dir,'bwa_ctl')}*{input_sample}*.bam")
            exp.sample_files[input_sample]['all_flagstat'] = glob_check(f"{cromwell_folder.format(exp_dir,'bwa_ctl')}*{input_sample}*.flagstat.qc")

            exp.sample_files[input_sample]['nodup_bam'] = glob_check(f"{cromwell_folder.format(exp_dir,'filter_ctl')}*{input_sample}*.nodup.bam")
            exp.sample_files[input_sample]['nodup_flagstat'] = glob_check(f"{cromwell_folder.format(exp_dir,'filter_ctl')}*{input_sample}*.nodup.flagstat.qc")

            # Remove Files
            glob_remove(f"{cromwell_folder.format(exp_dir, '*')}*fastq.gz")
            glob_remove(f"{cromwell_folder.format(exp_dir, '*')}*Align.gz")
            glob_remove(f"{non_shard_folder.format(exp_dir, '*')}*Align.gz")
            glob_remove(f"{cromwell_folder.format(exp_dir, '*_pr*')}*.bam")

        exp.sample_files[experiment] = {}
        exp.sample_files[experiment]['idr_optimal_peak'] = glob_check(f"{non_shard_folder.format(exp_dir,'reproducibility_idr')}optimal_peak.narrowPeak.gz")
        exp.sample_files[experiment]['idr_qc'] = glob_check(f"{non_shard_folder.format(exp_dir,'reproducibility_idr')}idr.reproducibility.qc")
        exp.sample_files[experiment]['overlap_peak'] = glob_check(f"{non_shard_folder.format(exp_dir,'reproducibility_overlap')}optimal_peak.narrowPeak.gz")
        exp.sample_files[experiment]['overlap_qc'] = glob_check(f"{non_shard_folder.format(exp_dir,'reproducibility_overlap')}overlap.reproducibility.qc")
        exp.sample_files[experiment]['qc_report'] = glob_check(f"{non_shard_folder.format(exp_dir,'qc_report')}qc.html")

    exp.tasks_complete.append('encode_results')

    return exp


def spike(exp):
    '''
    If calling from jupyter.  Change backend as needed.

    Align sequencing files to drosophila.
    '''
    import pandas as pd

    if len(exp.spike_samples) == 0:
        output('Not processing Spike-ins', log_file=exp.log_file, run_main=exp.run_main)
        exp.tasks_complete.append('Spike')
        return exp

    # Make QC folder
    spike_folder = make_folder(f'{exp.scratch}spike/')
    output('Processing samples with drosophila-spike in chromatin.', log_file=exp.log_file, run_main=exp.run_main)

    for sample in exp.spike_samples:
        bam = exp.sample_files[sample]['bam']
        spike_sample_folder = make_folder(f'{spike_folder}{sample}/')

        spike_command = [submission_prepend(),
                         f'cd {spike_sample_folder}'
                         f'samtools view -b -f 4 {bam} | samtools sort -n - | samtools fastq - > {spike_sample_folder}{sample}.bwa_unaligned.fastq',
                         f'bowtie2 -p 8 -x {exp.genome_indicies["spike_index"]} -U {spike_sample_folder}{sample}.bwa_unaligned.fastq -S {spike_sample_folder}{sample}.BDGP6.sam --very-sensitive-local -k 1 --no-unal',
                         f'samtools view -b -F 4 {spike_sample_folder}{sample}.BDGP6.sam | samtools sort - > {spike_sample_folder}{sample}.BDGP6.bam',
                         f'picard MarkDuplicates I={spike_sample_folder}{sample}.BDGP6.bam O={spike_sample_folder}{sample}.BDGP6.nodup.bam M={spike_sample_folder}{sample}.BDGP6.nodups.markdups.qc ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true',
                         f'samtools flagstat {spike_sample_folder}{sample}.BDGP6.nodup.bam > {spike_sample_folder}{sample}.unique_drosophila.flagstat.qc',
                         f'rm {spike_sample_folder}{sample}.BDGP6.sam {spike_sample_folder}{sample}.BDGP6.nodup.bam {spike_sample_folder}{sample}*.fastq'
                         ]

        exp.job_id.append(send_job(command_list=spike_command,
                                   job_name=f"{sample}_spike",
                                   job_log_folder=exp.job_folder,
                                   q='general',
                                   mem=10000,
                                   log_file=exp.log_file,
                                   project=exp.project,
                                   cores=2,
                                   run_main=exp.run_main
                                   ))

    # Wait for jobs to finish
    job_wait(exp.job_id, exp.log_file, exp.run_main)

    spike_reads = pd.DataFrame(index=['spike_reads', 'genome_reads'])

    for sample in exp.spike_samples:
        spike_sample_folder = f'{spike_folder}{sample}/'
        qc_file = f'{spike_sample_folder}{sample}.unique_drosophila.flagstat.qc'
        exp.sample_files[sample]['drosophila'] = qc_file

        with open(qc_file, 'r') as fp:
            spike_number = fp.read().split(' ')[0]

        with open(exp.sample_files[sample]['nodup_flagstat']) as fp:
            target_number = fp.read().split(' ')[0]

        spike_reads[sample] = [spike_number, target_number]

    exp.spike_reads = spike_reads.T
    condition_dict = pd.Series(exp.sample_df.Condition.values, index=exp.sample_df.Sample_Name).to_dict()

    exp.spike_reads['Replicate'] = [x.split('_')[-1] for x in exp.spike_reads.index.tolist()]
    exp.spike_reads['Condition'] = [condition_dict[x] for x in exp.spike_reads.index.tolist()]

    for name, spike_conditions in exp.spike_comparisons.items():
        out_dir = make_folder(f'{exp.scratch}spike/{name}')
        plot = spike_in_plot(exp.spike_reads, spike_conditions, name, out_dir)
        out_result(plot, f'{name.replace("_", " ")} Spike-In Comparison', run_main=exp.run_main)
        output(f'Spike-in comparison {name.replace("_", " ")} can be found here: {plot.replace(os.scratch, "")}')

    output(f'Spike-in counts:\n {spike_reads.T}', log_file=exp.log_file, run_main=exp.run_main)

    output('Spike-in alignment jobs finished.', log_file=exp.log_file, run_main=exp.run_main)

    # Generate one dataframe for all spike_counts

    output(f"Spike-in processing complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n", log_file=exp.log_file, run_main=exp.run_main)

    exp.tasks_complete.append('Spike')
    return exp
