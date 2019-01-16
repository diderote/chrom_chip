#!/usr/bin/env python3

import json
import shutil
import re
import os
from datetime import datetime

from chrome_chip.common import output, send_job, job_wait, job_pending, glob_check, make_folder, glob_remove
from chrome_chip.preprocess import stage


def encode3(exp):

    if 'Stage' not in exp.tasks_complete:
        output('Files not staged.\n', exp.log_file)
        exp = stage(exp)

    output('Running alignment and peak calling using ENCODE3 standards.', exp.log_file)
    output('ENCODE3 cromwell pipeline.', exp.log_file)

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
        UMI = True if UMI_list[0].lower() == 'yes' else False

        try:
            file_type = end_types[exp.sample_df[exp.sample_df.Condition == experiment]['Scratch_File1'].tolist()[0][-4:]]
        except KeyError:
            output(f"{exp.sample_df[exp.sample_df.Condition == experiment]['Scratch_File1'].tolist()[0]} not a valid file type for this pipeline.")

        genome = IPs[IPs.Condition == experiment]['Genome'].unique().tolist()
        if len(genome) > 1:
            raise IOError('Cannot align to more than one genome per condition.')

        chip_type = IPs[IPs.Condition == experiment]['ChIP Type'].unique().tolist()
        if len(chip_type) > 1:
            raise IOError('Cannot have more than one chip type (histone or TF) for a condition.')
        chip_type = 'Histone' if chip_type[0].lower() == 'histone' else 'tf'

        json_file = {}
        json_file['chip.pipeline_type'] = chip_type
        json_file['chip.paired_end'] = seq_type
        json_file['chip.genome_tsv'] = exp.genome_indicies['encode_tsv'][genome[0]]
        json_file['chip.bwa.mem_mb'] = 30000
        json_file['chip.macs2_mem_mb'] = 30000
        json_file['chip.peak_caller'] = 'macs2'
        json_file["chip.true_rep_only"] = False
        json_file["chip.dup_marker"] = "picard"
        json_file["chip.mapq_thresh"] = 30
        json_file["chip.regex_filter_reads"] = "chrM"
        json_file["chip.subsample_reads"] = 0
        json_file["chip.ctl_subsample_reads"] = 0
        json_file["chip.xcor_subsample_reads"] = 15000000
        json_file["chip.keep_irregular_chr_in_bfilt_peak"] = False
        json_file["chip.always_use_pooled_ctl"] = False
        json_file["chip.ctl_depth_ratio"] = 1.2
        json_file["chip.macs2_cap_num_peak"] = 500000
        json_file["chip.pval_thresh"] = 0.01
        json_file["chip.idr_thresh"] = 0.05
        json_file["chip.bwa_cpu"] = 4
        json_file["chip.bwa_mem_mb"] = 20000
        json_file["chip.bwa_time_hr"] = 48
        json_file["chip.filter_cpu"] = 2
        json_file["chip.filter_mem_mb"] = 20000
        json_file["chip.filter_time_hr"] = 24
        json_file["chip.bam2ta_cpu"] = 2
        json_file["chip.bam2ta_mem_mb"] = 10000
        json_file["chip.bam2ta_time_hr"] = 6
        json_file["chip.fingerprint_cpu"] = 2
        json_file["chip.fingerprint_mem_mb"] = 12000
        json_file["chip.fingerprint_time_hr"] = 6
        json_file["chip.xcor_cpu"] = 2
        json_file["chip.xcor_mem_mb"] = 16000
        json_file["chip.xcor_time_hr"] = 24
        json_file["chip.macs2_time_hr"] = 24
        json_file["chip.spr_mem_mb"] = 16000

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
        json_file['chip.align_only'] = True if final_stage == 'align' else False

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

        command_list = ['module rm python share-rpms65',
                        'source activate encode-chip-seq-pipeline',
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
                            cores=1
                            )

        exp.job_id.append(sent_job)
        job_pending(sent_job, exp.log_file)

    # Wait for jobs to finish
    job_wait(exp.job_id, exp.log_file)

    exp = encode_results(exp)

    exp.tasks_complete.append('ENCODE3')

    return exp


def UMI(exp):

    output('Deduplicating bam files using UMIs with UMI-tools.', exp.log_file)

    # exp.data_type = 'bam'

    out_dir = make_folder(f'{exp.scratch}UMI/')

    IPs = exp.IPs

    for experiment in IPs.Condition.unique().tolist():
        UMI = True if 'yes' in IPs[IPs.Condition == experiment]['UMI'].tolist() else False

        if not UMI:
            return exp

        else:
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

                command_list = ['module rm python share-rpms65',
                                'source activate chrome_chip',
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
                                           cores=1
                                           ))

                exp.sample_files[sample]['nodup_bam'] = nodup_bam
                exp.sample_files[input_sample]['nodup_bam'] = nodup_input

    job_wait(exp.job_id, exp.log_file)

    output('Dedplication complete.  Submitting deduplicated files for the remainder of processing.', exp.log_file)
    exp.tasks_complete.append('UMI')

    return encode3(exp)


def encode_results(exp):

    # exract file locations

    IPs = exp.IPs

    for experiment in IPs.Condition.unique().tolist():

        exp_dir = f'{exp.scratch}ENCODE3/{experiment}/'
        cromwell_folder = '{}cromwell-executions/chip/*/call-{}/shard-*/execution/'
        non_shard_folder = '{}cromwell-executions/chip/*/call-{}/execution/'
        output(f'Extracting results from {experiment}', exp.log_file)

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
        exp.sample_files[experiment]['idr_optimal_peak'] = glob_check(f"{non_shard_folder.format(exp_dir,'reproducibility_idr')}*{sample}optimal_peak.narrowPeak.gz")
        exp.sample_files[experiment]['idr_qc'] = glob_check(f"{non_shard_folder.format(exp_dir,'reproducibility_idr')}*{sample}idr.reproducibility.qc")
        exp.sample_files[experiment]['overlap_peak'] = glob_check(f"{non_shard_folder.format(exp_dir,'reproducibility_overlap')}optimal_peak.narrowPeak.gz")s
        exp.sample_files[experiment]['overlap_qc'] = glob_check(f"{non_shard_folder.format(exp_dir,'reproducibility_overlap')}overlap.reproducibility.qc")
        exp.sample_files[experiment]['qc_report'] = glob_check(f"{non_shard_folder.format(exp_dir,'qc_report')}qc.html")

    exp.tasks_complete.append('encode_results')

    return exp


def spike(exp):
    '''
    If calling from jupyter.  Change backend as needed.

    Align sequencing files to drosophila.
    '''

    output('Processing samples with drosophila-spike in chromatin.', exp.log_file)
    # Make QC folder
    spike_folder = make_folder(f'{exp.scratch}spike/')

    spike_list = [sample for sample in exp.samples if 'none' not in exp.IPs.loc[exp.IPs.Sample_Name == sample, 'Spike_Comparison'].tolist()]

    for sample in spike_list:
        bam = exp.sample_files[sample]['bam']

        spike_command = ['module rm python share-rpms65',
                         'source activate chrome_chip',
                         f'samtools view -b -f 4 {bam} | samtools sort -n - | samtools fastq - > {spike_folder}{sample}.bwa_unaligned.fastq',
                         f'bowtie2 -p 8 -x {exp.genome_indicies["spike_index"]} -U {spike_folder}{sample}.bwa_unaligned.fastq -S {spike_folder}{sample}.BDGP6.sam --very-sensitive-local -k 1 --no-unal',
                         f'samtools view -b -F 4 {spike_folder}{sample}.BDGP6.sam | samtools sort - > {spike_folder}{sample}.BDGP6.bam',
                         f'picard MarkDuplicates I={spike_folder}{sample}.BDGP6.bam O={spike_folder}{sample}.BDGP6.nodup.bam M={spike_folder}{sample}BDGP6.nodups.markdups.qc ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true',
                         f'samtools sort -n {spike_folder}{sample}.BDGP6.nodup.bam | samtools fastq - > {spike_folder}{sample}.BDGP6.nodup.fastq',
                         f'bowtie2 -p 4 -x exp.genome_indicies["genome_bt2_index"] -U {spike_folder}{sample}.BDGP6.nodup.fastq -S {spike_folder}{sample}.hg38.BDGP6.sam --very-sensitive-local -k 1',
                         f'samtools view -f 4 {spike_folder}{sample}.hg38.BDGP6.sam | samtools flagstat - > {spike_folder}{sample}.unique_drosophila.flagstat.qc',
                         f'rm {spike_folder}{sample}.BDGP6.sam {spike_folder}{sample}.BDGP6.bam {spike_folder}{sample}.BDGP6.nodup.bam {spike_folder}{sample}.BDGP6.nodup.fastq {spike_folder}{sample}.hg38.BDGP6.sam'
                         ]

        exp.job_id.append(send_job(command_list=spike_command,
                                   job_name=f"{sample}_spike",
                                   job_log_folder=exp.job_folder,
                                   q='general',
                                   mem=10000,
                                   log_file=exp.log_file,
                                   project=exp.project,
                                   cores=2
                                   ))

    # Wait for jobs to finish
    job_wait(exp.job_id, exp.log_file)

    for sample in exp.samples.values():
        exp.sample_files[sample]['drosophila'] = f'{spike_folder}{sample}.unique_drosophila.flagstat.qc'

    output('Spike-in alignment jobs finished.', exp.log_file)

    # Generate one dataframe for all spike_counts

    output(f"Spike-in processing complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n", exp.log_file)

    exp.tasks_complete.append('Spike')
    return exp
