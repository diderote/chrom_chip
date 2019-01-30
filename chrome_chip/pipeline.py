#!/usr/bin/env python3

from chrome_chip.parse import parse_config
from chrome_chip.common import validated_run
from chrome_chip.preprocess import stage, fastq_screen, trim, fastqc
from chrome_chip.encode_pipe import encode3, UMI, spike
from chrome_chip.qc import principal_component_analysis, preseq, final_qc, finish
from chrome_chip.overlaps import overlaps, annotation, heatmaps
# from chrom_chip.diff_binding import diff_binding


def pipeline(experimental_file):
        exp = parse_config(experimental_file, run_main=True)
        exp = validated_run('Stage', stage, exp)
        exp = validated_run('Fastq_screen', fastq_screen, exp)
        exp = validated_run('Trim', trim, exp)
        exp = validated_run('FastQC', fastqc, exp)
        exp = validated_run('ENCODE3', encode3, exp)
        exp = validated_run('UMI', UMI, exp)
        exp = validated_run('Spike', spike, exp)
        exp = validated_run('preseq', preseq, exp)
        exp = validated_run('PCA', principal_component_analysis, exp)
        exp = validated_run('Overlaps', overlaps, exp)
        # exp = validated_run('Diff_Bind', diff_binding, exp)
        exp = validated_run('Heatmaps', heatmaps, exp)
        exp = validated_run('Annotations', annotation, exp)
        exp = validated_run('MultiQC', final_qc, exp)
        exp = validated_run('Finished', finish, exp)
