#!/usr/bin/env python3

from chrom_chip.parse import parse_config
from chrom_chip.common import validated_run
from chrom_chip.preprocess import stage, fastq_screen, trim, fastqc
from chrom_chip.encode_pipe import encode3, encode_results, UMI, spike
from chrom_chip.qc import principal_component_analysis, preseq, final_qc, finish
from chrom_chip.overlaps import overlaps, annotation, heatmaps
# from chrom_chip.diff_binding import diff_binding


def pipeline(experimental_file):
        exp = parse_config(experimental_file)
        exp = validated_run('Stage', stage, exp)
        exp = validated_run('Fastq_screen', fastq_screen, exp)
        exp = validated_run('Trim', trim, exp)
        exp = validated_run('FastQC', fastqc, exp)
        exp = validated_run('ENCODE3', encode3, exp)
        exp = validated_run('UMI', UMI, exp)
        exp = validated_run('encode_results', encode_results, exp)
        exp = validated_run('Spike', spike, exp)
        exp = validated_run('preseq', preseq, exp)
        exp = validated_run('PCA', principal_component_analysis, exp)
        exp = validated_run('Overlaps', overlaps, exp)
        # exp = validated_run('Diff_Bind', diff_binding, exp)
        exp = validated_run('Heatmaps', heatmaps, exp)
        exp = validated_run('Annotations', annotation, exp)
        exp = validated_run('MultiQC', final_qc, exp)
        exp = validated_run('Finished', finish, exp)
