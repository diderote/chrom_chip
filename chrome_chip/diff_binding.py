#!/usr/bin/env python3


def diff_binding(exp):
    '''
    Differential Expression using

    Inputs
    ------
    exp.job_folder: '/path/to/job/log/folder'
    exp.log_file: 'log_file.txt'
    exp.scratch: out_dir folder
    exp.alignment_mode: 'gene' or 'transcript'
    exp.designs['colData']: pd.DataFrame()
    exp.designs['all_samples']: ['list','of','samples']
    exp.designs['design']: ex. '~ConditionA'
    exp.norm: 'ERCC','Median-Ratios','ERCC_Mixed','Empirical'
    exp.tasks_complete: []
    exp.count_matrix: pd.DataFrame()

    Optional
    --------
    exp.gc_norm: bool
    exp.gc_count_matrix
    exp.spike_counts
    exp.genome_indicies['ERCC_Mix']

    Outputs
    -------
    exp.de_results: {}
    saves file to exp.scratch/DESeq2_Results/...

    '''

    # output('Beginning DESeq2 differential expression analysis: {:%Y-%m-%d %H:%M:%S}\n'.format(datetime.now()), exp.log_file)

    return exp
