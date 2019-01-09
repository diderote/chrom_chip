#!/usr/bin/env python3

import reprlib

class Experiment:
    '''
    Experiment object for pipeline
    '''
    def __init__(self):
        self.tasks_complete = []
        self.job_id = []
        self.overlap_results = {}
        self.anno_results = {}
        self.overlaps = {}
        self.name = 'Experiment'
        self.genome_indicies = {'encode_tsv': {}}
        self.sample_files = {}

    def __repr__(self):
        exclude = ['job_id', 'name']
        experiment = f'{self.__dict__["name"]}'
        for key, value in self.__dict__.items():
            experiment += f'\n\t{key}: {reprlib.repr(value)}' if key not in exclude else ''
        return f'Experiment({experiment})'

__all__ = ['overlaps', 'common','preprocess','parse','pipeline','encode_pipe','diff_binding']
