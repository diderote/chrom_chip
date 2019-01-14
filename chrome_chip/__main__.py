#!/usr/bin/env python3

import argparse
import os
import sys

import papermill as pm

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from chrome_chip.pipeline import pipeline
from chrome_chip.common import send_job

parser = argparse.ArgumentParser()
parser.add_argument('--experimental_file', '-f', required=True, help='experimental yaml file', type=str)
parser.add_argument('--notebook', '-n', dest='notebook', help='run as jupyter notebook', action='store_true')
parser.add_argument('--no-notebook', dest='notebook', help='run as python script', action='store_false')
parser.add_argument('--template_notebook', '-t', required=False, help='location of template notebook', type=str)
parser.add_argument('--out_notebook', '-o', required=False, help='name of output notebook', type=str)
parser.add_argument('--submit', '-s', required=False, help='true will submit to LSF with conda intialization', action='store_true')
parser.add_argument('--project', '-p', required=False, help='LSF project name for submission', type=str)
parser.set_defaults(notebook=True, submit=False)
args = parser.parse_args()

if args.submit:
        cmd = f'python chrom_chip -f {args.experimental_file}'
        if not args.notebook:
            cmd += ' --no-notebook'
        else:
            cmd += f' -t {args.template_notebook}'
            if args.out_notebook:
                cmd += f' -o {args.out_notebook}'

        submission_header = ['module rm python share-rpms65', 'source activate chrom_chip', cmd]

        job_name = args.experimental_file.split('.')[0]

        send_job(command_list=submission_header,
                 job_name=job_name,
                 job_log_folder=f'{os.getcwd()}/',
                 q='general',
                 mem=3000,
                 log_file=f'bsub_{job_name}.log',
                 project=args.project,
                 cores=1,
                 submit=True)

else:
    if args.notebook:
        if (os.path.isfile(args.template_notebook) is False) or (args.template_notebook is False):
            raise IOError(f'Location of template notebook not found. Use -t option.')
        else:
            out_notebook = args.out_notebook if args.out_notebook else args.experimental_file.replace('yml', 'ipynb')
            pm.execute_notebook(args.template_notebook, out_notebook, parameters=dict(yaml_file=args.experimental_file), perpare_only=True)
            pm.execute_notebook(args.template_notebook, out_notebook, parameters=dict(yaml_file=args.experimental_file), log_output=True, report_mode=True)
    else:
        pipeline(args.experimental_file)
