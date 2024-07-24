#!/usr/bin/env python3

import argparse,os,pathlib,shutil,subprocess,math,sys,platform

def file_path(arg):
    arg = os.path.expanduser(arg)
    if os.path.isfile(arg): return arg
    else: raise RuntimeError(f'Not a file path: {arg}')

p = argparse.ArgumentParser(description='Parse benchmark results')
p.add_argument('out_file', type=file_path, help='Path to output file')

args = p.parse_args()

# Read file
parsed_gpu = 0
found = False
gpu = None
driver = None
cpu = None
cores = None
compiler = None
hostname = None
sequential = False
with open(args.out_file, 'r') as f:
    for l in f.readlines():
        if l.startswith('+'):
            continue
        elif parsed_gpu == 1:
            gpu = l.split(', ')[0].strip()
            driver = l.split(', ')[1].strip()
            parsed_gpu = 2
            continue
        elif l.startswith('Vendor'):
            continue
        elif l.startswith('Model name:'):
            cpu = l.split('Model name:')[1].strip()
            continue
        elif l.startswith('Core'):
            cores = l.split('Core(s) per socket:')[1].strip()
            continue
        elif l.startswith('name'):
            parsed_gpu = 1
            continue
        elif l.startswith('sequential'):
            sequential = True
            continue
        elif l.startswith('algorithm'):
            if found is True:
                continue
            else:
                found = True
                print(f'gpu,driver,cpu,#cores,seq,compiler,hostname,{l.strip()}')
                continue
        elif l.startswith('compiler'):
            compiler = l.split(':')[1].strip()
        elif l.startswith('hostname') or l.startswith('node'):
            hostname = l.split(':')[1].strip().replace('.nvidia.com', '')
        elif l.split(',')[0] in ['barnes-hut','all-pairs','all-pairs-collapsed','hilbert-tree']:
            print(f'{gpu},{driver},{cpu},{cores},{sequential},{compiler},{hostname},{l.strip()}')
