#!../venv/bin/python

import subprocess
import argparse
import os
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument('beamwidth')
parser.add_argument('sizes')
args = parser.parse_args()
[b_min, b_num_points, b_max] = args.beamwidth.split(':')
[s_min, s_num_points, s_max] = args.sizes.split(':')

dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)
script = os.path.join(dir_path, "start_sim")
os.environ["SIM_DIR"] = dir_path

points = np.linspace(float(b_min), float(b_max), num=int(b_num_points))
for point in points:
    proc = subprocess.Popen([script, '-bfng', '-p', f'150,{point}',
                             f'{s_min}:{s_num_points}:{s_max}'],
                            stdout=subprocess.PIPE)
    while True:
        line = proc.stdout.readline().decode()
        if 'END' in line:
            break
    proc.wait()

