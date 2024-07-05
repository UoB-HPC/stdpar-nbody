import argparse
import csv
import os
import sys

import numpy as np

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.append(parent_dir)
from plotter import read_points


def load_last_state(file_path):
    data = []

    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # skip the header row
        for row in reader:
            data.append([float(value) for value in row])

    positions = np.array(data)

    return positions


def load_barnes_hut_bin(file_path):
    points = read_points(file_path)

    # only load the last positions
    return points[-1, ...]


def get_state(file_path):
    if file_path.endswith('.bin'):
        return load_barnes_hut_bin(file_path)
    elif file_path.endswith('.csv'):
        return load_last_state(file_path)
    else:
        raise ValueError("Unknown file type")


def main(args):
    state1 = get_state(args.f1).swapaxes(0, 1)
    state2 = get_state(args.f2)

    error = np.abs(state1 - state2).sum()

    print(f'Total absolute difference is {error}')
    print(f'Absolute difference per body is {error / state1.shape[0]}')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=("Compare the output of different n-body simulators\n"
                     "It is assumed that\n"
                     "- .bin is Barnes-Hut simulation\n"
                     "- .csv is N_Body_Simulation\n"
                     )
    )

    parser.add_argument('f1', type=str, help="Path to a state file")
    parser.add_argument('f2', type=str, help="Path to a state file")

    main(parser.parse_args())
