import argparse
import csv
import struct

import numpy as np

def calc_constant():
    G_SI = 6.67428 * np.power(10.0, -11)

    # conversion factors
    meter_AU = 1.0 / (1.49597870691 * np.power(10.0, 11))
    second_days = 1.0 / 86400

    # scale G appropriately
    G_scaled = G_SI * (np.power(meter_AU, 3) / np.power(second_days, 2))

    return float(G_scaled)


def read_and_save(args):
    input_csv_path, output_binary_path = args.input_csv, args.output_bin
    skip_count = 0
    with open(input_csv_path, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)

        # skip header
        next(csv_reader)

        rows = []
        for row in csv_reader:
            mass = float(row[3])
            pos_x = float(row[4])
            pos_y = float(row[5])
            pos_z = float(row[6])
            vel_x = float(row[7])
            vel_y = float(row[8])
            vel_z = float(row[9])
            out_row = np.array((mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z))
            if not np.any(np.isnan(out_row)):
                rows.append(out_row)
            else:
                skip_count += 1

        size = len(rows)
        dimension = 3

    # save ouput
    total_mass = 0
    with open(output_binary_path, 'wb') as binary_file:
        binary_file.write(struct.pack('i', size))
        binary_file.write(struct.pack('i', dimension))
        # the simulation will be done in hours and AU
        # timestep
        binary_file.write(struct.pack('f', 1 / 24))
        binary_file.write(struct.pack('f', calc_constant()))

        # write body data
        row_count = int(args.prop * len(rows))
        for row in rows[:row_count]:
            total_mass += row[0]
            binary_file.write(struct.pack('f', row[0]))
            binary_file.write(struct.pack('f', row[1]))
            binary_file.write(struct.pack('f', row[2]))
            binary_file.write(struct.pack('f', row[3]))
            binary_file.write(struct.pack('f', row[4]))
            binary_file.write(struct.pack('f', row[5]))
            binary_file.write(struct.pack('f', row[6]))
    print(f'Saved {row_count} bodies')
    print(f'Total mass saved: {total_mass:.60g}')
    print(f'Skipped {skip_count} bodies')


def main():
    parser = argparse.ArgumentParser(description='Read n-body CSV and write to binary file.')
    parser.add_argument('input_csv', type=str, help='Input CSV file path')
    parser.add_argument('output_bin', type=str, help='Output binary file path')
    parser.add_argument('--prop', required=False, default=1, type=float, help='Proportion of dataset to use (default is 1)')

    args = parser.parse_args()

    read_and_save(args)


if __name__ == '__main__':
    main()

