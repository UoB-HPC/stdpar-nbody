import argparse
import csv
import struct


def read_and_save(input_csv_path, output_binary_path):
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
            rows.append((mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z))

        size = len(rows)
        dimension = 3

    # save ouput
    total_mass = 0
    with open(output_binary_path, 'wb') as binary_file:
        binary_file.write(struct.pack('i', size))
        binary_file.write(struct.pack('i', dimension))
        # timestep
        binary_file.write(struct.pack('f', 3600))
        # constant
        binary_file.write(struct.pack('f', 6.674e-11))

        # write body data
        for row in rows:
            total_mass += row[0]
            binary_file.write(struct.pack('f', row[0]))
            binary_file.write(struct.pack('f', row[1]))
            binary_file.write(struct.pack('f', row[2]))
            binary_file.write(struct.pack('f', row[3]))
            binary_file.write(struct.pack('f', row[4]))
            binary_file.write(struct.pack('f', row[5]))
            binary_file.write(struct.pack('f', row[6]))
    print(f'Total mass saved: {total_mass:.60g}')


def main():
    parser = argparse.ArgumentParser(description='Read n-body CSV and write to binary file.')
    parser.add_argument('input_csv', type=str, help='Input CSV file path')
    parser.add_argument('output_bin', type=str, help='Output binary file path')

    args = parser.parse_args()

    read_and_save(args.input_csv, args.output_bin)


if __name__ == '__main__':
    main()

