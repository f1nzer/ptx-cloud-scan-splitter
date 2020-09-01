#!/usr/bin/env python3
import os
import sys
import getopt

scan_index = 0

def get_scan(reader):
    rows_line = reader.readline()
    if rows_line == '':
        return False

    cols_line = reader.readline()
    rows = int(rows_line)
    cols = int(cols_line)
    print(f'Reading scan #{scan_index}: {rows}x{cols}')

    lines_count = rows * cols + 8
    output_file = os.path.join(output_folder, f'scan_{scan_index:06d}.ptx')
    with open(output_file, 'w') as writer:
        writer.write(rows_line)
        writer.write(cols_line)
        for _ in range(lines_count):
            writer.write(reader.readline())
    
    return True


cmd_help = 'split.py -i <input_file> -o <output_folder>'
try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:", ["ifile=", "ofolder="])
except getopt.GetoptError:
    print(cmd_help)
    exit(2)

for opt, arg in opts:
    if opt == '-h':
        print(cmd_help)
        exit()
    elif opt in ("-i", "--ifile"):
        input_file = arg
    elif opt in ("-o", "--ofolder"):
        output_folder = arg

with open(input_file, 'r') as reader:
    while get_scan(reader):
        scan_index += 1
        exit()
