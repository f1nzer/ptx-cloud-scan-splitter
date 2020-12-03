#!/usr/bin/env python3
import os

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from tempfile import NamedTemporaryFile
from glob import glob
from numpy import array, identity, matrix
from shutil import move
from multiprocessing import Pool
from itertools import islice
from functools import partial

def parse_cmd():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter,
                            usage='%(prog)s {s,split,b,bake} ptx_file')
    parser.add_argument('command',
                        choices=['s', 'split', 'b', 'bake', ],
                        help="Action to perform on input ptx file.")
    parser.add_argument('ptx_file', nargs='+')
    parser.add_argument('target_folder')

    args = parser.parse_args()

    if len(args.ptx_file) == 1:
        args.ptx_file = glob(args.ptx_file[0])
    args.ptx_file.sort()

    ptx_files = []
    for f in args.ptx_file:
        if os.path.splitext(f)[1] == ".ptx":
            ptx_files.append(f)
        else:
            print("Skiping invalid: %s" % f)
    args.ptx_file = ptx_files
    return args


def parse_ptx_header(header_txt):
    """Parse PTX header from list of strings, and return dictionary with header
    components.

    The returned dictionary also contains the value 'matrix' which is the
    header's transformation matrix transposed in to row major order."""
    header = {}
    header['columns'] = int(header_txt[0])
    header['rows'] = int(header_txt[1])
    header['scanner_position'] = [float(v) for v in header_txt[2].split()]
    header['scanner_axis_x'] = [float(v) for v in header_txt[3].split()]
    header['scanner_axis_y'] = [float(v) for v in header_txt[4].split()]
    header['scanner_axis_z'] = [float(v) for v in header_txt[5].split()]
    header['ptx_matrix'] = matrix(b';'.join(header_txt[6:]).decode('utf-8'))

    # Convenience members.
    header['matrix'] = header['ptx_matrix'].T
    header['n_points'] = header['rows'] * header['columns']
    return header


def list_ptx_header(header):
    """Read header dictionary and return header as list of string."""
    header_txt = [b'%i\n' % header['columns'],
                  b'%i\n' % header['rows'],
                  b'%.8f %.8f %.8f\n' % tuple(header['scanner_position']),
                  b'%.8f %.8f %.8f\n' % tuple(header['scanner_axis_x']),
                  b'%.8f %.8f %.8f\n' % tuple(header['scanner_axis_y']),
                  b'%.8f %.8f %.8f\n' % tuple(header['scanner_axis_z'])]
    mat = array(header['ptx_matrix'])
    for row in mat:
        header_txt.append(b'%.8f %.8f %.8f %.8f\n' % tuple(row))

    return header_txt


def copy_data(in_handle, out_handle, n_lines):
    "Copy data from open input file handle to open output file handle."
    for nn in range(n_lines):
        out_handle.write(in_handle.readline())


def create_header(mat, rows, cols):
    """create header using values from row major transformation matrix."""
    ptx_ar = array(mat).T

    return {'rows': rows,
            'columns': cols,
            'n_points': rows * cols,
            'matrix': mat,
            'ptx_matrix': mat.T,
            'scanner_axis_x': ptx_ar[0, 0:3].tolist(),
            'scanner_axis_y': ptx_ar[1, 0:3].tolist(),
            'scanner_axis_z': ptx_ar[2, 0:3].tolist(),
            'scanner_position': ptx_ar[3, 0:3].tolist()}


def split_file_by_scans(ptx_file: str):
    with open(ptx_file, 'rb') as ptx_handle:
        basename = os.path.splitext(os.path.split(ptx_file)[1])[0]
        scan_num = 0
        while True:
            header_txt = []
            while len(header_txt) < 10:
                header_txt.append(ptx_handle.readline())
            if not header_txt[-1]:
                break

            header = parse_ptx_header(header_txt)
            out_file = '%s_%0.3d.ptx' % (basename, scan_num)
            with open(out_file, 'wb') as out_handle:
                out_handle.writelines(list_ptx_header(header))
                copy_data(ptx_handle, out_handle, header['n_points'])

            scan_num += 1


def bake_line(header, raw_line: str):
    line = raw_line.split()
    coords = [float(v) for v in line[:3]]

    if (coords[0] == coords[1] == coords[2] == 0):
        return raw_line

    coords.append(1)
    coords = matrix(coords)
    coords = header['matrix'].dot(coords.T).T.tolist()
    return b'%s %s\n' % (b'%.5f %.5f %.5f' % tuple(coords[0][:3]), b' '.join(line[3:]))
    

def bake_file_by_scans(ptx_file: str, target_folder: str):
    print(ptx_file)
    with open(ptx_file, 'rb') as ptx_handle:
        basename = os.path.splitext(os.path.split(ptx_file)[1])[0]
        scan_num = 0
        while True:
            header_txt = []
            while len(header_txt) < 10:
                header_txt.append(ptx_handle.readline())
            if not header_txt[-1]:
                break

            header = parse_ptx_header(header_txt)
            func = partial(bake_line, header)

            identity_header = create_header(identity(4), header['rows'], header['columns'])
            with NamedTemporaryFile(delete=False) as tmp_handle:
                tmp_handle.writelines(list_ptx_header(identity_header))
                total_points = int(header['n_points'])
                done_points = 0
                while True:
                    progress = round(done_points / total_points * 100)
                    print(f'Scan {scan_num:03d}: {done_points} / {total_points} | {progress}%')

                    batch_size = min(10_000_000, total_points - done_points)
                    next_n_lines = list(islice(ptx_handle, batch_size))
                    if not next_n_lines:
                        break
                    
                    with Pool() as pool:
                        result = pool.map(func, next_n_lines)

                    tmp_handle.writelines(result)
                    done_points += len(next_n_lines)

                print()

            target_baked_file = '%s_%0.3d.ptx' % (basename, scan_num)
            target_baked_filename = os.path.join(target_folder, target_baked_file)
            move(tmp_handle.name, target_baked_filename)

            scan_num += 1


def handle_scan(command: str, ptx_file: str, target_folder: str):
    print("\nProcessing: %s" % ptx_file)
    if command in ['s', 'split']:
        split_file_by_scans(ptx_file, target_folder)
    elif command in ['b', 'bake']:
        bake_file_by_scans(ptx_file, target_folder)


def main():
    print('priv')
    args = parse_cmd()
    handle_scan(args.command, args.ptx_file[0], args.target_folder)


if __name__ == '__main__':
    main()
