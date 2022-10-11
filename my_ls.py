"""This script uses python's glob library to list contents of a directory of a given type. It is similar to ls command of unix shell."""

import argparse
import glob

def my_list(args):
    """Run the program"""
    dir = args.dir if args.dir[-1] == '/' else args.dir + '/'
    glob_input = dir + '*.' + args.suffix
    glob_output = sorted( glob.glob(glob_input) )
    for item in glob_output:
        print(item)
    


if __name__ == "__main__":
    USAGE = "my_ls.py [-h] dir suffix. List the file names in a given directory."
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dir', type=str, help= 'Input Directory Name')
    parser.add_argument('suffix', type=str, help="File suffix (eg. py, sh)")
    args = parser.parse_args()
    my_list(args)

