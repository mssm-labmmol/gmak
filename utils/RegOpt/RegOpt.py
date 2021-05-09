#!/usr/bin/env python3

from regopt.regopt import OptimizationDriver
import sys

def read_args(argv):
    """
    Returns
    -------
    args_dict : dictionary with the following key-value pairs
        'csv_files' : list of csv_files, one for each property
             'refs' : list of reference values for properties
          'weights' : list of weights for properties
    """
    args_dict = dict()
    args_dict['test_size'] = 0.25
    for i, arg in enumerate(argv):
        if (arg == '-f'):
            args_dict['csv_files'] = list()
            for subarg in argv[i+1:]:
                if (subarg[0] == '-'):
                    break
                else:
                    args_dict['csv_files'].append(subarg)
        if (arg == '-w'):
            args_dict['weights'] = list()
            for subarg in argv[i+1:]:
                if (subarg[0] == '-'):
                    break
                else:
                    args_dict['weights'].append(float(subarg))
        if (arg == '-r'):
            args_dict['refs'] = list()
            for subarg in argv[i+1:]:
                if (subarg[0] == '-'):
                    break
                else:
                    args_dict['refs'].append(float(subarg))
        if (arg == '-ts'):
            args_dict['test_size'] = float(argv[i+1])
    return args_dict

if __name__ == '__main__':

    if (len(sys.argv) == 1):
        print(f"usage: {sys.argv[0]} -f CSV_FILE... -r REF_VALUE... -w WEIGHT... [-ts TEST_SIZE]")
        exit(1)

    args_dict = read_args(sys.argv)

    driver = OptimizationDriver(test_size=args_dict['test_size'])

    for csv_file, ref, wei in zip(args_dict['csv_files'], args_dict['refs'], args_dict['weights']):
        driver.addPropertyFromCSV(csv_file, wei, ref)

    driver.fit()

    driver.optimize()

    driver.printResults()
