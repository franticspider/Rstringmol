"""oee_utils.py
   @author: Susan Stepney
   utils to find and read data files, for 
   "On the Open-Endendness of Detecting Open-Endedness"
   by Stepney and Hickinbotham, Artificial Life, 2023"""

import math
import pandas as pd

RUNDIR = 'rundata/runs/'    # location of run1 .. run8 data subdirs
TSTEP = 20000               # timestep
NRECORD = [101, 101, 101, 101, 101, 101, 101, 101,    # number of record points, runs 1-8
           9, 3, 6, 4, 7, 18, 5, 3, 2, 5, 7, 5           # runs 9-20
           ]


def get_reaction_filename(run, time):
    # filename = rundata/run<n|nn>/oeeRun<n>props<time>.csv
    # <run> = run number, 1--20
    # <time> = timestep, 0020000 - 2000000 in steps of 20000
    fn_str = RUNDIR + 'run' + str(run) + '/oeeRun' + str(run) + 'props{0:07d}'.format(time) + '.csv'
    return fn_str


def get_species_count_filename(run, time):
    # filename = rundata/run<n|nn>/sppcounts<time>.csv
    # <run> = run number, 1--20
    # <time> = timestep, 000000 - 980000, 1000000 - 2000000 in steps of 2000000
    if time < 1000000:
        fn_str = RUNDIR + 'run' + str(run) + '/sppcounts{0:06d}'.format(time) + '.csv'
    else:
        fn_str = RUNDIR + 'run' + str(run) + '/sppcounts{0:07d}'.format(time) + '.csv'
    return fn_str


def get_species_filename(run, time):
    # filename = rundata/run<n|nn>/species<nn>.csv
    # <run> = run number, 1--20
    fn_str = RUNDIR + 'run' + str(run) + '/species{0:02d}'.format(run) + '.csv'
    return fn_str


def get_selfself_filename(run):
    # filename = rundata/run<n|nn>/species<nn>.csv
    # <run> = run number, 1--20
    fn_str = RUNDIR + 'run' + str(run) + '/selfself{0:02d}'.format(run) + '.csv'
    return fn_str


def get_species_filename(run):
    # filename = rundata/run<n|nn>/species<nn>.csv
    # <run> = run number, 1--20
    fn_str = RUNDIR + 'run' + str(run) + '/species{0:02d}'.format(run) + '.csv'
    return fn_str


def get_qnn_filename(run):
    # filename = rundata/run<n|nn>/qnn<nn>.csv
    # <run> = run number, 1--20
    fn_str = RUNDIR + 'run' + str(run) + '/qnn{0:02d}'.format(run) + '.csv'
    return fn_str


def get_species_number_from_id(str):
    # format of str = 'sp<n>'
    return int(str[2:])  # remove 'sp', convert to number


def read_reaction_logfile(run, time, colnames):
    # run = run number, 1--20
    # time = timestep, 20000 - 2000000 in steps of 2000000
    # colnames = list of strings, column names in csv file to read

    fn_str = get_reaction_filename(run, time)

    # use skipinitialspace in case there are leading spaces in the column names in the csv file
    df = pd.read_csv(fn_str, skipinitialspace=True, usecols=colnames)

    return df


def read_qnn_logfile(run, colnames):
    # run = run number, 1--20
    # colnames = list of strings, column names to apply to unlabelled csv file

    fn_str = get_qnn_filename(run)

    # use skipinitialspace in case there are leading spaces in the column names in the csv file
    df = pd.read_csv(fn_str, skipinitialspace=True, names=colnames)

    return df


def read_species_logfile(run, colnames):
    # run = run number, 1--20
    # colnames = list of strings, column names in csv file to read

    fn_str = get_species_filename(run)

    # use skipinitialspace in case there are leading spaces in the column names in the csv file
    df = pd.read_csv(fn_str, skipinitialspace=True, usecols=colnames)

    return df


def read_species_count_logfile(run, time, colnames):
    # run = run number, 1--20
    # time = timestep, 000000 - 980000, 1000000 - 2000000 in steps of 2000000
    # colnames = list of strings, column names in csv file to read

    fn_str = get_species_count_filename(run, time)

    # use skipinitialspace in case there are leading spaces in the column names in the csv file
    df = pd.read_csv(fn_str, skipinitialspace=True, usecols=colnames)

    if 'n' in colnames:  # 'n' is the count of the number of obs -- occasionally zero
        df = df.loc[df['n'] != 0]  # select all the rows with non-zero counts

    return df


def read_selfself_logfile(run, colnames):
    # run = run number, 1--20
    # colnames = list of strings, column names in csv file to read

    fn_str = get_selfself_filename(run)

    # use skipinitialspace in case there are leading spaces in the column names in the csv file
    df = pd.read_csv(fn_str, skipinitialspace=True, usecols=colnames)

    return df


def read_species_logfile(run, colnames):
    # run = run number, 1--6
    # colnames = list of strings, column names in csv file to read

    fn_str = get_species_filename(run)

    # use skipinitialspace in case there are leading spaces in the column names in the csv file
    df = pd.read_csv(fn_str, skipinitialspace=True, usecols=colnames)

    return df


def read_number_of_species(run):
    # run = run number, 1--8

    fn_str = get_species_filename(run)

    # read last line of file
    with open(fn_str, 'r') as f:
        last_line = f.readlines()[-1].strip().split(',')
    n = get_species_number_from_id(last_line[0])  # first entry of line is 'sp<n>'

    return n

