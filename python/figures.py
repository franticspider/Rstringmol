#!/usr/bin/env python3

"""figure.py
   @author: Susan Stepney
   Produces all plots and table data for the paper 
   "On the Open-Endendness of Detecting Open-Endedness"
   by Stepney and Hickinbotham, Artificial Life, 2023"""

import math
import numpy as np

import oee_utils
import oee_plot


def qnn_plot(start, end):
    # QNN and cumulative QNN plots
    # start, end run number

    # qnn : no of each species at each time step
    # format of ID = 'spnnn'
    species_colnames = ['t', 'ID', 'qnn']
    all_ys = []
    allcum_ys = []

    for run in range(start, end + 1):
        T = oee_utils.NRECORD[run - 1]  # number of record points (including 0)
        ys = np.zeros(T)
        df = oee_utils.read_qnn_logfile(run, species_colnames)

        # populate qnn array with incremental qnn values
        for index, data_row in df.iterrows():  # each species
            t = data_row['t']  # execution timestep
            q = data_row['qnn']  # qnn of this species at this time
            ys[round(t / oee_utils.TSTEP)] += q  # convert execution time to (integer) logpoint time

        tot_qnn = sum(ys)
        print(run, 'total QNN =', tot_qnn)  # total qnn
        all_ys.append(ys)  # qnn, for small multiples grid plot
        allcum_ys.append(np.cumsum(ys))  # cumulative sum of qnn, for small multiples grid plot

    oee_plot.plot_grid(all_ys, 'all_qnn.pdf', ymax=0.35)
    oee_plot.plot_grid(allcum_ys, 'all_qnncum.pdf', ymax=2.9)


def pop_size():

    t0 = 0  # first time
    t1 = 2000000  # last time
    # t1 = 80000  # last time, for testing
    T = 101  # number of timesteps (including 0)

    # abundance : no of strings
    species_colnames = ['n', 'unb', 'act', 'pas']
    all_ys = []
    popln = []

    for run in range(1, 9):
        ys = np.zeros(T, dtype=int)

        # populate species array with number of times str is observed
        for i, time in enumerate(range(t0, t1 + 1, oee_utils.TSTEP)):
            # read in data about strings
            df = oee_utils.read_species_count_logfile(run, time, species_colnames)

            for index, data_row in df.iterrows():  # each species at this time
                ys[i] += data_row['n']  # all strings of species

        all_ys.append(ys)  # for small multiples grid plot
        ys = ys.tolist()
        print(run, '&', min(ys[2:51]), '&', min(ys[51:101]), '&', max(ys[2:51]),
              '&', max(ys[51:101]), '\\\\')  # min, max pop size over run
        popln.append(ys[2:])  # ignore t=0,1

    oee_plot.plot_grid(all_ys, 'all_pop_size.pdf', ymax=13000, yline=12500)
    oee_plot.plot_box(popln, 'all_pop_size_box.pdf', yticks=5)


def species_count():

    t0 = 0  # first time
    t1 = 2000000  # last time
    # t1 = 80000  # last time, for testing
    T = 101  # number of timesteps (including 0)

    # diversity : no of string types
    species_colnames = ['n']
    all_ys = []
    popln = []

    for run in range(1, 9):
        ys = np.zeros([2, T], dtype=int)  # 0 = total species; 1 = new species
        df_species = oee_utils.read_species_logfile(run, ['start'])

        for i, time in enumerate(range(t0, t1 + 1, oee_utils.TSTEP)):
            # read in data about strings: want number of (non-zero) rows
            df = oee_utils.read_species_count_logfile(run, time, species_colnames)
            ys[0, i] = len(df.index)

            # count number of rows in df_species that start at this time: new species
            ys[1, i] = df_species[df_species['start'] == time].count()

        all_ys.append(ys.T)  # for small multiples grid plot
        ys = ys[0, :].tolist()  # extract just the total species row
        print(run, '&', min(ys[2:51]), '&', min(ys[51:101]), '&', max(ys[2:51]),
              '&', max(ys[51:101]), '\\\\')  # min, max pop size over run
        popln.append(ys[2:])  # ignore t=0,1

    oee_plot.plot_grid(all_ys, 'all_species_count.pdf', ymax=1510)
    oee_plot.plot_box(popln, 'all_species_count_box.pdf', ymax=1510)


def species_12abundance():
    # 12 species with most individuals -- plots and tables

    t0 = 0  # first time
    t1 = 2000000  # last time
    # t1 = 80000  # last time, for testing
    T = 101  # number of timesteps (including 0)

    # abundance : no of each string type
    # format of ID = 'spnnn'
    species_colnames = ['ID', 'n']
    all_ys = []

    for run in range(1, 9):

        N = oee_utils.read_number_of_species(run)
        species = np.zeros([N, T], dtype=int)

        for i, time in enumerate(range(t0, t1 + 1, oee_utils.TSTEP)):
            # read in data about strings
            df = oee_utils.read_species_count_logfile(run, time, species_colnames)

            # populate species array with number of times str is observed
            for index, data_row in df.iterrows():  # each species
                id = oee_utils.get_species_number_from_id(data_row['ID'])  # species number
                n = data_row['n']  # number of observations of this species
                species[id - 1, i] += n  # species IDs start at 1

        print('run =', run, '; number of unique species =', species.shape[0])  # number of observed species
        # sort the species into descending order of max peak, to get consistent colours whenplotting thresholded species
        smax = np.amax(species, axis=1)  # maximum of each species
        perm = np.argsort(smax)  # permutation that sorts in ascending order
        perm = perm[::-1]  # reverse the permutation, for sorting in descending order
        species = species[perm, :]  # permute rows, so get same colours when plotting thresholded plots -- loses str ids

        # this stacked plot sloooooowly produces 2-4MB pdfs!
        oee_plot.plot_stacked(species, 'run' + str(run) + '_species_number_stacked.pdf', ymax=12500)

        thresh = 12  # number of peaks tp plots
        species = species[:thresh]  # top thresh peaks
        all_ys.append(species.T)
        species_12abundance_table(run, species, perm)

    oee_plot.plot_grid(all_ys, 'species_number_top_' + str(thresh) + '.pdf', ymax=6100)


def species_12abundance_table(run, species, perm):

    # extract more data about 12 highest peaks
    species_colnames = ['ID', 'seq']
    df_seq = oee_utils.read_species_logfile(run, species_colnames)

    peak_data = []
    for i, row in enumerate(species):  # iterate over rows of final array (12 peaks)
        peak = row.tolist()  # highest peak
        peak_nozero = [n for n in peak if n != 0]  # just the non-zero abundances
        start = peak.index(peak_nozero[0])  # index of start of peak (first non-zero occurrence)
        max_peak = max(peak_nozero)
        maxi = peak.index(max_peak)  # index of peak max
        peak.reverse()
        end = peak.index(peak_nozero[-1])  # index of end of reversed peak (last occurrence)
        end = len(peak) - end - 1
        id = perm[i] + 1    # species id (start at 1)
        seq_row = df_seq.loc[df_seq['ID'] == 'sp' + str(id)]   # row in species data for this ID
        seq_str = seq_row['seq'].values[0]
        peak_data.append([max_peak, start, maxi, end, seq_str, id])

    peak_data = np.array(peak_data)  # converts entries to strings...
    # sort the species into descending order of max peak, to get consistent colours whenplotting thresholded species
    start_times = peak_data[:, 1].astype(np.int)    # column of start times, converted back to ints
    perm = np.argsort(start_times)  # permutation that sorts in ascending order
    peak_data = peak_data[perm, :]  # permute rows into ascending start time order

    # print this data out in LaTeX table format
    for i, row in enumerate(peak_data):
        # print(row)
        print(row[0], '&', row[1], '&', row[2], '&', row[3], '& \\verb|' + row[4] + '| % sp' + str(row[5]))
        print('\\\\')


def length_lifetime(run):
    # extract data about species lifetimes - scatter plot

    T = 101  # number of timesteps (including 0)
    len_max = 150  # ignore longer strings
    len_min = 3  # the shortest string that can bind -- ignore shorter strings

    species_colnames = ['start', 'end', 'n', 'seq']
    df = oee_utils.read_species_logfile(run, species_colnames)
    xs = []   # species lifetime
    ys = []   # string len
    ss = []  # species abundance
    for index, data_row in df.iterrows():  # each species
        life = (data_row['end'] - data_row['start']) / oee_utils.TSTEP
        length = len(data_row['seq'])
        n = data_row['n']
        # long lived (=> n>1), and only sensible string lengths
        if life > 1 and length <= len_max and len_min <= length:
            xs += [life]
            ys += [length]
            ss += [math.sqrt(n)]
    oee_plot.plot_scatter(xs, ys, ss, 'run' + str(run) + '_length_life_scatter.pdf', xmax=T, ymax=len_max)


def string_length_diversity(run):

    t0 = 0  # first time
    t1 = 2000000  # last time
    # t1 = 80000  # last time, for testing
    T = 101  # number of timesteps (including 0)

    # string lengths : no of toggle
    species_colnames = species_colnames = ['n', 'seq']
    l_max = 150  # longest string length to consider
    L = l_max + 1

    tot_n = np.zeros([L + 1, T])

    for i, time in enumerate(range(t0, t1 + 1, oee_utils.TSTEP)):
        # read in data about strings
        df = oee_utils.read_species_count_logfile(run, time, species_colnames)

        # count number of times each length of string observed
        for index, data_row in df.iterrows():  # each string
            lstr = len(data_row['seq'])  # length of this string
            if lstr < l_max:  # ignore longer strings
                tot_n[lstr, i] += data_row['n']    # number of strings of this length

    # logs make small densities clearer; 'where' used to ignore zero values
    tot_n = np.log(tot_n, out=np.zeros_like(tot_n), where=(tot_n != 0))
    oee_plot.plot_density(tot_n, 'run' + str(run) + '_strlen_density_all.pdf', l_max)
    
    
def self_repl(run):

    tstep = 20000
    t0 = 0  # first time
    t1 = 2000000  # last time
    # t1 = 80000  # last time, for testing
    T = 101  # number of timesteps (including 0)

    # diversity : no of string types
    species_colnames = ['seq']  # species sequences present at each timestep
    self_colnames = ['actseq', 'pp_NoProduct', 'pp_SelfReplicator']  # self-self properties of all species

    ys = np.zeros([3, T], dtype=int)  # zeroed in case no string of that subtype

    df_self = oee_utils.read_selfself_logfile(run, self_colnames)

    # a list of seq that are selfRepl
    df_selfrepl = df_self.loc[df_self['pp_SelfReplicator']]  # select all the rows with True pp_SelfRepl
    self_repl_seq_lst = df_selfrepl['actseq'].tolist()   # get list of the corresponding sequences

    # a list of seq that are noProduct
    df_noprod = df_self.loc[df_self['pp_NoProduct']]  # select all the rows with True pp_NoProduct
    noprod_seq_lst = df_noprod['actseq'].tolist()   # get list of the corresponding sequences

    for i, time in enumerate(range(t0, t1 + 1, tstep)):
        # read in data about species
        df = oee_utils.read_species_count_logfile(run, time, species_colnames)
        species_count = len(df.index)  # total number of species

        # get species rows that are self_repl
        df0 = df[df['seq'].isin(self_repl_seq_lst)]
        selfrepl_count = len(df0.index)
        ys[0, i] = selfrepl_count  # count number of self repls

        # get species rows that are noprod
        df1 = df[df['seq'].isin(noprod_seq_lst)]
        noprod_count = len(df1.index)
        ys[1, i] = noprod_count  # count number of self repls

        ys[2, i] = species_count - selfrepl_count - noprod_count
        assert(ys[2, i] >= 0)

    # plot each measure on separate plot
    oee_plot.plot_stacked(ys, 'run' + str(run) + '_self_repl.pdf', ymax=1510)


def repl_reactions(run):

    t0 = oee_utils.TSTEP  # first time -- no reactions at time zero
    t1 = 2000000  # last time
    # t1 = 80000  # last time, for testing
    T = 101  # number of timesteps (including 0)

    # breaking down repl reactions into their networks subtypes
    reaction_colnames = ['nobs', 'np_Parasite1', 'np_Parasite2', 'np_Parasite',
                         'np_MutualRepl', 'np_Hypercycle', 'pp_biolRep']

    ys = np.zeros([4, T], dtype=int)  # zeroed in case no string of that subtype
    ys_repl = np.zeros([T], dtype=int)  # zeroed in case no string of that subtype
    ys_nonrepl = np.zeros([T], dtype=int)  # zeroed in case no string of that subtype

    # populate species array with number of times repl reaction is observed
    for i, time in enumerate(range(t0, t1 + 1, oee_utils.TSTEP)):
        # read in data about reactions
        df = oee_utils.read_reaction_logfile(run, time, reaction_colnames)

        for index, data_row in df.iterrows():  # each species at this time
            nobs = data_row['nobs']

            # Parasitic is "para1 or para2" : para <=> para1 \/ para2
            assert (data_row['np_Parasite'] == (data_row['np_Parasite1'] or data_row['np_Parasite2']))
            # a reaction cannot be both para1 and para2: not (para1 /\ para2)
            assert (not(data_row['np_Parasite1'] and data_row['np_Parasite2']))
            # every hypercycle is a mutual repl : hyper => mutual repl
            assert (not(data_row['np_Hypercycle']) or data_row['np_MutualRepl'])
            # repls are exactly parasite or mutual (incl hypercycle): repl <=> (para or mutual)
            assert (data_row['pp_biolRep'] == (data_row['np_Parasite'] or data_row['np_MutualRepl']))

            if data_row['pp_biolRep']:  # repls
                ys_repl[i + 1] += nobs

                if data_row['np_Parasite']:  # parasitic -- repl but not mutual
                    ys[0, i + 1] += nobs
                elif data_row['np_MutualRepl']:  # mutual
                    if data_row['np_Hypercycle']:  # hypercycle
                        ys[2, i + 1] += nobs  # plot last because often zero
                    else:  # mutual, but not hypercycle
                        ys[1, i + 1] += nobs
            else:
                ys[3,i + 1] += nobs # non-repl

    # plot with total line, to test separation done correctly
    # oee_plot.plot_stacked(ys, 'run' + str(run) + '_parasite_check.pdf', ymax=6250, ys_line=ys_repl)
    # plot without total line, for use in publ
    oee_plot.plot_stacked(ys, 'run' + str(run) + '_parasite.pdf', ymax=6250)
    oee_plot.plot_lines(ys.T, 'run' + str(run) + '_parasite_lines.pdf', ymax=6520, ylog=True)

    # ys_percent = ys / ys_repl
    # plot without total line, for use in publ -- no check line needed, as 100% clear
    # oee_plot.plot_stacked(ys_percent, 'run' + str(run) + '_parasite_pc.pdf')


def toggle_opcodes(run):

    t0 = oee_utils.TSTEP  # first time
    t1 = 2000000  # last time
    # t1 = 80000  # last time, for testing
    T = 101  # number of timesteps (including 0)

    # opcodes : no of toggle
    reaction_colnames = ['nobs', 'ctogg', 'nsteps']
    lt_max = 15  # largest number of toggle ops
    Lt = lt_max + 1

    tot_tog = np.zeros([Lt + 1, T], dtype=int)

    for i, time in enumerate(range(t0, t1 + 1, oee_utils.TSTEP)):  # t0+1 because nothing at t=0
        # read in data about reaction opcode use, steps
        df = oee_utils.read_reaction_logfile(run, time, reaction_colnames)

        # count number of times each opcode observed executing
        for index, data_row in df.iterrows():  # each reaction
            n_steps = data_row['nsteps']
            if n_steps < 1001:  # stringmol reactions terminated at 1001 steps
                nobs = data_row['nobs']  # number of observations of this reaction
                n = data_row['ctogg']  # number of toggle ops this reaction executes
                if n < Lt:
                    tot_tog[n, i + 1] += nobs

    # logs make small densities clearer; context to ignore zero element
    with np.errstate(divide='ignore'):
        oee_plot.plot_density(np.log(tot_tog), 'run' + str(run) + '_opcode_toggle_density.pdf', lt_max)

#================================================================================================
# the following functions print out the relevant plots in the paper.
# Plots not needed/wanted can be commented out

# consolidated plots of multiple runs

# print total QNN values, plots of QNN (figure 2), and plots of cumulative QNN (figure 3)
start, end = 1, 8     # full runs
# start, end = 9, 20    # early death runs
qnn_plot(start, end)

# print minimum and maximum population size per run (figure 5), 
# plots of population size and boxplots of ranges of populations size (figure 4)
pop_size()

# print minimum and maximum species counts in first and second half of runs (figure 7)
# plots of total and new species counts and boxplots of ranges of species counts (figure 6)
species_count()

# Print number of unique species per run (figure 7, unique column), 
# stacked plots of individuals per species (figure 8) and plots of 12 most abundant species (figure 9)
# prints data on 12 most abundant species (figures 18--25)
species_12abundance()

# individual plots per run
for run in range(start, end + 1):

    # scatter plot of string length against lifespan (figure 10)
    length_lifetime(run)
    
    # density plot of string length against time (figure 11)
    string_length_diversity(run)
    
    # stacked plot of self-replicator and other reaction types against time (figure 12)
    self_repl(run)

    # stacked plot (figure 13) and log plot (figure 14) of parasitic reaction types against time
    repl_reactions(run)
    
    # density plot of toggle opcode use against time, for a given run (figure 15) 
    toggle_opcodes(run)
    
    pass  # to allow commenting out all options
