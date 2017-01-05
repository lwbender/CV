'''
# veriipy.optimizers.py:
# optimization, fitting, and other related algorithms.
#
'''
import matplotlib as mpl
import pylab as plt
import numpy
#import matplotlib.dates as mpd
import os
import json
import multiprocessing as mpp
import time
import math

#from numba import jit
#
'''
def stitch_sequence_by_dydx(data_in, dy_factor=1.5, dy_dx_0=None, median_out=0.):
	# stitch together sequences by extreme fluctuations.
	# for now, force the input to be like [[x,y], ...]
	#
	dy_dx = [(y-data_in[j][1])/(x-data_in[j][0]) for j,(x,y) in enumerate([rw[0], [1] for rw in data_in])]
	dy_dx_0 = (dy_dx_0 or numpy.median(dy_dx) )
	#
	dy_break = dy_factor*dy_dx_0
	#
	dy = 0.
	output_data = [data_in[0]]
	for j, (x,y) in enumerate(data_in[1:]):
		dydx = dy_dx[j]
		#
		if dydx>dy_break: dy -= (y-data_in[j][1])
		output_data += [[x,y+dy]
	#
	return output_data
'''		
#
def get_sub_sequence_indices(data_in, dt_0=None, t_col=0, dt_0_factor=1.5):
	# a more general split_sequence_by_dt() like tool. find large discontinuities in data_in[t_col];
	# return a list if indices that mark the start of subsequences. append also len(data_in), so all subsequences are
	# defined by full_sequence[get_sub_sequence_indices()[j-1:j]]
	#
	# first, get intervals on X[t_col]:
	intervals = [rw[t_col] - data_in[j][t_col] for j,rw in enumerate(data_in[1:])]
	dt_0 = (dt_0 or numpy.median([x for x in intervals if x>0]))  # filter out the dt=0 data. these should not exist, but errors in firmware produce
	#															  # some junk data, in this case the same record is repeated over and over.
	#
	
	dt_break = dt_0_factor*dt_0
	#
	#sub_sequence_indices = [0] + [j+1 for dt,j in enumerate(intervals) if dt>dt_break] + [len(data_in)]
	#
	# note: it is possible that this will give us repeated values at the start/finish. the simplest solution is to use set(X).
	# is it faster to write a copy of the list and use if() to truncate the ends? particularly, by the time we re-sort the indices
	# (becasue set() does not preserve order). all of that said, i think if we get repeaded indices, it might actually be because
	# the sequences are actually broken, so let's just keep it for now.
	#
	return [0] + [j+1 for j,dt in enumerate(intervals) if dt>dt_break] + [len(data_in)]
#
def split_sequence_by_dt(data_in, dt_0=None, dt_0_factor=1.5, fignum=None):
    '''
    # break up sequence into sub-sequences separated by dt>dt_0.
    # after it works, port it to dosimeter_board or veriipy or something like that.
    # assume list-like and that t are in col 0 and x are in col 1. ([[t,x], ...])
    # ... and for now, assume list-like.
    # @ dt_0: expected interval; split for dt>dt_0
    # @ dt_0_factor: ... times this factor; split for dt>dt_0_factor*dt_0 (allow for a bit of variation)
    #
    # TODO: add an option to return the sequence start (and stop?) indices. probably the easiest/best way to
    # to this is to return a dict object like {j_start:[data], }
    '''
    #
    #print('data_in: ', data_in[0:10])
    dt_0 = (dt_0 or numpy.median([rw[0]-data_in[j][0] for j,rw in enumerate(data_in[1:])]))
    #print("dt_0: ", dt_0)
    #
    working_datas = [list(data_in[0]) + [0.]] + [[rw[0], rw[1], rw[0]-data_in[j][0]] for j,rw in enumerate(data_in[1:])]
    #print('working_datas: ', working_datas[0:10])
    #
    measurement_sequences = [[]]
    for rw in working_datas:
        #if rw[2]>dt_0_factor*dt_0 and not len(measurement_sequences[-1])==0:
        if rw[2]>dt_0_factor*dt_0:
            #print("new sequence...")
            measurement_sequences+=[[]]
            #
        measurement_sequences[-1] += [rw[0:2]]
        #
    #
    # fignum, for diagnostics:
    if not fignum is None:
    	plt.figure(fignum)
    	plt.clf()
    	for seq in measurement_squences:
    		plt.plot(*zip(*seq), ls='-', marker='.')
    #
    return measurement_sequences
#
def make_zeroed_working_copy(ary, cols=['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7']):
    new_ary = ary.copy()
    l_ary = len(ary)
    for col in cols: new_ary[col]=numpy.zeros(l_ary)
    return new_ary
#
