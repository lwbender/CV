#
import numpy
import scipy
import scipy.optimize as spo
#
import math
import pylab as plt
import datetime as dtm
import os
import sys
import random
import glob
#
import multiprocessing as mpp
#import  dosimeter_board as dbp
#Dosimeter = dbp.Board_Data
#
#import simple_filters as sfp
import veriipy.simple_filters as sfp
#
# catch some python2.x compatibility:
if sys.version_info.major>=3:
	# (of course, at this point, if we're somehow running python 2.x,just shoot me.)
	xrange = range	# catch python 2.x calls to xrange() (which python3 should have left alone).
#
import matplotlib.dates as mpd
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
#
# the new matplotlib uses a more complex iterator type class to handle colors, lines, etc. if we just want colors,
# we should just spell themout here.
#_colors =  mpl.rcParams['axes.color_cycle']		# ... and there's a simpler syntax for this somewhere...
_colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
colors_=_colors
#
def multi_plot(X, Ys, LSs='-', markers='.', LWs=1., lbls=None, fignum=0, ax=None, do_clf=True, figsize=None):
	if not hasattr(Ys[0],'__len__'): Ys=[Ys]
	Ys = [numpy.array(y) for y in Ys]
	if lbls=='': lbls=None
	#
	# note: this does not take into account an incomplete list... maybe a better way is to just use an iterator and/or modulus approach.
	if isinstance(LSs,str) or not hasattr(LSs,'__len__'): LSs=[LSs for _ in Ys]
	if isinstance(LWs,str) or not hasattr(LWs,'__len__'): LWs=[LWs for _ in Ys]
	if isinstance(markers,str) or not hasattr(markers,'__len__'): markers=[markers for _ in Ys]
	if lbls==None: lbls = list(range(len(Ys)))
	
	#print('stuff*: ', LSs, LWs, markers)
	#
	if ax==None:
		plt.figure(fignum, figsize=figsize)
		if do_clf: plt.clf()
		ax=plt.gca()
	#
	ax.plot(X,Ys[0],marker=markers[0],lw=LWs[0], ls=LSs[0], label=lbls[0])
	ax.plot([X[0], X[-1]], [Ys[0][0], Ys[0][0]], color='k', ls='-', alpha=.5, lw=1.5)
	#
	for j,Y in enumerate(Ys[1:]):
		#
		ax.plot(X,Y-(Y[0]-Ys[0][0]) ,marker=markers[j+1], lw=LWs[j+1], ls=LSs[j+1], label=lbls[j+1])
		#
	#

def plot_array_linear_polar(X, Ys, cols_to_plot=None, title_decoration='', verbose=0, periodicity=None, figsize=(6,4), subtract_lin_med=True, x_lbl='', y_lbl=''):
    #
    # ... a cool trick would be to guess periodicity by selecting a dominant Fourier mode...
    # for now, we could guess "1", but it's just a guess; let's assume that it's not periodic.
    #
    # ... and looking at it, this function hasn't been threshed out very well, so it probably needs some work.
    #
    # if not a dict, assume it's a row-array: [[a0,b0,c0], [a1,b1,c1],...]
    # but, so long as we are carefula bout the syntax we exploit, we can use a dict...
    # and we can morph this into a dict if it's just an array... but we'll work on that later.
    #
    if not hasattr(Ys, 'keys') and not hasattr(Ys, 'dtype'): Ys = {'col_{}'.format(j):col for j,col in enumerate(zip(*Ys))}
    if cols_to_plot is None: cols_to_plot=list(Ys.keys())
    periodicity = (periodicity or (max(X)-min(X)))
    #
    #for j, (key,val) in enumerate(Ys.items()):
    for j,key in enumerate(cols_to_plot):
        val = Ys[key]
        if verbose is not None and verbose>0:print('daig: ', j, key)
        #
        fg=plt.figure(figsize=figsize)
        ax1 = fg.add_axes([.1,.1,.4,.8])
        ax2 = fg.add_axes([.55, .1, .4,.8], projection='polar')
        #med_val = numpy.median(val)
        med_val_rad = numpy.min(val)
        med_val_lin = med_val_rad
        if not subtract_lin_med: med_val_lin=0
        #
        plt.suptitle('{}col: {}'.format(title_decoration,key))
        #ax=plt.gca()
        '''
<<<<<<< HEAD
        ax1.plot(X, val-med_val, color=colors_[j%len(colors_)], ls='-', marker='.')
        ax1.plot(X[0], val[0]-med_val, marker='s', color='r', ms=14, label='start')
        ax1.plot(X[-1], val[-1]-med_val, marker='o', color='m', ms=14, label='finish')
        #

        ax2.plot(numpy.array(X)*2.0*((numpy.pi/periodicity)), numpy.array(val)-med_val, color=colors_[j%len(colors_)], ls='-', marker='')
        ax2.plot([X[0]*2.0*((numpy.pi/periodicity))], val[0]-med_val, marker='s', color='r', ms=14, label='start')
        ax2.plot([X[-1]*2.0*((numpy.pi/periodicity))], val[-1]-med_val, marker='o', color='m', ms=14, label='finish')
        '''

        ax1.plot(X, val-med_val_lin, color=colors_[j%len(colors_)], ls='-', marker='.')
        ax1.plot(X[0], val[0]-med_val_lin, marker='s', color='r', ms=14, label='start')
        ax1.plot(X[-1], val[-1]-med_val_lin, marker='o', color='m', ms=14, label='finish')
        #
        ax2.plot(numpy.array(X)*2.0*((numpy.pi/periodicity)), numpy.array(val)-med_val_rad, color=colors_[j%len(colors_)], ls='-', marker='')
        ax2.plot([X[0]*2.0*((numpy.pi/periodicity))], val[0]-med_val_rad, marker='s', color='r', ms=14, label='start')
        ax2.plot([X[-1]*2.0*((numpy.pi/periodicity))], val[-1]-med_val_rad, marker='o', color='m', ms=14, label='finish')
        #

        ax2.legend(loc=0, numpoints=1)
        ax2.set_title('$R=Y_j - \\langle Y \\rangle$ \n')
        ax1.set_xlabel(x_lbl)
        ax1.set_ylabel(y_lbl)
	
if __name__=='__main__':
	pass
else:
	plt.ion()
	
