'''
# author: Mark R. Yoder, Ph.D.
# email: myoder@landauer.com
#		mark.yoder@gmail.com
#
# Code property of: Landauer Inc.
#
# pt100.py: class object and relevant code to model pt100 thermistors.
#
# author and development notes:
# let's have a go at doing this all in Python3, as it seems this is the direction we need to be moving.
#
# getting started:
# there are several ready-for-consumption scripts here. 
'''
#
# note: we should be able to configure auto-reload. to actually use reload(), first import reload:
from imp import reload	# (noting that we'll want to do this from the console).
import numpy
import scipy
import scipy.fftpack
import scipy.optimize as spo

import math
import pylab as plt
import datetime as dtm
import os
import sys
import random

import veriipy.simple_filters as sfp

# this might be problematically circular. we're doing this to facilitate auto-calibration, but maybe there's a better way...
#import dosimeter_board as dbp

#
# permit, but do not require multiprocessing:
try:
	import multiprocessing as mpp
	have_mpp=True
except:
	have_mpp=False
#
has_h5=True
try:
	import h5py
except:
	has_h5=False
	print("loading without HDF5 support.")
#
# catch some python2.x compatibility:
if sys.version_info.major>=3:
	xrange = range	# catch python 2.x calls to xrange() (which python3 should have left alone).
#
import matplotlib.dates as mpd
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
#
default_temp_data_file = 'TM_PCAPDatabase_20151020_temps.csv'
default_exp_data_file = 'TM_PCAPDatabase_20151020_exp.csv'
cap_cols_default = ['t_int', 't_ext_pt100', 'c39', 'c40', 'c41', 'c42', 'c43', 'c44', 'c46']
#
#default_thermistor_prams = [7.69465109e-01, -2.18659959e-01, 2.07580013e-02, -6.55228081e-04]		# note: this gives T in Kelvin. also, these from the script
default_thermistor_prams = [ 0.02548983,  0.01608673,  0.00397388,  0.00032487]		# these from Tim's x-ray training data, calibrated to Pt1000. but note Pt1000 is reading funny.
#																												# also, the calibration/fit needs to be better optimized to get a best fit.
#colors_ =  mpl.rcParams['axes.color_cycle']
colors_ = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
_colors=colors_
#
#
fs_label=18.
fs_title=18.
#
class Thermistor_exp(object):
	# a thermistor based on an exponential callibration.
	# this is formulated and calibrated for the (typically) C46 thermistor.
	#def __init__(self, therm_prams=None):
	def __init__(self, A=2.16763217451, B=-4.4007896636, C=-2.82696313623, D=1.6021460755):
		# ... so this isn't very versatile, and it's basically hard-coded to this 4 parameter format, but it's easy to follow and it's what we're using.
		self.A = A
		self.B = B
		self.C = C
		self.D = D
		#
	def get_T(self, ratio_value):
		# get temperature from ratio value... so this is pretty specific to picocap operations.
		# 
		return numpy.exp((self.A + self.B*ratio_value)/(1. + self.C*ratio_value + self.D*ratio_value**2.))
		#return = numpy.exp((2.16763217451 + (-4.4007896636*ratio_value))/(1 + (-2.82696313623*ratio_value) + (1.6021460755*ratio_value**2)))
	#
		
class Thermistor(object):
	def __init__(self, *prams, **kwargs):
		# standard thermistor device. well, standard_ish. a standard thermistor is expected to have a temperatur-resistance relationship like:
		# 1/T = a + b*ln(R) + d*[ln(R)]**3
		# however, in our application, we do not always get a direct measurement of resistance, we get a resistance ratio x=R/R0, where R0 might not be well known (it should
		# be, but sometimes it's just not). consequently, expanding out the x**3 term, we get a polynomial with non-zero (or even negligibly small) d=2 term:
		# A' = A + [ln(R0)]**3
		# B' = (B + 3*[ln(R0)]**2)
		# C' = 3*ln(R0)		## (does not exist in standard format)
		# D' = D
		# note that the 'proper' A,B,D, C=0 parameter can be solved for via a linear transformation A*(a,b,0,d) = (a',b',c',d')
		# in general, then, let's just assume a general polynomial expansion. we might even go nuts and do an n>3 expansion.
		#
		# for now, assume we're working on Celsius, so we'll do Kelvin conversions transparently.
		#
		# make a guess that if we get only 3 prams, we mean n=3 and c=0
		if prams==None or len(prams)==0 or list(prams)==[None]:
			prams = default_thermistor_prams
		if len(prams)==3:
			prams = prams[0:2] + [0.] + prams[2]
		self.prams=prams
		#
		kwargs['Tk']=kwargs.get('Tk', 273.15)	# kelvin conversion.
		# any others?
		self.__dict__.update(kwargs)
		#
		
	#
	def get_T(self, R, *prams):
		prams = (prams or self.prams)
		#
		Tk=self.Tk
		Tk=0.
		
		return (1.0/numpy.sum(p*(numpy.log(float(R)))**j for j,p in enumerate(prams))) - Tk
	#
	def calibrate_to_XT(self, XT, n=3, med_len=4, set_prams=True):
		# calibrate this thermistor to an RT dataset, akak [[X,T], ...]
		# note, X does not have to be properly corrected or nomralized; if it's not, we'll get a finite n=2 coefficient.
		# again, assume T are Celsius... or anyway, that we use Tk as if they were; we can get all Kelvin by setting self.Tk=0.
		#
		# med_len: apply a median filter.
		#
		X_therm, T_calib = zip(*XT)
		if med_len>0:
			X_therm = sfp.filter_med_N(X_therm, med_len)
			T_calib = sfp.filter_med_N(T_calib, med_len)
		#
		A = [[x**d for d in range(n+1)] for x in numpy.log(X_therm)]
		b = [1./(t+self.Tk) for t in T_calib]
		#
		lsqs = numpy.linalg.lstsq(A,b)[0]
		#
		if set_prams: self.prams=lsqs
		#
		return lsqs
	
	'''
	# we might need to find a better, less ciruclar way to do this. pt100 objects are meant to be inherited by or members of Board_Data objects, not the other way around...
	#
	def calibrate_default(data_in='data/march23/exported-dosimeter-data-2016-03-24-092116.csv', col_therm='C46', col_T='T_Int', T0=-12., therm_factor=1.0, med_len=5, fignum=None):
		B=dbp.Board_Data(data_in=data_in)
		#
		# 1) at this time, temperature calculations are messed up; we need to recalibrate the pt1000 measurements (aka, what is the reference capacitor?).
		# 2) something screwy with the board/firmware is messing up the lower T domain.
		#
		X,T = zip(*[[rw[col_therm]/therm_factor, rw[col_T]] for rw in B.board_data if rw[col_T]>=T0])
		XT = list(zip(simple_filters.filter_med_N(X,med_len), simple_filters.filter_med_N(T, med_len)))
		#
		prams =  self.calibrate_thermistor(XT,n=3)
		#
		#f_poly = lambda x,P:numpy.sum([p*x**j for j,p in enumerate(P)])
		#f_T = lambda x,P: 1./numpy.sum([p*numpy.log(x)**j for j,p in enumerate(P)])
		f_T = self.get_T
		#
		if fignum!=None:
			plt.figure(fignum)
			plt.clf()
			print("prams: ", prams)
			#
			plt.plot([t for x,t in XT], [f_T(x,prams)-273.15 for x,t in XT], 'o')
		self.prams=prams
		return prams
	'''

class PT100(object):
	# class object for PT100 thermister. keep various callibration/conversion constants and funcions here.
	#
	# we need to add R_reff, which is the reference capacitor in the PicoCap circuit. we've been using R0 for this, somewhat accidentally.
	def __init__(self, R0=None, R100=None, R_reff=None, A=None, B=None, C=None, T_min=None, T_max=None, verbose=0):
		'''
		#if not R0==None and R100==None and A==None and B==None and C==None:
		# 	print("Warning: setting non-default values for key PT_100 parameters")
		# T_min, T_max are expected min/max values. the default (presently) values are (or might be) -200 < T < 850.
		# which represents the operational range of the PT_100 device, according to manufacturer white papers.
		# we might use these values to help decide which are the correct roots when we estimate temperature from
		# the solution of a 3rd or 4th order polynomial.
		'''
		#
		for key,val in locals().items():
			if key!='self' and val!=None:
				#print("Warning: setting key parameter to non-default value:: %s:%s" % (str(key), str(val)))
				if verbose: print("Setting key parameter to non-default value:: %s:%s" % (str(key), str(val)))
				
		#
		# note: these constants come from the pt100 spec sheet.
		R0    = (R0 or 100.)
		#R100  = (R100 or 138.51)
		R100  = (R100 or 1.3851*R0)
		#
		# if we don't know otherwise, guess that R_reff = R0
		R_reff = (R_reff or R0)
		#R_reff = 100.
		#
		A     = (A or 3.9803e-3)			# in C^-1
		B     = (B or -5.775e-7)			# in C^-2
		C     = (C or -4.183e-12)			# in C^-4
		T_min = (T_min or -200.)
		T_max = (T_max or 850.)
		#
		self.T_hi_lo = 0.		# temperature where we switch between high/low temp.formula.
		#
		# note: we'll be solving a polynomial to get temperature from measured resistance.
		# there are two solutions for -200C < T < 0C and for 0C < T < 850C. of course, we have to first solve both of these
		# to know which one to use. these are, respectively:
		# (1 - R_T/R_0) + A*T + B*T**2 -100*C*T**3 + C*T**4 = 0
		# (1 - R_T/R_0) + A*T + B*T**2
		self.__dict__.update(locals())
		#
	def get_temperature_from_ratio(self, R_T, R_reff=None, low_high_return=None):
		'''
		# get temperature from raw, resistance ratio measurement, R_T=R/R_0
		# (R_T is a bit misleading; it's not a resistance; it's the measured resistance ratio).
		# R0 defaults to class scope, then to 100.0
		'''
		if R_T==None: return None
		#
		R_reff = (R_reff or self.R_reff)
		R_reff = (R_reff or 100.0)
		#
		return self.get_temperature(R_T*R_reff, low_high_return=low_high_return)
		
		
	def get_temperature(self, R_T, low_high_return=None):
		'''
		# first, get temperature tuples; then guess which is the correct measurement.
		# i'm guessing the correct measurement must be real (with minimal imaginary component, but allowing for FP error there)
		# and preliminary  analyses suggest the root with minimum absolute value.
		# preliminary analysis also suggests that it's always the last root in the list (smallest real)
		# this behavior seems to be very consistent (and probably makes sense; eigen-values should come back ranked)
		# if this is not sufficient, down the road, develop criteria based on minimizing the Im. component (as close to 0 as possible) and
		# min. absolute value. initially excluding values outside the white paper specifications (see below) can also be a good start.
		# note: min(numpy.abs([T_values})) probably returns the correct item (with a bit of coersion).
		#
		# note the minimum operating value for R_T, R_min ~ 18.52 corresponds to -200C; R_max ~ 390.48.
		# we can use these criteria to chuck bogus readings and also to find the correct root solution (correct temperature T).
		#
		# low_high_return: 'low': return low value, 'high': return high value; any other: continue and do the hi/lo analysis (return the most likely correct number).
		# 
		'''
		#
		#D_low  = self._Temps_low(R_T)
		#D_high = self._Temps_high(R_T)
		#
		# do we just assume everything is real?
		T_low  = numpy.real(self._Temps_low(R_T)[-1])
		T_high = numpy.real(self._Temps_high(R_T)[-1])
		T0 = self.T_hi_lo
		#
		#  now 4 possibilities: one of each T_low, T_high is above/below 0; both are above/below.
		# ... weird behavior sometimes this returns an array type, not a float. so we should trap this. i think returning float(x) will suffice.
		if T_low >= T0 and T_high >= T0: return float(T_high)
		if T_low <  T0 and T_high <  T0: return float(T_low)
		if (T_low >= T0 and T_high <  T0) or (T_low < T0 and T_high >= T0):
			# return a weighted average like ((T1-T0)*T1 + (T2-T0)*T2)/abs(T2-T1)
			#
			# for these cases, the two temperatures will almost certainly be very close to one another, so an unweighted average would probably suffice.
			# but this is probably more correct.
			return float(T_low*abs(T_low-T0) + T_high*abs(T_high-T0))/abs(T_high-T_low)
		
	#
	def _temp_tuple_test(self, N=100, R_low=84.27, R_high=119.4):
		# temp tuple testing script.
		# R_low, R_high correspond to approximately -40C, 50C.
		#
		dR = R_high-R_low
		R=random.Random()
		#return {r:self.get_temperature_tuples(rtype='list', R_T=R_low + (r*dR)) for k,r in enumerate([R.random() for j in range(N)])}
		#return {r:self.get_temperature_tuples(rtype='list', R_T=R_low + (r*dR)) for r in [R.random() for j in range(N)]}
		#
		T_roots = []
		#R_Ts = []
		for j in range(N):
			#
			this_R = R_low + (R.random()*dR)
			rt = self.get_temperature_tuples(rtype='dict', R_T = this_R)
			T_roots += [[this_R] + [x for x in rt['low']] + [x for x in rt['high']]]
			#R_Ts += [this_R]
		#
		T_roots.sort(key=lambda rw: rw[0])
		#
		plt.figure(0)
		plt.clf()
		plt.title('temperature tuple test.\n(last col from each solution (3,5) should be correct T)')
		plt.ylabel('Temperature $T$')
		plt.xlabel('$R_T$')
		#
		# ... and plot it:
		# colors_
		R_Ts = [rw[0] for rw in T_roots]
		for j,col in enumerate(list(zip(*T_roots))[1:]):
			#
			zo = 5
			mkr='.'
			lw=1.5
			if j>3:
				zo=4
				mkr='^'
				lw=2.
			clr = colors_[j%len(colors_)]
			plt.plot(R_Ts, numpy.real(col), ls='-', lw=lw, marker=mkr, color=clr, label='re_col: %d' % j)
			plt.plot(R_Ts, numpy.imag(col), ls='--', marker='s', color=clr, zorder=1)
		plt.legend(loc=0, numpoints=1)
	#	
	def get_temperature_tuples(self, R_T, rtype='dict'):
		# calculate temperature (tuple?)
		# for now, force the use of class member values. solve the manufacturer polynomial to get temperature.
		#
		# use numpy.roots(p) where p = (highest_rank_coef, ..., 0-rank coef)
		# let's further modularize:
		Ts_low  = self._Temps_low(R_T)
		Ts_high = self._Temps_high(R_T)
		#
		if rtype=='dict':
			return {'low':Ts_low, 'high':Ts_high}
		elif rtype == 'tuple':
			return (Ts_low, Ts_high)
		else:
			return [Ts_low, Ts_high]
	def _Temps_low(self, R_T):
		return numpy.roots((self.C, -100.*self.C, self.B, self.A, (1.-R_T/self.R0)))	# note: from pt100 white paper, this factor of 100 is correct; it is *not* supposed to be
		#																										# the reference capacitor, accidentally hard-coded. it is, in fact, in reference to T=100oC
	def _Temps_high(self, R_T):
		return numpy.roots((self.B, self.A, (1.-R_T/self.R0)))
	#
#
class PT1000(PT100):
	# PT1000 subclass of PT100. the only difference (i think) is that the reference resistance R0=1000, not 100.
	#
	def __init__(self, *args, **kwargs):
		# change R0 to 1000 if it's None. we could force this subclass to always have R0=1000, but we want these objects to be seamless and versatile.
		# aka, if we grow accustomed to using PT1000, and PT100 fades into the sunset, we want it to remain versatile, and we don't want to need to know
		# about its parent class(es).
		#
		R0=1000
		if len(args)>0:
			#if args[0]==None: args[0] = R0
			args[0] = (args[0] or R0)
			#
		#elif 'R0' in kwargs.keys():
		else:
			#kwargs['R0'] = (kwargs['R0'] or R0)
			kwargs['R0'] = (kwargs.get('R0',R0) or R0)
			#
		super(PT1000, self).__init__(*args, **kwargs)
#
