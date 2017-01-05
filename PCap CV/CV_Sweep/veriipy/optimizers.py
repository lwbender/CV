'''
# veriipy.optimizers.py:
# optimization, fitting, and other related algorithms.
#
# i think the thing to do, eventually, is to integrate these functions and classes (this module) into
# a more modular, generalized, portable "optimizers" library. these are really "auto_optimizers"; we can 
# put these together with ROC tools and other optimization algorithms.
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

from numba import jit
#
#
def auto_optimize_cols(A,Bs,col_names=None, min_len_factor=.5, step_size=1, chi_expon=1.5, n_cpu=None):
	# a procedural wrapper to use Auto_Optimizer and return a
	return None
	

class Auto_Optimizer(object):
	# would be benefit from interiting mpp.Pool(), or just contain it, and maybe manually map some of the functions?
	# we could build this like we built Global_ETAS, where we have a base Optimizer class, variations of that, and then an MPP handler that inherits from Optimizer
	# and tranparelty handles the mpp parts. or, we just create a single class that recognizes the spp/mpp implementations.
	#
	# earlier versions of this parallelized by distributing multiple optimization jobs to processors. we'll actually break down each optimizatio job.
	# so we'll parallelize by splitting over the fit-range of each optimization job,
	# and we'll distirubt the start-points so that each processor gets approximately the same load, aka:
	#  - step_size_mpp = step_size*n_cpu
	#  - start_indices = [j for j in range(n_cpu)]
	#
	# note: as per pickling and piping overhead, for ~1000 rows, the spp solution is much faster than full mpp, so we should put in some sort of logic to guess the best
	# cpu-count (of course, if it's not mpp.cpu_count(), then it's only a second or two, so who cares).
	#
	def __init__(self, A,B, col_name='mycol', min_len_factor=.5, step_size=1, chi_expon=1.5, n_cpu=None, auto_run=True):
		#
		self.__dict__.update(locals())
		#self.A=A
		#self.B=B		# result vector(s), transposed, so b1, b1, b3 --> B = [[b10, b20, b30], [b11,b21,b31], ...]
		#
		n_cpu=(n_cpu or mpp.cpu_count())
		if n_cpu<1: n_cpu = mpp.cpu_count()-n_cpu
		#
		if len(A)<1500: n_cpu=1
		#
		#
		if auto_run:
			self.best_prams = self.optimize()
		#
	#
	def optimize(self, A=None, B=None, n_cpu=None):
		n_cpu=(n_cpu or self.n_cpu)
		n_cpu=(n_cpu or mpp.cpu_count())
		#
		if A is None: A=self.A
		if B is None: B=self.B
		#
		#print('n_cpu: ', n_cpu)
		if n_cpu==1:
			self.best_prams = self.auto_optimize_array(A,B)
		else:
			self.best_prams = self.auto_optimize_mpp_handler(A,B, n_cpu)
		#
		return self.best_prams
		#
	#
	def auto_optimize_mpp_handler(self,A,B, n_cpu=None):
		if A is None: A=self.A
		if B is None: B=self.B
		n_cpu = (n_cpu or self.n_cpu)
		#
		mpp_step_size = self.step_size*n_cpu
		#
		P = mpp.Pool(n_cpu)
		resultses = [P.apply_async(self.auto_optimize_array, args=(A,B), kwds={'start_index':j, 'step_size':mpp_step_size}) for j in range(n_cpu)]
		P.close()
		P.join()
		#
		best_fits = [res.get() for res in resultses]
		#print('best_fits: ', best_fits)
		best_fit = sorted(best_fits, key=lambda rw: (None if len(rw[2])==0 else numpy.sqrt(numpy.mean(rw[2]))/max(1, (len(A)-rw[0]))**self.chi_expon))[0]
		#
		return best_fit
	#
	# auto_optimize_array can fundamentally be a procedural function, outside this class. we'll just write a connecter to it here.
	def auto_optimize_array(self, A=None,B=None, min_len_factor=None, step_size=None, start_index=0, chi_expon=None,fignum=None, do_clf=True, ax=None, show_all_lsqs=True, do_aggregate=True):
		#
		# @A: matrix of model coefficients (aka, T,zeta, etc. series).
		# @B: matrix (transposed vectors) of outputs, so:
		#  A \dot p' = b, or B = [b0,b1,b2...]^T, and A \dot P' = B, and p are the parameters for one vector, P are the parameters for multiple vectors.
		# using more than one output vector b amounts to solving for a set of training parameters from multiple sources.
		# this approach, however, might be problematic. the residual is minimized for each vector b_j, but the starting point (and therefore the selected set of fits)
		# is chosen by the combined minimum chi-square, so it is probably more rigorous to simply fit each column completely independently.
		#
		if A is None: A = self.A
		if B is None: B = self.B
		if min_len_factor == None: min_len_factor=self.min_len_factor
		if step_size==None: step_size=self.step_size
		chi_expon = (chi_expon or self.chi_expon)
		#
		return auto_optimize_array(A=A,B=B,min_len_factor=min_len_factor,step_size=step_size,start_index=start_index,chi_expon=chi_expon, fignum=fignum, do_clf=do_clf, ax=ax,show_all_lsqs=show_all_lsqs,do_aggregate=do_aggregate)
		'''
		# walk through starting positions 0,1,... and to a least-squares fit. pick the subset that gives the best fit.
		# this is largely ripped out from get_optimal_training_sets(), and in fact the matrix part of get_optimal_training_sets() should be ripped out and reference this function.
		#
		# self-optimize Ab' = b. for now, assume our concern is the leading parts of A,b, so we'll start with the full sequence and walk down the starting point until we get an optimal chi-square.
		# for now, force b --> vector.
		#
		#	
		len_A = len(A)
		len_rw = len(A[0])
		#
		#print('lens: ', len(A), len(b))
		#
		max_starting_index = max(1, int(min_len_factor*len(A)))
		#
		#lsqs = [[j] + list(numpy.linalg.lstsq(A[j:], b[j:])[0:2]) for j in range(int(len(A)-1-len_rw))]
		# this works, but only shows the subset. do we want to see the full distribution of lsqs?
		if show_all_lsqs:
			#lsqs = [[j] + list(numpy.linalg.lstsq(A[j:], B[j:])[0:2]) for j in range(start_index, int(len(A)-1-len_rw), step_size)]
			lsqs = [[j] + list(numpy.linalg.lstsq(A[j:], B[j:])[0:2]) for j in range(start_index, min(max_starting_index, int(len(A)-1-len_rw)), step_size)]
		else:
			lsqs = [[j] + list(numpy.linalg.lstsq(A[j:], B[j:])[0:2]) for j in range(start_index, min(max_starting_index, int(len(A)-1-len_rw)), step_size)]
		#
		#print('len(lsqs): ', len(lsqs), start_index, step_size)
		# sometimes we get a weird outcome, so trap cases where lsqs() returns junk.
		#chis = [[x[0],(None if len(x[2])==0 else numpy.sqrt(x[2][0])/max(1, (len(A)-x[0]))**chi_expon)] for x in lsqs] 	# aka, reduced_chi-square/sqrt(n_dof) like quantity.
		chis = [[k, x[0],(None if len(x[2])==0 else numpy.sqrt(numpy.mean(x[2]))/max(1, (len(A)-x[0]))**chi_expon)] for k,x in enumerate(lsqs)] 	# aka, reduced_chi-square/sqrt(n_dof) like quantity.
		#
		#print('chis: ', chis[0:5])
		
		# now, get the optimal chi:
		chis.sort(key=lambda rw: rw[-1])
		k_chi, j_chi, chi_0 = chis[0]		# for spp (start_index=0, step_size=1), k_chi==j_chi... right?
		#
		if fignum!=None:
			plt.figure(fignum)
			if do_clf: plt.clf()
			if ax==None: ax=plt.gca()
			#
			clr = colors_[len(plt.gca().lines)%len(colors_)]
			ax.plot(*zip(*[(x,y) for x,y in chis]), ls='-', color=clr, lw=2.5)
			ax.plot([j_chi], [chi_0], marker='o', color=clr)
		#
		# do we want to necessarily aggregate here? why not. we can change it later, or make it optional.
		# lsqs returns like [[a0,a1,a2,...], [b0,b1,b2,...], ]
		# if we pass only a single b vector, this is still true, but each parameter is an l=1 vector, so the operation does not change.
		# same idea with the second value/array returned by lsqs.
		#
		if do_aggregate:
			#
			return [lsqs[k_chi][0], numpy.array([numpy.mean(X) for X in lsqs[k_chi][1]]), [numpy.mean(lsqs[k_chi][2])]]
		else:
			#
			return lsqs[k_chi]
		#
		'''
#
def auto_optimize_multiple_arrays(A=None,bs=None, min_len_factor=.5, step_size=1, start_index=0, chi_expon=1.5,fignum=None, do_clf=True, ax=None, show_all_lsqs=True, do_aggregate=True,n_cpu=None):
	'''
	# (this appears to work.. but needs more work. see **prms and params sent in mpp mode as well.).
	# auto_optimize; parallelize by each column. we see in Auto_Optimize() that for relatively small sets (certainly len(A)<=1000 or so), the mpp overhead is more costly than the benefits.
	# so for many operations, it will be faster to just parallelize by full sequence, so there's (presumably) less piping.
	#
	# A: matrix of model coefficients
	# bs: array of output vectors. note this is different than some of the ther auto_optimize_array() inputs; for this, we're going to force b -> a vector, and forego
	# simultaneous fits.
	#   SO, A --> a fitting matrix, [[a0,c0, d0,...], [a1,c1, d1, ...], ...], bs -> [b1, b2, b3...] ], where b_j are vectors of observables (aka, MOScap measurement seq.).
	'''
	print('running auto-optimize for multiple vectors.')
	original_inputs = {key:val for key,val in locals().items() if not key in ['n_cpu']}
	if n_cpu==1:
		prms = {key:val for key,val in original_inputs.items() if key not in ['bs']}
		prms['n_cpu']=1
		#
		# note: we might need touse this weird format for b/B, basically to facilitate arrays of multiple b vectors... or we rewrite the optimize function.
		#return [auto_optimize_array(B=[[x] for x in b], **XX) for b in bs]
		return [auto_optimize_array(A=A,B=b, **prms) for b in bs]
	#
	n_cpu = (n_cpu or mpp.cpu_count())
	P = mpp.Pool(n_cpu)
	#
	print('mpp auto-optimizing: ', A.shape, bs.shape)
	#resultses = [P.apply_async(auto_optimize_array, args=(A,b), kwds={'start_index':start_index, 'step_size':step_size}) for b in bs]
	opt_params = {'start_index':start_index, 'step_size':step_size, 'min_len_factor':min_len_factor, 'chi_expon':chi_expon, 'do_aggregate':do_aggregate}
	#
	# do these necessarioly come back in the correct order? i think they do (or backwards) because they are loaded into their respective arrays in order.
	# in other words, we can label them from the calling side, as opposed to inforcing a data structure like:
	# resultses = {col_name:P.apply_async(auto_optimize_array, args=(A,b), kwds=opt_params) for col_name,b in bs.items()]
	#  ... and then b_f = {key:res.get() for key,res in resultes.items()} ... or we can do it as a 2D list.
	print('auto-fitting. len(bs): {}, {}'.format(len(list(bs)), numpy.shape(bs)))
	resultses = [P.apply_async(auto_optimize_array, args=(A,b), kwds=opt_params) for b in bs]
	#
	P.close()
	P.join()
	#
	best_fits = [res.get() for res in resultses]
	#print('best_fits: ', best_fits)
	#best_fit = sorted(best_fits, key=lambda rw: (None if len(rw[2])==0 else numpy.sqrt(numpy.mean(rw[2]))/max(1, (len(A)-rw[0]))**chi_expon))[0]
	#
	return best_fits	
#
def auto_optimize_mpp_handler(self,A,B, step_size=1, n_cpu=None):
	'''
	# an mpp handler to auto-fit a single sequence. split up the job by splitting up the starting points. for short sequences, the overhead seems to overwhelm the
	# benefit of mpp, and it is probably better to use either a spp function or auto_optimize_multiple_arrays(). nominally, this handler is adept for optimizing one or two
	# really long sequences.
	'''
	#if A is None: A=self.A
	#if B is None: B=self.B
	#n_cpu = (n_cpu or self.n_cpu)
	#
	n_cpu = (n_cpu or mpp.cpu_count())
	#
	mpp_step_size = step_size*n_cpu
	#
	P = mpp.Pool(n_cpu)
	resultses = [P.apply_async(auto_optimize_array, args=(A,B), kwds={'start_index':j, 'step_size':mpp_step_size}) for j in range(n_cpu)]
	P.close()
	P.join()
	#
	best_fits = [res.get() for res in resultses]
	#print('best_fits: ', best_fits)
	best_fit = sorted(best_fits, key=lambda rw: (None if len(rw[2])==0 else numpy.sqrt(numpy.mean(rw[2]))/max(1, (len(A)-rw[0]))**self.chi_expon))[0]
	#
	return best_fit
#	
def auto_optimize_array(A=None,B=None, min_len_factor=.5, step_size=1, start_index=0, chi_expon=2,fignum=None, do_clf=True, ax=None, show_all_lsqs=True, do_aggregate=True):
	'''
	# self-optimize Ab' = b. for now, assume our concern is the leading parts of A,b, so we'll start with the full sequence and walk down the starting point until we get an optimal chi-square.
	# for now, force b --> vector.
	# note: this can be used as a SPP worker for an MPP handler like auto_optimize_multiple_arrays(). for smaller arrays, this approach is much faster than MPP models that break
	# up each data set for MPP handling. for 2 or 3 (ana, n<n_cpu) very large data sets (don't know what that means yet), the latter approach is probably optimal, and i think it is coded
	# up somewhere.
	#
	# @A: matrix of model coefficients (aka, T,zeta, etc. series).
	# @B: matrix (transposed vectors) of outputs, so:
	#  A \dot p' = b, or B = [b0,b1,b2...]^T, and A \dot P' = B, and p are the parameters for one vector, P are the parameters for multiple vectors.
	# using more than one output vector b amounts to solving for a set of training parameters from multiple sources.
	# this approach, however, might be problematic. the residual is minimized for each vector b_j, but the starting point (and therefore the selected set of fits)
	# is chosen by the combined minimum chi-square, so it is probably more rigorous to simply fit each column completely independently.
	#
	# walk through starting positions 0,1,... and to a least-squares fit. pick the subset that gives the best fit.
	# this is largely ripped out from get_optimal_training_sets(), and in fact the matrix part of get_optimal_training_sets() should be ripped out and reference this function.
	#
	'''
	#
	#	
	len_A = len(A)
	len_rw = len(A[0])
	#
	#print('lens: ', len(A), len(b))
	#
	max_starting_index = max(1, int(min_len_factor*len(A)))
	#
	#lsqs = [[j] + list(numpy.linalg.lstsq(A[j:], b[j:])[0:2]) for j in range(int(len(A)-1-len_rw))]
	# this works, but only shows the subset. do we want to see the full distribution of lsqs?
	if show_all_lsqs:
		#lsqs = [[j] + list(numpy.linalg.lstsq(A[j:], B[j:])[0:2]) for j in range(start_index, int(len(A)-1-len_rw), step_size)]
		lsqs = [[j] + list(numpy.linalg.lstsq(A[j:], B[j:])[0:2]) for j in range(start_index, min(max_starting_index, int(len(A)-1-len_rw)), step_size)]
	else:
		lsqs = [[j] + list(numpy.linalg.lstsq(A[j:], B[j:])[0:2]) for j in range(start_index, min(max_starting_index, int(len(A)-1-len_rw)), step_size)]
	#
	#print('len(lsqs): ', len(lsqs), start_index, step_size)
	# sometimes we get a weird outcome, so trap cases where lsqs() returns junk.
	# lsqs are like [[start_index, params, chi_sqrs], ...]
	# where chi_sqrs are an array (in general) because we might pass a matrix, not a simple vector, to the lstsq() operation, and so we average them (mean(x[2])  below)
	# for our cost-function.
	#
	#chis = [[x[0],(None if len(x[2])==0 else numpy.sqrt(x[2][0])/max(1, (len(A)-x[0]))**chi_expon)] for x in lsqs] 	# aka, reduced_chi-square/sqrt(n_dof) like quantity.
	chis = [[k, x[0],(None if len(x[2])==0 else numpy.sqrt(numpy.mean(x[2]))/max(1, (len(A)-x[0]))**chi_expon)] for k,x in enumerate(lsqs)] 	# aka, reduced_chi-square/sqrt(n_dof) like quantity.
	#
	#print('chis: ', chis[0:5])
	
	# now, get the optimal chi:
	chis.sort(key=lambda rw: rw[-1])
	k_chi, j_chi, chi_0 = chis[0]		# for spp (start_index=0, step_size=1), k_chi==j_chi... right?
	#
	if fignum!=None:
		plt.figure(fignum)
		if do_clf: plt.clf()
		if ax==None: ax=plt.gca()
		#
		clr = colors_[len(plt.gca().lines)%len(colors_)]
		ax.plot(*zip(*[(x,y) for x,y in chis]), ls='-', color=clr, lw=2.5)
		ax.plot([j_chi], [chi_0], marker='o', color=clr)
	#
	# do we want to necessarily aggregate here? why not. we can change it later, or make it optional.
	# lsqs returns like [[a0,a1,a2,...], [b0,b1,b2,...], ]
	# if we pass only a single b vector, this is still true, but each parameter is an l=1 vector, so the operation does not change.
	# same idea with the second value/array returned by lsqs.
	#
	if do_aggregate:
		#
		return [lsqs[k_chi][0], numpy.array([numpy.mean(X) for X in lsqs[k_chi][1]]), [numpy.mean(lsqs[k_chi][2])]]
	else:
		#
		return lsqs[k_chi]
		#

#
def optimizer_test(data_file=os.path.join( os.path.split(os.path.abspath(__file__))[0],'data/test_optimization_data.json'), n_cpu=None):
	# define some tests. some reasonable unit tests would also be to modify the b column values (aka, *=1.5 or something) and see that 1) the outputs change,
	# and 2) if we modify the single column and the 2-of-same column, we get the same outputs.
	# ... and other mods too...
	#
	AB = json.load(open(data_file,'r'))
	#
	print('load an optimizer object, then run some unit diagnostics.')
	print('All fitting solutions should produce the same output. eventually, we will build the check-sum for this into the unit test itself.')
	opter = Auto_Optimizer(AB['A'],numpy.array(AB['B']), n_cpu=1, auto_run=False)
	#
	print('run a simple, spp optimization')
	best_prams_spp1 = opter.optimize(n_cpu=1)
	print('best prams, spp1.: ', best_prams_spp1)
	#
	# (and a little bit of a different calling signature, but amounting to the same thing)
	# fit two of the same column:
	print('fit to copy of target columns. output should be the same as spp1')
	opter = Auto_Optimizer(AB['A'],numpy.array([[b[0],b[0]] for b in AB['B']]), n_cpu=1)
	print('best prams2: ', opter.best_prams)
	#
	mpp_pramses=[]
	#
	for k in range(1,mpp.cpu_count()+1):
		# now, try with multiple processes:
		print('run a simple, mpp[{}] optimization'.format(k))
		t0=time.time()
		mpp_pramses += [opter.optimize(n_cpu=k)]
		print('best prams[dt={}], mpp{}.: {}'.format(time.time()-t0, k, mpp_pramses[-1]))
	#
	print('******')
	xx = opter.auto_optimize_array(start_index=0,step_size=1)
	print('xx: ', xx)
	xx = opter.auto_optimize_array(start_index=1,step_size=2)
	print('xx: ', xx)
	xx = opter.auto_optimize_array(start_index=0,step_size=2)
	print('xx: ', xx)
	#
#
def jit_fact_test(N,n):
	# script to test @jit compiler.
	obj = Jit_tester()
	obj.jit_fact_test(5,5)
	print('*******************\n\n')
	obj.jit_fact_test(N,n)
	#
#
@jit
def n_fact_jit(x):
	# calc x! with a for-loop
	# this is just a test, so for a fraction, just take the int(x) and chuck a warning.
	if x%1!=0:
		print('WARNING: non-integer value submitted: %f. calculating, instead, %d!' % (x,int(x)))
		x=int(x)
	#
	x=float(x)
	#xx=1
	#for j in range(1,x+1): xx*=j
	for j in range(1,int(x)): x*=j
	#
	return x
class Jit_tester(object):
	# can we make numbas.jit work?
	def jit_fact_test(self,N,n=100):
		#
		# test jit with a factorial test.
		t0=time.time()
		for j in range(n): x=self.n_fact_jit(N)
		print('sJIT: N! = {}'.format(x))
		t1=time.time()
		print('finished: {}/{}'.format(t1-t0, math.log((t1-t0))))
		#
		# test jit with a factorial test.
		t0=time.time()
		for j in range(n): x=n_fact_jit(N)
		print('xJIT: N! = {}'.format(x))
		t1=time.time()
		print('finished: {}/{}'.format(t1-t0, math.log((t1-t0))))
		#
		t0=time.time()
		for j in range(n): x=self.n_fact(N)
		print('pypy: N! = {}'.format(x))
		t1=time.time()
		print('finished: {}/{}'.format(t1-t0, math.log((t1-t0))))
		#
		t0=time.time()
		for j in range(n): x=math.factorial(N)
		print('math: N! = {}'.format(x))
		t1=time.time()
		print('finished: {}/{}'.format(t1-t0, math.log((t1-t0))))
	#
	@jit
	def n_fact_jit(self,x):
		# calc x! with a for-loop
		# this is just a test, so for a fraction, just take the int(x) and chuck a warning.
		if x%1!=0:
			print('WARNING: non-integer value submitted: %f. calculating, instead, %d!' % (x,int(x)))
			x=int(x)
		#
		x=float(x)
		#xx=1
		#for j in range(1,x+1): xx*=j
		for j in range(1,int(x)): x*=j
		#
		return x
	#
	def n_fact(self,x):
		# calc x! with a for-loop
		# this is just a test, so for a fraction, just take the int(x) and chuck a warning.
		if x%1!=0:
			print('WARNING: non-integer value submitted: %f. calculating, instead, %d!' % (x,int(x)))
			x=int(x)
		#
		x=float(x)
		#xx=1
		#for j in range(1,x+1): xx*=j
		for j in range(1,int(x)): x*=j
		#
		return x
	
	

	
