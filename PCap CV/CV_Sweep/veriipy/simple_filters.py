
import numpy
import math
#
def filter_extreme_dy(data_in, dy0=None, dy_factor=10.0):
	# like filter_large_dy_dx, but i think simpler...
	# remove obviously huge spikes that are clearly artifacts. replace with a linear interpolation. since we're only taking
	# the y-values, just use an average (assume even spacing).
	# also assume quasi-stationary...  ??? for now, just assume the first dp is ok.
	#
	data_in=data_in.copy()
	if dy0==None: dy0 = dy_factor*numpy.sqrt(numpy.median([(x-data_in[j])**2. for j,x in enumerate(data_in[1:])]))
	#
	# not much to say about the first data point; hopefully it's not a total disaster
	for j,x in enumerate(data_in[1:]):
		if abs(x-data_in[j])>dy0: data_in[j+1]=data_in[j]
	return data_in

def filter_Null(x, *args, **kwargs):
	return x

def filter_extreme_dy_dx(data_in, dy_dx_0=None, dy_dx_factor=10.0):
	# an extreme value filter based on increasing rate, dy/dx (like filter_extreme_dy, but consider dx too...
	# for think in terms of deg C/second?? dunno.
	#
	# data_in: should be like [[x,y], ...]
	#
	if dy_dx_0==None: 
		dy_dx = [(y-data_in[j][1])/(x-data_in[j][0]) for j, (x,y) in enumerate(data_in[1:])]
		dy_dx_0 = dy_dx_factor*numpy.sqrt(numpy.mean([z*z for z in dy_dx]))
	dy_dx_0 = abs(dy_dx_0)
	#
	# note we want to work through this list as it evoloves, so we need to recalculate dy/dx as we 'correct' it.
	r_data = [list(rw) for rw in data_in]	#[data_in[0]]
	for j,(x,y) in enumerate(r_data[1:]):
		if abs((y-r_data[j][1])/(x-r_data[j][0]))>dy_dx_0: r_data[j+1][1]=r_data[j][1]
		#if abs(dy_dx[j])>dy_dx_0: r_data[j+1][1]=r_data[j][1]
		
	#return  list([data_in[0]]) + [rw if abs(dy_dx[j])<dy_dx_0 else data_in[j] for j,rw in enumerate(data_in[1:])]
	#return  list([data_in[0]]) + [rw if abs(dy_dx[j])<dy_dx_0 else [rw[0], numpy.mean(data_in[j][1] + data_in[min(j+2,len(data_in)-1)][1])] for j,rw in enumerate(data_in[1:])]
	#return  list([data_in[0]]) + [rw if abs(dy_dx[j])<dy_dx_0 else [rw[0], [rw[0], data_in[j][1]] for j,rw in enumerate(data_in[1:])]
	return r_data
	
	
def filter_large_dy_dx(data_in, dy0=None, dy_percent=.8):
	# filter out rows with large change dy/dx.
	# note: data_in are like [[x,y],...]
	# note: this just removes the offending data. nominally, we want to interpolate... so write a hybrid of this filter + filter_extreme_dy()
	#
	# add an index k to re-sort when we're done.
	xydy = [[x,y,(y-data_in[k][1])/(x-data_in[k][0]), k] for k, (x,y) in enumerate(data_in[1:])]
	xydy.sort(key=lambda rw: rw[2])
	#
	j_thresh = int(len(xydy)*dy_percent)
	dy_thresh = xydy[j_thresh]
	#
	#use k to re-sort and return:
	#xydy = sorted(xydy[0:j_thresh], key=lambda rw: rw[3])
	#return [[rw[0], rw[1]] for rw in xydy]
	return [[rw[0], rw[1]] for rw in sorted(xydy[0:j_thresh], key=lambda rw: rw[3])]
#
def filter_med_N(data_in, med_len=5):
	# median filter:
	# ... we can simplify this syntax.
	#n=len(data_in)
	return numpy.array([numpy.median(data_in[max(j+1-med_len,0):j+1]) for j,x in enumerate(data_in)])
#
#
def compound_filter(f_filter, data_in, n=2, *args, **kwargs):
	# apply a filter n times.
	#
	r_val = [x for x in data_in]
	#
	for j in range(n):
		r_val = f_filter(r_val, *args, **kwargs)
	#
	return numpy.array(r_val)
	#
#
def filter_hi_lo(data_in, f_len=5, lo=.2,hi=.8):
	# a "high-low" or "80-20" type filter. identify the hi, lo percentile threshold. truncate values at that threshold.
	#
	# - i think there is nothing for this but to make a copy of the data, so it might be a bit cumbersome for
	# long sequences.
	# - also, we introduce a lot of overhead to optimally handle the first f_len of the sequence. for long sequences, it 
	# might be better to introduce a compromise, like maybe a median filter for n<f_len, then the hi-lo filter. otherwise,
	# we have to re-calc and copy a bunch of variables at each step.
	#
	lo=min(hi,lo)
	hi=max(hi,lo)
	#
	data_out=[]
	k_0 = int(numpy.floor(lo*(f_len-1)))
	k_1 = int(numpy.ceil(hi*(f_len-1)))
	#print('range: ', k_0, k_1)
	#for j_end in range(1,len(data_in)+1):
	for j_end,x in enumerate(data_in):
		#
		j_end += 1
		#
		sub_seq = [xx for xx in data_in[max(0,j_end-f_len):j_end]]
		#sub_seq_0 = sub_seq[:]
		# this will be cumbersome for long sequences. but these values become fixed for n>f_len, so let's just drop an "if" in here:
		#if len(sub_seq)<f_len:
		if j_end<=f_len:
			l_sq = len(sub_seq)-1
			k_0 = int(numpy.floor(lo*l_sq))
			k_1 = int(numpy.ceil(hi*l_sq))
		
		#
		sub_seq.sort()
		x_min = sub_seq[k_0]
		x_max = sub_seq[k_1]
		#
		#x_out=x
		#x_out = min(x_max,x_out)
		#x_out = max(x_min,x_out)
		#data_out += [x_out]
		#
		data_out += [max(x_min,min(x_max,x))]
		#
		#if data_out[-1]!=x or True: print('[%d] sq, k_0, k_1: ' % j_end, sub_seq, k_0, k_1, x_min, x_max, x, data_out[-1], len(sub_seq))
	#
	#print('final j_end: ', j_end, lo, hi, f_len)
	return numpy.array(data_out)

# 
def filter_mean_N(data_in, med_len=5):
	# median filter:
	# @med_len: length for averaging; using "med_len" instead of "mean_len", or something, to be consistent with filter_med_N()
	#return [numpy.mean(data_in[max(0,j-med_len):j]) for j in range(len(data_in))]
	#n=len(data_in)
	#return [data_in[0]] + [numpy.median(data_in[max(j-med_len,0):j]) for j in range(1,n)]
	return [numpy.mean(data_in[max(j+1-med_len,0):j+1]) for j,x in enumerate(data_in)]
	
def filter_start_end(data_in, start_index=0, stop_index=None):
	'''
	# a filter function. the idea is that a list or dict (or container in general) of filters can be passed to a function or class object; the function then
	# runs each filter in order to augment the data. this is a crayon simple. we'll have to (maybe) refine the input parameter and calling standard, but
	# note that this model is pretty modular; an execution list can look like FP={'filter':filter_function, 'params':{dict of params}}
	# the calling script would look like for filter, params in FP.items(): X = filter(X, **params)
	# and of course X can be included in FP like {'data_in':X}. we can also pass as a list: [[filter_function, {params}], ...]
	'''
	#
	# note that X[0:None] should return the whole list -- at least in Python 2.7 and 3.4
	start_index = min(len(data_in)-1, start_index)
	return data_in[start_index:stop_index]
