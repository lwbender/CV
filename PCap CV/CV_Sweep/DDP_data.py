# this should now be independent of dosimeter_board; dosimeter_board can (presumably) be depricated.
#
#from dosimeter_board import  *
from pt100 import PT1000
from pt100 import PT100
from pt100 import Thermistor_exp
#
# standard imports:
import numpy
import scipy

import math
import pylab as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

import os
import sys
import random
#
import multiprocessing as mpp
#
import datetime as dtm
import pytz
import matplotlib.dates as mpd
# custom stuff:
import pt100
from pt100 import PT1000
from pt100 import PT100
from pt100 import Thermistor_exp
# note: veriipy is a "package", likely included in the repository as a git "submodule"
import veriipy
import veriipy.simple_filters as sfp
import veriipy.data_type_wranglers as dtwp
import veriipy.sql_tools as sqp
#
#colors_ =  mpl.rcParams['axes.color_cycle']
colors_ = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
#
# temporary default assignments:
default_data_v05 = 'data/consolidated/F518E3BAC5EA_ddp_export_p05.csv'
# older format (and there are others older than this):
default_data_v00 = 'data/consolidated/9C194F09_ddp_export_0615_inv_disc.csv'
#
class DDP_Data_Base(object):
	# a toy class to pretend to be the DDP_data object.
	# eventually, we'll integreate these components into Board_Data(). for now, i think, we'll just 
	#use it to make a workflow demo.
	#
	# note also: there's so much legacy code in Board_Data that it's probably faster to just port that code into 
	# this new code-base, as needed
	#
	# we're going to, for the time being, work off of a presumed configuration:
	# PC0 is the refcap value
	# PC1-6 are MOScaps (or equivalent)
	# PC7 is a thermistor, though it will probably be going away... and replaced by a MOScap in inversion?
	#
	def __init__(self, MOScap_cols=None, pt_100_R0=None, pt_100_R100=None, pt_100_A=None, pt_100_B=None, pt_100_C=None, 
				 pt_100_T_min=None, pt_100_T_max=None, data_start_index=0, data_stop_index=None, 
				 raw_conversion_factor=2**-21, zetas_filter_len=24, C0=220000., cols_fitting_model=None, cols_therm={}, **kwargs):
		#
		#self.__dict__.update(locals())
		#
		self.__dict__.update(locals())
		if MOScap_cols is None: MOScap_cols = ['PC%d' % k for k in range(1,7)]
		self.MOScap_cols = ['PC%d' % k for k in range(1,7)]
		#
		# TODO: set default filer parameters from expected measurement interval; we'll typically want a 24 hr. interval.
		self.default_hi_lo_params = {'f_len':25, 'hi':.9, 'lo':.1}
		#
		# note: this default will probably change; we expect to have a fitting device in one of the MOScap slots (maybe PC1), but we don't really know which one yet.
		# ... or if we'll have a T_Int...
		if cols_fitting_model is None: cols_fitting_model=['T_Int', 'zetas', 'PC1', 'time']
		self.cols_fitting_model=cols_fitting_model
		#
		self.therm_cols = ['PC7']
		self.PT1000_cols = ['T_Int']
		self.numerical_cols = ['T_Int', 'T_Ext', 'PC0','PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7']
		self.reff_cap_col = 'PC0'
		#
		# PT1000 parameters (None's will use the default values in the PT100(0) class):
		self.pt_100_parameters = {'R0':pt_100_R0, 'R100':pt_100_R100, 'A':pt_100_A, 'B':pt_100_B, 'C':pt_100_C, 'T_min':pt_100_T_min, 'T_max':pt_100_T_max}
		self.pt1000 = PT1000(**self.pt_100_parameters)
		self.pt100 = self.pt1000		# legacy nomenclature...
		self.Thermistor = Thermistor_exp()		# ...and for now, just use default values.
		#
		# this will call the load_data_to_raw() for the various subclasses. we need a default subclass; maybe it will inherit from Board_Data()
		# probably all the inputs, inits, etc. need to be reconsidered.
		# ... or can we stop calling this board_data, and call it something not dumb? ... later. it will be easier to do a global
		# replace of board_data than datas...
		self.board_data = self.load_data_to_raw(r_type='recarray')
		#
		# did we get anything? if we don't pass a data source (a filename or sql query data), 
		# we'll return None, [], or something like that.
		# TODO: (this can probably be one "if" statement... later)
		if not self.board_data is None and not len(self.board_data)==0: self.process_raw_data()
		if not self.board_data is None and not len(self.board_data)==0: self.calc_shortcuts()

	def process_raw_data(self):
		self.raw_to_ratios()
		self.ratios_to_physical()
		#
		# calculate "zetas", liner term from MOScap-refcap analysis
		# we might not always want to do this, but we'll figure that out later...
		# now, any remaining pre-calculated values:
		self.zetas = self.get_zetas(col_zeta=self.reff_cap_col, filter_lo=.1, filter_hi=.9, filter_len=self.zetas_filter_len, mean_start=.33, mean_end=1.)
		#
		# TODO:
		# evaluate how to best implement extra smoothing for zetas.
		# and zetas are noisy; let's add a level of median filtering.
		#self.zetas = sfp.filter_med_N(self.zetas, self.zetas_f_len)
		#
	def calc_shortcuts(self):
		# elapsed time and intervals shortcuts:
		# if hasattr(self, 'dtype'):
		#
		self.min_f_time=min([rw[-1] for rw in self.board_data])
		#self.min_f_time = min(self['f_date'])
		self.elapsed_time = numpy.array([0.] + numpy.array([rw[-1]-self.min_f_time for rw in self.board_data]))
		self.intervals = numpy.array([t-self.elapsed_time[j] for j,t in enumerate(self.elapsed_time[1:])])
		#
	#
	@property 
	def intervals_seconds(self):
		return self.intervals*3600*24
	@property
	def elapsed_time_seconds(self):
		return numpy.array(self.elapsed_time)*3600.*24.
	@property
	def elapsed_time_min(self):
		return numpy.array(self.elapsed_time)*60.*24.
	#
	def raw_to_ratios(self,data_in=None, raw_conversion_factor=None, C0=None):
		# for now, expect a recarray...
		# just
		#
		if data_in is None: data_in=self.board_data
		#
		raw_conversion_factor = (raw_conversion_factor or self.raw_conversion_factor)
		raw_conversion_factor = (raw_conversion_factor or 2.**-21)
		#
		for col in self.numerical_cols: data_in[col]*=(raw_conversion_factor)
		#	
	#
	def ratios_to_physical(self, data_in=None, C0=None, MOScap_cols=None):
		#
		# TODO: rename T_Int to T_pt (or something like that)... or maybe the PT1000 is going away entirely?
		#
		C0 = (C0 or self.C0)
		C0 = (C0 or 220000.)
		#
		if MOScap_cols is None: MOScap_cols=self.MOScap_cols
		#
		# first, convert moscap cols to capacitance:
		for col in MOScap_cols: self[col]*=C0
		#
		# now, PT1000:
		for col in self.PT1000_cols:
			# this will almost always be just 1 col, 'T_Int', but making this a list defined in __init__ gives us some flexibility.
			for j,x in enumerate(self[col]): self[col][j] = self.pt1000.get_temperature_from_ratio(x)
		#
		for col in self.therm_cols:
			for j,x in enumerate(self[col]):self[col][j] = self.Thermistor.get_T(x)
		#
	#
	def get_zetas(self, col_zeta = 'PC0', filter_lo=.1, filter_hi=.9, filter_len=None, mean_start=.33, mean_end=1.):
		# helper function to get "zetas", or the transformed reference capacitor values.
		# filter_lo/hi: hi,lo for filter_hi_lo()
		# mean_start/end: domain over which to take average to guess RC value.
		#
		# ... but we might not have reff-cap data:
		if not col_zeta in self.keys():
			#
			return None
		#
		filter_len = (filter_len or self.filter_len)
		filter_len = (filter_len or 15)
		#
		# get filtered reff_cap intervals (note we reserve the option for a mean-filter as well as the hi_lo filter).
		C_rs = sfp.filter_mean_N(sfp.filter_hi_lo(self[col_zeta], f_len=filter_len, lo=filter_lo, hi=filter_hi),1)
		#
		# now, the adjustment for linear capacity, C = C0 + \beta V
		# mean <RC>: assume the system "settles", but we can't assume the settling function, so let's
		# just average the last 1/3 of the sequence or so.
		RC = numpy.mean(C_rs[int(mean_start*len(C_rs)):int(mean_end*len(C_rs))])
		return numpy.array([(1.-numpy.exp(-x/RC))*RC/x for x in C_rs])
	#
	#
	def time_to_index(self, t_seconds=0):
		# find index of some time t.
		# future: if rev: search backwards.
		j=0
		for j,t in enumerate(self.elapsed_time_seconds):
			if t>=t_seconds: break
		#
		return j
	#
	#def get_fitting_matrix(self, cols_model=['T_Int', 'PC1', 'zetas', 'time'], col_b=None, j_start=0, j_end=None, f_filt=sfp.filter_hi_lo, filter_params={'f_len':10, 'hi':.9, 'lo':.1}):
	def get_fitting_matrix(self, cols_model=['T_Int', 'PC1', 'zetas', 'time'], col_b=None, j_start=0, j_end=None, f_filt=sfp.filter_Null, filter_params={}):
		#
		# TODO: handle default filter_params from class scope variables (self.filter_hi_lo_params, etc.). maybe we also need to handle a default length, or
		# pass variables as lists (*args, as opposed to **kwargs), but to do that, we need to get better at querying the target function for an input signature.		
		#
		j_end = (j_end or len(self.board_data))
		#
		# we might get some custom columns:
		special_cols_map = {'z':'zetas', 'zeta':'zetas', 'zetas':'zetas', 't':'time', 'time':'time'}
		# (later, we might specify elapsed_time_min/seconds, ''->days.
		#
		cols_basic   = [col for col in cols_model if col in self.dtype.names]
		cols_special = [special_cols_map.get(col,None) for col in cols_model if not col in cols_basic and not special_cols_map.get(col,None) is None]
		#
		# get the basic columns; then add special cols. we might also treat some columns specially here, like T_int.
		#A = zip(*[[1 for _ in range(j_start, j_end)]]+ [f_filt(self.board_data[cl][j_start:j_end], **filt_params) for cl in cols_model] + therm_cols)
		A = [[1 for _ in range(j_start, j_end)]]+ [f_filt(self.board_data[cl][j_start:j_end], **filter_params) for cl in cols_basic]
		if 'zetas' in cols_special: A+=[self.zetas[j_start:j_end]]
		if 'time'  in cols_special: A+=[self.elapsed_time_seconds[j_start:j_end]]
		#
		# any others?
		A = list(zip(*A))
		#
		if col_b==None:
			return list(A)
		else:
			# if it's just one column, passed as just a column name, not a list:
			if isinstance(col_b, str):
				return list(A),f_filt(self.board_data[col_b][j_start:j_end])
			else:
				# it's a list, presumably of column names.
				# note, you can do a simultaneous leastsquares fit of A,b for b --> [[C00, C01, C02...], [C10, C11, C12,...],...]
				return list(A), list(zip(*[f_filt(self.board_data[col][j_start:j_end]) for col in col_b]))
				
			
		#
	def get_autofit_col_parameters(self, cols=None, cols_model=['T_Int', 'zetas', 'PC1', 'time'], start_index=0, min_len_factor=.5, step_size=1,
	 chi_expon=1.5, fignum=None, do_clf=None, ax=None, show_all_lsqs=False, n_cpu=None, f_filter=sfp.filter_hi_lo, filter_params={'f_len':1, 'lo':0., 'hi':1.}, j_start=0, j_end=None):
		'''
		# run an auto-fitter to get best fit parameters for each column.
		# @cols: cols to fit. usually all of them.
		# @cols_model: model data, aka C_model = numpy.dot(cols_model,params). this can include special fields, like "zetas" (from "Measurement_Reference_Cap"), various thermistors etc.
		#              (most of which have gone away) and maybe we'll handle the inverse sensor here as well.
		'''
		#
		if cols_model is None: cols_model = self.cols_fitting_model
		if cols is None: cols=self.MOScap_cols
		#
		if f_filter is None:
			f_filter=sfp.filter_Null
			filter_params = {}
		#
		#print('cols: ', cols)
		# ... and we'll use veriipy.optimizer auto_optimizer...
		A,B = self.get_fitting_matrix(cols_model=cols_model, col_b=cols, f_filt=f_filter, filter_params=filter_params )
		#
		# save most recent A,B (is there a smarter way to handle this?)
		#self.fitting_A = A
		#self.fitting_B = B
		#
		print('B-data: ', len(B), numpy.shape(B))
		# B returns in a transposed matrix (for output vectors b0,b1,b2..., we get B = [[b00,b10,b20,...],[b01,b11,b21], ...]
		# if we pass B to the fitting algorithm, it will give an aggregated fit for all coluns. we want to pass each (transposed) vector.
		#best_fits=None
		best_fits = veriipy.optimizers.auto_optimize_multiple_arrays(A=A,bs=list(zip(*B)), min_len_factor=min_len_factor, step_size=step_size, start_index=start_index, chi_expon=chi_expon,
		fignum=fignum, do_clf=do_clf, ax=ax, show_all_lsqs=show_all_lsqs, do_aggregate=True, n_cpu=n_cpu)
		#
		return {col:bf for col,bf in zip(cols,best_fits)}
	#
	def export_to_csv(self, fout='dosimeter_data.csv'):
		# export to CSV...
		with open(fout,'w') as f:
			f.write('#Dosimeter data export\n#data_file:\t{}\n'.format(self.data_in))
			f.write('#!%s\n' % '\t'.join(list(self.board_data.dtype.names) + ['elapsed_time_sec', 'zetas']))
			#[f.write('%s\n' % str(list(self.elapsed_time_seconds)[j])) for j,rw in enumerate(self.board_data))]
			# '\t'.join([str(x) for x in rw]+
			[f.write('%s\n' % ('\t'.join([str(x.decode() if hasattr(x,'decode') else x) for x in rw] + [str(self.elapsed_time_seconds[j]), str(self.zetas[j])]) )) for j,rw in enumerate(self.board_data)]
		#
		return None
	#

	#
	# trailer stuff: general class behavior definitions and redefinitions of built-in functions.
	# make this class act like its data structure member:
	def __iter__(self): return iter(self.board_data)
	def __getitem__(self,x): return self.board_data.__getitem__(x)
	def __setitem__(self,j,x): return self.board_data.__setitem__(j,x)
	def __len__(self): return self.board_data.__len__()
	def keys(self): return self.board_data.dtype.names
	def items(self): return list(zip(*self.board_data.tolist()))		# might need to write this as: list(board_data), for compatibility.
	@property
	def dtype(self): return self.board_data.dtype
	#

	#
#
class DDP_Data_sql(DDP_Data_Base):
	def __init__(self, ble_id='F518E3BAC5EA', start_date=None, end_date=None, sql_and='', col_map=None,
				 sql_select=None, sql_con=None, *args, **kwargs):
		# do stuff...
		#
		#self.default_col_map = {'tempint': 'T_Int', 'tempext':'T_Ext', 'refcap':'PC0', 'c39':'PC1_raw', 'c40':'PC2_raw', 'c41':'PC3_raw', 'c42':'PC4_raw', 'c43':'PC5_raw', 'c44':'PC6_raw', 'c46':'PC7_raw', 'measurementtimestamp':'timestamp', 'servertimestamp':'creation_date'}
		self.default_col_map = {'tempint': 'T_Int', 'tempext':'T_Ext', 'refcap':'PC0', 'c39':'PC1', 'c40':'PC2', 'c41':'PC3', 'c42':'PC4', 'c43':'PC5', 'c44':'PC6', 'c46':'PC7', 'measurementtimestamp':'timestamp', 'servertimestamp':'creation_date'}
		if col_map is None: col_map=self.default_col_map
		#
		# we get some bogus timestamps of Unix_day_one (1969-1-1), so let's set a min_date to at least the start of this project.
		# what's the correct format???
		#start_date=(start_date or dtm.datetime(2016,1,1, tzinfo=pytz.timezone('UTC')))
		#
		self.__dict__.update(locals())
		#
		super(DDP_Data_sql, self).__init__(*args, **kwargs)
	#
	def load_data_to_raw(self, ble_id=None, start_date=None, end_date=None, sql_and='', col_map=None, sql_select=None, r_type='dict', sql_con=None):
		#
		# allow a call to retun no data... though maybe we make a special subclass for that purpose later...
		if ble_id == '': return None
		#
		ble_id = (ble_id or self.ble_id)
		start_date = (start_date or self.start_date)
		end_date = (end_date or self.end_date)
		sql_and = (sql_and or self.sql_and)
		if col_map is None: col_map = self.col_map
		if not hasattr(col_map, 'keys'): col_map = {x:x for x in col_map}
		col_map_inv = {val:key for key,val in col_map.items()}
		#
		sql_select = (sql_select or self.sql_select)
		r_type=(r_type or self.r_type)
		sql_con = (sql_con or self.sql_con)
		#
		# now, construct a sql query:
		# ... we could just get the raw json data... and we eventually will probably. but for now, let's just get the tabulated bits.
		# and using structured sql can be tricky, so in some cases, just substitute a string... and i'm not sure how well structured sql works
		# with pymssql, so this might be sloppy...
		#
		if sql_select is None:
			cols_str = ', '.join(['%s as %s' % (x,y) for x,y in col_map.items()])
			#
			sql_select = 'select {cols} from fieldtrialmeasurements where dosimeterid=\'{ble}\' '.format(cols=cols_str, ble=ble_id)
			if not start_date is None: sql_select = sql_select + ' and measurementtimestamp>=\'{}\' '.format(start_date)
			if not end_date   is None: sql_select = sql_select + ' and measurementtimestamp<=\'{}\' '.format(end_date)
			sql_select = sql_select + ' ' + sql_and
		#print("sql_select: ", sql_select)
		if sql_con==None:
			sql_con = sqp.get_mssql_connection()
		#
		csr=sql_con.cursor()
		csr.execute(sql_select)
		#
		# let's just assume we have plenty of memory...
		datas=csr.fetchall()
		#
		#print('dtas:\n', datas[0:5])
		#
		# this is sloppy, but it works:
		for j,rw in enumerate(datas):
			datas[j]=list(rw)
			for k,x in enumerate(rw):
				datas[j][k] = dtwp.guess_type(x, int_to_float=True)
		#print('dtas2:\n', datas[0:5])
		#
		#print('\n\n')
		col_names = [rw[0] for rw in csr.description]
		#date_index = col_names.index(col_map_inv['timestamp'])
		date_index = col_names.index('timestamp')
		#
		datas = list(zip(*zip(*datas), [mpd.date2num(rw[date_index]) for rw in datas]))
		#for j,rw in enumerate(datas):
		#	datas[j] += [mpd.date2num(rw[date_index])]
		col_names += ['f_date']
		#
		datas.sort(key=lambda rw: rw[-1])
		#
		if r_type=='dict':
			return {col:data for col,data in zip(col_names, zip(*datas))}
		elif r_type=='recarray':
			# for now, assume all floats:
			# ... except the timestamp/datetime field(s) which must be converted to numpy64.
			d_formats = []
			for j,x in enumerate(datas[0]):
				if isinstance(x, str): d_formats += ['S%d' % max(1,max([len(rw[j]) for rw in datas]))]
				elif isinstance(x, dtm.datetime): d_formats += ['M8[us]']
				else: d_formats += [type(x).__name__]
			#
			#print("col_names: ", col_names)
			return numpy.core.records.fromarrays(list(zip(*datas)), names=col_names, formats=d_formats)
		else:
			return [col_names] + datas	   
#
#
class DDP_Data_recary(DDP_Data_Base):
	# ... given a recarray:
	# (note; for this version, given a recarray in which the pre-calculations (aka, raw to physical, etc.) are done, aka, a copy of the object during the workflow.
	# TODO: this sub-class is not yet tested... at all.
	def __init__(self, data_ary, *args, **kwargs):
		#
		self.board_data = data_ary.copy()
		super(DDP_Data_recary, self).__init__(*args, **kwargs)
	#
	def load_data_to_raw(self, *args, **kwargs):
		# write a null load_data function
		return self.board_data
	def process_raw_data(self, *args, **kwargs):
		# overwrite "process raw"; raw data have been processed. maybe we put o condition on this so we can alternately
		# pass an array of raw vals?
		self.zetas = self.get_zetas(col_zeta=self.reff_cap_col, filter_lo=.1, filter_hi=.9, filter_len=self.zetas_filter_len, mean_start=.33, mean_end=1.)
		return None
	#
class DDP_Data_csv(DDP_Data_Base):
	def __init__(self, data_csv=default_data_v05, csv_delim=',', col_map=None, *args, **kwargs):
		# ... and the remaining parameters are common to all Board_Data (or we'll rename it? maybe DDP_data,
		# DDP_handler?) objects. we can either cut-paste call signature(s) or use *args, **kwargs...
		# do stuff...
		self.data_csv=data_csv
		self.csv_delim=csv_delim
		# "_raw" denotes the format of the input columns; we're going to map them to not "_raw"
		#if col_map is None or len(col_map)==0: col_map = ['T_Int', 'T_Ext', 'PC0_raw','PC1_raw', 'PC2_raw', 'PC3_raw', 'PC4_raw', 'PC5_raw', 'PC6_raw',
		#							   'PC7_raw', 'timestamp', 'creation_date']
		#
		self.col_map=col_map
		#
		super(DDP_Data_csv, self).__init__(*args, **kwargs)
	#
	def load_data_to_raw(self, data_csv=None, delim=None, j_start=None, j_stop=None, col_map=[], r_type='dict'):
		# you know... sqlLite would probably be a pretty good way to have done this...
		#
		# we'll want these cols:
		# T_Int, T_Ext, PC0_raw,PC1 _raw, PC2_raw, PC3_raw, PC4_raw, PC5_raw, PC6_raw, PC7_raw
		# @ col_map: map (dict) of columns to extract from .csv and their exported names. we'll provide defaults, and if
		#   a list is passed, then convert it to an identity-map:
		#
		# cols: do we want the BLE Serial number on a row-by row basis? is is a non-numerical type; we could keep one copy
		#		but we also want to be sure we have the right one... for sql queries, we'll almost certainly filter on that val,
		#		so for now, let's assume it's been pre-filtered and skip it. if necessary, we'll add a filter here.
		# col_map is like: {'this_col_name':'standard_col_name', ...}
		#
		# TODO: catch the case where a recarray is passed, and just return. also write an subclass that accepts a list
		# of data + cols.
		#
		
		#
		# yoder: move this map assignment to post-read-col names. this way, we can load up some default col_maps for compatibility.
		#col_map = (col_map or self.col_map)
		#if col_map is None or len(col_map)==0: col_map = ['T_Int', 'T_Ext', 'PC0_raw','PC1_raw', 'PC2_raw', 'PC3_raw', 'PC4_raw', 'PC5_raw', 'PC6_raw',
		#							   'PC7_raw', 'timestamp', 'creation_date']
		## is col_map a list (or not-dict)?
		#if not hasattr(col_map, 'keys'): col_map = {x:x.replace('_raw', '') for x in col_map}
		#col_map_inv = {val:key for key,val in col_map.items()}   # an inverse map, so we can look up an input column name.
		#
		delim	= (delim or self.csv_delim)
		data_csv = (data_csv or self.data_csv)		
		#
		# allow a null data set (then we can load the data semi-manually)
		if data_csv =='' or data_csv is None or not os.path.isfile((data_csv or '')): return None
		#
		# start/stop indices:
		# first, have they been provided externally?
		j_start = (j_start or self.data_start_index)
		j_stop  = (j_stop or self.data_stop_index)
		#
		# now, be sure they make sense:
		j_start = (j_start or 0)
		# j_stop: we don't yet know how big the file is, so let's just handle this in the loop.
		
		self.data_input_src = data_csv
		#
		# open the csv file:
		# first row is column headers
		# second row is header-data (Discharge REsistor, cytle Time, Number Average, Time Between Measurements, VDD, REf_C Rated Value).
		# for now. let's keep the long names and replace ' ' with '_'.
		#
		# notes: the first row is like: data,data,data,,,,,,,	# where ,,,,, is a bunch of empty sites.
		#	subsequent rows are like: ,,,,, data,data,data...	# (repeat above).
		#
		board_data = []	# later we'll cast as a recarray. 
		#
		# open the file and load data...
		# guess the delimiter? ... but let's integrate this with the general reading of the file.
		#if delim==None:
		with open(data_csv, 'r') as f_in:
			#################
			# first, gather meta-data, header-data, and other overhead things:
			# for now, assume the format is like: comments, col_names (not marked with #! or anyting, headers, datas.
			this_rw='#'   #we'll be spinning through the input, skipping comments.
			while this_rw[0] in ('#', chr(10), chr(13), chr(9), chr(32)):
				this_rw=f_in.readline()
			#
			# end of comments. this shouild be the col names:
			if delim==None:
				for d in (chr(9), ',', chr(32)):
					# tab, ",", {space} -- in that order.
					if d in this_rw:
						delim = d
						#print('delimter chosen: %s/%d' % (d, ord(d)))
						break
			#
			# now, we know the delimiter and we should be looking at the column names.
			#
			# in older versions, we have 2 cols named T_Int and T_Ext. the first is the "raw" values; keep those.. while preserving the sequence
			col_names = [cn.strip() for cn in this_rw.strip().split(delim)]
			#
			# if a column name is repeaded, keep the first one:
			#col_names_index = {cn:j for j,cn in enumerate(col_names)}	   # use this to pull data from rows by col.
			col_names_index = {}
			for j,cn in enumerate(col_names):
				if not cn in col_names_index.keys(): col_names_index[cn]=j
			#
			# define a col_map, mapping these col_names to standard col names.
			col_map = (col_map or self.col_map)
			if col_map is None or len(col_map)==0:
				if 'C39_raw' in col_names:
					# older format:
					#col_map = {'T_Int':'T_Int', 'T_Ext':'T_Ext', 'Measurement Reference Cap':'PC0_raw', 'C39_raw':'PC1_raw', 'C40_raw':'PC2_raw',
					# 'C41_raw':'PC3_raw', 'C42_raw':'PC4_raw', 'C43_raw':'PC5_raw', 'C44_raw':'PC6_raw', 'C46_raw':'PC7_raw', 'timestamp':'timestamp', 'creation_date':'creation_date'}
					col_map = {'T_Int':'T_Int', 'T_Ext':'T_Ext', 'Measurement Reference Cap':'PC0', 'C39_raw':'PC1', 'C40_raw':'PC2',
					 'C41_raw':'PC3', 'C42_raw':'PC4', 'C43_raw':'PC5', 'C44_raw':'PC6', 'C46_raw':'PC7', 'timestamp':'timestamp', 'creation_date':'creation_date'}
				else:
					# newer format:
					col_map = ['T_Int', 'T_Ext', 'PC0_raw','PC1_raw', 'PC2_raw', 'PC3_raw', 'PC4_raw', 'PC5_raw', 'PC6_raw',
										   'PC7_raw', 'timestamp', 'creation_date']
			# is col_map a list (or not-dict)?
			#print('col_map: ', col_map)
			if not hasattr(col_map, 'keys'): col_map = {x:x.replace('_raw', '') for x in col_map}
			col_map_inv = {val:key for key,val in col_map.items()}   # an inverse map, so we can look up an input column name.			
			#
			# what are the indices of the columns we'll be using?
			active_cols, active_col_indices = list(zip(*sorted([[key,col_index] for key, col_index in col_names_index.items()if key in col_map.keys()], key=lambda rw:rw[1])))
			#
			#print('active_cols: ', active_cols)
			#print('active_col_indices: ', active_col_indices)
			#
			# for some reason, the management of these headers is constantly in flux, so be prepared to manage this 
			#headers/data parsing bit on a regular basis.
			# ... and honestly, we do ****-all with these, so the most important thing to do here is skip the line.
			header_vals = [dtwp.guess_type(x) for x in f_in.readline().strip().split(delim) if not x in ('', None)]	  # and maybe other white-space...
			#header_len = len(header_vals)
			# (alternatively, take this from the number of leading blanks in the data rows -- which seems to be more dependable...)
			header_len=6
			#
			while len(header_vals)<header_len: header_vals+=[0.]
			# in older versions, we tack C_reff onto the end of this, but it will be smarter to just nest it into the parent class.
			# we won't deal with it here.
			#
			header_cols = [cl.strip().replace(chr(32), '_') for cl in col_names[0:header_len]]
			data_cols   = [cl.strip().replace(chr(32), '_') for cl in col_names[header_len:]]
			#
			# sometimes we get repeated column names:
			for c_n in set(data_cols):
				#
				# generalize for multiple redundancies:
				col_count=0		 # counter for number of redundant columns.
				#if data_cols.count(c_n)>1:
				while data_cols.count(c_n)>1:
					#data_cols[data_cols.index(c_n)] = data_cols[data_cols.index(c_n)] + '_raw'
					data_cols[data_cols.index(c_n)] = data_cols[data_cols.index(c_n)] + '%d' % col_count
					col_count+=1
				#
			#
			# a little trick to retrieve the data like X.col1:
			# (except that this is going to be an intermediate object, so we can skip it for now)
			#self.__dict__.update({col:val for col,val in zip(header_cols, header_vals)})
			#
			# we need the column index of the measurement date... whose name keeps changing.
			# this is going to be a sloppy mess for the time being...
			#
			# let's assume that, moving forward, we'll be dealing only with 'timestamp' col and timestamp type.
			# for backwards compatibility in dates, see original dosmiter_board.py code. search for "if record_data_type"
			#
			######################
			# Process Data #######
			######################
			# so now, we're to the data, and we start counting rows.
			for j_rw, rw in enumerate(f_in):
				#
				# we're breaking on header_len, so we need to detect headers in each row.
				#print('rw: ', '**', delim, '**', rw, ' ** ', header_len)
				#
				if rw[0] in (chr(32),chr(9), chr(13), chr(10), '#'): continue
				if rw[header_len+2]==',': continue
				#				
				l0 = len(rw.split(delim))
				#
				# so going forward, do we want all the columns, or just the cols of interest?
				# it's going to be a hassel to handle the versioning. for now, let's just gather everything. we'll end up
				# constructing a reduced set in the parent class.
				#board_data += [[dtwp.guess_type(x.strip(), int_to_float=True) for x in rw.strip().split(delim)[header_len:]]]
				#
				# .. but we only want "active" indices:
				rws = rw.strip().split(delim)[0:]
				#
				# ... and actually, we get back plain string values from fieldtrials as well; it might have been smarter to just load the col.
				# data and ship off to one place for type conversion...
				board_data += [[dtwp.guess_type(rws[k].strip(), int_to_float=True) for k in active_col_indices]
							   + [mpd.date2num(dtm.datetime.utcfromtimestamp(float(rws[col_names_index[col_map_inv['timestamp']]])))]]
				#
				# now, get the measurement data and add a float_date
				#record_dtm = rws[col_map[col_map_inv['timestamp']]]   # look up what the 'timestamp' input col name was; get that index.
				#board_data[-1]+=[mpd.date2num(datetime.datetime.utcfromtimestamp(record_dtm))]
				#
				# we might get an empty row (for some reason)
				while len(board_data[-1])==0:
					#board_data.pop()
					board_data = board_data[:-1]
					continue
				#
			#
		#, now sort and chop the data. then massage things that need massaging here (like date types):
		# now, massage the date (and other early-massaging) and calculate an integet "f_date" field... which is the same_ish as timestamp, except in day-units.
		#
		# for backwards compatibility, see dosimeter_board.py, Board_Data.load_data_csv() or .load_data(). hopefully, however, we can leave this behind us.
		board_data.sort(key=lambda rw:rw[-1])
		board_data = board_data[j_start:j_stop]
		#
		# now, save/return the data as a dict object: {col_name: col_data, ...}
		# the main idea is that, in this format, it will be easy to manipulate, add, remove columns... easier than a recarray...
		# it might be time to start using PANDAS.
		#out_cols = list(active_cols) + ['f_date']
		out_cols = [col_map[col] for col in active_cols] + ['f_date']
		if r_type=='dict':
			return_data = {col:data for col,data in zip(out_cols, zip(*board_data))}
		elif r_type=='recarray':
			# for now, assume all floats:
			d_formats = []
			for j,x in enumerate(board_data[0]):
				if isinstance(x, str): d_formats += ['S%d' % max(1,max([len(rw[j]) for rw in board_data]))]
				elif isinstance(x, dtm.datetime): d_formats += ['M8[us]']
				else: d_formats += [type(x).__name__]
			#
			#return numpy.core.records.fromarrays(list(zip(*board_data)), names=active_cols, formats=d_formats)
			return_data = numpy.core.records.fromarrays(list(zip(*board_data)), names=out_cols, formats=d_formats)
			#
		else:
			return_data = [out_cols] + board_data
		#
		return return_data
#
# now some general procedual functions that don't necessarily require class scope (though we might integrate these later for performance).

						
  
