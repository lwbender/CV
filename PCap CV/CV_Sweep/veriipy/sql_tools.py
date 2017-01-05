'''
This code is property of LANDAUER Inc. and may not be replicated without permission.
#
# primary author:
#  mark yoder: myoder@landauer.com
#
# synopsis:
# some sample/training code and/or tools for connecting python to MS SQL and MySQL server(s).
#
# Getting pymssql for Anaconda Python:
#  pymssql: python microsoft SQL (mssql) driver
#  mssql is not as commonly used with Python and other non MS systems as MySQL and other OpenSource projects
#  that prioritize cross-platform compatibility. pymssql is not available through the direct, official "conda"
#  management system, there are -- however, several pymssql packages available. here's how to get it:
#  (how it should work...)
#  1) hey, maybe something changed. just try to install it:
#    > conda install pymssql
#  2) if that doesn't work (and it probably won't), follow conda's good counsel:
#    > anaconda search -t conda pymssql
#  3) you'll be given some third party options; you can try to install one like:
#    > conda install -c yikelu pymssql=2.0.1
#   ("-c" means "channel", yikelu is somebody on GitHub, pymssql is the project and the "=2.0.1" is the version,
#    but it's optional.
#  4) BUT, these all seem to be Python2 packages, and so far all of these proper approaches have failed, so i just did this:
#  ***** and this works: *****
#  5) pip install pymssql
#   you can also download the pymssql package (google it; there's a pymssql.org site) and install the downloaded package
#   using pip or pip3 (there is a syntax to use pip to install from a local folder).
#
# mysql.connector:
# install mysql.connector like:
#  > conda install mysql-connector-python
'''
#
import getpass
import pymssql
#import mysql.connector
import json

import os
import datetime as dtm
import pytz

def get_mssql_connection(user_id=None, db_name=None, db_server=None, bb_name=None, db_domain='LANDAUER', do_test=True, verbose=False):
	'''
	# get a connection to an mssql database; user provides inputs as needed. start with the veriify fieldtrials service as defaults...
	# we can change the default behavior as desired, but for now, let's:
	#  - assume that if a parameter is provided (not None), we know it's what we wanted, so don't prompt.
	#  - it might be smart to always prompt for password, so we dont' go typing our passwords into viewable screens, but we can leave that to operator discretion.
	#  - what else?
	#
	'''
	#
	db_name_default   = 'vedb_t01'
	db_server_default = 'loqdb240.int.landauerinc.com'
	db_domain_default = 'LANDAUER'
	user_id_default = '{}\\{}'.format((db_domain or db_domain_default), getpass.getuser())
	#
	db_server    = db_server
	db_server = (db_server or db_server_default)
	db_name      = db_name
	db_name   = (db_name or db_name_default)
	db_domain    = db_domain 
	db_domain = (db_domain or db_domain_default)
	user_id      = user_id   
	user_id = (user_id or user_id_default)
	#
	passwd = getpass.getpass('Password for ({}): '.format(user_id))
	#
	connect_data = {'user':user_id, 'server':db_server, 'database':db_name, 'password':passwd}
	if verbose: print('***')
	#for key,val in connect_data.items(): print('{}: {}'.format(key,val)
	for key,val in connect_data.items():
		print('{}: {}'.format(key,(val if not key in ['password'] else '******')) )
	#
	if verbose: print('creating pymssql conection object: ')
	#
	con_mssql = pymssql.connect(**connect_data)
	#
	# test conection:
	if do_test:
		csr=con_mssql.cursor()
		csr.execute('select count(*) from fieldtrialmeasurements')
		print('sql_test: ', [rw for rw in csr])
		print('successful test?')
	#
	return con_mssql
#
def test_mssql_connection(con, sql_query=None):
	sql_query = (sql_query or 'select count(*) from fieldtrialmeasurements')
	csr = con.cursor()
	#
	try:
		csr.execute(sql_query)
		for rw in csr:
			print('sql_test successful; row_count={}'.format(rw))
		#
	except:
		print('sql connect failed on query = \'{}\''.format(sql_query))
	#
#
	
#
