#!/usr/bin/env python

import sys
from glob import glob
def ArgvToDict(argv_list=sys.argv, required=[], optional={}, verbose=True):
	assert type(argv_list) == list and type(required) == list and type(optional) == dict
	argv_dict = {}
	if '-h' not in argv_list and '--help' not in argv_list:
		if len(argv_list) > 1:
			argv_list = ' '.join(argv_list).replace(' --', ' -').split(' -')[1:] # already deleted sys.argv[0]
			for opt_list in argv_list:
				opt_list = opt_list.split(' ', 1)
				if ',' in opt_list[-1]:
					opt_list = [opt_list[0],] + opt_list[-1].split(',')
				elif '*' in opt_list[-1] or '?' in opt_list[-1]:
					opt_list = [opt_list[0],] + glob(opt_list)
				try:
					opt_list.remove('')
				except ValueError:
					pass
				opt_list_len = len(opt_list)
				# "-option value":
				# -option
				if len(opt_list[0]) == 1:
					opt_list[0] = '-' + opt_list[0]
				elif len(opt_list[0]) > 1:
					opt_list[0] = '--' + opt_list[0]
				# value
				if opt_list_len == 1:
					argv_dict.update({opt_list[0]: True})
				elif opt_list_len == 2:
					argv_dict.update({opt_list[0]: opt_list[-1]})
				elif opt_list_len >= 3:
					argv_dict.update({opt_list[0]: opt_list[1:]})
				# To add optional arguments:
				for opt_argv in optional.iteritems():
					if opt_argv[0] not in argv_dict.keys():
						argv_dict.update({opt_argv[0]: opt_argv[-1]})
					else:
						pass
			if verbose:
#				# To check the arguments:
				all_argv = required + optional.keys()
				unknown_argv_len, lack_argv_len = 0, 0
				# 1. unrecognized
				if all_argv != []:
					unknown_argv = []
					for key in argv_dict.iterkeys():
						if key not in all_argv:
							unknown_argv.append('\"%s\"' % key)
					unknown_argv_len = len(unknown_argv)
					if unknown_argv_len > 0:
						print 'Unrecognized argument' + 's'*int(round((unknown_argv_len-1.0)/unknown_argv_len)) + ': '\
							  + ', '.join(unknown_argv) + '.'
				# 2. required
				if required != []:
					lack_argv = []
					for key in required:
						if key not in argv_dict.keys():
							lack_argv.append('\"%s\"' % key)
					lack_argv_len = len(lack_argv)
					if lack_argv_len > 0:
						print 'Require the argument' + 's'*int(round((lack_argv_len-1.0)/lack_argv_len)) + ': '\
							  + ', '.join(lack_argv) + '.'
				if unknown_argv_len + lack_argv_len > 0:
					print 'Please try again.'
					argv_dict = {}
				# 3. check all args roughly
				else:
					assert len(argv_dict.keys()) == len(all_argv)
			else:
				assert len(argv_dict.keys()) == len(all_argv)
		else:
			pass
	else:
		argv_dict = {'-h': True}
	return argv_dict
	