#!/usr/bin/env pypy

import sys, re, os, gc, commands
from copy import deepcopy
from Argv import ArgvToDict as ATD
from multiprocessing import Pool
from subprocess import check_call as CC
from datetime import datetime as DT

# TODO: trim reads in SAM according to their quality

def ExtractDict(sam_file_chunk):
	#try:
		#chunk = int(sam_file_chunk.split('.')[-1])
	#except ValueError:
		#chunk = 1
	with open(sam_file_chunk, 'r') as inf:
		sam_dict_chunk = {}
		for line in inf:
			# For some unknown resons, Longranger could generate SAM files with all reads
			# prefixed "@". If so, use the second condition which may not support reads' 
			# name only contain 3 charaters.
			# To filer out SAM headers:
			#if line[0] != '@':
			if line[3] != '\t':
			#if line[0] != '@' or line[3] != '@': # More precise and stringent, but slower.
				line = line.strip().split('\t', 1)
				BC = re.findall(BC_pattern, line[-1])
				BX = re.findall(BX_pattern, line[-1])
				if len(BC) == 1 and len(BX) == 1:
					BCBX = '%s-%s' % (BC[0], BX[0])
					sam_dict_chunk_value = sam_dict_chunk.get(BCBX, (0, []))
					count = sam_dict_chunk_value[0] + 1
					seqname_list = deepcopy(sam_dict_chunk_value[-1])
					seqname_list.append(line[0])
					sam_dict_chunk.update({BCBX: (count, deepcopy(seqname_list))})
					#if sam_dict_chunk_value != None:
						#count = sam_dict_chunk_value[0] + 1
						#seqname_list = deepcopy(sam_dict_chunk_value[-1])
						#seqname_list.append(seqname)
						#sam_dict_chunk.update({BCBX: [count, deepcopy(seqname_list)]})
					#else:
						#sam_dict_chunk.update({BCBX: [1, [seqname,]]})
		# debug
		for v in sam_dict_chunk.itervalues():
			assert v[0] == len(v[-1])
		#
	return sam_dict_chunk

def MergeDict(dicts):
	#assert type(dicts) != str and len(dicts) >= 2
	d_merged = deepcopy(dicts[0])
	dicts[0] = ''; gc.collect()
	i = 0
	for d in dicts[1:]:
		assert type(d) == dict
		i += 1
		for item in d.iteritems():
			d_merged_value = d_merged.get(item[0], (0, []))
			count = int(d_merged_value[0]) + int(item[-1][0])
			seqname_list = deepcopy(d_merged_value[-1]) + deepcopy(item[-1][-1])
			d_merged.update({item[0]: (count, deepcopy(seqname_list))})
			#if d_merged_value == None:
				#d_merged.update({item[0]: item[-1]})
			#else:
				#count = int(d_merged_value[0]) + int(item[-1][0])
				#seqname_list = deepcopy(d_merged_value[-1]) + deepcopy(item[-1][-1])
				#assert count == len(deepcopy(seqname_list))
				#d_merged.update({item[0]: [count, deepcopy(seqname_list)]})
		d, dicts[i] = '', ''; gc.collect()
	return d_merged

def File2Dict(in_file):
	d = {}
	with open(in_file, 'r') as inf:
		for line in inf:
			if line[0] != '#':
				line = line.strip().split('\t')
				value = d.get(line[0])
				if value == None:
					count = int(line[1])
					seqname_list = deepcopy(line[2].split(','))
				else:
					count = int(value[0]) + int(line[1])
					seqname_list = deepcopy(value[-1]) + deepcopy(line[2].split(','))
				#assert count == len(seqname_list) # It is wrong, because seqname has duplicates.
				d.update({line[0]: [count, deepcopy(seqname_list)]})
	return d

def HelpMsg():
	print '''
Commands:
    pypy 10xdownsample.py -i SAM/dict [-c cutoff -p threads_num]
	'''	
if __name__ == '__main__':
	args = ATD(argv_list = sys.argv, required = ['-i'], optional = {'-c':5000, '-p':1, '-h':False})
	#print args
	if args.get('-h') or args == {}:
		HelpMsg()
		quit()
	in_file, cutoff, threads = args.get('-i'), args.get('-c'), int(args.get('-p'))
	suffix = in_file.split('.')[-1]
	print '%s' % in_file
	if suffix == 'sam':
		sam_file = in_file
		BC_pattern = re.compile(r'BC:Z:([ATGC]{8})')
		BX_pattern = re.compile(r'BX:Z:([ATGC]{16})-1')	
		if threads > 1:
			print 'Parsing SAM file using %s threads...' % threads
			tp1 = DT.now()
			print '    Assigning jobs...',
			sys.stdout.flush()	
			jobs = Pool(processes = threads)
			chunk_files = commands.getoutput(r'find %s -name "%s.chunk.*"' % \
			                (os.path.split(in_file)[0], os.path.split(in_file)[-1])).split('\n')
			if len(chunk_files) != threads:				
				CC('rm -f %s.chunk.* && split -d -n l/%s %s %s.chunk.' % \
				   (sam_file, str(threads), sam_file, sam_file), shell=True)
				print 'Done.'
			else:
				print 'Chunk files exsit.'
			sam_file_list = []
			for index in xrange(threads):
				sam_file_list.append('%s.chunk.%02d' % (sam_file, index))
			print '    Extracting info...',
			sys.stdout.flush()
			dicts = jobs.map(ExtractDict, sam_file_list)
			print 'Done.'
			print '    Merging...',
			sys.stdout.flush()
			sam_dict = MergeDict(dicts)
			print 'Done.'
			del dicts; gc.collect()
			tp2 = DT.now()
			print '    Done in %ss.' % (tp2-tp1).seconds
		else:
			print 'Parsing SAM file...',
			sys.stdout.flush()
			tp1 = DT.now()
			sam_dict = ExtractDict(sam_file)
			tp2 = DT.now()
			print 'Done.\n    Done in %ss.' % (tp2-tp1).seconds
	elif suffix == 'dict':
		tp1 = DT.now()
		print 'Parsing DICT file...',
		sys.stdout.flush()
		sam_dict = File2Dict(in_file)
		tp2 = DT.now()
		print 'Done.\n    Done in %ss.' % (tp2-tp1).seconds
	else:
		print 'Not supported format: %s' % suffix.upper()
		quit()
	print 'Outputting...',
	sys.stdout.flush()
	tp1 = DT.now()
	# initialize for output
	m = 0
	i = 0
	total_seq = 0
	filter_count = {}
	CC('rm -f %s.dict' % sam_file, shell=True)
	if type(cutoff) == str:
		cutoff = [int(cutoff),]
	for c in cutoff:
		cutoff[i] = int(cutoff[i])
		CC('rm -rf ./%s_trimmed_seqname_list_cutoff%s' % (sam_file, c), shell=True)
		CC('mkdir -p ./%s_trimmed_seqname_list_cutoff%s' % (sam_file, c), shell=True)
		filter_count.update({int(c): (0, 0)})
		i += 1
	#
	for bcbx in sam_dict.iteritems():
		m += 1
		count = int(bcbx[-1][0])
		seqname_list = deepcopy(bcbx[-1][-1])
		# double-check
		if suffix == 'sam':
			assert count == len(seqname_list)
		seqname_list = {}.fromkeys(seqname_list).keys() # delete duplicates
		total_seq = total_seq + len(seqname_list)
		# 1. yield dict file
		if suffix == 'sam':
			with open('%s.dict' % sam_file, 'a') as dictoutf:
				dictoutf.write('%s\t%s\t%s\n' % (bcbx[0], count, ','.join(seqname_list)))
		# 2. yield seqname_query file
		bc = bcbx[0].split('-')[0]
		for c in cutoff:
			assert type(c) == int
			if count >= c:
				filter_count[c] = ((filter_count.get(c)[0] + 1),\
				                   (filter_count.get(c)[-1] + len(seqname_list)))
				with open('./%s_trimmed_seqname_list_cutoff%s/%s' % (sam_file, str(c), bc), 'a')\
				     as seqnameoutf:
					seqnameoutf.write('\n'.join(seqname_list) + '\n')
	print 'Done.'
	# release memory
	del sam_dict; gc.collect()
	for fc in filter_count.iteritems():
		print '    cutoff=%s -> ' % fc[0] + \
		      '%.2f' % (100-fc[-1][0]*100.0/m) + '% of 10XG barcodes ' + \
			  'along with %.2f' % (100-fc[-1][-1]*100.0/total_seq) + '% of reads were filtered out.' 
	tp2 = DT.now()
	print '    Done in %ss.\nAll set.' % (tp2-tp1).seconds