#!/usr/bin/env pypy

import sys, re, os, gc, commands
from copy import deepcopy
from Argv import ArgvToDict as ATD
from multiprocessing import Pool
from subprocess import check_call as CC
from datetime import datetime as DT

# TODO: trim reads in SAM according to their quality

def ExtractDict(sam_file_chunk):
	try:
		chunk = int(sam_file_chunk.split('.')[-1])
	except ValueError:
		chunk = 1
	with open(sam_file_chunk, 'r') as inf:
		sam_dict_chunk = {}
		for line in inf:
			if line[0] != '@':
				line = line.strip().split('\t', 1)
				BC = re.findall(BC_pattern, line[-1])
				BX = re.findall(BX_pattern, line[-1])
				if len(BC) == 1 and len(BX) == 1:
					BCBX = '%s-%s' % (BC[0], BX[0])
					seqname = line[0]
					sam_dict_chunk_value = sam_dict_chunk.get(BCBX)
					if sam_dict_chunk_value != None:
						count = sam_dict_chunk_value[0] + 1
						seqname_list = deepcopy(sam_dict_chunk_value[-1])
						seqname_list.append(seqname)
						sam_dict_chunk.update({BCBX: [count, deepcopy(seqname_list)]})
					else:
						sam_dict_chunk.update({BCBX: [1, [seqname,]]})
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
			d_merged_value = d_merged.get(item[0])
			if d_merged_value == None:
				d_merged.update({item[0]: item[-1]})
			else:
				count = int(d_merged_value[0]) + int(item[-1][0])
				seqname_list = deepcopy(d_merged_value[-1]) + deepcopy(item[-1][-1])
				assert count == len(deepcopy(seqname_list))
				d_merged.update({item[0]: [count, deepcopy(seqname_list)]})
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
	if args.get('-h') or args == {}:
		HelpMsg()
		quit()
	in_file, cutoff, threads = args.get('-i'), int(args.get('-c')), int(args.get('-p'))
	suffix = in_file.split('.')[-1]
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
				CC('rm -f *.sam.chunk.* && split -d -n l/%s %s %s.chunk.' % \
				   (str(threads), sam_file, sam_file), shell=True)
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
			#CC('rm -f %s.chunk.*' % sam_file, shell=True)
		else:
			print 'Parsing SAM file...',
			sys.stdout.flush()
			tp1 = DT.now()
			sam_dict = ExtractDict(ExtractDict)
			tp2 = DT.now()
			print 'Done.\n    Done in %ss.' % (tp2-tp1).seconds
		CC('mkdir -p ./trimmed_seqname_list_cutoff%s' % str(cutoff), shell=True)
	elif suffix == 'dict':
		tp1 = DT.now()
		print 'Parsing DICT file...',
		sys.stdout.flush()
		sam_dict = File2Dict(in_file)
		tp2 = DT.now()
		print 'Done.\n    Done in %ss.' % (tp2-tp1).seconds
	else:
		print 'Not supported format: %s' % suffix
		quit()
	print 'Outputting to "./trimmed_seqname_list_cutoff%s/"...' % str(cutoff),
	sys.stdout.flush()
	tp1 = DT.now()
	n, m = 0, 0
	for bcbx in sam_dict.iteritems():
		m += 1
		count = int(bcbx[-1][0])
		seqname_list = deepcopy(bcbx[-1][-1])
		assert count == len(seqname_list)
		seqname_list = {}.fromkeys(seqname_list).keys() # delete duplicates
		# 1. yield dict file
		if suffix == 'sam':
			with open('./trimmed_seqname_list_cutoff%s/%s.dict' % (str(cutoff), sam_file), 'a') as dictoutf:
				dictoutf.write('%s\t%s\t%s\n' % (bcbx[0], count, ','.join(seqname_list)))
		# 2. yield seqname_query file
		if count >= cutoff:
			bc = bcbx[0].split('-')[0]
			n += 1
			with open('./trimmed_seqname_list_cutoff%s/%s' % (str(cutoff), bc), 'a') as seqnameoutf:
				seqnameoutf.write('\n'.join(seqname_list) + '\n')
	print 'Done.'
	# release memory
	del sam_dict; gc.collect()
	print '    %.2f' % (100-n*100.0/m) + '% of 10XG barcodes were filtered out.' 
	tp2 = DT.now()
	print '    Done in %ss.\nAll done.' % (tp2-tp1).seconds