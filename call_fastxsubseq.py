#!/usr/bin/env pypy

import sys, os, commands, re, gc
from Argv import ArgvToDict as ATD
from copy import deepcopy
sys.path.append('/home/wuyang/softwares/my_scripts/')
from fastx_subseq.fastx_subseq import Fastx

def InputToFileList(inp):
	file_list = []
	if type(inp) == list:
		for f in inp:
			if os.path.isfile(f):
				file_list.append(f)
			elif os.path.isdir(f):
				file_list += commands.getoutput('find %s -maxdepth 1 -name "*" -type f' % f).split('\n')
	else:
		if os.path.isfile(inp):
			file_list = [inp, ]
		elif os.path.isdir(inp):
			file_list = commands.getoutput('find %s -maxdepth 1 -name "*" -type f' % inp).split('\n')
			# "-maxdepth 1" limits 'find' command not to access subfolders
	return file_list

def HelpMsg():
	print '''
Commands:
pypy call_fastxsubseq.py --fastqs FASTA/Q --lists query_list -v
'''

if __name__ == '__main__':
	args = ATD(argv_list=sys.argv, required=['--fastqs', '--lists'], optional={'-h':False, '-v':False, '--indices':[]}, verbose=True)
	print args
	if args.get('-h') or args == {}:
		HelpMsg()
		quit()	
	fastqs, lists = InputToFileList(args.get('--fastqs')), InputToFileList(args.get('--lists'))
	indices = args.get('--indices')
	verbose = args.get('-v')
	bc_p = re.compile(r'([ATGCN]{8})')
	lists_dir, indices_tmp = [], []
	if type(indices) == str:
		indices = [indices, ]
	for l in lists:
		bc = re.findall(bc_p, l)
		if len(bc) == 1:
			if indices == []:
				indices_tmp.append(bc[0])
				lists_dir.append(os.path.split(l)[0])
			else:
				if bc[0] in indices:
					lists_dir.append(os.path.split(l)[0])
				else:
					continue
		else:
			continue
	indices = {}.fromkeys(indices_tmp).keys()
	lists_dir = {}.fromkeys(lists_dir).keys()
	#assert len(indices) > 0 and indices[0] != ''
	fastqs_dict = {}
	for f in fastqs:
		bc = re.findall(bc_p, f)
		if len(bc) == 1 and bc[0] in indices:
			files = fastqs_dict.get(bc[0])
			if files == None:
				fastqs_dict.update({bc[0]: [f,]})
			else:
				assert type(files) == list
				files.append(f)
				fastqs_dict.update({bc[0]: deepcopy(files)})
	fastqs = []
	if verbose:
		n = 0
		for f in fastqs_dict.itervalues():
			fastqs += f
		fastqs_num = len(fastqs)
	# Not going to use multiprocess because imported function consumes memory a lot
	for index in fastqs_dict.iterkeys():
		for fastq_file in fastqs_dict.get(index):
			#print fastq_file
			if verbose:
				n += 1
				print ('(%s/%s) Extracting from "%s":' % (str(n), str(fastqs_num), fastq_file))
			fastq_info = Fastx(fastq_file, verbose)
			fastq_info.ExrtactInfo()
			for indir in lists_dir:
				seqname_listfile = '%s/%s' % (indir, index)
				if os.path.isfile(seqname_listfile):
					outdir = indir + '/extracted_sequences/'
					fastq_info.FetchSeq(seqname_listfile, outdir)
			fastq_info.ReleaseMemory()