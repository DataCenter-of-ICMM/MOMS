#!/usr/bin/env python
'''
python program for hybrid denovo assembly using multiple-channel bionano optical map data and NGS data

Prerequisite Softwares:
1. OMBlast (https://academic.oup.com/bioinformatics/article/33/3/311/2584477)

Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences
MOMS is licensed under the Mulan PSL v1.
You can use this software according to the terms and conditions of the Mulan PSL v1.
You may obtain a copy of Mulan PSL v1 at:
    http://license.coscl.org.cn/MulanPSL
THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND, EITHER EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT, MERCHANTABILITY OR FIT FOR A PARTICULAR
PURPOSE.
See the Mulan PSL v1 for more details.
'''

import os
import sys
import re
import subprocess
import time
from time import  strftime
from configobj import ConfigObj
import argparse

__version__ = "0.1.%s"%(filter(str.isdigit, "$Revision: 54$"))
DESCRIPTION = "MOMS (Multiple-channel Optical Map Scaffolder) -- version %s"%(__version__)
COPYRIGHT = "Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences. All Rights Reserved."

SGESetting = subprocess.check_output('if [ -f "$SGE_ROOT/default/common/settings.sh" ]; then echo -n "$SGE_ROOT/default/common/settings.sh"; else echo -n ""; fi', shell=True)

singleThreadModule = {
	'main'			: 1,
	'inputor'		: 1,
	'cmapconvertor'	: 1, 
	'cmapaligner'	: 0,
	'cmapresolver'	: 1,
	'scaffolder'	: 0,
	'sandwichscaff'	: 1,
	'bngaligner'	: 0,
	'ngsaligner'	: 0,
	'reporter'		: 1
}

class Message:
	'''
	Class for message processing
	'''
#set (https://misc.flogisoft.com/bash/tip_colors_and_formatting)
	reset_all = 0;
	bold = 1;
	dim = 2;
	underlined = 4;
	blink = 5;
	inverted = 7;
# foreground (text) color
	fg_default = 39;
	fg_black = 30;
	fg_red = 31;
	fg_green = 32;
	fg_yellow = 33;
	fg_blue = 34;
	fg_magenta = 35;
	fg_cyan = 36;
	fg_white = 97;
# background color
	bg_default = 49;
	bg_black = 40;
	bg_red = 41;
	bg_green = 42;
	bg_yellow = 43;
	bg_blue = 44;
	bg_magenta = 45;
	bg_cyan = 46;

	@staticmethod
	def title(msg):
		print "\033[%d;%d;%dm%s\033[%dm"%(Message.bold, Message.fg_white, Message.bg_black, msg, Message.reset_all);
		sys.stdout.flush();

	@staticmethod
	def subtitle(msg):
		print "\033[%d;%d;%dm%s\033[%dm"%(Message.bold, Message.fg_blue, Message.bg_black, msg, Message.reset_all);
		sys.stdout.flush();

	@staticmethod
	def info(msg):
		print "\n\033[%d;%d;%dm%s\033[%dm"%(Message.bold, Message.fg_green, Message.bg_black, msg, Message.reset_all);
		sys.stdout.flush();

	@staticmethod	
	def error(msg):
		print "\033[%d;%d;%dm%s\033[%dm"%(Message.bold, Message.fg_red, Message.bg_black, msg, Message.reset_all);
		sys.stdout.flush();

	@staticmethod
	def run_info(cmd):
		print "\033[%d;%d;%dm[%s] Run Command:\033[%dm %s"%(Message.bold, Message.fg_magenta, Message.bg_black, strftime("%H:%M:%S"), Message.reset_all, cmd);
		sys.stdout.flush();

	@staticmethod	
	def time():
		print "[\033[%d;%d;%dm%s\033[%dm]"%(Message.reset_all, Message.fg_yellow, Message.bg_black, strftime("%Y-%m-%d %I:%M:%S %p"), Message.reset_all);
		sys.stdout.flush();

#### Functions
def processArgs():
	parser = argparse.ArgumentParser(description=DESCRIPTION)
	parser.add_argument('-i', dest='fasta', help='--input file resulting from NGS assembly [.fasta]', required=True)
	parser.add_argument('-b', dest='cmaps', nargs='+', help='--BNG cmap file for a enzyme [.cmap]', required=True)
	parser.add_argument('-o', dest='output', help='--output folder', required=True)
	parser.add_argument('-f', dest='force', default='no', help='--force, to override the output folder (default: no)', required=False)
	parser.add_argument('-s', dest='stage', default=None, help='--stage of scaffolding to run/skip', required=False)
	parser.add_argument('-m', dest='mancuts', nargs='*', default=None, help='--manual cut file', required=False)
	parser.add_argument('-c', dest='cluster', default='yes', help='--run in a SGE cluster (default: yes)', required=False)
	parser.add_argument('-t', dest='num_threads', type=int, default=8, help='--number of threads (default: 8)', required=False)
	parser.add_argument('--conf', dest='conf', default=None, help='--designated configuration file', required=False)
	parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))
	args = parser.parse_args()
	inputs = [args.fasta] + args.cmaps;
	bInputOK = True
	for fl in inputs:
		if not os.path.exists(fl):
			Message.error("Input file \"%s\" does not exist"%(fl));
			bInputOK = False
	if not bInputOK:
		sys.exit()

	return args;

# Set global variables
args = processArgs();

def showWelcomeInfo(args):
	# Initialize display parameters
	dicts = {}
	dicts['FASTA file (-i)'] = args.fasta
	for i, bng in zip(range(len(args.cmaps)), args.cmaps):
		key = "BNG cmap file for enzyme%d (-b)"%i
		dicts[key] = bng
	dicts['Output folder (-o)'] = args.output
	dicts['Num of threads (-t)'] = args.num_threads
	if args.mancuts != None:
		for i, mancut in zip(range(len(args.mancuts)), args.mancuts):
			key = "Manual cut file %d (-m)"%i
			dicts[key] = mancut
	if args.force.lower() == 'yes':
		dicts['Override output (-f)'] = 'yes'
	if args.cluster.lower() == 'yes':
		dicts['Run in a SGE cluster (-c)'] = 'yes'
	if args.stage != None:
		dicts['Stage of scaffolding (-s)'] = args.stage
	if args.conf != None:
		dicts['Configuration (--conf)'] = args.conf

	intLength = 0
	for k in dicts.keys():
		if len(k) > intLength:
			intLength = len(k)
	intLength += 1

	# display
	Message.subtitle(DESCRIPTION);
	Message.title(COPYRIGHT);
	print "\n Parameter used:"
	for k,v in dicts.iteritems():
		print "%s: %s"%(k.rjust(intLength), v)
	sys.stdout.flush();

def writeQueueScript(queueName, dictVars, submit=True):
	'''
	Generate scripts for publishing SGE jobs
	'''

	cfg = Queue(queueName)
	workPath = os.path.abspath(dictVars['WORKDIR']) if dictVars.has_key('WORKDIR') else args.output
	tempPath = "%s/scripts/queue"%os.path.dirname(os.path.abspath(sys.argv[0]))
	tempName = "%s.sh"%queueName
	inFile = "%s/%s"%(tempPath, tempName)
	if not os.path.exists(inFile):
		inFile = "%s/default.sh"%tempPath
	outFile = "%s/%s.sh"%(workPath, cfg.sname)

	pubVars = {
		'WORKDIR'	: workPath,
		'SCRNAME'	: cfg.sname,
		'QUEUE'		: cfg.queue,
		'JOBNAME'	: cfg.jname,
		'PENAME'	: cfg.pe,
		'SLOTS'		: 1 if singleThreadModule.has_key(queueName) and singleThreadModule[queueName]==1 else cfg.slots,
		'MEM'		: cfg.mem,
		'PRIORITY'	: cfg.priority
	}

	if queueName=='main':
		pubVars['JOBNAME'] = "%s-%s%s"%(cfg.jname, strftime("%y%m%d.%H%M%S"), workPath.replace('/', '+'))

	with open(inFile, 'r') as f:
		temp = f.read()

	for k,v in pubVars.iteritems():
		temp = temp.replace(k, str(v))
	if dictVars != None:
		for k,v in dictVars.iteritems():
			temp = temp.replace(k, str(v))

	if not os.path.isdir(workPath):
		os.system("mkdir -p %s"%(workPath))

	with open(outFile, 'w') as f:
		f.write(temp)

	if submit:
		os.system("source %s; qsub -e /dev/null %s"%(Config.sgesetting, outFile));

#### Classes
class Util:
	'''
	The class for Utilities
	'''
	@staticmethod
	def run(cmd):
		Message.run_info(cmd);
		retcode = os.system(cmd);
		if retcode != 0:
			Message.error("Error occurs!");
			sys.exit()

	@staticmethod
	def mkdir(path):
		return Util.run("mkdir -p %s"%(path))

	@staticmethod
	def rmdir(path):
		return Util.run("rm -rf %s"%(path))

	@staticmethod
	def change_ext(path, ori, new):
		return os.path.dirname(os.path.abspath(path)) + "/" + re.sub("\.%s$"%ori, ".%s"%new, os.path.basename(path));

	@staticmethod
	def abs_dirname(path):
		return os.path.dirname(os.path.abspath(path))

	@staticmethod
	def basename(path):
		return os.path.basename(path)

class Config:
	'''
	The class for configuration
	'''
	PROGRAM = sys.argv[0]
	PATH = Util.abs_dirname(PROGRAM)
	config_file = Util.change_ext(PROGRAM, "py", "conf") if args.conf == None else args.conf;
	config = ConfigObj(config_file)
	parameters = config['parameters']
	paths = config['paths']
	spath = PATH + "/" + paths['scripts.dir'];
	apath = PATH + "/" + paths['aligner.dir'];
	xmldir = paths['xml.dir'] if 'xml.dir' in paths else "";

	infasta = os.path.abspath(args.fasta)
	incmaps = args.cmaps
	for i in range(len(incmaps)):
		incmaps[i] = os.path.abspath(incmaps[i])

	nthreads = args.num_threads
	bforce = int(args.force.lower() == 'yes')
	bcluster = int(args.cluster.lower() == 'yes')
	bqueue = int(args.cluster.lower() == 'auto')
	outpath = os.path.abspath(args.output)
	fasta_prefix = re.sub("\..*$", "", os.path.basename(infasta))
	assembly_prefix = "assembly"
	bng_suffix = "adjusted_cut"
	ngs_suffix = "cut"
	scaffold_suffix = "hybrid"

	config_file = "%s/queue.conf"%PATH
	config = ConfigObj(config_file)
	paths = config['paths']
	sgesetting = subprocess.check_output('if [ -f "%s" ]; then echo -n "%s"; else echo -n ""; fi'%(paths['sge.setting'], paths['sge.setting']), shell=True)

class Queue:
	'''
	The class for queue configuration
	'''
	def __init__(self, typeName):
		config_file = "%s/queue.conf"%Util.abs_dirname(sys.argv[0])
		config = ConfigObj(config_file)
		cfg = config[typeName]
		self.sname		= cfg['script_name'] if cfg.has_key('script_name') else 'start%s'%typeName.upper()
		self.queue		= cfg['queue_name'] if cfg.has_key('queue_name') else ''
		self.jname		= cfg['job_name'] if cfg.has_key('job_name') else typeName
		self.pe			= cfg['pe'] if cfg.has_key('pe') else ''
		self.slots		= Config.nthreads if Config.nthreads>0 else cfg['slots']
		self.mem		= cfg['mem'] if cfg.has_key('mem') else '1g'
		self.priority	= cfg['priority'] if cfg.has_key('priority') and cfg.has_key('priority')>-1024 and cfg.has_key('priority')<=1024 else 0

class Program:
	'''
	The base class for a program
	'''
	steps = []
	cnt = 0
	def __init__(self, name, desc, description, depends):
		Program.steps.append(self)
		Program.cnt += 1
		self.no = Program.cnt
		self.program = Config.spath + "/" + name
		self.desc = desc
		self.description = description
		self.depends = depends
		assert(Config.outpath != "")
		self.outdir = "%s/Step-%02d_%s"%(Config.outpath, self.no, desc.replace(" ", "_"))
		self.status = "%s/status.txt"%self.outdir
		self.config()

	def config(self):
		print "invoked config() from Class %s" % self.__class__.__name__

	def preprocess(self):
		pass

	def waiting(self):
		if Config.bqueue:
			while True:
				status = int(subprocess.check_output("grep -iE 'done|error' %s 2>/dev/null | wc -l"%(self.status), shell=True))
				if status > 0:
					break
				time.sleep(1)

	def isDone(self):
		print "invoked isDone() from Class %s" % self.__class__.__name__

	def process(self):
		print "invoked process() from Class %s" % self.__class__.__name__

	def run(self):
		Message.info("Step %d: %s"%(self.no, self.description))
		self.preprocess()
		if not self.isDone() or Config.bforce:
			if Config.bforce and os.path.exists(self.outdir):
				Util.rmdir(self.outdir)
			Util.mkdir(self.outdir)
			self.process()
			self.waiting()
			if not self.isDone():
				Message.error("Not finished");
				sys.exit()
			Message.time()

class Inputor(Program):
	'''
	The base class for input files preparation
	'''
	def config(self):
		self.fasta = Config.infasta
		self.cmaps = ','.join(Config.incmaps)
	
	def isDone(self):
		num1 = int(subprocess.check_output("ls -1 %s/contig-length.txt %s 2>/dev/null | wc -l"%(self.outdir, self.fasta), shell=True));
		num2 = int(subprocess.check_output("ls -1 %s/%s_*.cmap 2>/dev/null | wc -l"%(self.outdir, Config.assembly_prefix), shell=True));
		num3 = int(subprocess.check_output("grep ':' %s/enzymes.txt 2>/dev/null | wc -l"%(self.outdir), shell=True));
		return (num1 == 2) and (num2 == len(Config.incmaps) and (num2 == num3))
	
	def process(self):
		cmd=("%s %s %s %s")%(self.program, self.fasta, self.cmaps, self.outdir);
		if Config.bqueue:
			writeQueueScript('inputor', {'WORKDIR': self.outdir, 'COMMAND': "bash %s"%cmd})
		else:
			Util.run(cmd)

class CmapConvertor(Program):
	'''
	The base class for FASTA to CMAP convertor
	'''
	def config(self):
		self.indir = self.depends[0].outdir
		self.fasta = Util.basename(self.depends[0].fasta)
		self.minnicks = int(Config.parameters['fa2cmap.minLabels']);
		self.minklen = int(Config.parameters['fa2cmap.minLength']);

	def preprocess(self):
		enzyme_file = ("%s/enzymes.txt")%(self.indir);
		enzymes = [line.rstrip('\n') for line in open(enzyme_file)];
		self.enzymes = ','.join(enzymes);

	def isDone(self):
		num = int(subprocess.check_output("grep -m 1 '^100%% done$' %s/status.txt 2>/dev/null | wc -l"%(self.outdir), shell=True));
		return (num == 1)

	def process(self):
		cmd=("%s %s/%s %s %s")%(self.program, self.indir, self.fasta, self.enzymes, self.outdir)
		if (self.minnicks > 0) or (self.minklen > 0):
			cmd += " %d"%self.minnicks
			if self.minklen > 0:
				cmd += " %d"%self.minklen
		
		if Config.bqueue:
			writeQueueScript('cmapconvertor', {'WORKDIR': self.outdir, 'COMMAND': "bash %s"%cmd})
		else:
			Util.run(cmd)

class CmapRescaler(Program):
	'''
	The base class for rescaling BNG cmap
	'''
	def config(self):
		self.refdir = self.depends[0].outdir; # cmapconvertor
		self.refpre = Config.fasta_prefix;
		self.qrydir = self.depends[1].outdir; # inputor
		self.qrypre = Config.assembly_prefix

	def preprocess(self):
		self.nfiles = int(subprocess.check_output("ls -1 %s/%s_*.cmap 2>/dev/null | wc -l"%(self.qrydir, self.qrypre), shell=True));

	def isDone(self):
		num = int(subprocess.check_output("ls -1 %s/*.cmap 2>/dev/null | wc -l"%(self.outdir), shell=True));
		return (num == self.nfiles)

	def process(self):
		cmd = ("%s %s/%s %s/%s %s")%(self.program, self.refdir, self.refpre, self.qrydir, self.qrypre, self.outdir);
		if Config.nthreads != 0:
			cmd += " %d"%Config.nthreads
			if Config.xmldir:
				cmd += " %s"%Config.xmldir
		else:
			if Config.xmldir:
				cmd += " 0 %s"%Config.xmldir

		if Config.bqueue:
			writeQueueScript('cmapaligner', {'WORKDIR': self.outdir, 'COMMAND': "bash %s"%cmd})
		else:
			Util.run(cmd)

class ChimeraResolver(Program):
	'''
	The base class for Chimaeral cmap resolver
	'''
	def config(self):
		self.refdir = self.depends[0].outdir; # cmapconvertor
		self.refpre = Config.fasta_prefix
		self.qrydir = self.depends[1].outdir; # cmaprescaler
		self.qrypre = Config.assembly_prefix

	def preprocess(self):
		self.nfiles = int(subprocess.check_output("ls -1 %s/*.cmap 2>/dev/null | wc -l"%(self.qrydir), shell=True));

	def isDone(self):
		num1 = int(subprocess.check_output("ls -1 %s/%s_*_adjusted_cut.cmap 2>/dev/null | wc -l"%(self.outdir, self.qrypre), shell=True));
		num2 = int(subprocess.check_output("ls -1 %s/%s_*_cut.cmap 2>/dev/null | wc -l"%(self.outdir, self.refpre), shell=True));
		return (num1 == self.nfiles) and (num2 == self.nfiles)

	def process(self):
		cmd = ("%s %s/%s %s/%s %s")%(self.program, self.refdir, self.refpre, self.qrydir, self.qrypre, self.outdir);
		if Config.nthreads != 0:
			cmd += " %d"%Config.nthreads
			if Config.xmldir:
				cmd += " %s"%Config.xmldir
		else:
			if Config.xmldir:
				cmd += " 0 %s"%Config.xmldir

		if Config.bqueue:
			writeQueueScript('cmapresolver', {'WORKDIR': self.outdir, 'COMMAND': "bash %s"%cmd})
		else:
			Util.run(cmd)

class HybridScaffolder(Program):
	'''
	The base class for Hybrid Scaffolder
	'''
	def config(self):
		self.indir = self.depends[0].indir; # cmapconvertor
		self.fasta = self.depends[0].fasta; # cmapconvertor
		self.encdir = self.depends[0].outdir; # cmapconvertor
		self.encpre = Config.fasta_prefix
		self.resdir = self.depends[1].outdir; # cmapresolver
		self.ngspre = Config.fasta_prefix
		self.bngpre = Config.assembly_prefix
		self.bfillgap = int(Config.parameters['scaffold.single.fillgap']);

	def isDone(self):
		num1 = int(subprocess.check_output("ls -1 %s/*.%s.cmap 2>/dev/null | wc -l"%(self.outdir, Config.scaffold_suffix), shell=True));
		num2 = int(subprocess.check_output("ls -1 %s/*.agp 2>/dev/null | wc -l"%(self.outdir), shell=True));
		num3 = int(subprocess.check_output("ls -1 %s/*.fasta 2>/dev/null | wc -l"%(self.outdir), shell=True));
		return (num1 == len(Config.incmaps)) and (num2 == len(Config.incmaps)) and (num3 == len(Config.incmaps))

	def process(self):
		cmd = ("%s %s/%s %s/%s %s/%s %s/%s %s %d")%(self.program, self.resdir, self.ngspre, self.resdir, self.bngpre, self.indir, self.fasta, self.encdir, self.encpre, self.outdir, self.bfillgap);
		if Config.nthreads != 0:
			cmd += " %d"%Config.nthreads
			if Config.xmldir:
				cmd += " %s"%Config.xmldir
		else:
			if Config.xmldir:
				cmd += " 0 %s"%Config.xmldir

		if Config.bqueue:
			writeQueueScript('scaffolder', {'WORKDIR': self.outdir, 'COMMAND': "bash %s"%cmd})
		else:
			Util.run(cmd)

class SandwichScaffolder(Program):
	'''
	The base class for Sandwich Scaffolder
	'''
	def config(self):
		self.scfdir = self.depends[0].outdir; # scaffolder
		self.scfsuf = Config.scaffold_suffix
		self.resdir = self.depends[1].outdir; # cmapresolver
		self.ngspre = Config.fasta_prefix

	def isDone(self):
		num = int(subprocess.check_output("ls -1 %s/multicolors.cmap 2>/dev/null | wc -l"%(self.outdir), shell=True));
		return (num == 1)

	def process(self):
		cmd = ("%s %s %s %s/%s %s")%(self.program, self.scfdir, self.scfsuf, self.resdir, self.ngspre, self.outdir);

		if Config.bqueue:
			writeQueueScript('sandwichscaff', {'WORKDIR': self.outdir, 'COMMAND': "bash %s"%cmd})
		else:
			Util.run(cmd)

class BNGAligner(Program):
	'''
	The base class for BNG data to scaffolds aligner
	'''
	def config(self):
		self.refdir = ("%s")%(self.depends[0].outdir); # sandwichScaff
		self.qrydir = self.depends[1].outdir; # cmapresolver
		self.qrypre = Config.assembly_prefix
		self.qrysuf = Config.bng_suffix
		self.bfillgap = int(Config.parameters['scaffold.mono.fillgap']);

	def preprocess(self):
		self.nfiles = int(subprocess.check_output("ls -1 %s/mono/*.cmap 2>/dev/null | wc -l"%(self.refdir), shell=True));

	def isDone(self):
		num = int(subprocess.check_output("ls -1 %s/single/*.xmap 2>/dev/null | wc -l"%(self.outdir), shell=True));
		if (num != self.nfiles) or (num != len(Config.incmaps)):
			return False;
		num = int(subprocess.check_output("ls -1 %s/{multicolors,unified}.cmap 2>/dev/null | wc -l"%(self.outdir), shell=True));
		return (num == 2);

	def process(self):
		cmd = ("%s %s %s/%s %s %s %d")%(self.program, self.refdir, self.qrydir, self.qrypre, self.qrysuf, self.outdir, self.bfillgap);
		if Config.nthreads != 0:
			cmd += " %d"%Config.nthreads
			if Config.xmldir:
				cmd += " %s"%Config.xmldir
		else:
			if Config.xmldir:
				cmd += " 0 %s"%Config.xmldir

		if Config.bqueue:
			writeQueueScript('bngaligner', {'WORKDIR': self.outdir, 'COMMAND': "bash %s"%cmd})
		else:
			Util.run(cmd)

class NGSAlignerUni(Program):
	'''
	The base class for NGS sequences to scaffolds aligner
	'''
	def config(self):
		self.refdir = ("%s")%(self.depends[0].outdir); # bngaligner
		self.qrydir = self.depends[1].outdir; # chimeraresolver
		self.qrypre = Config.fasta_prefix
		self.qrysuf = Config.ngs_suffix
		self.encdir = self.depends[2].outdir; # cmapconvertor
		self.indir = self.depends[3].outdir; # inputor
		self.fasta = self.depends[3].fasta # inputor

	def preprocess(self):
		enzyme_file = ("%s/enzymes.txt")%(self.indir);
		enzymes = [line.rstrip('\n') for line in open(enzyme_file)];
		self.enzymes = ','.join(enzymes);

	def isDone(self):
		num = int(subprocess.check_output("ls -1 %s/{combined.xmap,combined_NGS_coord_translation.txt} 2>/dev/null | wc -l"%(self.outdir), shell=True));
		return (num == 2);

	def process(self):
		cmd = ("%s %s %s/%s %s %s %s %s %s")%(self.program, self.refdir, self.qrydir, self.qrypre, self.qrysuf, self.fasta, self.enzymes, self.encdir, self.outdir);
		if Config.nthreads != 0:
			cmd += " %d"%Config.nthreads
			if Config.xmldir:
				cmd += " %s"%Config.xmldir
		else:
			if Config.xmldir:
				cmd += " 0 %s"%Config.xmldir

		if Config.bqueue:
			writeQueueScript('ngsaligner', {'WORKDIR': self.outdir, 'COMMAND': "bash %s"%cmd})
		else:
			Util.run(cmd)

class NGSAligner(Program):
	'''
	The base class for NGS sequences to scaffolds aligner
	'''
	def config(self):
		self.refdir = ("%s")%(self.depends[0].outdir); # bngaligner
		self.qrydir = self.depends[1].outdir; # cmapresolver
		self.qrypre = Config.fasta_prefix
		self.qrysuf = Config.ngs_suffix

	def isDone(self):
		num = int(subprocess.check_output("ls -1 %s/{combined.xmap,combined_NGS_coord_translation.txt} 2>/dev/null | wc -l"%(self.outdir), shell=True));
		return (num == 2);

	def process(self):
		cmd = ("%s %s %s/%s %s %s")%(self.program, self.refdir, self.qrydir, self.qrypre, self.qrysuf, self.outdir);
		if Config.nthreads != 0:
			cmd += " %d"%Config.nthreads
			if Config.xmldir:
				cmd += " %s"%Config.xmldir
		else:
			if Config.xmldir:
				cmd += " 0 %s"%Config.xmldir

		if Config.bqueue:
			writeQueueScript('ngsaligner', {'WORKDIR': self.outdir, 'COMMAND': "bash %s"%cmd})
		else:
			Util.run(cmd)

class SeqReporter(Program):
	'''
	The base class for final scaffolds reporter
	'''
	def config(self):
		self.alndir = self.depends[0].outdir; # ngsaligner
		self.alnfile = "combined.xmap";
		self.cmapdir = self.depends[1].outdir; # sandwichScaff
		self.cmapfile = "multicolors.cmap";
		self.fadir = self.depends[2].outdir; # inputor
		self.fafile = Util.basename(self.depends[2].fasta); # inputor
		self.encdir = self.depends[3].outdir; # cmapconvertor
		self.encpre = Config.fasta_prefix;
		self.outpre = Config.fasta_prefix

	def isDone(self):
		num = int(subprocess.check_output("ls -1 %s/%s{.agp,.gap,_NCBI.fasta,_NOT_SCAFFOLDED.fasta}  2>/dev/null | wc -l"%(self.outdir, self.outpre), shell=True));
		return (num == 4);

	def process(self):
		cmd = ("%s %s/%s %s/%s %s/%s %s/%s %s/%s")%(self.program, self.alndir, self.alnfile, self.fadir, self.fafile,
				 self.cmapdir, self.cmapfile, self.encdir, self.encpre, self.outdir, self.outpre);

		if Config.bqueue:
			writeQueueScript('reporter', {'WORKDIR': self.outdir, 'COMMAND': "bash %s"%cmd})
		else:
			Util.run(cmd)

# MAIN ENTRY POINT
def main():
	'''
	The main function.
	'''
	if(not os.path.exists(Config.apath + "/RefAligner")):
		cmd=("%s/setup.sh")%(Config.spath);
		Util.run(cmd);
		if(not os.path.exists(Config.apath + "/RefAligner")):
			Message.error("Not finished");
			sys.exit();

# start running the whole program
	showWelcomeInfo(args);
	if(not os.path.exists(Config.outpath)):
		Util.mkdir(Config.outpath)

# prepare modules for running
	inputor = Inputor("prepare-input.sh", "Input", "Prepare input files", [])
	cmapconvertor = CmapConvertor("fa2cmap.sh", "NGS CMAP encoding", "Encode NGS contigs to cmap file(s)", [inputor])
	cmaprescaler = CmapRescaler("rescale-cmaps.sh", "BNG CMAP rescaling", "Rescale BNG cmap file(s)", [cmapconvertor, inputor])
	chimeraresolver = ChimeraResolver("resolve-chimeras.sh", "Chimeras resolution", "Detect and resolve chimeral cmaps", [cmapconvertor, cmaprescaler])
	scaffolder = HybridScaffolder("hybrid-scaffold.sh", "Pre-scaffold", "Perform single-enzyme hybrid scaffolding using chimera-resolved cmaps", [cmapconvertor, chimeraresolver])
	if len(args.cmaps)>1:
		sandwichScaff = SandwichScaffolder("sandwich-scaffold.sh", "Final Scaffold", "Perform multi-enzyme scaffolding mediated by NGS contigs", [scaffolder, chimeraresolver])
		bngaligner = BNGAligner("align-final-bng.sh", "BNG anchoring", "Align BNG data to scaffolds", [sandwichScaff, chimeraresolver])
		if int(Config.parameters['anchor.unichannel'])>0:
			ngsaligner = NGSAlignerUni("align-final-ngs-uni.sh", "NGS anchoring", "Align NGS contigs to scaffolds", [bngaligner, chimeraresolver, cmapconvertor, inputor])
		else:
			ngsaligner = NGSAligner("align-final-ngs.sh", "NGS anchoring", "Align NGS contigs to scaffolds", [bngaligner, chimeraresolver])
		reporter = SeqReporter("report.sh", "Report", "Report the final scaffolds in AGP/FASTA format", [ngsaligner, sandwichScaff, inputor, cmapconvertor])

# start running individual steps
	for prog in Program.steps:
		prog.run()

# for direct script invoking
if __name__ == "__main__":
	if Config.bcluster:
		if Config.sgesetting == '':
			Message.error('Error: SGE initialization script "settings.sh" doesn\'t exist, please check the file "queue.conf".')
			sys.exit()
		cmd = "python %s -i %s -b %s -o %s -t %d -c auto"%(os.path.abspath(sys.argv[0]), Config.infasta, ' '.join(Config.incmaps), Config.outpath, Config.nthreads)
		writeQueueScript('main', {'WORKDIR': Config.outpath, 'COMMAND': cmd})
		sys.exit();
	main()
