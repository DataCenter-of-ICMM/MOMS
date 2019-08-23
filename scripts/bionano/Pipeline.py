import os, sys, shutil
#import traceback
import xml.etree.ElementTree as ET
import argparse, copy, time

import PairwiseModule as pm
import AssemblyModule as am
import RefinementModule as refmod
import GroupedRefinementModule as grprefmod
import CharacterizeModule as cm
import SampleCharModule as scm
import SVModule as svm
import Multithreading as mthread
import utilities as util
import AlignModule as alignmod
import subprocess

"""@package Pipeline 
Main routine for defining and running De Novo Assembly
"""


util.setVersion("$Id: Pipeline.py 6590 2017-06-29 17:32:11Z wandrews $")


class varsPipeline():
    """Pipeline Variables Class.
    
        This class parses the command line and input files.  It stores the 
        runtime parameters and determines the process flow for the other
        modules.  Along with the process flow, the global file system is
        created, and binary dependencies are checked.
    """

    def __init__(self):
        self.infoReportVersion = "2.3.0" #see self.prerunLog
        #see self.getScriptPaths for use of these scripts dirs members
        self.scriptsTopDir = "Analysis" #this is a subdir of Pipeline which contains all external scripts
        self.scriptsSubdirs = ["SV", "RUtil"] #subdirs of scriptsTopDir
        self.scriptsSVdirs = ["InDelConfidence", "MergeInversions", "CopyNumberProfiles"] #sub-dirs of scriptsTopDir/SV
        self.scriptsIndelConf_path = "InDelConfidence.r"
        self.scriptsMergeInv_linux = "Merge_Inversions.pl"
        self.scriptsMergeInv_win   = "MergeInversions.exe"
        self.scriptsMergeInv_outsuf = "_inversions.smap"
        self.scriptsMergeInv_ver = 3567 #minimum refaligner version required by this script
        self.scriptsCopyNum_top = "000_COPY_NUMBER_WRAPPER.r"
        self.scriptsCopyNum_data = "DATA" #for windows, I can't put paths together manually
        self.scriptsCopyNum_hg19dir = "hg19" #subdir of data
        self.scriptsCopyNum_hg38dir = "hg38" #subdir of data
        self.scriptsCopyNum_hg19 = "hg19_chromosome_bspq.cmap" #in hg19dir
        self.scriptsCopyNum_hg19bed = "hg19_gaps.bed" #in hg19dir
        self.scriptsCopyNum_hg38 = "hg38_bspq.cmap" #in hg38dir
        self.scriptsCopyNum_hg38bed = "hg38_gaps_5kb.bed" #in hg38dir
        self.scriptsCompress    = "compressPipeline.sh" #in main Pipeline dir
        self.scriptsCompressName = "pipeline_results"
        self.localRoot = None # Top folder for analysis results
        self.toolsFolder = None # Target of executable files
        self.nThreads = 0       # Meaningful for local execution, how many concurrent jobs to run.
        self.nPairwiseJobs = 0
        self.nExtensionIter = 0 # Number of extension-merge iterations after assembly
        self.extensionCount = 0 # Running count of extension merge iterations
        self.expID = 'expDev'   # Prefix for all resulting output files
        self.ref = None         # Reference genome as cmap file (None if no reference)
        self.refDeresed = None  # Reference cmap after referenceProcess (CharacterizeModule)
        self.bypassCount = 0    # Number of stages to bypass (command line: -B)
        self.currentStageCount = 0   # Current Stage number, to be compared with bypass count
        self.executeCurrentStage = True  # Result of bypass count / currentStageCount comparison
        self.curCharacterizeTargets = [] 
        ## @var curCharacterizeTargets
        # Latest merged cmaps for characterization
        self.curCharacterizeCmaps = [] # Latest Cmaps fo characterization
        self.numCharacterizeJobs = 1 # Load balance the characterize work into numCharacterizeJobs jobs
        self.totAssemblyLenMb = 0 # characterize
        self.imgFile = None
        self.bnxFileIn = None
        self.bnxFile = None
        self.bnxStatsFile = None
        self.bnxNmol     = 0
        self.bnxTotLenbp = 0
        self.subsampleMb = 0
        self.doNoise = False
        self.noiseOnly = False 
        self.time=False
        self.perf=False
        self.startTime = 0 #time.time() #pipeline start time
        self.cleanup = 0 
        self.genomeSize = 100 
        self.onCluster = False
        self.clusterArgumentsFileIn = None
        self.clusterLogDir = None
        self.wipe = False
        self.idFile = ''
        self.sorted_file = ''
        self.error = 0
        self.warning = 0
        self.message = ''
        self.optArgumentsFile = None
        #list blocks which are _always_ required in optArguments.xml (the -a parameter)
        #this will be added to depending on command line parameters
        #also all of these are required in clusterArguments also
        self.requiredArgGroups = [#'imgDetection', #required only for -I
                                  #'sampleChar', #required only for -n
                                  #'bnx_sort', #required because triangle mode is always on -- moved to requiredOptArgGroups
                                  #'noise0', #moved to requiredOptArgGroups
                                  'pairwise',
                                  'assembly',
                                  'refineA',
                                  'refineB',
                                  #'refineNGS', #required only for -f
                                  #'refineFinal', #NOT required for -u, otherwise required
                                  'extension',
                                  'merge',
                                  'characterizeDefault',
                                  #'mapLambda', #required only for -L
                                  #'lambdaFilter', #required only for -L
                                  ]
        #these two were in requiredArgGroups, but they're not required for clusterArgs
        self.requiredOptArgGroups = ['bnx_sort',
                                     'noise0']
        self.requiredClusterArgGroups = ['splitting', #for splitBNX (PairwiseModule)
                                         'cmapMerge']
        #these groups are copied from clusterArgs to optArgs if present in former (but should be present in latter for runs on stand-alone node): values will be replaced by the argument list
        self.memoryArgGroups = {'largeNodeMem':None, 'tinyNodeMem':None, 'splitNodeMem':None}

        #self.shortArgs = ['-T','-j','-N','-M','-i','-I','-b','-l','-t','-B','-e','-r','-H','-P','-s','-n','-x','-c','-g','-C','-w','-a','-L','-p','-d','-f','-F','-u','-U','-D','-S','-V','-W','-A','-G'] #this is used in celery implementation only (keep up to date)--no -v for celery implementation -- obsolete
        self.argData = {}
        self.clusterArgData = {}
        self.stageComplete = None
        #self.curContigPrefix = None
        self.inputContigPrefix = ''
        self.outputContigPrefix = ''
        self.inputContigFolder = None
        self.outputContigFolder = None
        self.curContigCount = 0
        self.latestMergedCmap = None
        self.prefixUsed = []
        #self.logstdouterr = False #not used
        #self.logdir = '' #not used
        self.lambdaRef = None
        self.ngsInDir = None
        self.skipFinalRefine = None
        self.doAlignMol = True
        self.characterizeDirName = 'alignref'
        self.refineFinal1 = "refineFinal1" #this is used in self.copyToRefineFinal and DNPipeline.
        self.alignMolDir = "" #set in AlignModule and used in SVModule
        self.alignMolvrefName = 'alignmolvref'
        self.alignMolvrefMergeName = "merge" #see AlignModule.mergeRcmaps
        self.bedFile = "" #optional, only if -G used
        self.refaligner_version = 0 #see prerunLog
        self.autoNoise = False #-y command line arg
        self.scanScaleArg = "-ScanScaling" #the RefAligner option
        self.scanScaleSection = "autoNoise1" #the section (Pipeline module) in which this option is used
        self.doScanScale = True #False #see def checkScanScaling: default is to use scan scaling
        self.rescaleSuffix = "_rescaled" #the file will be <path+prefix><this>.bnx
        self.groupSV = 0 #now a command line arg, -s
        self.pipeReportFile = None
        self.memoryLogpath = ""
        self.stageFolder={}
        self.largemem = "largeNodeMem" #use this optargs section to move jobs from mic to host (SVModule and GroupedRefinementModule)
        self.groupRefineHostFr = 0
        self.NumHostJobs = 0 # If 0, run fixed fraction self.groupRefineHostFr of number of jobs on Host, otherwise run this many jobs on Host (Set during Stage 0 in Refine.groupContigs)
        self.groupRefineHostCA = "mediumHostJob" #if not present, fall back on 'autoNoise0'
        self.stdoutRecheck = True # redo CheckStdout check after all jobs have completed (can be slow, but needed to get a complete report of errors)
        self.GroupSize = 1 # Multiply group sizes by this factor to reduce number of jobs (Only needed for Genomes larger than Human)
        self.usecolorArg = 0

    def toExecuteCurrentStage(self, quiet=False):
        """Check if the current stage is to be excecuted, determined by the 
        bypass parameter (command line: -B).
        """
        self.currentStageCount += 1
        if self.currentStageCount > self.bypassCount:
            self.executeCurrentStage = True
            if not quiet : outstr = "Executing"
        else:
            self.executeCurrentStage = False
            if not quiet : outstr = "Skipping"
        if not quiet :
            outstr += " stage number %i" % self.currentStageCount #, ":", self.stageComplete #problem with printing stageComplete is that it's the previous stage and the next stage is not known
            #print outstr
            self.updatePipeReport(outstr+"\n\n")
        return self.executeCurrentStage
    

    def prerunChecks(self):
        """Check inputs and run time parameters before starting assembly
        
        Parse the command line, check the integrity of binary and input files, 
        check existance of required files, check contridictions in flow based 
        on inputs. Finally assign pipeline variables and construct output file
        structure.
        """
        self.parseCommandLine()
        self.checkFlow() #this is more or less just command line option consistency checking--do before input checking
        if self.error: #if bad command line args, do not do anything else
            return
        self.assignFileTargets()
        if self.error: #if localRoot is not set
            return
        self.checkInputs()
        self.checkToolPaths()
        if self.error: 
            return
        self.setupFileStructure() #this is where bnx is copied and dirs are created
        self.checkCluster()

        self.RefAlignerBinOrig = self.RefAlignerBin
        self.AssemblerBinOrig = self.AssemblerBin
        self.DMstaticBinOrig = self.DMstaticBin
        if self.onCluster:
                self.RefAlignerBin += "${BINARY_SUFFIX:=}"
                self.AssemblerBin += "${BINARY_SUFFIX:=}"
                self.DMstaticBin += "${BINARY_SUFFIX:=}"
        if self.onCluster:
            self.parseArguments(readingClusterFile=True)
        self.parseArguments()
        if not self.error: #if bad xml this is redundant
            self.checkRequiredArguments() #move this to here from parseArguments bc call once for both optArgs and clusterArgs
        self.checkDependencies()
        if self.error:
            return
        #self.setCoverageBasedParams() #move this to DNPipeline.constructData so that it's after bnx_sort

    
    def parseCommandLine(self):
        parser = argparse.ArgumentParser(description='Pipeline for de novo assembly - BioNano Genomics')
        parser.add_argument('-T', help='Available threads per Node [default 1]', default=1, type=int)
        parser.add_argument('-j', help='Max Threads per job [default -T]', dest='maxthreads', default=0, type=int)
        parser.add_argument('-jp', help='Max Threads per pairwise or stage0 job [default -T arg]', dest = 'maxthreadsPW', default = 0,type=int)
        parser.add_argument('-N', help='Number of split bnx files; number of pairwise jobs is N*(N-1)/2 (optional, default 2)', default=2,type=int)
        parser.add_argument('-G', help="Bed file for gaps, used in structural variation (SV) detection to check for SV overlap with reference gaps", dest='bed', default="")
        parser.add_argument('-i', help='Number of extension and merge iterations (default=1, must be in range [0,20], use 0 to skip)', dest='iter', default=1, type=int)
        #parser.add_argument('-I', help='File with listed paths for image processing; no longer supported--do not use (use without .bnx)', dest='img', default=None) #deprecate
        parser.add_argument('-b', help='Input molecule (.bnx) file, required', dest = 'bnx', type=str)
        parser.add_argument('-l', help='Location of output files root directory, required, will be created if does not exist; if does exist, will overwrite contents (may be error-prone)', dest='local', type=str)
        parser.add_argument('-t', help='Location of executable files (RefAligner and Assembler, required)', dest='tools', type=str)
        bypassHelp = 'Skip steps, using previous result. <= 0:None, 1:ImgDetect, 2:NoiseChar/Subsample, 3:Pairwise, 4:Assembly, 5:RefineA, 6:RefineB, 7:merge0, 8+(i-1)*2:Ext(i), 9+(i-1)*2:Mrg(i), N+1:alignmol'
        parser.add_argument('-B', help=bypassHelp,dest='bypass', default = '0', type=int)
        parser.add_argument('-e', help='Output file prefix (optional, default = exp)', dest='exp', default = 'exp')
        parser.add_argument('-r', help='Reference file (must be .cmap), to compare resulting contigs (optional)',dest='ref', default=None)
        #parser.add_argument('-s', help='Make a subsampling of the input set in Mb. (optional)',dest='sub', default = 0) #retire this option
        #parser.add_argument('-s', help='SV jobs configuration: 0 = single job (required for correct haplotype calls), 1 = single job per contig (not recommended), 2 = grouped (default 0; optional)', dest='groupsv', type=int, default=0) #retire
        #parser.add_argument('-n', help='Evaluate single molecule noise characterization', dest='noise', action='store_true') #retire
        parser.add_argument('-x', help='Exit after auto noise (noise characterization), do not preform de novo assembly', action='store_true')
        parser.add_argument('-c', help='Remove contig results (0 - keep all (default), 1 - remove intermediate files, 2 - store in sqlite, 3 - store in sqlite and remove)', dest='cleanup', default=0, type=int)
        #parser.add_argument('-g', help='Organism genome size estimate in megabases, used for tuning assembly parameters [optional, if > 0, will modify parameters, if == 0, ignored, must be float]', dest='size', default=0, type=float) #retire this option
        parser.add_argument('-C', help='Run on cluster, read XML file for submission arguments (optional--will not use cluster submission if absent)', dest='cxml', default=None)
        parser.add_argument('-w', help='Wipe clean previous contig results', dest='wipe', action='store_true')
        parser.add_argument('-a', help='Read XML file for parameters (required)', dest='xml', default=None)
        #parser.add_argument('-L', help='Lambda phage is spiked in, used for molecule scaling (only used with -I input)', dest='lambdaRef', default=None) #deprecate
        parser.add_argument('-p', help='Log performance in pipelineReport 0=None, 1=time, 2=perf, 3=time&perf (default=1)', dest='perf', default=1,type=int)
        parser.add_argument('-d', help='Retired option (contig subdirectories), always enabled.',dest='dirs', action='store_true')
        #parser.add_argument('-f', help='Directory which contains NGS contigs as cmaps for use in refineNGS (no longer supported).', dest='ngsarg', default=None) #retire
        #parser.add_argument('-F', help='When using -f, disable min contigs check for all stages prior to refineNGS (no longer supported).', dest='ngsskiparg', action='store_true') #retire
        parser.add_argument('-u', help='Do not perform final refinement (not recommended).', dest='nofinalref', action='store_true')
        #parser.add_argument('-U', help='Group contigs in refinement and extension stages', dest='groupcontigs', action='store_true', default=0)
        #specify no required args with 'nargs="?"': if this is omitted, default is used; if no arg, const is used. In this case, I want both to be 1.
        parser.add_argument('-U', help='Group contigs in refinement and extension stages [default ON, use 0 to disable]', dest='groupcontigs', nargs="?", type=int, default=1, const=1)
        #parser.add_argument('-D', help='Distribute contig characterization into parts. (optional, default=1)',dest='nCharJobs', default = 1, type=int) #retire this option
        parser.add_argument('-v', help='Print version; exit if argument > 1 supplied.', nargs="?", dest='version', const=1, default=0, type=int)
        #parser.add_argument('-S', help='Log stdout/stderr of all RefAligner calls in files with same names as output files and suffix .stdout (default on--use this flag to turn off).', dest='stdoutlog', default=1, action='store_false')
        parser.add_argument('-V', help='Detect structural variations. Default: only after final stage (normally refineFinal); if argument 2, also after refineB; if argument 0, disable.', dest='runsv', default=1, type=int)
        #parser.add_argument('-W', help='Reserved for celery implementation (do not use).',dest='workorder')
        parser.add_argument('-A', help='Align molcules to final contigs (ON by default, use this to turn off).', default=True, action='store_false')
        parser.add_argument('-y', help='Automatically determine noise parameters (requires reference; optional, default off)', default=False, action='store_true')
        #parser.add_argument('-Y', help='Disable scan scaling in auto noise (default on with reference and -y)', default=True, action='store_false') #I don't see why we need this...
        parser.add_argument('-m', help='Disable molecule vs reference alignments (default on with reference)', default=True, action='store_false')
        parser.add_argument('-H', help='Use HG19 (human genome) as reference, loaded from Analysis/SV/CopyNumberProfiles/DATA. Overrides -r argument. Use HG38 if argument 2 is supplied. [Default OFF]', nargs="?", type=int, default=0, const=1)
        parser.add_argument('-f', help='Run this fraction of grouped jobs on host (0.2 if no arg) [default 0]', nargs="?", const=0.2, default=0, type=float)
        parser.add_argument('-J', help='Number of threads on host for grouped jobs (has no effect without -f) [default 48]', default=48, type=int)
        parser.add_argument('-z', help='Zip pipeline results (default ON, use this to turn off).', default=True, action='store_false')
        parser.add_argument('-E', help='ReCheck stdout completeness for completed jobs (default ON, use this to turn off).', default=True, action='store_false')
        parser.add_argument('-W', help='Multiply group sizes by this factor to reduce number of jobs (for Genomes larger than Human)', default=1, type=int)
        parser.add_argument('-F', help='Color channel: replace -usecolor X in optArgs with this, must be either 1 or 2 [default OFF]', default=0, type=int)
        
        result = parser.parse_args()

        #if result.version : #old behavior was to only print on -v 1; new is to always print and log, and exit on > 1
        # need to wait until setupFileStructure though to have the pipereport available
        if result.version > 1 :
            print( "Pipeline Version: %s" % util.getMaxVersion() )
            sys.exit(0)
        
        #print "result.local:", result.local
        self.localRoot = (os.path.abspath(result.local) if result.local else "")
        #print "self.localRoot:", self.localRoot
        self.toolsFolder = (result.tools if result.tools else "")
        self.nThreads = result.T
        self.nPairwiseJobs = result.N #if result.N != -999 else result.T)
        self.pairwiseTriangleMode = True
        self.groupContigs=result.groupcontigs
        
        self.maxthreads   = (result.maxthreads if result.maxthreads > 0 else result.T)
        self.maxthreadsPW = (result.maxthreadsPW if result.maxthreadsPW > 0 else result.T )
        #print "maxthreadsPW= ", self.maxthreadsPW # DEBUG
        self.maxthreadsHost = result.J
        self.nExtensionIter = result.iter
        self.imgFile = None #result.img #deprecate
        self.expID = result.exp
        self.ref = result.ref
        self.autoNoise=result.y
        self.bypassCount = result.bypass
        self.bnxFileIn = result.bnx
        #self.subsampleMb = float(result.sub)
        self.noiseOnly = result.x
        if self.noiseOnly and not self.autoNoise :
            self.warning += 1
            self.message += "  WARNING: exiting after auto noise (-x) but auto noise is not enabled: enabling auto noise\n"
            self.autoNoise = True
        self.doNoise = False #result.noise #now _no longer_ assume doNoise if noiseOnly; retire this option
        self.cleanup = result.cleanup
        if self.cleanup < 0 or self.cleanup > 3 :
            self.error += 1
            self.message += "  ERROR: -c (cleanup) must be 0-3; %s supplied\n" % str(self.cleanup)
        self.genomeSize = 0 #float(result.size) #retire this option
        if result.cxml:
            self.onCluster = True
            self.clusterArgumentsFileIn = result.cxml
        self.optArgumentsFileIn = result.xml
        self.wipe = result.wipe
        self.lambdaRef = None #result.lambdaRef #deprecate
        self.time = False
        self.perf = False
        if result.perf == 1 or result.perf == 3 :
            self.time = True
        if result.perf == 2 or result.perf == 3 :
            self.perf = True
        if result.perf > 3 or result.perf < 0:
            self.error += 1
            self.message += '  ERROR:   -p must be between 0-3 input: %d\n' % result.perf
        self.workorder = False #result.workorder 
        self.contigSubDirectories = True #result.dirs #always True now
        #self.ngsInDir = result.ngsarg #retire
        self.ngsBypass = 0 #if not using, 0 is sufficient
        #if -F (ngsskiparg), must know how many stages to skip. The number is 6, ie, up to and including refineB
        #if result.ngsskiparg :
        #    self.ngsBypass = 7
        self.skipFinalRefine = result.nofinalref
        #self.numCharacterizeJobs = result.nCharJobs #retired
        self.stdoutlog = True #result.stdoutlog #always use this
        if result.runsv < 0 :
            self.message += '  ERROR:   argument -V must be integer >= 0: %i\n' % result.runsv
            self.error += 1
        self.runSV = result.runsv
        self.doAlignMol = result.A
        self.doAlignMolvRef = result.m
        self.bedFile = result.bed
        self.load_hg19 = False
        self.load_hg38 = False
        if result.H < 0 or result.H > 2 : #0 = no load, 1 = hg19, 2 = hg38
            self.message += '  ERROR: invalid -H argument: must be 0 (no loading), 1 (hg19), or 2 (hg38)\n'
            self.error += 1
        elif result.H == 1 :
            self.load_hg19 = True
        elif result.H == 2 :
            self.load_hg38 = True
        #self.doScanScale = result.Y
        self.groupSV = 0 #result.groupsv #retire option: always use single job
        self.groupRefineHostFr = result.f
        if self.groupRefineHostFr < 0 or self.groupRefineHostFr > 1. :
            self.message += "  ERROR: argument -h must be >= 0 AND <= 1. (%f)" % self.groupRefineHostFr
            self.error += 1
        self.argCompress = result.z
        self.stdoutRecheck = result.E
        if not self.stdoutRecheck :
            print "\t Not Re-Checking Stdout completeness\n"
        self.GroupSize = result.W
        if self.GroupSize != 1 :
            print "\t Multiplying Contig Group sizes by " , self.GroupSize
        if result.F < 0 or result.F > 2 :
            self.message += "  ERROR: argument -F must be >= 0 AND <= 2. (%f)" % result.F
            self.error += 1
        self.usecolorArg = result.F
        
    def setupFileStructure(self):
        aboveLocalRoot = os.path.split(self.localRoot.rstrip("/"))[0] #if localRoot ends in "/", need the rstrip
        if not util.checkDir(aboveLocalRoot, makeIfNotExist=False) :
            self.message += '  ERROR:   Parent Directory Non-existant: %s\n' % aboveLocalRoot
            self.error += 1
            return
        if not util.checkDir(self.localRoot) : #this will make if not exists
            self.message += '  ERROR:   Could not create %s\n' % self.localRoot
            self.error += 1
            return
    
        if self.wipe and os.path.exists(self.contigFolder):
            shutil.rmtree(self.contigFolder)
            
        if not self.noiseOnly :
            if not util.checkDir(self.contigFolder) : #this will make if not exists
                self.message += '  ERROR:   Could not create %s\n' % self.contigFolder
                self.error += 1
                return
            if not util.checkFile(self.idFile) :
                f1 = open(self.idFile, 'w')
                f1.close()
                if not util.checkFile(self.idFile) :
                    self.message += '  ERROR:   Bad ID file %s\n' % self.idFile
                    self.error += 1

        #if contigSubDirectories, then each one has its own alignref which is made inside CharacterizeModule
        if not(os.path.exists(self.contigAlignTarget)) and self.ref and not self.contigSubDirectories :
            os.mkdir(self.contigAlignTarget)
        #also wipe the align folder
        if self.wipe and os.path.exists(self.alignFolder):
            shutil.rmtree(self.alignFolder)
        if not self.noiseOnly :
            if not util.checkDir(self.alignFolder) :
                self.message += '  ERROR:   Could not create %s\n' % self.alignFolder
                self.error += 1
                return

        if self.onCluster:
            if not(os.path.exists(self.clusterLogDir)):
                os.mkdir(self.clusterLogDir)
        if self.bnxFileIn:
			if self.bnxFileIn == "/dev/stdin":
				f = open(self.bnxFile, 'w')
				for line in sys.stdin:
					f.write(line)
				f.close()
			else:
				self.bnxFileIn = os.path.abspath(self.bnxFileIn)
				if not util.checkFile(self.bnxFileIn, ".bnx") :
					self.message += '  ERROR:   bnx file does not end in ".bnx" or is not readable: %s\n' % self.bnxFileIn
					self.error += 1
					return
				#this file is not used after bnx_sort, which is the first stage, so no point copying if bypass
				if not(self.bnxFileIn == self.bnxFile) and self.bypassCount < 1 :
					shutil.copy(self.bnxFileIn, self.bnxFile)
        if self.optArgumentsFileIn :
            self.optArgumentsFileIn = os.path.abspath(self.optArgumentsFileIn)
            if not(self.optArgumentsFileIn == self.optArgumentsFile) and os.path.exists(self.optArgumentsFileIn) :
                shutil.copy(self.optArgumentsFileIn, self.optArgumentsFile)
        if self.clusterArgumentsFileIn:
            self.clusterArgumentsFileIn = os.path.abspath(self.clusterArgumentsFileIn)
            if not(self.clusterArgumentsFileIn == self.clusterArgumentsFile):
                shutil.copy(self.clusterArgumentsFileIn, self.clusterArgumentsFile)

        if self.ref :
            if not util.checkDir(self.refFolder) :
                self.updatePipeReport( "  ERROR: could not make ref dir %s\n" % self.refFolder )
                return
            newref = os.path.join(self.refFolder, os.path.basename(self.ref))
            if self.ref != newref : #if this is true, shutil.copy will raise an 'Error'
                shutil.copy(self.ref, newref)
            self.ref = newref
        
        util.InitStatus(os.path.join(self.localRoot, "status.xml"))
    
        pipeReport = ' '.join(sys.argv) + '\n\n'
        infoReport = ' '.join(sys.argv) + '\n\n'
        self.updatePipeReport(pipeReport, printalso=False)
        self.updateInfoReport(infoReport)
    #end setupFileStructure
        
    def checkInputs(self):
        if self.load_hg19 or self.load_hg38 : #if this parameter is used, it will override self.ref
            if self.ref : #warn user their reference is being ignored
                self.message += '  WARNING: ignoring reference (%s) because -H used--using HG19/38 instead\n' % self.ref
                self.warning += 1
            #some of this manipulation is also done in getScriptPaths, but I want to use the ref checking here
            refdir = os.path.join( os.path.dirname(os.path.realpath(__file__)), self.scriptsTopDir ) #Pipeline/Analysis
            refdir = os.path.join( refdir, self.scriptsSubdirs[0] ) #Analysis/SV
            refdir = os.path.join( refdir, self.scriptsSVdirs[2] ) #Analysis/SV/CopyNumberProfiles
            refdir = os.path.join( refdir, self.scriptsCopyNum_data ) #Analysis/SV/CopyNumberProfiles/Data
            if self.load_hg19 :
                refdir   = os.path.join( refdir, self.scriptsCopyNum_hg19dir ) #Analysis/SV/CopyNumberProfiles/Data/hg19
                self.ref = os.path.join( refdir, self.scriptsCopyNum_hg19 ) #Analysis/SV/CopyNumberProfiles/Data/hg19/hg19_chromosome_bspq.cmap
                self.bedFile= os.path.join( refdir, self.scriptsCopyNum_hg19bed ) #Analysis/SV/CopyNumberProfiles/Data/hg19/hg19_gaps.bed
            elif self.load_hg38 :
                refdir   = os.path.join( refdir, self.scriptsCopyNum_hg38dir ) #Analysis/SV/CopyNumberProfiles/Data/hg38
                self.ref = os.path.join( refdir, self.scriptsCopyNum_hg38 ) #Analysis/SV/CopyNumberProfiles/Data/hg38/...
                self.bedFile= os.path.join( refdir, self.scriptsCopyNum_hg38bed ) #Analysis/SV/CopyNumberProfiles/Data/hg38/...
            #if above paths are not present, below will catch it
        #now require ref is a cmap or spots file
        if self.ref:
            if not util.checkFile(self.ref, ".cmap") and not util.checkFile(self.ref, ".spots") :
                self.message += '  ERROR:   Reference not found or not readable or invalid file type (must be .spots or .cmap): %s\n' % self.ref
                self.error += 1
                self.ref = None
            else :
                self.ref = os.path.abspath(self.ref)
        if not self.ref and self.doNoise : #used to be or self.noiseOnly, now allow image processing with no ref
            self.message += '  ERROR: Cannot run sample characterization experiment (-n) with no reference (-r)\n'
            self.error += 1
        if not self.ref and self.autoNoise :
            if self.noiseOnly :
                self.message += '  ERROR: exiting after auto noise (-x) but no reference supplied, or invalid: please supply reference (-r) OR remove -x.\n'
                self.error += 1
            else :
                self.message += '  WARNING: Disabling auto noise detection (-y, default ON) because no reference supplied, or invalid.\n'
                self.warning += 1
                self.autoNoise = False
        if not self.ref and self.runSV : #require ref for runSV
            self.message += '  WARNING: Disabling SV detection (-V, default ON) because no reference supplied, or invalid.\n'
            self.warning += 1
            self.runSV = False
        #also require it is a cmap, because spots doesn't work with -i input which is required for SVModule
        elif self.ref and self.ref.endswith(".spots") :
            self.message += '  WARNING: Disabling SV detection (-V, default ON) because reference supplied is spots file: use .cmap to enable SV detection.\n'
            self.warning += 1
            self.runSV = False
        #ref also required for alignmolvref
        if not self.ref and self.doAlignMolvRef :
            self.message += '  WARNING: Disabling alignment of molecules vs reference (-m, default ON) because no reference supplied, or invalid.\n'
            self.warning += 1
            self.doAlignMolvRef = False
        #bed file check -- reference is already checked
        if self.bedFile and not self.ref :
            self.message += '  WARNING: bed file supplied (%s) but no reference supplied (-r), ignoring bed file\n' % self.bedFile
            self.warning += 1
            self.bedFile = ""
        elif self.bedFile :
            if not util.checkFile(self.bedFile, ".bed") :
                err = ': bed file not found or invalid file type (must be .bed): %s\n' % self.bedFile
                if self.load_hg19 or self.load_hg38 : #if this parameter is used, it will override self.ref
                    self.message += '  WARNING'+err
                    self.warning += 1 #change to warning in case auto-loading breaks
                else :
                    self.message += '  ERROR'+err
                    self.error += 1 #error if not load_hg19
                self.bedFile = ""
            else :
                self.bedFile = os.path.abspath(self.bedFile)
        #check bnx/image inputs
        if not self.bnxFileIn and not self.imgFile :
            self.message += '  ERROR:   Must provide either -I or -b as input\n'
            self.error += 1
        elif self.bnxFileIn and self.imgFile:
            self.message += '  ERROR:   Cant provide two experiment sources use -img or -bnx exclusive\n'
            self.error += 1
        if self.bnxFileIn:
            if not self.bnxFileIn == "/dev/stdin" and not util.checkFile(self.bnxFileIn, ".bnx") :
                self.message += '  ERROR:   Bnx File: %s not found or invalid type.\n' % self.bnxFileIn
                self.error += 1
        if self.imgFile:
            if not util.checkFile(self.imgFile) : #just a text file, no type checking
                self.message += '  ERROR:   Img File: %s  not found\n' % self.imgFile
                self.error += 1
            if self.bypassCount :
                self.message += '  ERROR:   Img File supplied, but bypassed (%i) : cannot bypass image processing\n' % self.bypassCount
                self.error += 1
                
        if self.optArgumentsFileIn:
            if not util.checkFile(self.optArgumentsFileIn) :
                self.message += '  ERROR:   XML File not found or not readable: %s\n' % self.optArgumentsFileIn
                self.error += 1
        else:
            self.message += '  ERROR:   -a XML File required\n' 
            self.error += 1
        
        if self.onCluster:
            if not os.path.exists(self.clusterArgumentsFileIn) :
                self.message += '  ERROR:   XML Cluster File not found %s\n' % self.clusterArgumentsFileIn
                self.error += 1
        
        if self.onCluster and self.imgFile:
            self.message +=  '  ERROR:   -C and -I, cant run img detection on cluster (I/O issue)\n' % self.imgFile
            self.error += 1
        
        if self.lambdaRef:
            if not os.path.exists(self.lambdaRef) :
                self.message +=  '  ERROR:   -L, Lambda Reference %s not found\n' % self.lambdaRef
                self.error += 1

        if self.ngsInDir :
            self.checkNgsDir() #for ngs dir, need various checks
            
        if self.numCharacterizeJobs > 1 and not(self.ref):
            self.message += '  WARNING:  Distributing characterize jobs (-D %d) without a reference (-r)\n' % self.numCharacterizeJobs
            self.message += '       Setting -D 1 and continuing\n'
            self.numCharacterizeJobs = 1
            self.warning += 1

        #check that -F (ngsBypass) isn't used without -f (ngsInDir), otherwise exitTestMinimumContigs is broken
        if self.ngsBypass and not self.ngsInDir :
            self.message +=  '  ERROR:   -F (skip minimum contigs check: %s) can only be used with -f (ngs contigs dir: %s) ' % (str(self.ngsBypass), self.ngsInDir)
            self.error += 1


    #setup ngs data members based on ngsInDir, which is a command line arg
    def checkNgsDir(self) :
        #check dir exists
        if not util.checkDir(self.ngsInDir, checkWritable=False, makeIfNotExist=False) :
            self.message += '  ERROR:   NGS dir does not exist %s\n' % self.ngsInDir
            self.error += 1
            return #abort if dir not found

        #check it has cmaps, set data member for it
        hascmap = False #guilty until proven innocent
        self.ngsContigPrefix = ""
        for ngsf in os.listdir(self.ngsInDir) :
            #If no "_contig", no valid prefix (this format is required by findContigs()), so no valid cmap
            if util.checkCmap( os.path.join(self.ngsInDir, ngsf) ) and util.hasSubStr(ngsf, "_contig") :
                hascmap = True #true after just one
                #the prefix is everything before "_contig". 
                newprefix = ngsf[:ngsf.find("_contig")]
                #check that prefixes are consistent -- require at least two (ie, confirm previous)
                if self.ngsContigPrefix and self.ngsContigPrefix == newprefix :
                    break #have confirmation, so done getting prefix
                else : #first time found
                    self.ngsContigPrefix = newprefix 
        if not hascmap :
            self.message += '  ERROR:   NGS dir does not contain any cmaps %s\n' % self.ngsInDir
            self.error += 1
        #print "self.ngsContigPrefix = ", self.ngsContigPrefix #debug

                
    def checkFlow(self):
        '''Check if command line arguments make sense--call immediately after parseCommandLine.'''
        if self.nPairwiseJobs <= 0:
            self.message += '  ERROR:   Number of Pairwise Jobs (-N) must be >= 1\n'
            self.error += 1
        if self.nThreads <= 0:
            self.message += '  ERROR:   Maximum number of threads (-T) must be >= 1\n'
            self.error += 1
        if self.maxthreads <= 0:
            self.message += '  ERROR:   number of threads per job (-j) must be >= 1\n'
            self.error += 1
        if self.maxthreads > self.nThreads:
            self.message += '  ERROR:   Maximum number of threads (-T) must be >= threads per job (-j)\n'
            self.error += 1
        if self.bypassCount > 4 and self.wipe:
            self.message += '  ERROR:   -B bypassing assembly and -w wiping previous contigs: remove either -B or -w\n'
            self.error += 1
        if self.subsampleMb < 0:
            self.message += '  ERROR:   -s subsample must be positive number\n'
            self.error += 1
        if self.nExtensionIter < 0 or self.nExtensionIter > 20 :
            self.message += '  ERROR:   -i (number of extension-merge stages) must be 0 or positive number <= 20\n'
            self.error += 1
        if self.genomeSize < 0 :
            self.message += '  ERROR:   -g (genome size) must be > 0\n'
            self.error += 1
        if self.groupSV < 0 or self.groupSV > 2 :
            self.message += '  ERROR:   -s (grouped SV) must be 0, 1, or 2\n'
            self.error += 1
            

    def assignFileTargets(self):
        
        if not self.localRoot :
            self.message += '  ERROR:   invalid output dir (-l argument): %s\n' % self.localRoot
            self.error += 1
            return            
        self.localRoot = os.path.abspath(self.localRoot)
        self.contigFolder = os.path.join(self.localRoot, 'contigs')
        self.alignFolder = os.path.join(self.localRoot, 'align')
        self.refFolder = os.path.join(self.localRoot, 'ref')
        if self.onCluster:
            self.clusterLogDir = os.path.join(self.localRoot,'ClusterLogs')
        self.alignMolvrefDir = os.path.join(self.contigFolder, self.alignMolvrefName)
        
        self.stageFolder["Contigs"] = self.contigFolder
        self.stageFolder["Align"] = self.alignFolder

        self.toolsFolder = os.path.abspath(self.toolsFolder)
        self.RefAlignerBin = os.path.join(self.toolsFolder, 'RefAligner')
        self.AssemblerBin = os.path.join(self.toolsFolder, 'Assembler')
        self.DMstaticBin = os.path.join(self.toolsFolder, 'DM-static')
        self.bnxFile = os.path.join(self.localRoot, 'all.bnx')
        self.bnxStatsFile = os.path.join(self.localRoot, 'molecule_stats.txt')
        self.bnxTarget = os.path.join(self.localRoot, 'bnxOut') #only used in ImageProcessingModule
        self.alignTarget = os.path.join(self.alignFolder, 'alignOut')
        self.contigPathTxtFile = os.path.join(self.contigFolder, 'curContigs')
        self.idFile = os.path.join(self.contigFolder, 'ID')
        self.contigAlignTarget = os.path.join(self.contigFolder, 'alignref')
        
        self.pipeReportFile = os.path.join(self.localRoot, self.expID + '_pipelineReport.txt')
        self.infoReportFile = os.path.join(self.localRoot, self.expID + '_informaticsReport.txt')
        self.infoReportFileSimple = os.path.join(self.localRoot, self.expID + '_informaticsReportSimple.txt')
        self.optArgumentsFile = os.path.join(self.localRoot, self.expID + '_optArguments.xml')
        self.clusterArgumentsFile = os.path.join(self.localRoot, self.expID + '_clusterArguments.xml')
        if os.path.exists(self.pipeReportFile):
            fileCt = 1
            while(True):
                testFileTag = '_%02d.txt' % fileCt
                testFileName = self.pipeReportFile.replace('.txt', testFileTag)
                if fileCt >= 99 or not os.path.exists(testFileName) :
                    self.pipeReportFile = testFileName
                    self.infoReportFile = self.infoReportFile.replace('.txt', testFileTag)
                    #self.infoReportFileSimple = self.infoReportFileSimple.replace('.txt', testFileTag)
                    self.optArgumentsFile = self.optArgumentsFile.replace('.xml', '_%02d.xml' % fileCt)
                    self.clusterArgumentsFile = self.clusterArgumentsFile.replace('.xml', '_%02d.xml' % fileCt)
                    break
                fileCt += 1
        self.memoryLogpath = os.path.join(self.localRoot, "memory_log.txt")
    #end assignFileTargets

    def checkToolPaths(self):
        if not util.checkExecutable(self.RefAlignerBin):
            self.message += '  ERROR: RefAligner not found or not executable, exiting\n'
            self.error += 1
        if self.bypassCount < 5 and not util.checkExecutable(self.AssemblerBin) :
            self.message += '  ERROR: Assembler not found or not executable, exiting\n'
            self.error += 1
        if not self.bnxFileIn and not util.checkExecutable(self.DMstaticBin) :
            self.message += '  ERROR: DM-static not found or not executable, exiting.\n'
            self.error += 1
        

    def checkCluster(self):
        if self.onCluster:
            try:
                import drmaa
            except:
                self.message += '  ERROR:  option -C:  Error importing drmaa library\n'
                self.message += '  ERROR:   try: sudo apt-get install python-drmaa\n'
                self.error += 1
            #try:
            #    sgeRoot = os.environ['SGE_ROOT']
            #except:
            #    self.message += ' option -C:  Requires enviroment variable "SGE_ROOT" to be set\n'
            #    self.error += 1
            #else: #only if no exception
            #    if not util.checkDir(sgeRoot, checkWritable=False, makeIfNotExist=False) :
            #        self.message += ' option -C:  Enviroment variable "SGE_ROOT" is not a valid dir: %s\n' % sgeRoot
            #        self.error += 1
    #end checkCluster

    def getAssemblyParams(self,mbInput):
        """Sets the pvalue thresholds based on genome size input parameter.
        """
    
        Tval = 10/(self.genomeSize*1e6) # Sets the pvalue threshold for matches based on genome size
        TEval = Tval/10
        T = str(Tval)
        TE = str(TEval)
        aveCov = mbInput / self.genomeSize
        #changing the minimum coverage is dangerous--below 5 maps, results shouldn't be trusted, so don't change these
        #minCov = str(min(20, max(2,int(aveCov/5))))#*
        #minAvCov = str(min(10, max(1,int(aveCov/10))))#*
        #new format is list--use small utility fn
        self.replaceParam('pairwise' , '-T', T)
        self.replaceParam('assembly' , '-T', T)
        self.replaceParam('refineA'  , '-T', T)
        self.replaceParam('refineB'  , '-T', TE)
        self.replaceParam('extension', '-T', TE)
        self.replaceParam('extension', '-TE', TE)

        infoReport = '  Average Coverage %3.2fx, of a %3.1f Mb Genome (user input)\n' % (aveCov, self.genomeSize)
        if aveCov < 30:
            infoReport += '  WARNING: Average Coverage %1.1fx, is below Min. Recommended (30x)\n' % aveCov
        self.updateInfoReport(infoReport, printalso=True)


    def replaceParam(self, stage, param, newval) :
        if not param in self.argData[stage] : #if it's not there, add it
            self.argData[stage].extend( [param, newval] )
            return
        idx = self.argData[stage].index(param) #if it is there, replace it
        self.argData[stage][idx+1] = newval #replace the next element
        

    def parseArguments(self, readingClusterFile=False):
        """Read XML file for runtime arguments
        
        """
        if not(readingClusterFile):
            targetFile = self.optArgumentsFileIn
        else:
            targetFile = self.clusterArgumentsFileIn

        try:
            tree = ET.parse(targetFile)
        except Exception,e: #the actual type is 'xml.etree.ElementTree.ParseError', but lets just catch 'em all
            self.error += 1
            self.message += '  ERROR: Invalid xml file %s: %s\n' % (targetFile, e)
            return
        # try to get xml.etree.ElementTree version to check compatibility
        requiredxmlver = "1.3.0" #from /usr/lib/python2.7/xml/etree/ElementTree.py
        xmlver = "" #try to get this from ET (xml.etree.ElementTree)
        try:
            #print "XML PARSER VERSION =", ET.VERSION 
            xmlver = ET.VERSION
        except :
            self.warning += 1
            self.message += '  WARNING: version of xml.etree.ElementTree not found\n'
        if xmlver != requiredxmlver and not readingClusterFile : #no need to do both files
            self.warning += 1
            self.message += '  WARNING: version of xml.etree.ElementTree is not recommended version (%s)\n' % requiredxmlver
        
        #print '  Using %s for arguments' % targetFile
        allArguments = {}
        for argGroup in tree.getroot().getchildren() : #loop on modules (pairwise, assembly, etc)
            curArgs = [] #use list to preserve order
            for aflag in argGroup.getchildren() : #each arg for this module
		if sys.version_info < (2,7,0):
			elements=aflag.getiterator()
		else:
			elements=aflag.iter()
                for pair in elements : #iterate in order: each element in module
                    include = False
                    if pair.tag == "include" :
                        include = True
                    elif pair.get("attr") : #skip comment lines; no attr in include tag
                        curArgs.append( pair.get("attr") ) #ie, -T
                    for key in sorted(pair.attrib.keys()) : #you have to do this for 'valN' to be in order
                        if key.startswith("val") :
                            if include and allArguments.has_key(pair.get(key)) :
                                include = False #use this to check to make sure something got included
                                curArgs.extend( allArguments[pair.get(key)] )
                            else :
                                curArgs.append( pair.get(key) )
                    if include : #if still True, nothing is included--this is an error
                        self.error += 1
                        self.message += "  ERROR: optArguments include value \"%s\" not found for group %s\n" % (pair.get(key), argGroup.tag)
            #special code for -contigsplit:
            # for common sections, this will get added again--check if it's already there
            if "-contigsplit" in curArgs and not curArgs[curArgs.index("-contigsplit")+3].endswith("ID") :
                #RefAligner requires two floats following -contigsplit, then the count file--assume two floats are present
                curArgs.insert(curArgs.index("-contigsplit")+3, self.idFile)
            if self.usecolorArg > 0 and '-usecolor' in curArgs :
                t = curArgs.index("-usecolor")
                if len(curArgs) > t+1 and not curArgs[t+1].startswith("-") : #replace
                    curArgs[ t+1 ] = str(self.usecolorArg)
                else :
                    curArgs.insert(t+1, str(self.usecolorArg))

	    ascii=1
	    for arg in curArgs:
		try:
			arg.decode('ascii')
		except:
			print("**** ERROR in xml file \"%s\": argument \"%s\" is not plain ASCII (UTF-8: \"%s\")\n" % (targetFile, arg.encode('ascii', 'replace'), arg.encode('UTF-8')))
			ascii=0

	    if not ascii:
		sys.exit(1)

            allArguments[argGroup.tag] = curArgs
            #load data into memoryArgGroups if clusterArgs
            if readingClusterFile and argGroup.tag in self.memoryArgGroups.keys() : 
                self.memoryArgGroups[argGroup.tag] = curArgs
            #load data from memoryArgGroups if optArgs and data is present for this stage
            elif ( not readingClusterFile and argGroup.tag in self.memoryArgGroups.keys() and
                   self.memoryArgGroups[argGroup.tag] ) :
                #print "Replacing section", argGroup.tag, ":", allArguments[argGroup.tag], ":", self.memoryArgGroups[argGroup.tag]
                allArguments[argGroup.tag] = self.memoryArgGroups[argGroup.tag]
        #end loop on xml
        if not(readingClusterFile):
            self.argData = allArguments
        else:
            self.clusterArgData = allArguments


    def checkRequiredArguments(self):
        """check that the xml file has all of the argument sets which will be 
        used by any pipeline module.
        """
        #argDataKeys = self.argData.keys() #all of the arg groups read in from the xml
        #eles of requiredArgGroups are always required
        for arg in self.requiredArgGroups :
            self.checkOptArgs    (arg, "required stage")
            self.checkClusterArgs(arg, "required stage")

        for arg in self.requiredOptArgGroups :
            self.checkOptArgs    (arg, "required stage")

        for arg in self.requiredClusterArgGroups :
            self.checkClusterArgs(arg, "required stage")

        #check for required arg sets based on command line flags
        if self.imgFile :
            self.checkOptArgs    ("imgDetection", "but image file %s supplied" % self.imgFile)
            self.checkClusterArgs("imgDetection", "but image file %s supplied" % self.imgFile)

        if self.doNoise :
            self.checkOptArgs    ("sampleChar", "but noise characterization (-n) specified")
            self.checkClusterArgs("sampleChar", "but noise characterization (-n) specified")

        if self.lambdaRef :
            self.checkOptArgs    ("mapLambda", "but lambda ref %s supplied" % self.lambdaRef)
            self.checkClusterArgs("mapLambda", "but lambda ref %s supplied" % self.lambdaRef)
            self.checkOptArgs    ("lambdaFilter", "but lambda ref %s supplied" % self.lambdaRef)
            self.checkClusterArgs("lambdaFilter", "but lambda ref %s supplied" % self.lambdaRef)

        if self.ngsInDir :
            self.checkOptArgs    ("refineNGS", "but ngs contig dir %s supplied" % self.ngsInDir)
            self.checkClusterArgs("refineNGS", "but ngs contig dir %s supplied" % self.ngsInDir)

        if not self.skipFinalRefine :
            self.checkOptArgs    ("refineFinal", "but not skipping final refinement (-u)")
            self.checkClusterArgs("refineFinal", "but not skipping final refinement (-u)")

        if self.runSV :
            self.checkOptArgs    ("svdetect", "but detecting SVs (-V, default ON)")
            self.checkClusterArgs("svdetect", "but detecting SVs (-V, default ON)")

        if self.groupContigs :
            groupstages = ["refineB0", "refineB1", "extension0", "extension1", "refineFinal0", "refineFinal1"]
            for stg in groupstages :
                self.checkOptArgs    (stg, "running grouped refine/extension (-U)")
                self.checkClusterArgs(stg, "running grouped refine/extension (-U)")

        if self.doAlignMol or self.doAlignMolvRef :
            self.checkOptArgs    ("alignmol", "but aligning molcules to final contigs (-V, default ON) or reference (-m, default ON)")
        if self.doAlignMol :
            self.checkClusterArgs("alignmol", "but aligning molcules to final contigs (-V, default ON)")
        if self.doAlignMolvRef :
            self.checkClusterArgs("alignmolvref", "but aligning molcules to reference (-m, default ON)")

        #print "autoNoise:", self.autoNoise, "checkoptargs:", self.checkOptArgs("autoNoise0", soft=True), self.checkOptArgs("autoNoise1", soft=True)
        if self.autoNoise and not (self.checkOptArgs("autoNoise0", soft=True) or self.checkOptArgs("autoNoise1", soft=True) ) :
            self.warning += 1
            err = "AutoNoise turned off because autonoise section (autoNoise0 or autoNoise1) is not present in optArguments.xml"
            self.message += "  WARNING: "+err+"\n"
            #self.updatePipeReport( err+"\n" )
            util.LogError("warning", err)
            self.autoNoise=False
            
        if self.autoNoise :	    
            self.checkOptArgs    ("autoNoise0", "but autoNoise enabled (-y)")
            self.checkOptArgs    ("autoNoise1", "but autoNoise enabled (-y)")
            self.checkClusterArgs("autoNoise0", "but autoNoise enabled (-y)")
            self.checkClusterArgs("autoNoise1", "but autoNoise enabled (-y)")
    #end checkRequiredArguments

    def checkOptArgs(self, stageName, errstr="", soft=False) :
        '''Utility method to simplify checkRequiredArguments.'''
        if not self.argData.has_key(stageName) or not len(self.argData[stageName]) :
            if soft:
		    return False
            self.error += 1
            self.message += '  ERROR: no parameter definitions for %s in %s' % (stageName, self.optArgumentsFileIn)
            if errstr :
                self.message += ', '+errstr
            self.message += '\n'
            return False
        return True

    def checkClusterArgs(self, stageName, errstr="") :
        '''Utility method to simplify checkRequiredArguments.'''
        if self.onCluster and not self.clusterArgData.has_key(stageName) : #if not on cluster, never error
            self.error += 1
            self.message += '  ERROR: no cluster parameter definitions for %s in %s' % (stageName, self.clusterArgumentsFileIn)
            if errstr :
                self.message += ', '+errstr
            self.message += '\n'


    def checkScanScaling(self) :
        """Check if -ScanScaling is in "autoNoise1" section of opt arguments file
        (argData). If so, and option to disable scan scaling is used, remove it."""
        #if not using autoNoise, silently disable scan scaling
        if not self.autoNoise or not self.argData.has_key(self.scanScaleSection) : #argsListed will raise exception if False
            self.doScanScale = False 
            return
        autonoise = self.argsListed(self.scanScaleSection) # "autoNoise1"
        havess = (self.scanScaleArg in autonoise) # "-ScanScaling" #the RefAligner option
        if not havess :
            if self.doScanScale : #print warning since intention was to run but won't
                err = "Warning: option for scan scaling not found in optArguments.xml file; not performing scan scaling"
                self.updatePipeReport( err+"\n" )
                util.LogError("warning", err)
            self.doScanScale = False #won't do it if it's not there
        elif not self.doScanScale : #if doScanScale is True, nothing else necessary
            idx = autonoise.index(self.scanScaleArg)
            nrm = 1
            for ele in autonoise[idx+1:] :
                if ele[0] != "-" : #not a new argument
                    nrm += 1
                else : #new arg
                    break
            del autonoise[idx:idx+nrm]
            self.doScanScale = False
            #print warning message that scanscaling is disabled
            err = "Warning: scan scaling in optArguments.xml will be removed because -Y was specified; not performing scan scaling"
            self.updatePipeReport( err+"\n\n" )
            util.LogError("warning", err)
        

    def setCoverageBasedParams(self):
        """Set coverage threshold values based on input coverage depth."""
        if self.genomeSize:
            print('  Reading coverage depth from bnx')
            #dset = molecule.moleculeDataset(500, PitchNm = 700.) #old code
            #dset.readBnxFile(self.bnxFile, mbOnly=True)
            #mbInput = dset.MegabasesDetected
            bnx = util.bnxfile(self.sorted_file+".bnx", thresh=[0])
            self.getAssemblyParams(bnx.getTotLenMb())
            
    #this feature is retired
    #def subSampleBnx(self):
    #    """Subselect molecule from bnx file randomly for reduced coverage input. """
    #    if self.subsampleMb:
    #        self.error = molecule.subsampleBnx(self, bypass = not(self.executeCurrentStage))
            
    def argsListed(self, moduleName):
        """Get arguments for stage moduleName from self.argData
        which were read from optArgumentsFileIn xml file in self.parseArguments."""
        return self.argData[moduleName]


    def getClusterArgs(self, moduleName, category=None):
        """Get cluster arguments for stage moduleName from self.clusterArgData
        which were read from clusterArgs xml file in self.parseArguments."""
        if not self.onCluster :
            return None
        #new behavior: if key is present, proceed properly
        # if not, print an error
        if category:
		mname=moduleName+category
	else:
		mname=moduleName
        if self.clusterArgData.has_key(mname) :
            return ' '.join(self.clusterArgData[mname])
        if category and moduleName!="default":
		return self.getClusterArgs("default", category)
        print '  Failed to read cluster Arguments for module %s section %s' % (moduleName, mname)
        return ''

    
    def prerunLog(self):
        """At start of DNPipeline, call this method to log information before any
        processing starts. So far, only Pipeline version and start time.
        """
        #record version here in pipereport only
        self.updatePipeReport( "  Pipeline Version: %s\n\n" % util.getMaxVersion() )
        raver = util.getRefAlignerVersion(self.RefAlignerBinOrig) #now return int
        if raver :
            self.updatePipeReport( "  RefAligner Version: %i\n\n" % raver )
            self.refaligner_version = raver
        else : #None or False means unable to get version
            raver = ""
            err = "Unable to determine RefAligner version"
            self.updatePipeReport( err+"\n\b" )
            util.LogError("warning", err)
        #print "self.refaligner_version =", self.refaligner_version #debug
        #requirement from Software: introduce a version string in the informatics report to aid parsing
        self.updateInfoReport("Informatics Report Version: %s\n\n" % self.infoReportVersion)
        #put time here because it looks better after the Prerun Tests in pipelineCL.py
        self.startTime = time.time() #time since Epoch
        self.updatePipeReport( "  Pipeline start time: %s\n\n" % time.ctime(self.startTime) )
        util.initMemoryLog(self.memoryLogpath)
        self.checkScanScaling()

    def printMessage(self):
        """This is called from pipelineCL.py only if there is an error or warning.
        Log and print message, and also version info if error.
        """
        if self.pipeReportFile and util.checkFile(self.pipeReportFile) : #if error occurs in checkFlow, this is necessary
            self.updatePipeReport(self.message+"\n")
        else :
            print self.message
        if self.error :
            print( "Pipeline Version: %s\n" % util.getMaxVersion() )
        sys.stdout.flush()

		
    def FinalCleanup(self):
        """Remove intermediate cmap files for simpler resulting file structure
        Feature requested by software group (-c command line option)
        """
        if self.cleanup == 0 or self.cleanup == 2 :
            return
        #finalMergedPrefix = self.prefixUsed[-1]
        #src = os.path.join(self.contigFolder, finalMergedPrefix + '.cmap')
        #finalCmap = self.expID + '.cmap'
        #dst = os.path.join(self.contigFolder, finalCmap)
        #shutil.move(src, dst)
        #finalPrefix = finalMergedPrefix.lower()
        #allFiles = os.listdir(self.contigFolder)
        #for curFile in allFiles:
        #    if curFile.startswith(finalPrefix):
        #        continue
        #    if curFile == finalCmap:
        #        continue
        #    pathname = os.path.join(self.contigFolder, curFile)
        #    if os.path.isdir(pathname):
        #        continue
        #    try:
        #        os.remove(pathname)
        #    except:
        #        pass
        #above old code
        util.LogStatus("progress", "stage_start", "FinalCleanup")
        #remove these dirs: all are sub-dirs of contigs, do not remove refineFinal{,_sv}, auto_noise, alignmolvref
        remove = ["unrefined", "refineA", "refineB", "extension", "mrg"]
        for qfile in os.listdir(self.contigFolder) :
            fullpath = os.path.join(self.contigFolder, qfile)
            if not os.path.isdir(fullpath) :
                continue
            for rstr in remove :
                if qfile.find(rstr) == -1 :
                    continue
                try :
                    shutil.rmtree(fullpath)
                    self.updatePipeReport("FinalCleanup: removed dir %s\n" % fullpath)
                except:
                    pass
        util.LogStatus("progress", "stage_complete", "FinalCleanup")
    #end FinalCleanup
	
    def ImportSQLite(self, folder=None):
        """Use irys_db_tool.R to create sqlite files with data in contig sub-dirs.
        Folder argument is key of self.stageFolder, or None if using self.contigFolder.
        Only call script if cleanup >= 2 supplied on command line. Requires Rscript."""
        if self.cleanup < 2 :
            return
        if not self.Rscript_path :
            return

        if folder == None :
            folder = self.contigFolder
        elif self.stageFolder.has_key(folder) :
            folder = self.stageFolder[folder]
        else :
            self.updatePipeReport("Warning in ImportSQLite: folder \"%s\" missing\n" % str(folder))
            return

        outarg = folder+".sqlite"
        
        if os.path.exists(outarg): #this will crash anyway
            self.updatePipeReport("Warning in ImportSQLite: output file \"%s\" already exists\n" % str(outarg))
            #traceback.print_exc(file=sys.stdout)
            return
        
        args = [self.Rscript_path, os.path.join(self.pipeline_dir, "irys_db_tool.R") , "load", outarg,  folder]
        logstr = ("Loading folder %s into %s\n") % (folder, outarg)
        self.updatePipeReport(logstr) #, printalso=False)
        #print(args)
        so=open(folder+".sqlite.stdout", 'w')
        se=open(folder+".sqlite.stderr", 'w') #todo: combine these two
        subprocess.Popen(args, stdout=so, stderr=se).wait()
        so.close()
        se.close()
    #end ImportSQLite

        
    def prepareContigIO(self, contigPrefix, name=None, updateDirs=True):
        ''' Set the input folder and prefix, as well as the output 
        '''
        self.inputContigPrefix = self.outputContigPrefix
        self.outputContigPrefix = contigPrefix
        if self.contigSubDirectories:
            if updateDirs:
                self.inputContigFolder = self.outputContigFolder
                self.outputContigFolder = os.path.join(self.contigFolder, contigPrefix)
                #will make if not exist, only returns False if already exists or can't make; don't make dir if errors
                if not self.error and not util.checkDir(self.outputContigFolder) : 
                    self.updatePipeReport("ERROR in varsPipeline.prepareContigIO: bad dir: "+self.outputContigFolder+"\n")
        else:
            self.inputContigFolder = self.contigFolder
            self.outputContigFolder = self.contigFolder
        self.prefixUsed.append(self.outputContigPrefix)
        if name :
		self.stageFolder[name]=self.outputContigFolder


    def mergeIntoSingleCmap(self, outdir=""):
        """Merge contigs from current stage into single cmap file for subsequent processing.
        Default output dir is self.outputContigFolder--override this with outdir if supplied"""
        debug = False #True
        #log = False #not used
        if debug :
            print "\ncontigFolder:", self.outputContigFolder #debug #this is really the input folder for this merge
            print "outputContigPrefix:", self.outputContigPrefix #debug #the file prefix of contig cmaps
            print "contigPathTxtFile:", self.contigPathTxtFile #debug #write all contig cmap paths to this text file
        
        # new variable : self.numCharacterizeJobs - distribute characterization job (by size)
        contigFiles, contigIDs = self.findContigs(self.outputContigFolder, self.outputContigPrefix, txtOutput=self.contigPathTxtFile)
        self.curContigCount = len(contigFiles)
        if debug :
            print "curContigCount:", self.curContigCount, "\n"
        if len(contigFiles) == 0:
            self.latestMergedCmap = None
        if len(contigFiles) == 1:
            self.latestMergedCmap = os.path.join(self.outputContigFolder, contigFiles[0])

        mergeTargets = [self.contigPathTxtFile]
        if self.numCharacterizeJobs > 1:
            mergeTargets += self.curCharacterizeTargets
        
        stageName = 'Cmap_Merge_' + self.outputContigPrefix #self.curContigPrefix
        mergeJobSet = mthread.jobWrapper(self,stageName,clusterArgs=self.getClusterArgs('cmapMerge'))
        self.curCharacterizeCmaps = []
        for i,mergeTarget in enumerate(mergeTargets):
            cargs = [self.RefAlignerBin]
            cargs += ['-if', mergeTarget]
            if i == 0:
                mergeFileName = util.uniquifyContigName(self.outputContigPrefix)
                self.prefixUsed.append(mergeFileName)
                if outdir :
                    outputTarget = os.path.join(outdir, mergeFileName)
                else :
                    outputTarget = os.path.join(self.outputContigFolder, mergeFileName)
                expectedResultFile = outputTarget + '.cmap'
                keyResult = expectedResultFile
                if self.numCharacterizeJobs == 1:
                    self.curCharacterizeCmaps.append(expectedResultFile)
            else:
                outputTarget = os.path.join(self.outputContigFolder, mergeFileName + '_partial_%d' % (i))
                expectedResultFile = outputTarget + '.cmap'
                self.curCharacterizeCmaps.append(outputTarget + '.cmap')
            cargs += ['-merge', '-o', outputTarget, '-f']
            if self.stdoutlog :
                cargs.extend( ['-stdout', '-stderr'] )
            jobName = 'mrg_%s_%d' % (self.outputContigPrefix,i)
            #if log : #not used
            #    outpre = os.path.join(self.logdir, "log_merge_"+self.stageComplete)
            #    sJob = mthread.singleJob(cargs, jobName, expectedResultFile, jobName, stdOutFile=outpre+"stdout.log", stdErrOutFile=outpre+"stderr.log", clusterLogDir=self.clusterLogDir)
            #else :
            sJob = mthread.singleJob(cargs, jobName, expectedResultFile, jobName, clusterLogDir=self.clusterLogDir)            
            mergeJobSet.addJob(sJob)
        if self.executeCurrentStage:
            mergeJobSet.multiThreadRunJobs(1, callLogStatus=False) # run locally; no need for status of these
            mergeJobSet.doAllPipeReport()
        if os.path.exists(keyResult):
            self.latestMergedCmap = keyResult
        else:
            self.latestMergedCmap = None
            #self.curContigCount = 0
        
    def findContigs(self,contigFolder, searchString, txtOutput=None):
        """Using string match, find phase-specific contigs
        
        Output is list of load balanced cmap list
        Optional output is load balanced cmap files or single cmap files
        """
        
        fileList = os.listdir(contigFolder)
        searchString += '_contig'
        contigFiles = []
        contigIDs = []
        fileSizes = []
        for fileName in fileList:
            match = fileName.find(searchString)
            if match == -1 or fileName.endswith(".tmp"):
                continue
            prefixEnd = match + len(searchString)
            targetFile = os.path.join(contigFolder, fileName)
            statRes = os.stat(targetFile)
            if statRes.st_size == 0:
                continue
            suffixStart = fileName.find('.cmap')
            cmapString = '.cmap'
            if fileName.endswith(cmapString):
                suffixStart = len(fileName) - len(cmapString)
                #cmapRefinedString = '_refined.cmap'
                #if fileName.endswith(cmapRefinedString):
                #    suffixStart = fileName.__len__() - cmapRefinedString.__len__()
            else:
                continue
            match = fileName.find('condensed')
            if match != -1:
                continue
            contigNumber = fileName[prefixEnd:suffixStart]
            if any([x.isalpha() for x in contigNumber]):
                #print "skipping contig %s, %s" % (fileName, contigNumber) #debug
                continue
            contigFiles.append(targetFile)
            contigIDs.append(contigNumber)
            fileSizes.append(statRes.st_size)
        
        if not(txtOutput):
            return contigFiles, contigIDs

        f1 = open(txtOutput, 'wb')
        for contigFile in contigFiles:
            f1.write(contigFile + '\n')
        self.curCharacterizeTargets = [txtOutput]
        
        # Split up the Cmaps for Characterization
        nBins = min(self.numCharacterizeJobs, contigFiles.__len__())
        if self.numCharacterizeJobs > 1 and self.ref:
            self.curCharacterizeTargets = []
            binnedContigFiles = BinContigsBySize(contigFiles, fileSizes, nBins)
            for i,ContigGroup in enumerate(binnedContigFiles):
                partialTxtOutput = txtOutput + '_%d' % (i+1)
                f1 = open(partialTxtOutput, 'wb')
                for contigFile in ContigGroup:
                    f1.write(contigFile + '\n')
                f1.close()
                self.curCharacterizeTargets.append(partialTxtOutput)
        return contigFiles, contigIDs
    #end findContigs


    def copyToRefineFinal(self, verbose=False) :
        """In the case when refineFinal doesn't run, copy the last
        stage's results (should always be refineB) to the expected
        location of the refineFinal contigs so that IrysView can
        find them.
        For now, the sv folder is not copied.
        """
        source = self.copyRFdir # this should be where the contigs are, ie, should be refineB -- use new data member
        pref   = os.path.split(source)[0]
        target = os.path.join(pref, self.expID + "_" + self.refineFinal1) # ie, 'exp_refineFinal1'
        if not util.checkDir(target) : #will make if doesn't exist; if can't make or exists and not a dir, report error and return
            err = "Error making refineFinal dir %s, or exists and not dir; cannot copy from %s" % (target, source)
            self.updatePipeReport( err+"\n" )
            util.LogError("error", err)
            return
        self.updatePipeReport( "Copying contigs from " + source + "\nto " + target + "\n\n" )
        replace = ['rB', 'refineB1'] #replacing--replace with refineFinal1
        #loop recursively in source using os.walk
        for root, dirs, files in os.walk(source) :
            #if 'alignmol' in dirs : 
            #    dirs.remove('alignmol') #on second thought, keep it
            if root == source :
                usetarget = target
            else : #here, root is sub-dir of source, ie, alignref
                usetarget = os.path.join(target, os.path.split(root)[1])
                util.checkDir(usetarget) #make the alignref subdir
            for ifile in files :
                path = os.path.join(root, ifile)
                if verbose :
                    print "copying", path, "to", usetarget
                shutil.copy(path, usetarget)
                if path.endswith(".xmap") : #for xmaps, copy twice, the second time changing the name
                    for rep in replace :
                        copy = False
                        pathbase,name = os.path.split(path)
                        if name.find(rep) != -1 :
                            copy = True
                            outname = os.path.join(usetarget, name.replace(rep, self.refineFinal1))
                        elif name.find(rep.upper()) != -1 :
                            copy = True
                            outname = os.path.join(usetarget, name.replace(rep.upper(), self.refineFinal1.upper()))
                        if copy :
                            copytarget = os.path.join(usetarget, outname)
                            if verbose :
                                print "copying", path, "to", copytarget
                            shutil.copy(path, copytarget)
    #end copyToRefineFinal

        
    def updatePipeReport(self, content, printalso=True):
        if printalso :
            print content, #no extra newline
        #if content.startswith("ERROR") : #this is normally redundant or not necessary
        #    util.LogError("error", content)
        f1 = open(self.pipeReportFile, 'a')
        f1.write(content) 
        f1.close()

    #default printalso False, unlike PipeReport
    def updateInfoReport(self, content, printalso=False, simple=False, simplenl=False): 
        if printalso :
            print content, #no extra newline
        f1 = open(self.infoReportFile, 'a')
        f1.write(content) 
        f1.close()
        if simple :
            f1 = open(self.infoReportFileSimple, 'a')
            f1.write(content)
            if simplenl :
                f1.write("\n")
            f1.close()

    #if not ismerge, we want to exit. If it is, we don't, ie, for bacteria which merge to single contig.
    #return 0 for no error; return 1 for 0 contigs (no further processing); return 2 for 1 contig (skip to final refine)
    def exitTestMinimumContigs(self, minCount, ismerge=False):
        '''Set pipeline exit condition if contigs too few to continue
        
        '''
        #for ngsBypass, in addition to skipping the stages themselves (in toExecuteCurrentStage),
        # also skip this check (min contigs) because there need not be any results of assembly and refineA/B
        if self.ngsBypass and self.currentStageCount <= self.ngsBypass :
            return 0
        if self.curContigCount <= minCount:
            if not ismerge or self.curContigCount <= 0 : #if you have 0 contigs, always an error
                exitString = '\nERROR: Contig count %d <= %s minimum %d, exiting\n' % (self.curContigCount, self.stageComplete, minCount)
                retval = 1
                util.LogError("critical", "stage %s did not produce minimum number of contigs" % self.stageComplete)
            else :
                exitString = '\nContig count %d <= %s minimum %d\n' % (self.curContigCount, self.stageComplete, minCount)
                util.LogError("warning", "single contig after %s" % self.stageComplete)
                retval = 2
            self.updatePipeReport(exitString)
            return retval
        #print "exitTestMinimumContigs:", self.stageComplete, ":", self.curContigCount #debug
        return 0


    def runJobs(self, module, name="") :
        '''varsPipeline.runJobs: wrapper to module.runJobs:
        will just use the jobWrapper class in Multithreading.py.
        Note: all modules _must_ have the runJobs method implemented.
        '''
        if not self.workorder : #default (SGE, IrysView, local)
            module.runJobs()
        #else : #celery ONLY -- legacy code no longer used
        #    from opensgi.cloud.api.tasks import runJobGroup
        #    runJobGroup(self.workorder, name, module) # Celery Specific


    # targetFile is a list of bnxfiles which are in outputList
    # this is used here for SampleCharModule, and in PairwiseModule
    # this was previously a standalone fn called writeIntermediate
    def writeListToFile(self, outputList, targetFile):
        f1 = open(targetFile, 'w')
        for val in outputList: 
            writeVal = os.path.abspath(val)
            f1.write(writeVal + '\n')
        f1.close()   

    def checkDependencies(self) :
        '''Check whether R and perl are installed, record results.'''
        self.R_path = None
        self.Rscript_path = None
        self.perl_path = None
        self.platform = util.getPlatform() #1 for linux, 2 for windows, None for error
        if not self.platform :
            errstr = "Unable to detect platform--calls to external scripts disabled"
            util.LogError("warning", errstr)
            self.updatePipeReport(perl[1]+"\n\n")
            return

        if True : #now use for copy number (SVModule)
            R = util.checkExternalDependencies('R', islinux=(self.platform == 1))
            self.R_path = R[0] #will be None if not found
            if R[1] : #if not None, error occurred
                util.LogError("warning", R[1])
                self.updatePipeReport(R[1]+"\n\n")
            # print R #test
        if True : #now use for copy number (SVModule)
            Rscript = util.checkExternalDependencies('Rscript', islinux=(self.platform == 1))
            self.Rscript_path = Rscript[0] #will be None if not found
            if Rscript[1] : #if not None, error occurred
                util.LogError("warning", Rscript[1])
                self.updatePipeReport(Rscript[1]+"\n\n")
            # print R #test
            
        if self.platform == 1 : #linux only--for windows we will use executable, so no need to check for perl
            perl = util.checkExternalDependencies('perl', islinux=(self.platform == 1))
            self.perl_path = perl[0] #will be None if not found
            if perl[1] : #if not None, error occurred
                util.LogError("warning", perl[1])
                self.updatePipeReport(perl[1]+"\n\n")
            # print perl #test
            # print util.runJob([self.perl_path, '-e' ,'print "Hello\n"'], returnstdout=True, printfail=True) #example one-liner
        self.getScriptPaths() #restore me
        #check perf iff perf will be used, ie, this condition:
        if self.perf > 1 :
            perf = checkExternalDependencies('perf', islinux=(self.platform == 1))
            if perf[0] == None :
                self.error += 1
                self.message += '  ERROR: perf is enabled (-p >= 2) but not found or not executable (%s)\n' % (perf[1])
    #end checkDependencies

    def getScriptPaths(self) :
        '''Load data members for external scripts.'''

        self.pipeline_dir   = os.path.dirname(os.path.realpath(__file__))
        script_dir     = os.path.join(self.pipeline_dir, self.scriptsTopDir) #Analysis
        sv_dir         = os.path.join(script_dir  , self.scriptsSubdirs[0]) #Analysis/SV
        self.rutil_dir = os.path.join(script_dir  , self.scriptsSubdirs[1]) #Analysis/RUtil
        indelcon_dir   = os.path.join(sv_dir      , self.scriptsSVdirs[0]) #Analysis/SV/InDelConfidence
        mergeinv_dir   = os.path.join(sv_dir      , self.scriptsSVdirs[1]) #Analysis/SV/MergeInversions
        self.scriptsIndelConf   = os.path.join(indelcon_dir, self.scriptsIndelConf_path) #just script path
        self.scriptsCopyNum_dir = os.path.join(sv_dir      , self.scriptsSVdirs[2]) #Analysis/SV/CopyNumberProfiles
        #For windows, we use an executable, but for linux, use perl_path and script name.
        # Make this a list for easier use in SVModule
        if self.platform == 1 : #linux
            self.scriptsMergeInv = [self.perl_path, os.path.join(mergeinv_dir, self.scriptsMergeInv_linux)]
        else : #windows
            self.perl_path = 1 #Hack to enable executable for windows to avoid installing perl--careful if more perl scripts are introduced which are not executable
            self.scriptsMergeInv = [os.path.join(mergeinv_dir, self.scriptsMergeInv_win)]

        #more to add in the future...

    def doCompression(self) :
        '''Call bash script to tar/gzip results'''

        if not self.argCompress :
            return

        stageName = "compress"
        scriptpath = os.path.join(self.pipeline_dir, self.scriptsCompress)
        if not util.checkFile(scriptpath) :
            err = "missing Pipeline compression script: %s" % scriptpath
            util.LogError("error", err)
            self.updatePipeReport("Error:"+err+"\n")
            return

        util.LogStatus("progress", "stage_start", stageName) #after above bc check if bypass (executeCurrentStage)

        #Note: in Multithreading.py, /bin/bash is hardcoded as well;
        #if problem, replace with call to checkExternalDependencies (and update for bash)
        script =  ["/bin/bash", scriptpath]
        script += ["--targetFolder", self.localRoot]
        script += ["--outputFolder", self.localRoot]
        script += ["--prefix", self.scriptsCompressName, "--gzip", "--verbose"]
        #print script #debug
        error = False
        stdout = joberr = ""
        retc = 1
        self.updatePipeReport("Starting compression script:\n"+" ".join(script)+"\n")
        try :
            retc, stdout, joberr = util.runJob(script, returnstdout=True) #bc returnstdout is True, return tuple
        except OSError : 
            err = "compress job failed: %s" % joberr
            util.LogError("error", err)
            self.updatePipeReport("Error:"+err+"\n")
            error = True

        if stdout :
            log = os.path.join(self.localRoot, self.scriptsCompressName+".log")
            self.updatePipeReport("Creating compress job log: "+log+"\n")
            f = open(log, 'w') #save to log file
            f.write(stdout)
            f.close()
        result = os.path.join(self.localRoot,self.scriptsCompressName+".tar.gz")
        if not error and (retc or not util.checkFile(result)) : #retc = 0 is success
            err = "compress failed; missing compression result: %s" % result
            util.LogError("error", err)
            self.updatePipeReport("Error:"+err+"\n")

        util.LogStatus("progress", "stage_complete", stageName)
    #end doCompression

#end class varsPipeline



## Execution of De Novo assembly workflow

class DNPipeline():
    """ Execution of De Novo assembly workflow.    
    Flow is determined by varsPipeline class.
    This class is essentially a wrapper for all of the Pipeline modules.
    Called from pipelineCL.py.
    """

    def __init__(self):
        self.lastCharacterize = None #used to store results of characterize module calls
    
            
    def run(self, varsP):
        """Main pipeline method.        
        Executes pipeline operations prescribed by varsPipeline class.
        """ 
        os.putenv('O64_OMP_SET_AFFINITY', 'false')  
        
        varsP.prerunLog() #general information in log

        if self.constructData(varsP): #returns 0 on successful completion, 1 on failure or exit after sample char
            self.finalizePipeline(varsP)
            return #end of Pipeline
        
        #varsP.subSampleBnx() #this feature is retired

        #pairwise, assembly, refineA/B, merge0, iterative extend-merge
        try:
            skiplevel = self.runAllStages(varsP) 
        except RuntimeError :
            self.finalizePipeline(varsP)
            return #end of Pipeline

        #this is for the case of merging to a single contig after refine B (another refinement isn't necessary)
        copyrb = False #copy refineB to refineFinal, but after runFinalStages in order to wait for threads
        if not skiplevel : #0 or None
            try:
                skiplevel = self.runRefineFinal(varsP)
            except RuntimeError :
                self.finalizePipeline(varsP)
                return #end of Pipeline
        elif skiplevel != 1 : #1 is error, so no copy
            varsP.copyRFdir = varsP.outputContigFolder #copy this bc it gets reset in SVModule below
            copyrb = True

        #alignmol and SV: as long as you have at least one contig, these should run
        self.runFinalStages(varsP, not skiplevel or skiplevel == 2) #2nd arg will determine whether to run 

        if copyrb :
            varsP.copyToRefineFinal() #for IrysView

        print
        ers = util.SummarizeErrors(varsP=varsP)
        #don't cleanup if error so we can debug; if critical, wouldn't have reached this point anyway
        if not ers.has_key('error') : 
            varsP.FinalCleanup()

        self.finalizePipeline(varsP, ers)
    #end run
            

    def constructData(self, varsP):
        """Build bnx files from raw image data, or prepare input bnx for stages to follow.
        Optionally run single molecule noise characterization.
        """
        # Perform Image Processing:
        varsP.toExecuteCurrentStage(quiet=True)
        if not varsP.bnxFileIn :
            #new return of this is 1 for failure, False or None for success
            import ImageProcessingModule as ipm #this is only used here
            if ipm.performImageAnalysis(varsP, bypass=not(varsP.executeCurrentStage)) :
                return 1
        
        alignmod.getAlignStats(varsP, None, bnxpath=varsP.bnxFile) #pre-sort (ie, pre-filter) bnx stats
        
        if not varsP.noiseOnly: #if noiseOnly, no need for sort/split
            try:
                pm.sortBNX(varsP) #raise RuntimeError if job fails
            except RuntimeError :
                return 1

        varsP.setCoverageBasedParams() #move this to DNPipeline.constructData so that it's after bnx_sort
              
        if varsP.autoNoise:
            try :
                scm.autoNoise(varsP)
            except RuntimeError :
                return 1

        if not varsP.noiseOnly: #if noiseOnly, no need for sort/split
            try:
                pm.splitBNX(varsP) #raise RuntimeError if job fails
            except RuntimeError :
                return 1

        # Sample Characterization
        if varsP.toExecuteCurrentStage(quiet=True) and varsP.doNoise:
            util.LogStatus("progress", "stage_start", "Noise characterization")
            stageSC = scm.SampleChar(varsP)
            varsP.runJobs(stageSC, "SampleChar")
            stageSC.checkResults() 
            util.LogStatus("progress", "stage_complete", "Noise characterization")

        if varsP.noiseOnly:
            varsP.updatePipeReport('Finished noise evaluation -noiseonly, exiting\n') #print and log sucessful completion message
            return 1 #like failure in that it causes pipeline to exit
        
        return 0
    #end constructData

    def DoMolVRef(self, arefmod, varsP):
	util.LogStatus("progress", "stage_start"   , arefmod.stageName)
	varsP.runJobs(arefmod, "alignmolvref")
	arefmod.checkResults()
	util.LogStatus("progress", "stage_complete", arefmod.stageName)
        arefmod.getAlignStats()

    def runAllStages(self, varsP) :
        """Run stages: Pairwise, Assembly, RefineA, RefineB, merge0,
         and iterative extend-merge stages (but not refineFinal, SV, or alignmol).

         Stage execution entails the following:
         - Instantiate Stage Class (from included module)
         - Populate Stage Jobs - hardcoded + varsPipeline
         - Run Stage Jobs
         - Log Stage Data: computation data and informatics data
        """

        # Perform Pairwise Comparison:
        util.LogStatus("progress", "stage_start", "Pairwise")
        stagePW = pm.Pairwise(varsP)
        if varsP.toExecuteCurrentStage():
            #stagePW.runJobs()               # Single Node or Grid Engine
            varsP.runJobs(stagePW, "Pairwise")
        util.LogStatus("progress", "stage_complete", "Pairwise")
        if stagePW.checkResults() : # Update report here -- returns 1 if all align files are missing
            return 1 #1 is failure (a-la bash)
        
        util.LogStatus("progress", "stage_start", "Assembly")
        # Denovo Assembly
        stageASSMBL = am.Assemble(varsP)
        if varsP.toExecuteCurrentStage():
            #stageASSMBL.runJobs()
            varsP.runJobs(stageASSMBL, "Assembly")
        stageASSMBL.checkResults()
        util.LogStatus("progress", "stage_complete", "Assembly")
        if varsP.exitTestMinimumContigs(0):
            return 1 #1 is failure (a-la bash)
        self.runCharacterize(varsP)

        if varsP.doAlignMolvRef : #move this here instead of before pairwise because if these jobs are run on host, they can hold up the assembly job (also run on host)
            #arefmod = alignmod.AlignRefModule(copy.deepcopy(varsP)) #constructor will call generateJobList--use copy in case background
            #no more AlignRefModule--instead, AlignModule
            arefmod = alignmod.AlignModule(copy.deepcopy(varsP), doref=True) #constructor will call generateJobList--use copy in case background
            if varsP.executeCurrentStage : #this is not an enumerated stage, but if bypassing, skip it
                if varsP.onCluster : #only background on cluster
                    mthread.start_new_thread(self.DoMolVRef, (arefmod, varsP))
                else :
		    self.DoMolVRef(arefmod, varsP)
            else :
                arefmod.checkResults() #want this here even if bypassing so that report has mol stats in it
            #if not arefmod.checkRequiredCopyNumberArgs() : #required for copy number
            #    varsP.doAlignMolvRef = False #disable copy number
            #no longer disable copy number
            #arefmod.checkRequiredCopyNumberArgs() #required for copy number
        elif not varsP.doAlignMol : #need to get mol stats for all_sorted.bnx if not doing alignmol
            alignmod.getAlignStats(varsP, None, bnxpath=varsP.sorted_file+".bnx") 
        
        # Contig Refine A
        if varsP.groupContigs:
		stageRFNA = grprefmod.Refine("refineA", varsP) #new class requires first arg string
	else:
		stageRFNA = refmod.Refine("refineA", varsP) #new class requires first arg string
        if varsP.toExecuteCurrentStage():
            #stageRFNA.runJobs()
            varsP.runJobs(stageRFNA, "refineA")
            stageRFNA.checkResults()
        else :
             stageRFNA.endStage()
        if varsP.exitTestMinimumContigs(0):
            return 1 #1 is failure (a-la bash)
    
 	varsP.ImportSQLite("refineA") #no characterize for refineA
        
        # Contig Refine B
        if varsP.groupContigs:
		stageRFNB = grprefmod.Refine("refineB0", varsP) #new class requires first arg string
		if varsP.toExecuteCurrentStage():
                    varsP.runJobs(stageRFNB, "refineB0")
                    stageRFNB.checkResults()
                else :
                    stageRFNB.endStage()

		stageRFNB = grprefmod.Refine("refineB1", varsP) #new class requires first arg string
		if varsP.toExecuteCurrentStage():
                    varsP.runJobs(stageRFNB, "refineB1")
	else :
		stageRFNB = refmod.Refine("refineB", varsP) #new class requires first arg string
		if varsP.toExecuteCurrentStage():
                    varsP.runJobs(stageRFNB, "refineB")
        stageRFNB.checkResults() #need to call this even if not executing bc otherwise curContigCount isn't set in mergeIntoSingleCmap, so exitTest fails when it shouldn't (sometimes)
        skiplevel = varsP.exitTestMinimumContigs(1, True) #treat refineB like a merge, because 1 contig means skip following
        if skiplevel == 0 or skiplevel == 2 :
            self.lastCharacterize = self.runCharacterize(varsP)
        if skiplevel == 1 : #0 is success: no longer exit on single contig
            return skiplevel #1 for no contigs (exit), 2 for single contig (skip to SV, alignmol)

        if varsP.runSV > 1 : #runSV: for now, anything > 1 is refineB
	    util.LogStatus("progress", "stage_start", "SVrefineB")
            stageSV = svm.SVdetect(varsP) 
            #stageSV.runJobs()
            varsP.runJobs(stageSV, "SVrefineB")
            stageSV.checkResults()
	    util.LogStatus("progress", "stage_complete", "SVrefineB")

        # NGS Contig Refine
        #because toExecuteCurrentStage is inside if ngsInDir, the current counting of stages is preserved when not using this
        if varsP.ngsInDir :
	    util.LogStatus("progress", "stage_start", "refineNGS")
            stageRFNngs = refmod.Refine("refineNGS", varsP) #new class requires first arg string
            if varsP.toExecuteCurrentStage():
                varsP.runJobs(stageRFNngs, "refineNGS")
            stageRFNngs.checkResults()
            #not sure about exitTestMinimumContigs -- may want to continue if no ngs but have denovo
            self.runCharacterize(varsP) #see how your refined contigs compare to unrefined vis-a-vis reference
	    util.LogStatus("progress", "stage_complete", "refineNGS")

	varsP.toExecuteCurrentStage()
	stageMRG = refmod.Merge(varsP, "Merge0")
	for j in range(26):
		stageMRG.generateJobList()
		if varsP.executeCurrentStage:
			varsP.runJobs(stageMRG, "merge_0")
		if stageMRG.mergeComplete():
			break
	self.lastCharacterize = self.runCharacterize(varsP)
        skiplevel = varsP.exitTestMinimumContigs(1, ismerge=True) 
        if skiplevel : #if this job failed, you get 1, if single contig get 2
            #don't return 2 bc that skips refineFinal--for all subsequent stages, return 0 or 1 only
            return (skiplevel if skiplevel == 1 else 0) 

        # Iterative Extension and Merging
        for i in range(varsP.nExtensionIter):
	    if varsP.groupContigs:
		stageEXT = grprefmod.Refine("extension0", varsP)
		if varsP.toExecuteCurrentStage():
                    #stageEXT.runJobs()
                    varsP.runJobs(stageEXT) #, "extension0_%i"%i) #arg name not used
                    stageEXT.checkResults()
		else :
                    stageEXT.endStage()
	
		stageEXT = grprefmod.Refine("extension1", varsP)
		if varsP.toExecuteCurrentStage():
                    #stageEXT.runJobs()
                    varsP.runJobs(stageEXT) #, "extension1_%i"%i) #arg name not used
                    stageEXT.checkResults("_"+str(i+1))
		else :
                    stageEXT.endStage()

	    else :
		stageEXT = refmod.Extension(varsP)
		if varsP.toExecuteCurrentStage():
                    #stageEXT.runJobs()
                    varsP.runJobs(stageEXT, "extension_%i"%i)
                    stageEXT.checkResults()
		else :
                    stageEXT.endStage()
            self.lastCharacterize = self.runCharacterize(varsP)
            skiplevel = varsP.exitTestMinimumContigs(1, ismerge=True)
            if skiplevel : #see comment at merge0
                return (skiplevel if skiplevel == 1 else 0) 
			
            varsP.toExecuteCurrentStage()
            #stageMRG = refmod.Merge(varsP) #instead of calling constructor again, use new method
            stageMRG.resetStage( "Merge%d" % (i+1) )
            for j in range(26):
                stageMRG.generateJobList()
                if varsP.executeCurrentStage:
                    #stageMRG.runJobs()
                    varsP.runJobs(stageMRG, stageMRG.stageName)
                if stageMRG.mergeComplete():
                    break
            self.lastCharacterize = self.runCharacterize(varsP)
            skiplevel = varsP.exitTestMinimumContigs(1, ismerge=True)
            if skiplevel : #see comment at merge0
                return (skiplevel if skiplevel == 1 else 0) 
    #end runAllStages


    def runRefineFinal(self, varsP):
        """Run refineFinal only.
        """

        # done with extend - merge; do final refinement
        #  if no extend-merge, previous stage was refineB, so no further refinement is necessary
        if not varsP.skipFinalRefine and varsP.nExtensionIter > 0 :
	    if varsP.groupContigs:	
		stageRFNfinal = grprefmod.Refine("refineFinal0", varsP) #new class requires first arg string
		if varsP.toExecuteCurrentStage():
                    #stageRFNfinal.runJobs()
                    varsP.runJobs(stageRFNfinal, "refineFinal0")
                    stageRFNfinal.checkResults()
		else :
                    stageRFNfinal.endStage()
		
		stageRFNfinal = grprefmod.Refine(varsP.refineFinal1, varsP) #new class requires first arg string (see varsP.__init__)
		if varsP.toExecuteCurrentStage():
                    #stageRFNfinal.runJobs()
                    varsP.runJobs(stageRFNfinal, varsP.refineFinal1)
                    stageRFNfinal.checkResults()
		else :
                    stageRFNfinal.endStage()
	    else :
		stageRFNfinal = refmod.Refine("refineFinal", varsP) #new class requires first arg string
		if varsP.toExecuteCurrentStage():
                    #stageRFNfinal.runJobs()
                    varsP.runJobs(stageRFNfinal, "refineFinal")
                    stageRFNfinal.checkResults()
		else :
                    stageRFNfinal.endStage()
            if varsP.exitTestMinimumContigs(0):
                return 1 #if 1, skip final stages
            self.lastCharacterize = self.runCharacterize(varsP)
            self.runCharacterizeFinal(varsP) #run this for both grouped and not
    #end runRefineFinal


    def runFinalStages(self, varsP, execute) :
        """Run alignmol and SV.
        Bool execute determines whether stages are run. If false, just call wait_all_threads.
        """

        #align molecules to contigs
        if execute and varsP.doAlignMol : 
            util.LogStatus("progress", "stage_start", "alignmol")
            amod = alignmod.AlignModule(varsP) #constructor will call generateJobList -- call so varsP.alignMolDir is set
            if varsP.toExecuteCurrentStage(): #alignmol can be bypassed
                varsP.runJobs(amod, "alignmol")
                amod.checkResults() #calls doAllPipeReport but _not_ getAlignStats
            util.LogStatus("progress", "stage_complete", "alignmol")

	util.LogStatus("progress", "stage_start", "AsyncWait")
	mthread.wait_all_threads() #ensure start_new_thread jobs (Characterize) finished
	util.LogStatus("progress", "stage_complete", "AsyncWait")

        #runSV: only if previous stage was not SVdetect
        if ( execute and varsP.runSV and self.lastCharacterize and
             (varsP.ngsInDir or varsP.nExtensionIter or not varsP.skipFinalRefine) ) :
            plist = self.lastCharacterize.curCharacterizeFileRoots #this is a list whose only element is output path+prefix of last job (only zero if didn't run)
            if plist :
                charnoise = scm.readNoiseParameters(plist[0]) #empty string will return {}
                stageSV = svm.SVdetect(varsP, charnoise)
                #if varsP.refDeresed : #set in SVdetect constructor (really in CharacterizeModule.referenceProcess)
                #if is not necessary: if this fails, no jobs are created, so runJobs will do nothing
                varsP.runJobs(stageSV, "SVfinal") 
                stageSV.checkResults()
                varsP.ImportSQLite("svdetect")

        util.logMemory(varsP.memoryLogpath, varsP.startTime, "align_stats") #add before and after amod.getAlignStats
        #wait until here to call getAlignStats because it needs result of last characterize
        #this relies on no calls to toExecuteCurrentStage between here and alignmol above
        if execute and varsP.doAlignMol and varsP.executeCurrentStage : 
            #Characterize (actually MapClassesRev) stores assembly size in totAssemblyLenMb -- copy this to local varsP
            if self.lastCharacterize : # None if never ran any characterize
                varsP.totAssemblyLenMb = self.lastCharacterize.varsP.totAssemblyLenMb
            amod.getAlignStats()
            
        #print(varsP.stageFolder)    
        #move the index of stageFolder into this method to check it
 	varsP.ImportSQLite("refineFinal1") #must be after async wait for last characterize
    
        util.logMemory(varsP.memoryLogpath, varsP.startTime, "align_stats") #add before and after amod.getAlignStats

        varsP.doCompression()
    #end runFinalStages(self, varsP)


    def finalizePipeline(self, varsP, ers=None):
        """Final Pipeline logging on exit."""

        elp = time.time() - varsP.startTime #elapsed time in seconds
        varsP.updatePipeReport( "  Pipeline end time: %s\n  Elapsed time: %.2fm; %.2fh; %.2fd\n\n" % (time.ctime(), elp/60., elp/3600., elp/3600./24) )

        if ers == None :
            ers = util.SummarizeErrors(varsP=varsP)
        if ers.has_key('critical') : #Pipeline did not complete
		varsP.updatePipeReport("Pipeline has failed\n") 
		util.LogStatus("progress", "pipeline", "failure")
        elif ers.has_key('error') : #different status log message in this case
		varsP.updatePipeReport("Pipeline has completed with errors\n") 
		util.LogStatus("progress", "pipeline", "complete_with_errors")
	else : #only warnings is success
		varsP.updatePipeReport("Pipeline has successfully completed\n")
		util.LogStatus("progress", "pipeline", "success")

    
    def runCharacterize(self, varsP):
        """Run the characterization module.
        """
        if not varsP.executeCurrentStage: #bypass
            dc = cm.dummyCharacterize(varsP) #just get a path for SV noise parameters
            if dc.curCharacterizeFileRoots : #if characterize failed before bypass, this will be empty, so characterize will still run
                return
            else :
                err = "running characterize for %s" % varsP.stageName
                varsP.updatePipeReport( "warning : "+err+"\n" )
                util.LogError("warning", err)

        stageCHAR = cm.Characterize(copy.deepcopy(varsP))
        mthread.start_new_thread(self.doRunCharacterize, (stageCHAR, varsP, copy.deepcopy(varsP.stageName)))
        return stageCHAR #for getting its data members (really its varsP data member)

    #same as run characterize except not required; run only if argData (xml) has section characterizeFinal
    def runCharacterizeFinal(self, varsP):
        """Run the characterization module using the flag for characterizeFinal.
        """
        charf = "characterizeFinal"
        if varsP.ref and varsP.executeCurrentStage and varsP.argData.has_key(charf) and varsP.argData[charf] :
            stageCHAR = cm.Characterize(copy.deepcopy(varsP), 1) #second arg is 1 for characterizeFinal arguments
            mthread.start_new_thread(self.doRunCharacterize, (stageCHAR, varsP, varsP.stageName))
            return stageCHAR #for getting its data members (really its varsP data member)
	return
	    
    
    def doRunCharacterize(self, stageCHAR, varsP, stageName):
        """Utility method in order to multithread the pipeline code so that
        processing can continue while characterization runs.
        """
        #stageCHAR.runJobs()
        varsP.runJobs(stageCHAR, "Characterize") #here, celery won't know what stage this is, but that's ok
        stageCHAR.checkResults()
        #if not stageName in varsP.stageFolder: #this is moved inside ImportSQLite
	#	print("Missing " + stageName)
	#	print(varsP.stageFolder)
	if not stageName in ["refineFinal1"]: #wait for refineFinal because of CharacterizeFinal
		varsP.ImportSQLite(stageName)

#end class DNPipeline

    
def writeToFile(fileName, content):
    f1 = open(fileName, 'w')
    f1.write(content) 
    f1.close()


def BinContigsBySize(ContigFiles, FileSizes, nBins):
    sortDict = {}
    for i,ContigFile in enumerate(ContigFiles):
        FileSize = FileSizes[i]
        if sortDict.has_key(FileSize):
            sortDict[FileSize].append(ContigFile)
        else:
            sortDict[FileSize] = [ContigFile]
    SortedFileSizes = sorted(sortDict.keys(), reverse = True)
    
    binnedContigFiles = []
    ct = 0
    for FileSize in SortedFileSizes:
        CurContigFiles = sortDict[FileSize]
        for CurContigFile in CurContigFiles:
            if ct/nBins == 0:
                binnedContigFiles.append([CurContigFile])
            else:
                binnedContigFiles[ct%nBins].append(CurContigFile)
            ct += 1
    return binnedContigFiles
    
