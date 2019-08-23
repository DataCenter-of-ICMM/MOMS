
description = """Standalone script for running alignmol and/or merging alignmol output files from Pipeline."""

import os, sys
import argparse
import time

def runAlignMol() :    
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-q', dest='queryDir', help='Path to merged cmap to align molecules (-b) to OR alignmol dir from Pipeline for merge (if latter, no alignments are performed), required', type=str)
    parser.add_argument('-b', dest='bnx', help='Input molecule (.bnx) file, required if aligning molecules', type=str)
    #parser.add_argument('-b', dest='bnx', help='Input molecule (.bnx) file OR path to dir containing split bnx pieces, required if aligning molecules', type=str) #I should add the split feature; for now, just do single bnx
    parser.add_argument('-a', dest='optArguments', help='Path to optArguments.xml (optional, default optArguments_human.xml in Pipeline dir if found, otherwise required)', default="", type=str)
    parser.add_argument('-r', help='If this flag is used, alignmolvref arguments are used, otherwise alignmol arguments are used (default alignmol; optional)', dest='ref', action='store_true')
    parser.add_argument('-o', dest='outputDir', help='output dir (optional, defaults to sub-dir of input map dir called "alignmol")', default="", type=str)
    parser.add_argument('-t', dest='RefAligner', help='Path to RefAligner or dir containing it (required)', type=str) 
    parser.add_argument('-T', dest='numThreads', help='Total number of threads (cores) to use (optional, default 4)', default=4, type=int)
    parser.add_argument('-j', dest='maxthreads', help='Threads per Job, -maxthreads (non-cluster only;optional, default 4)', default=4, type=int)
    parser.add_argument('-e', dest='errFile', help='.err file to use for noise parameters--will supersede noise parameters in the optArgument supplied (but that file must still be supplied for non-noise parameters)--should be from autoNoise', default="", type=str)
    parser.add_argument('-E', dest='errbinFile', help='.errbin file to use for noise parameters--will supersede noise parameters in the optArgument supplied (but that file must still be supplied for non-noise parameters)--should be from autoNoise', default="", type=str)
    parser.add_argument('-F', help='Color channel: replace -usecolor X in optArgs with this, must be either 1 or 2 [default OFF]', default=0, type=int)
    parser.add_argument('-p', dest='pipelineDir', help='Pipeline dir (optional, defaults to script dir, or current directory)', default="", type=str)
    parser.add_argument('-cd', help='Control data directory for copy number profiles [optional]', dest='controldir', type=str, default="")
    parser.add_argument('-pd', help='Parameters directory for copy number profiles [optional]', dest='paramdir', type=str, default="")
    parser.add_argument('-op', help='Outlier Probability for copy number profile (optional)', dest='outlierp', default=0.001, type=float) # 12202017 lzhang 0.95-->0.001
    result = parser.parse_args()

    #outprefix = "exp_refineFinal1" #this is the default; assume for now -- not used

    #check all Pipeline dependencies
    if result.pipelineDir :
        cwd = result.pipelineDir
    else :
        cwd = os.path.split(os.path.realpath(__file__))[0] #this is path of this script
        if not os.path.isfile(os.path.join(cwd,"utilities.py")) : #if still not here, last try is actual cwd
            cwd = os.getcwd() #still check this below

    #this is the only one imported here and in runCharacterize
    if not os.path.isfile(os.path.join(cwd,"utilities.py")):
        print "ERROR: utilities.py missing in dir", cwd, "check -p argument, or run this script in Pipeline dir"
        sys.exit(1)
    import utilities as util

    if not os.path.isfile(os.path.join(cwd,"AlignModule.py")):
        print "ERROR: AlignModule.py missing in dir", cwd, "check -p argument, or run this script in Pipeline dir"
        sys.exit(1)
    import AlignModule as alignmod

    if not util.checkFile(os.path.join(cwd,"Pipeline.py")):
        print "ERROR: Pipeline.py missing in dir", cwd, "check -p argument, or run this script in Pipeline dir"
        sys.exit(1)
    import Pipeline

    if not util.checkFile(os.path.join(cwd,"mapClasses.py")):
        print "ERROR: mapClasses.py missing in dir", cwd, "check -p argument, or run this script in Pipeline dir"
        sys.exit(1)
    import mapClasses as mc

    #input dir
    if not result.queryDir :
        print "ERROR: Query (-q) argument not supplied."
        sys.exit(1)
    qrypath = os.path.realpath(result.queryDir)
    print "qrypath= ", qrypath # DEBUG

    if util.checkDir(qrypath, checkWritable=False, makeIfNotExist=False) : #output elsewhere so not writeable is ok
        runaligns = False
    elif util.checkCmap(qrypath) :
        runaligns = True
    else :
        print "ERROR: Query argument ("+qrypath+") not found or not a dir or cmap. Check -q argument."
        sys.exit(1)

    #this check isn't really necessary...make it a warning -- left over from runAlignMerge.py
    #if not os.path.split(qrypath)[1].endswith("alignmol") :
    #    print "Warning: Query dir ("+qrypath+") does not end with 'alignmol'; please be sure this is a Pipeline alignmol dir\n"

    #RefAligner -- check for either path to RefAligner, or dir containing it, depending on cluster args
    rabin = "" #need empty string for generateJobList even though no jobs are run
    if runaligns :
        rabin = result.RefAligner
        #replicate Pipeline behavior: RefAligner is always required
        if os.path.isdir(rabin) :
            rabin = os.path.join(rabin, "RefAligner")
        if not util.checkExecutable(rabin):
            print "ERROR: RefAligner not found or not executable at", rabin, "\nPlease supply RefAligner dir or full path as -t arg."
            sys.exit(1)

    #optargs file
    optargs = None
    if runaligns and result.optArguments : #supplied on command line
        optargs = result.optArguments
        if not util.checkFile(optargs, ".xml") :
            print "optArguments path is supplied ("+optargs+") but not found or doesn't end in .xml, check -a argument."
            sys.exit(1)
    elif runaligns : #load from Pipeline dir if running alignments
        optargs = os.path.join(cwd,"optArguments_human.xml")
        if not util.checkFile(optargs):
            print "optArguments.xml missing in Pipeline directory ("+cwd+"). Try supplying path explicitly using -a."
            sys.exit(1)

    #output dir
    if not result.outputDir :
        outdir = os.path.join(qrypath, "merge") #should be same as in AlignModule
    else :
        outdir = os.path.realpath(result.outputDir)
    if os.path.isdir(outdir) :
        if not util.checkDir(outdir) : #check writeable
            print "\nERROR: Output dir is not writeable:\n", outdir, "\n"                
            sys.exit(1)
        #this is ok here
        #elif outdir == contigdir :
        #    print "\nERROR: Output dir cannot be same as input dir:\n", outdir, "\n"                
        #    sys.exit(1)                
        print "\nWARNING: Output dir already exists, results will be overwritten:\n", outdir, "\n"
    elif not util.checkDir(outdir) : #does not exist, make, if False, can't make or not writeable
        print "\nERROR: Output dir cannot be created or is not writeable:\n", outdir, "\n"
        sys.exit(1)
    
    #bnx file
    bnxfile = result.bnx
    if bnxfile : #must check for empty string BEFORE you do realpath, or it returns cwd
        bnxfile = os.path.realpath(bnxfile)
        if not util.checkFile(bnxfile, ".bnx") :
            print "ERROR: bnx file supplied but not found or incorrect suffix:", bnxfile
            sys.exit(1)
    elif runaligns :
        print "ERROR: bnx file not supplied but running alignments; please supply bnx file as -b argument"
        sys.exit(1)

    #nthreads
    nthreads = result.numThreads
    if nthreads <= 0 :
        print "ERROR: Number of threads value invalid (must be > 0): %i" % nthreads
        sys.exit(1)

    #maxthreads
    maxthreads = result.maxthreads
    if maxthreads <= 0 :
        print "ERROR: Max threads value invalid (must be > 0): %i" % maxthreads
        sys.exit(1)
    elif nthreads < maxthreads :
        print "Warning: num threads (-T: %i) < max threads (-j: %i): increasing num threads to equal max threads\n" % (nthreads, maxthreads)
        nthreads = maxthreads

    #.errbin file
    errbinfile = result.errbinFile
    if errbinfile :
        errbinfile = os.path.realpath(result.errbinFile)
        if not util.checkFile(errbinfile, ".errbin") :
            print "ERROR: errbin file supplied but not found or incorrect suffix:", errbinfile
            sys.exit(1)

    #.err file
    errfile = result.errFile
    if errfile and errbinfile :
        print "Warning: .err and .errbin arguments supplied; ignoring .err file"
        errfile = ""
    elif errfile :
        errfile = os.path.realpath(result.errFile)
        if not util.checkFile(errfile, ".err") :
            print "err file supplied but not found or incorrect suffix:", errfile
            sys.exit(1)

    if errfile and not util.checkFile(os.path.join(cwd,"SampleCharModule.py")):
        print "SampleCharModule.py missing in dir", cwd, "check -p argument, or run this script in Pipeline dir"
        sys.exit(1)
    elif errfile :
        import SampleCharModule as scm

    if result.F < 0 or result.F > 2 :
        print "Error: argument -F must be >= 0 AND <= 2. (%i)" % result.F
        sys.exit(1)

    doref = result.ref
    if result.controldir and not util.checkDir(result.controldir, checkWritable=False, makeIfNotExist=False) :
        print '  ERROR: control dir for CNV not found: %s\n' % (result.controldir)
        sys.exit(1)
    cnv_controldir = result.controldir
    if result.paramdir and not util.checkDir(result.paramdir, checkWritable=False, makeIfNotExist=False) :
        print '  ERROR: param dir for CNV not found: %s\n' % (result.paramdir)
        sys.exit(1)
    cnv_paramdir = result.paramdir
    if result.outlierp < 0 or result.outlierp > 1. :
        print '  ERROR: outlier probability for CNV must be between 0 and 1: %f\n' % (result.outlierp)
        sys.exit(1)
    cnv_outlierp = result.outlierp

    #DONE checking arguments

    print "Using output dir", outdir
    if runaligns :
        print "Aligning", bnxfile, "\nTo", qrypath, "\n"
    else :
        print "Merging", qrypath, "\n"

    startTime = time.time() #time since Epoch
    memory_log = os.path.join(outdir, "memory_log.txt")
    util.initMemoryLog(memory_log)

    varsP = Pipeline.varsPipeline()
    varsP.RefAlignerBin        = rabin
    varsP.contigFolder         = "" #not used but needs to be an attr
    varsP.outputContigFolder   = "" #not used but needs to be a string attr
    varsP.pipeReportFile = os.path.join(outdir, "alignmol_jobs_log.txt")
    varsP.infoReportFile = os.path.join(outdir, "alignmol_log.txt")
    varsP.infoReportFileSimple = os.path.join(outdir, "alignmol_log_simple.txt")
    varsP.doAlignMolvRef = False #only used in AlignModule.getAlignStats for whether to skip stats: False for no skip
    varsP.usecolorArg = result.F
    varsP.expID = "exp" #this is default for normal Pipeline; hardcode here
    varsP.cnv_controldir = cnv_controldir
    varsP.cnv_paramdir = cnv_paramdir
    varsP.cnv_outlierp = cnv_outlierp
    util.InitStatus( os.path.join(outdir, "status.xml") )
    varsP.checkDependencies()

    if runaligns :
        varsP.optArgumentsFileIn   = optargs
        varsP.latestMergedCmap     = qrypath #if !doref, need this one
        varsP.ref                  = qrypath #and if doref, need this one
        varsP.nThreads             = nthreads #necessary otherwise job won't start -- max threads per node
        varsP.maxthreads           = maxthreads #threads per job
        p = os.path.split(qrypath)[1]
        varsP.outputContigPrefix   = p[:p.rfind(".")] #filename prefix
        varsP.stdoutlog    = True #use -stdout -stderr
        varsP.sorted_file = bnxfile[:bnxfile.rfind(".")] #enables the mol fraction align in AlignModule.getAlignStats
        if qrypath.endswith(".cmap") : #enable the mol stats
            varsP.totAssemblyLenMb = mc.multiCmap(qrypath, lengthonly=True).totalLength / 1e6

        varsP.memoryLogpath  = os.path.join(outdir, "memory_log.txt")
        varsP.parseArguments() #parses optArgumentsFile
        #varsP.checkDependencies() #move up for Rscript_path (for CNV)
        varsP.RefAlignerBinOrig = rabin
        varsP.updatePipeReport( ' '.join(sys.argv) + '\n\n' , printalso=False)
        varsP.prerunLog() #general information in log -- needed for refaligner_version

        noisep = {}
        if errbinfile :
            noisep = {"readparameters": errbinfile}
            #print "Using noise parameters from "+errbinfile+"\n" #move below
        elif errfile :
            noisep = scm.readNoiseParameters(errfile.replace(".err",""))
            if noisep.has_key('readparameters') : #remove this because it's redundant, and it can cause problems with RefAligner compatibility
                del noisep['readparameters']
            if not noisep : #readNoiseParameters returns empty dict on failure
                print "ERROR reading noise parameters, check .err file:", errfile
                sys.exit(1)
            #redundant with below?
            print "Using noise parameters from "+errfile+":\n" + " ".join(["-"+str(k)+" "+str(v) for k,v in noisep.iteritems()])+"\n"

        #some code from SampleCharModule to load args into noise0
        infoReport="Loaded noise parameters:\n"
        klist = ["FP", "FN", "sf", "sd", "sr", "bpp", "readparameters"] #hardcoding parameters is kind of bad, but it fixes the order without using OrderedDict.
        #noiseargs = self.varsP.argsListed('noise0') #not necessary
        for v in klist :
            if not noisep.has_key(v) :
                continue
            param=str(noisep[v])
            util.LogStatus("parameter", "auto_"+v, param)
            infoReport+=v+":"+param+"\n"
            varsP.replaceParam("noise0", "-"+v, param)
        varsP.updateInfoReport(infoReport + '\n', printalso=True)

    else :
        print "Getting file list from", qrypath
        outFileList = getOutFileList(util, qrypath)
        if not outFileList :
            print "ERROR: Query dir ("+qrypath+") does not contain alignmol data. Check -q argument."
            sys.exit(1)
        else :
            print "Found", len(outFileList), "alignment results"
    #end if runaligns

    amod = alignmod.AlignModule(varsP, doref, outdir, bnxfile) #constructor will call generateJobList

    if runaligns :
        amod.runJobs()
	amod.checkResults()
    else :
        amod.outFileList = outFileList
        p = os.path.split(outFileList[0])[1]
        if p.count("_") > 1 : #expect something like "EXP_REFINEFINAL1_4"
            #p = p[:p.rfind("_")+1] #remove integer suffix
            p = p[:p.rfind("_")] #remove integer suffix (and underscore)
        #else :
        #    p += "_" #because mrgstr is appended
        varsP.outputContigPrefix = p

    if not runaligns or len(amod.jobList) > 0 :
        amod.getAlignStats()

    if runaligns :
        print
        util.SummarizeErrors(varsP=varsP) #below changed, just call the fn and ignore return
        #copy from Pipeline.py
        #if util.SummarizeErrors(varsP=varsP)==0:
        #    varsP.updatePipeReport("Pipeline has successfully completed\n") 
        #    util.LogStatus("progress", "pipeline", "success")
        #else:
        #    varsP.updatePipeReport("Pipeline has completed with errors\n") 
        #    util.LogStatus("progress", "pipeline", "failure")

    #BELOW OLD CODE

    return

    #in Pipeline, this is called first
    #print "Calling getAlignStats:" #but it won't work without varsP atm; skip it
    #getAlignStats(self.varsP, self.outFileList, self.varsP.totAssemblyLenMb, isref=False, mergepath=self.mergedir)
    #getAlignStats(self.varsP, self.outFileList, self.varsP.totAssemblyLenMb, isref=False, mergepath=self.mergedir)

    print "Calling mergeMap"
    print outFileList[0] #, "\n", outputdir #moved above
    util.logMemory(memory_log, startTime, "mergeMap_start")
    #mergeMap(self.varsP, self.outFileList, mergepath=self.outputdir) #varsP is optional
    alignmod.mergeMap(None, outFileList, outputdir) 
    util.logMemory(memory_log, startTime, "mergeMap_end")

    print "Calling mergeRcmaps"
    util.logMemory(memory_log, startTime, "mergeRcmaps_start")
    #mergeRcmaps(outFileList, outdir, varsP=None, splitByContig=None, stageName="alignmol") :
    alignmod.mergeRcmaps(outFileList, outputdir, splitByContig=True, stageName=outprefix) 
    util.logMemory(memory_log, startTime, "mergeRcmaps_end")

    print "Calling split_XMap_byContig" #split_XMapQcmap_byContig"
    util.logMemory(memory_log, startTime, "split_XMap_byContig_start")
    #xmapdict = alignmod.split_XMap_byContig(outFileList, outputdir, stageName=outprefix) #old
    xmapdict = alignmod.split_XMap_byContig_new(outFileList, outputdir, stageName=outprefix)
    util.logMemory(memory_log, startTime, "split_XMap_byContig_end")

    print "Calling split_Qcmap_byContig" 
    util.logMemory(memory_log, startTime, "split_Qcmap_byContig_start")
    #alignmod.split_Qcmap_byContig(outFileList, outputdir, xmapdict) #old
    alignmod.split_Qcmap_byContig_new(outFileList, outputdir, xmapdict, stageName=outprefix) #new: better performance
    util.logMemory(memory_log, startTime, "split_Qcmap_byContig_end")

    print "AlignMerge successfully completed"
#end getArgs


def getOutFileList(util, inpath) :
    olist = util.getListOfFilesFromDir(inpath, "xmap") #get all xmaps
    return [x.replace(".xmap","") for x in olist]


if __name__ == "__main__" :
    args = runAlignMol()

