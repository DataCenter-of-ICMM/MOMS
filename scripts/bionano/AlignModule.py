import os
import shutil #for rmtree
import Multithreading as mthread
import utilities as util
import mapClasses as mc

"""
@package AlignModule
Align molecules to final contigs, typically refine Final contigs.
"""

util.setVersion("$Id: AlignModule.py 6374 2017-05-16 17:13:18Z wandrews $")


class AlignModule(mthread.jobWrapper):
    """Class for alignment of molecules to final contigs.
    """

    def __init__(self, varsP, doref=False, outputdir=None, bnxin=None):
        """doref determines parameter set from optargs.
        outputdir not needed for Pipeline, but used in runAlignMol.py.
        If bnxin supplied, will run single job with it.
        """
        self.varsP = varsP
        self.doref = doref
        self.singleJob = not doref #if this is 'True', then submit single job: enable by default for alignmol only
        self.bnxin = bnxin #see generateJobList

        self.argStageName = 'alignmol' #use arguments from alignmol (optArgs, not clusterArgs)
        if not doref :
            self.stageName = 'alignmol' #also name of dir which is sub-dir of varsP.outputContigFolder
            self.alignTarget = os.path.join(varsP.outputContigFolder, self.stageName) #output dir
            self.varsP.alignMolDir = self.alignTarget #store in varsP for subsequent processing
        else :
            self.stageName = self.varsP.alignMolvrefName #also name of dir which is sub-dir of localRoot
            self.alignTarget = os.path.join(self.varsP.contigFolder, self.stageName) #output dir
        if outputdir :
            self.alignTarget = outputdir

        util.checkDir(self.alignTarget) #will make if doesn't exist
        self.mergedir = os.path.join(self.alignTarget, self.varsP.alignMolvrefMergeName) #copy from AlignRefModule

        ca = self.varsP.getClusterArgs(self.stageName if not self.singleJob else "autoNoise0")
        super(AlignModule, self).__init__(self.varsP, self.stageName, clusterArgs=ca) 

        self.outFileList = []
        self.generateJobList()
        self.logArguments()


    def generateJobList(self):
        """AlignModule.generateJobList: create RefAligner jobs for aligning molecules to contigs.
        """
        #for runAlignMol, this method is called but not used: exit if RefAlignerBin is empty
        if not self.varsP.RefAlignerBin :
            return

        #the contigs are obtained from varsP.latestMergedCmap--check its validity, a return will mean no jobs, and no jobs is now handled in multiThreadRunJobs.
        if not self.doref and ( not self.varsP.latestMergedCmap or
                                not util.checkCmap(self.varsP.latestMergedCmap) ) :
            err = "Error in AlignModule.generateJobList: varsP.latestMergedCmap is not set or not valid cmap; skipping %s" % self.stageName
            self.varsP.updatePipeReport(err+"\n")
            util.LogError("error", err)
            return

        #Note: noise parameters should be fixed becuase when bnx is split, -M
        # would find different parameters for different contigs. Use noise0.

        baseargs = [self.varsP.RefAlignerBin]
        if not self.doref :
            baseargs += ['-ref', self.varsP.latestMergedCmap] #reference is latest merged cmap
            mappref = os.path.split(self.varsP.latestMergedCmap)[1]
            mappref = mappref[:mappref.find(".")]
        else :
            baseargs += ['-ref', self.varsP.ref] 
            mappref = self.stageName #use stageName also for output filename

        noiseargs = self.varsP.argsListed('noise0')
        haverefargs = False
        try : #argsListed does not check key
            refargs = self.varsP.argsListed(self.stageName) #'alignmolvref'
            haverefargs = True
        except KeyError : #this is same as old behavior
            #refargs = self.varsP.argsListed('noise0') + self.varsP.argsListed(self.argStageName) #old
            refargs = self.varsP.argsListed(self.argStageName) #new
        #refargs = noiseargs + refargs

        if haverefargs :
            self.jobargs = refargs

        #single job with bnxin (constructor)
        if self.bnxin or self.singleJob :
            #the 'mapped-unsplit' option below means no further file manipulation is necessary: result goes in merge dir
            util.checkDir(self.mergedir)
            outarg = os.path.join(self.mergedir, self.varsP.expID + "_" + self.varsP.refineFinal1) # ie, 'exp_refineFinal1'
            self.outFileList.append( outarg ) #file prefixes
            jobargs = baseargs + ['-o', outarg]
            if self.bnxin :
                jobargs += ['-i', self.bnxin]
            else :
                jobargs += ['-i', self.varsP.sorted_file+".bnx"] #the rescaled bnx (but missing suffix)

            stdoutf = None
            if self.varsP.stdoutlog : #remember, these must be after -o
                jobargs.extend( ['-f', '-stdout', '-stderr'] )
                stdoutf = outarg+".stdout"
            jobargs += ['-maxthreads', str(self.varsP.nThreads)] # NOTE : always use all threads for alignmol jobs
            #add noise0 before alignmol (stageName) so that the latter can override the former
            jobargs += noiseargs
            vetostr = '(intervals.txt|\.map)$'
            if self.varsP.onCluster : #if cluster, need to enclose in quotes (same line in SVModule.py)
                vetostr = "\'"+vetostr+"\'"
            jobargs.extend( ['-output-veto-filter', vetostr] )
            if self.singleJob : #new RefAligner feature for alignmol: produce separate xmap/qcmap for each contig (ref)
                jobargs.extend( ['-mapped-unsplit', '0'] )
            #fix maxmem: the optargs will be overridden in this module only, for group == 0 only, if largemem section exists
            maxm = '-maxmem'
            if maxm in refargs and self.varsP.argData.has_key(self.varsP.largemem):
                largs = self.varsP.argsListed(self.varsP.largemem) #should be just ['-maxmem', '128']
                if len(largs) >= 2 and util.getIntFromString(largs[1]) != None :
                    self.varsP.updatePipeReport("AlignModule: changing -maxmem from "+refargs[refargs.index(maxm)+1]+" to "+largs[1]+" due to host job\n")
                    refargs[refargs.index(maxm)+1] = largs[1]
            jobargs += refargs

            s1Job = mthread.singleJob(jobargs, self.stageName, outarg+".err", self.stageName, maxThreads=self.varsP.nThreads, clusterLogDir=self.varsP.clusterLogDir, expectedStdoutFile=stdoutf)
            self.addJob(s1Job)
            return #and this is the only job

        #loop over the split bnxs, make one job per bnx
        for idx in range(1,self.varsP.nPairwiseJobs+1) :

            outarg = os.path.join(self.alignTarget, mappref+"_"+str(idx))
            self.outFileList.append( outarg ) #file prefixes
            jobargs = baseargs + ['-o', outarg]
            idxstr = "_%s_of_%s" % (idx, self.varsP.nPairwiseJobs)
            jobargs += ['-i', self.varsP.bnxFile.replace(".bnx", idxstr+".bnx")]

            stdoutf = None
            if self.varsP.stdoutlog : #remember, these must be after -o
                jobargs.extend( ['-f', '-stdout', '-stderr'] )
                stdoutf = outarg+".stdout"
            jobargs += ['-maxthreads', str(self.varsP.nThreads)]  # NOTE : always use all threads for alignmol jobs
            #add noise0 before alignmol (stageName) so that the latter can override the former
            jobargs += noiseargs
            #if idx != 1 : #keep _r for first job only -- copied from SVModule
            #    jobargs.extend( ['-output-veto-filter', '_r.cmap$'] ) #need this for copy number; do NOT veto
            jobargs.extend( ['-output-veto-filter', 'intervals.txt$'] ) #this feature not in old RefAligner
            jobargs += refargs

            s1Job = mthread.singleJob(jobargs, self.stageName+idxstr, outarg+".xmap", self.stageName+idxstr, maxThreads=self.varsP.nThreads, clusterLogDir=self.varsP.clusterLogDir, expectedStdoutFile=stdoutf)
            self.addJob(s1Job)

    #end generateJobList


    def runJobs(self):
        self.multiThreadRunJobs(self.varsP.nThreads)


    def checkResults(self):
        self.doAllPipeReport()


    def checkRequiredCopyNumberArgs(self):
        """Arguments required for copy number are -minlen 150 and -T 9.
        Note that there are many others that are required, but not checked.
        Return True for correct values, False for incorrect
        """
        minlenarg = '-minlen' #the RefAligner arg
        minlenval = '150' #the required value
        pvarg     = '-T' #the RefAligner arg
        pvval     = '1e-9' #the required value--assume the format is correct also...
        if not self.jobargs : #this means the 'alignmolvref' section is not in optargs
            err = "Warning in AlignModule.checkRequiredCopyNumberArgs: missing arguments for %s" % self.stageName
            util.LogError("warning", err)
            self.varsP.updatePipeReport(err+"\n")
            return False
        else : #have it, make sure it has -T 1e-9
            if not pvarg in self.jobargs :
                err = "Warning in AlignModule.checkRequiredCopyNumberArgs: missing argument %s in arg group %s" % (pvarg, self.stageName)
                util.LogError("warning", err)
                self.varsP.updatePipeReport(err+"\n")
                return False
            pv = self.jobargs[self.jobargs.index(pvarg)+1]
            if pv != pvval :
                err = "Warning in AlignModule.checkRequiredCopyNumberArgs: argument %s in arg group %s has value %s but required value is %s" % (pvarg, self.stageName, pv, pvval)
                util.LogError("warning", err)
                self.varsP.updatePipeReport(err+"\n")
                return False

        sortgroup = 'bnx_sort'
        sortargs = self.varsP.argsListed(sortgroup) #this is a required stage; make sure it has -minlen 150
        if not minlenarg in sortargs :
            err = "Warning in AlignModule.checkRequiredCopyNumberArgs: missing argument %s in arg group %s" % (minlenarg, sortgroup)
            util.LogError("warning", err)
            self.varsP.updatePipeReport(err+"\n")
            return False
        minlen = sortargs[sortargs.index(minlenarg)+1]
        if minlen != minlenval :
            err = "Warning in AlignModule.checkRequiredCopyNumberArgs: argument %s in arg group %s has value %s but required value is %s" % (minlenarg, sortgroup, minlen, minlenval)
            util.LogError("warning", err)
            self.varsP.updatePipeReport(err+"\n")
            return False

        return True
    #end checkRequiredCopyNumberArgs


    def getAlignStats(self):
        """Open output files of alignment jobs and report on statistics.
        """
        #MapClassesRev stores totAssemblyLenMb
        self.varsP.updatePipeReport("Starting AlignModule Align Stats stage for %s\n" % self.stageName, printalso=True)
	util.LogStatus("progress", "stage_start", "%s_stats" % self.stageName)
        if self.doref :
            reflen = mc.multiCmap(self.varsP.ref, lengthonly=True).totalLength / 1e6
        else :
            reflen = self.varsP.totAssemblyLenMb
        #in Pipeline, bnx stats are redundant on second call (which is alignmol) if reference alignments 
        skip = self.varsP.doAlignMolvRef and not self.doref
        #if util.checkDir(self.mergedir, makeIfNotExist=False) :
        #    self.varsP.updatePipeReport("AlignModule merge dir exists; removing to regnerate files: %s\n" % self.mergedir, printalso=True)
        #    shutil.rmtree(self.mergedir)
        getAlignStats(self.varsP, self.outFileList, reflen, isref=self.doref, mergepath=self.mergedir, skipbnx=skip)
        #single RefAligner job with new arg means none of this is necessary
        if not (self.bnxin or self.singleJob) :
            mergeMap(self.varsP, self.outFileList, mergepath=self.mergedir) 
            splitByContig = (2 if self.doref else 0) #see mergeRcmaps
            stageName = (self.varsP.alignMolvrefName if self.doref else "")
            mergeRcmaps(self.outFileList, self.mergedir, self.varsP, splitByContig, stageName)
            xmapDict = split_XMap_byContig_new( self.outFileList, self.mergedir, self.varsP, stageName)
            split_Qcmap_byContig_new(self.outFileList, self.mergedir, xmapDict, self.varsP, stageName)
        self.varsP.updatePipeReport("Finished AlignModule Align Stats stage for %s\n\n" % self.stageName, printalso=True)
	util.LogStatus("progress", "stage_complete", "%s_stats" % self.stageName)
    #end getAlignStats

#end class AlignModule



def getAlignStats(varsP, outFileList, reflen=0, isref=False, mergepath="", bnxpath=None, skipbnx=False) :
    '''Standalone fn for alignment statistics for both AlignModule and AlignRefModule.
    reflen should be in Mb. If mergepath supplied, put merged .err there.
    If bnxpath == None, assume varsP.sorted_file; otherwise, just report stats of this
    file and ignore outFileList.
    '''

    statonly = False #bnx stats only
    #skipbnx = False #.err file processing only, but still report raw covrage : move to arg
    simple = True #arg to updateInfoReport: put all in this file (but noise parameters are always False)
    if bnxpath == None :
        if not varsP.sorted_file : #for runAlignMol, this is empty: nothing to do in this case
            skipbnx = True
        else :
            bnxpath = varsP.sorted_file+".bnx" #set in PairwiseModule.sort_BNX even if bypassed, but needs suffix
    else : #if bnxpath != None :
        statonly = True
    if not skipbnx and not util.checkFile(bnxpath) :
        varsP.updatePipeReport("Warning in AlignModule.getAlignStats: bnxpath supplied but not found: %s\n" % bnxpath)
        return

    #find the minlen used for bnx_sort, which is a required arg set
    sortargs = []
    if varsP.argData.has_key('bnx_sort') : #for runAlignMol.py
        sortargs = varsP.argsListed('bnx_sort')
    minlen = 0
    validminlen = False
    if "-minlen" in sortargs :
        minlen = sortargs[sortargs.index("-minlen")+1] #next ele should be the len, if next ele isn't in list, the sort job will fail
        minlen = util.getIntFromString(minlen) #returns None if can't cast to int
        if minlen :
            validminlen = True

    if not validminlen and bnxpath == None and sortargs :
        varsP.updatePipeReport("Warning in AlignModule.getAlignStats: unable to obtain minlen from bnx_sort arguments; defaulting to 0\n")
    if bnxpath != None : #if bnxpath, ignore minlen
        minlen = 0

    mapn = ("reference" if isref else "assembly")
    nmol = 0 #total n mol above minlen
    totlen = 0 #total mol len above minlen
    if util.checkFile(bnxpath) :
        #the bnxfile class is very wasteful. replace with below
        #bnx = util.bnxfile(bnxpath, [minlen]) #second arg are minlen thresholds, just use one for now
        head = "Molecule Stats (%s):\n" % bnxpath
        outstr = head
        moldict = util.simpleBnxStats(bnxpath, minlen)
        nmol = moldict["nmol"]
        totlen = moldict["totlen"]
        #if isref : #this is the same for isref or not, but just print twice bc no easy way to tell if was printed previously
        outstr +=  "Total number of molecules: %6i\n" % nmol
        outstr += ("Total length (Mbp)       : %10.3f\n") % totlen
        outstr += ("Average length (kbp)     : %10.3f\n") % moldict["avglen"]
        outstr += ("Molecule N50 (kbp)       : %10.3f\n") % moldict["n50"]
        outstr += ("Label density (/100kb)   : %10.3f\n") % moldict["labdensity"]
        #    if reflen : #disable the "Genome Cov" line bc its redundant with Ref Cov below
        #        bnx.molstats[minlen].genomesizemb = 0 
        #    outstr += str(bnx.molstats[minlen]) 
        #nmol = bnx.molstats[minlen].nmol
        #totlen = bnx.molstats[minlen].totlen

        if reflen : 
            cov = totlen / reflen #totlen is in Mb
            outcov = ("Raw coverage of %s (X): %10.3f\n") % (mapn, cov)
        if skipbnx : 
            varsP.updateInfoReport(head + outcov + "\n", printalso=True, simple=simple)
        elif isref or reflen or statonly : #if neither, nothing to print
            if reflen :
                outstr += outcov
            varsP.updateInfoReport(outstr + "\n", printalso=True, simple=simple)
    elif not skipbnx :
        varsP.updatePipeReport("Warning in AlignModule.getAlignStats: missing bnx path:"+bnxpath+"\n")

    if statonly or outFileList == None : #no align files to process
        return

    #lastly, load .xmaps and .errs from alignmol jobs and report on stats
    totmaplen = 0 #sum of lengths of mapped portions of all molecules, on reference
    totmapqrylen = 0 #sum of lengths of mapped portions of all molecules, on query
    totconf = 0 #sum of confidence of all alignments
    nalign = 0 #total number of alignments
    fplist = [] #lists for error rates
    fprlist = []
    fnlist = []
    bpplist = []
    nmaplist = [] #from .err
    gmaplist = [] #from .err
    llrmlist  = []; llrgmlist = []; bppsdlist = []
    sflist = []; sdlist = []; srlist = []; reslist = []; resdlist = []
    header = ""
    err = None #will be the alignParams object if any .err files are found
    mappref = ""
    mergeerr = True #disable for single job (.err already exists in that case)
    if len(outFileList) > 0 :
        mappref = getMergeFilename(outFileList[0]) #make function to unify with same convention in mergeMap
    if mergepath and mergepath != os.path.dirname(outFileList[0]) : #debug
        print "Warning in getAlignStats: mergepath =", mergepath, "outFile dir =", os.path.dirname(outFileList[0])
    if len(outFileList) == 1 and not util.checkFile(outFileList[0]) : #this means it is a prefix: get all the files
        outFileList = util.getListOfFilesFromDir(mergepath,".xmap")
        outFileList = map(lambda x: x.rstrip(".xmap"), outFileList)
    for outpath in outFileList : #these are file prefixes
        if util.checkFile(outpath+".xmap") :
            xmap = mc.xmap(outpath+".xmap")
            nalign += len(xmap.xmapLookup)
            totmaplen += xmap.getSumMappedRefLen() #in kb
            totmapqrylen += xmap.getSumMappedQryLen() #in kb
            totconf += sum([x.Confidence for x in xmap.xmapLookup.values()])
        else :
            varsP.updatePipeReport("Warning in AlignModule.getAlignStats: missing xmap:"+outpath+".xmap"+"\n")
        if not util.checkFile(outpath+".err") and outpath == outFileList[0] and mappref :
            test = os.path.join(mergepath,mappref[:-1])+".err" #getMergeFilename leaves trailing "_"
            if util.checkFile(test) :
                outpath = os.path.join(mergepath,mappref[:-1])
                mergeerr = False
            else : #debug
                print "Warning in getAlignStats: missing .err file:", test
        if util.checkFile(outpath+".err") :
            err = mc.alignParams(outpath+".err")
            if not header :
                header = err.header
            fplist.append(err.fp)
            fprlist.append(err.fprate)
            fnlist.append(err.fn)
            bpplist.append(err.bpp)
            reslist.append(err.res)
            nmaplist.append(err.nmaps)
            gmaplist.append(err.goodmaps)
            llrmlist.append(err.llrm)
            llrgmlist.append(err.llrgm)
            bppsdlist.append(err.bppsd)
            sflist.append(err.sf)
            sdlist.append(err.sd)
            srlist.append(err.sr)
            resdlist.append(err.ressd)
        #else : #debug
        #    print "Warning in getAlignStats: missing .err file:", outpath+".err"

    #nalign from xmap should be the same as goodmaps from .err
    sumgoodmaps = sum(gmaplist)
    #if sumgoodmaps != nalign : #suppress this warning bc single job does not have separate .err files
    #varsP.updateInfoReport("Warning in getAlignStats: n mol align differ in .err files (%i) and .xmaps (%i)\n" % (sumgoodmaps, nalign), printalso=True)
    if totmaplen or totconf or nalign :
        outstr =  "Molecules aligned to %s:\n" % mapn
        outstr += "Total number of aligned molecules  : %9i\n" % nalign
        outstr += "Fraction of aligned molecules      : %13.3f\n" % (float(nalign)/nmol if nmol else 0)
        outstr += "Total molecule align length (Mbp)  : %11.1f\n" % (totmapqrylen / 1e3) #Mb
        outstr += "Total %9s align length (Mbp) : %11.1f\n" % (mapn, (totmaplen / 1e3)) #Mb
        if reflen > 0 : 
            outstr += ("Effective coverage of %9s (X): %13.3f\n") % (mapn, (totmaplen / 1e3 / reflen)) #totlen is in kb
        outstr += "Average aligned length (kbp)       : %11.1f\n" % (totmapqrylen/nalign if nalign else 0)
        outstr += "Fraction aligned length            : %13.3f\n" % (totmapqrylen/1e3/totlen if totlen else 0) #totmapqrylen is in kb, totlen is in mb
        #outstr += "Tot confidence    : %11.1f\n" % totconf
        outstr += "Average confidence                 : %11.1f\n" % (totconf/nalign if nalign else 0)
        varsP.updateInfoReport(outstr, printalso=True, simple=simple, simplenl=True) #simple new line bc extra newline below missing
    avgfp  = (sum(fplist)/len(fplist)   if len(fplist) else 0)
    avgfpr = (sum(fprlist)/len(fprlist) if len(fprlist) else 0)
    avgfn  = (sum(fnlist)/len(fnlist)   if len(fnlist) else 0)
    avgbpp = (sum(bpplist)/len(bpplist) if len(bpplist) else 0)
    avgres = (sum(reslist)/len(reslist) if len(reslist) else 0)
    avgllr = (sum(llrmlist)/len(llrmlist) if len(llrmlist) else 0)
    avgllg = (sum(llrgmlist)/len(llrgmlist) if len(llrgmlist) else 0)
    avgbps = (sum(bppsdlist)/len(bppsdlist) if len(bppsdlist) else 0)
    avgsf  = (sum(sflist)/len(sflist) if len(sflist) else 0)
    avgsd  = (sum(sdlist)/len(sdlist) if len(sdlist) else 0)
    avgsr  = (sum(srlist)/len(srlist) if len(srlist) else 0)
    avgrsd = (sum(resdlist)/len(resdlist) if len(resdlist) else 0)
    if avgfp or avgfn or avgbpp :
        outstr =  "Avg FP(/100kb)    : %12.2f\n" % avgfp
        outstr += "Avg FP ratio      : %13.3f\n" % avgfpr
        outstr += "Avg FN ratio      : %13.3f\n" % avgfn
        outstr += "Avg bpp           : %11.1f\n" % avgbpp
        outstr += "Avg sf            : %13.3f\n" % avgsf
        outstr += "Avg sd            : %13.3f\n" % avgsd
        outstr += "Avg sr            : %13.3f\n" % avgsr
        varsP.updateInfoReport(outstr + "\n", printalso=True, simple=False) #do not put these in simple
    if err and mergeerr and mergepath : #have an error file (alignParams) object
        util.checkDir(mergepath)
        mrgstr = (varsP.alignMolvrefMergeName if varsP else "merge")
        outpath = os.path.join(mergepath, mappref+mrgstr+".err")
        err.fp = avgfp
        err.fn = avgfn
        err.sf = avgsf
        err.sd = avgsd
        err.bpp = avgbpp
        err.res = avgres
        err.nmaps = sum(nmaplist)
        err.llrm  = avgllr
        err.goodmaps = sumgoodmaps
        err.llrgm = avgllg
        err.bppsd = avgbps
        err.fprate = avgfpr
        err.sr = avgsr
        err.ressd = avgrsd
        err.writeToFile(outpath)

#end getAlignStats



def mergeMap(varsP, outFileList, mergepath) :
    """outFileList is list of path+prefixes--each should have a .map file:
    merge them to a merged .map file in dir mergepath."""

    outFileList.sort() #sort to ensure reproducibility (order of entries)
    maplist = []
    for outpath in outFileList : #these are file prefixes
        if util.checkFile(outpath+".map") :
            maplist.append(outpath+".map")
        elif varsP :
            varsP.updatePipeReport("Warning in AlignModule.mergeMap: missing map: "+outpath+".map"+"\n")
        else :
            print "Warning in AlignModule.mergeMap: missing map: "+outpath+".map"+"\n"

    if not len(maplist) : #nothing to merge
        return

    if not util.checkDir(mergepath) :
        varsP.updatePipeReport("Warning in AlignModule.mergeMap: merge path invalid: "+mergepath+"\n")
        return

    headstart = ["#", "S", "M"] #last two lines of header start with "Software" and "MappedMoleculeId"
    #header = ""
    headerdone = False
    #data = ""
    lineno = 1 #can't just append: need to change index in first column
    sep = "\t"
    mappref = getMergeFilename(outFileList[0]) #also in getAlignStats
    mrgstr  = (varsP.alignMolvrefMergeName if varsP else "merge") #same for vref and not
    outpath = os.path.join(mergepath, mappref+mrgstr+".map")
    f1 = open(outpath, 'w')
    for path in maplist :
        f = open(path)
        for line in f :
            if line[0] in headstart and not headerdone :
                #header += line
                f1.write(line)
            elif line[0] not in headstart :
                tokens = line.split()
                tokens[0] = str(lineno)
                #data += sep.join(tokens)+"\n" #newline was stripped by split
                f1.write(sep.join(tokens)+"\n")
                lineno += 1
        headerdone = True
        f.close()

    #f1.write(header+data) 
    f1.close()
#end mergeMap


def getMergeFilename(mappref) :
    mappref = os.path.split(mappref)[1] #this is path + filename prefix, but possibly with integer suffix
    #if mappref.count("_") > 1 : #rather than checking multiple, check after any _
    pos = mappref.rfind("_")
    if pos != -1 :
        suf = mappref[pos+1:] #check if suffix is integer
        if util.getIntFromString(suf) != None :
            mappref = mappref[:pos+1] #remove integer suffix
        elif pos < len(mappref)-1 : #same as no suffix below, don't make "__"
            mappref += "_" #because mrgstr is appended
    else :
        mappref += "_" #because mrgstr is appended
    return mappref
#end getMergeFilename

    
#def split_XMap_byContig(outFileList, mergepath, varsP=None, stageName="alignmol") :
#    """outFileList is list of path+prefixes--each should have a .xmap and _q.cmap file:
#    split into one per contig."""
#    logOrPrintError("Start split_XMapQcmap_byContig", varsP, warn=False)
#    xmapFilelist = []
#    for outpath in outFileList : #these are file prefixes
#        if util.checkFile(outpath+".xmap") :
#            xmapFilelist.append(outpath+".xmap")
#        else :
#            err_msg = "Warning in AlignModule.split_XMapQcmap_byContig: missing xmap: "+outpath+".xmap"
#            logOrPrintError(err_msg, varsP)
#
#    if not len(xmapFilelist) : #nothing to merge
#        err_msg = "Error in AlignModule.split_XMapQcmap_byContig: no xmaps found"
#        logOrPrintError(err_msg, varsP)
#        return
#    
#    outFileList.sort() #sort to ensure reproducibility (order of entries)
#    xmapDict={}
#    aHeaderTemplate = ''
#    for path in xmapFilelist :
#        aXmap = mc.xmap(sourceFile=path)
#        if aHeaderTemplate == '' : 
#           aHeaderTemplate = aXmap.header
#        for xentry in aXmap.xmapLookup.values() :
#            refC = xentry.contigRef
#            if not xmapDict.has_key(refC):
#                xmapDict[refC] = mc.xmap() 
#            cIndex = len(xmapDict[refC].xmapLookup) + 1
#            xmapDict[refC].xmapLookup[cIndex] = xentry
#        # end of aXmap entry
#    # end of maplist
#    filepref = (varsP.outputContigPrefix if varsP else stageName) #same as line in mergeRcmaps
#    # read header and then write out the xmap file:      
#    for refContig in xmapDict.keys() :
#        xmapDict[refContig].header = aHeaderTemplate
#        outpref = os.path.join(mergepath, filepref+'_contig'+str(refContig))
#        #print "Writing", filepref+".xmap"
#        xmapDict[refContig].editHeaderQueryMaps(outpref+'_q.cmap')         
#        xmapDict[refContig].editHeaderMaps(outpref+'_r.cmap', query=False)  
#        xmapDict[refContig].writeToFile(outpref+'.xmap')
#    # end of header 
#    logOrPrintError("split_XMapQcmap_byContig: wrote %i xmaps" % len(xmapDict), varsP, warn=False)
#    return xmapDict
#end split_XMap_byContig

def split_XMap_byContig_new(outFileList, mergepath, varsP=None, stageName="") :
    """outFileList is list of path+prefixes--each should have a .xmap and _q.cmap file:
    split into one per contig."""
    logOrPrintError("Start split_XMapQcmap_byContig", varsP, warn=False)
    xmapFilelist = []
    for outpath in outFileList : #these are file prefixes
        if util.checkFile(outpath+".xmap") :
            xmapFilelist.append(outpath+".xmap")
        else :
            err_msg = "Warning in AlignModule.split_XMapQcmap_byContig: missing xmap: "+outpath+".xmap"
            logOrPrintError(err_msg, varsP)

    if not len(xmapFilelist) : #nothing to merge
        err_msg = "Error in AlignModule.split_XMapQcmap_byContig: no xmaps found"
        logOrPrintError(err_msg, varsP)
        return

    #filepref = (varsP.outputContigPrefix if varsP else stageName) #same as line in mergeRcmaps
    filepref = (varsP.outputContigPrefix if varsP and stageName == "" else stageName) #same as line in mergeRcmaps

    outFileList.sort() #sort to ensure reproducibility (order of entries)
    xmapLineDict = {} #if you store the number of lines here, you can avoid counting every time it's opened
    xmapMolDict = {} #store the molecule IDs here for use in split_Qcmap_byContig
    newxmaplist = [] #store paths of output xmaps to fix their headers
    header = "" #get header of first file
    with open(outFileList[0]+".xmap") as f1 :
        for line in f1 : #no readline--that will iterate over each char in line instead of line itself
            if line[0] == '#':
                header += line
            else :
                break
    for path in xmapFilelist :
        f = open(path)
        for line in f :
            if line[0] == "#" : #get header separately above
                continue
            #I don't think there's any way to avoid split, except looping over chars, but that's probably just as slow
            tokens = line.split()
            try:
                qryid = int(tokens[1])
                refid = int(tokens[2])
            except:
                continue

            outpref = os.path.join(mergepath, filepref+'_contig'+str(refid)) 
            if not outpref in newxmaplist : #this loop is every line; don't duplicate
                newxmaplist.append(outpref) #prefixes
            outf = open(outpref+".xmap", "a+") #make a new file if not exists; if does, points to end of file
            if not refid in xmapLineDict :
                xmapLineDict[refid] = 1
                xmapMolDict[refid] = [qryid]
                outf.write(header) #write header to disk
            else :
                xmapLineDict[refid] += 1
                #because xmapMolDict is used to make the _q.cmap, its entries should be unique
                #assert xmapMolDict[refid].count(qryid) == 0, ("dup molid %i, path %s" % (qryid, path))
                if not qryid in xmapMolDict[refid] :
                    xmapMolDict[refid].append(qryid)
            #outf.write("\t".join([str(xmapLineDict[refid])]+tokens[1:])+"\n")
            tokens[0] = str(xmapLineDict[refid])
            outf.write("\t".join(tokens)+"\n")
            outf.close() #avoid keeping too many file handles open at the expense of re-open many times
        #end for line in f
        f.close()
    #end for xmapFilelist

    #need to fix headers still, ie, the editHeaderMaps/QueryMaps: must re-read and -write files
    for path in newxmaplist :
        with open(path+".xmap", "r") as f :
            lines = f.readlines()
        with open(path+".xmap", "w") as f :
            for line in lines :
                if line.find("Query Maps") != -1 :
                    line = line.split(":")[0] + ":\t" + path + "_q.cmap" + "\n"
                elif line.find("Reference Maps") != -1 :
                    line = line.split(":")[0] + ":\t" + path + "_r.cmap" + "\n"
                f.write(line)
    logOrPrintError("split_XMapQcmap_byContig: wrote %i xmaps" % len(xmapMolDict), varsP, warn=False) #reproduce original fn
    if 0 :
        bad = False 
        print "DEBUG:"
        for xl in xmapMolDict.values() : #list of mols
            for i in xl :
                if xl.count(i) > 1 :
                    bad = True
                    print i
        if bad :
            print xmapMolDict
        print "DEBUG\n"
    return(xmapMolDict)
#end split_Xmap_byContig_new


#deprecate in favor of split_Qcmap_byContig_new
#def split_Qcmap_byContig(outFileList, mergepath, xmapDict, varsP=None, stageName="alignmol") :
#    """outFileList is list of path+prefixes--each should have a .xmap and _q.cmap file:
#    split into one per contig."""
#
#    # readin all _q.cmap:
#    qcmapFilelist = []
#    for outpath in outFileList : #these are file prefixes
#        if util.checkFile(outpath+"_q.cmap") :
#            qcmapFilelist.append(outpath+"_q.cmap")
#        else :
#            err_msg = "Warning in AlignModule.split_XMapQcmap_byContig: missing _q.cmap: "+outpath+"_q.cmap"
#            logOrPrintError(err_msg, varsP)
#
#    if not len(qcmapFilelist) : #nothing to merge
#        err_msg = "Error in AlignModule.split_XMapQcmap_byContig: no _q.cmaps found"
#        logOrPrintError(err_msg, varsP)
#        return     
#    mergedQCMap = mc.multiCmap() 
#    cmap_header = ''  
#    for pathIndex in range(0, len(qcmapFilelist)) :
#        tmpMCmap = mc.multiCmap(infile=qcmapFilelist[pathIndex])
#        #this will work without BestRef because each entry in different files should be identical...
#        mergedQCMap.cmapdict.update(tmpMCmap.cmapdict)
#        mergedQCMap.totalLength += tmpMCmap.totalLength
#        if pathIndex == 0: 
#            cmap_header = tmpMCmap.header
#    filepref = (varsP.outputContigPrefix if varsP else stageName) #same as line in mergeRcmaps
#    # write out _q.cmap:
#    for refContig in xmapDict.keys() :
#        outQmapFile = os.path.join(mergepath, filepref+'_contig'+str(refContig)+'_q.cmap')
#        aNewMCMap = mc.multiCmap() 
#        aNewMCMap.header = cmap_header
#        for xmapEntry in xmapDict[refContig].xmapLookup.values() :
#            qContigID = int(xmapEntry.contigQry)
#            if not aNewMCMap.cmapdict.has_key(qContigID) : #...and nothing prevents the same cmap being used in two different files
#                aNewMCMap.cmapdict[qContigID] = mergedQCMap.cmapdict[qContigID]
#                aNewMCMap.totalLength += mergedQCMap.cmapdict[qContigID].length
#        aNewMCMap.writeToFile(outQmapFile) 
#    # end of write out _q.cmap          
#    logOrPrintError("split_XMapQcmap_byContig: wrote %i _q.cmaps" % len(xmapDict), varsP, warn=False) 
#end split_Qcmap_byContig


def split_Qcmap_byContig_new(inFileList, mergepath, xmapDict, varsP=None, stageName="") :
    # readin all _q.cmap:
    qcmapFilelist = []
    for outpath in sorted(inFileList) : #these are file prefixes--sort to ensure reproducibility
        if util.checkFile(outpath+"_q.cmap") :
            qcmapFilelist.append(outpath+"_q.cmap")
        else :
            err_msg = "Warning in AlignModule.split_XMapQcmap_byContig: missing _q.cmap: "+outpath+"_q.cmap"
            logOrPrintError(err_msg, varsP)

    if not len(qcmapFilelist) : #nothing to merge
        err_msg = "Error in AlignModule.split_XMapQcmap_byContig: no _q.cmaps found"
        logOrPrintError(err_msg, varsP)
        return     

    header = "" #get header of first qcmap
    with open(qcmapFilelist[0]) as f1 :
        for line in f1 : #no readline--that will iterate over each char in line instead of line itself
            if line[0] == '#':
                header += line
            else :
                break
    #create all output files, header only
    #filepref = (varsP.outputContigPrefix if varsP else stageName) #same as line in mergeRcmaps
    filepref = (varsP.outputContigPrefix if varsP and stageName == "" else stageName) #same as line in mergeRcmaps
    for contigid in xmapDict.keys() :
        outQmapFile = os.path.join(mergepath, filepref+'_contig'+str(contigid)+'_q.cmap')
        f1 = open(outQmapFile, "w")
        f1.write(header)
        f1.close()
    #convert xmapDict to a molDict: keys are molids, and values are contig ids -- this should speed up below
    molDict = {}

    #for cid,xmap in xmapDict.iteritems() :
        #for xmapentry in xmap.xmapLookup.values() :
            #if not molDict.has_key(xmapentry.contigQry) : #new mol
            #    molDict[xmapentry.contigQry] = [xmapentry.contigRef]
            #else :
            #    molDict[xmapentry.contigQry].append(xmapentry.contigRef)
    #old xmapDict was contigid:"xmap object"; new one is contigid:"list of mol ids"
    for cid,molids in xmapDict.iteritems() :
        for molid in molids :
            if not molDict.has_key(molid) : #new mol
                molDict[molid] = [cid]
            else :
                molDict[molid].append(cid)
    #print "DEBUG:\n", molDict, "DEBUG\n" #debug
    #read input files, find all contigs to which each molecule aligns, write to that qcmap
    nmol = 0
    molid = -1
    for qcmap in qcmapFilelist :
        previd = 0 #molecule id from _q.cmap, int to compare with xmap.contigQry
        molstr = "" #all the lines in the _q.cmap for this molecule
        f1 = open(qcmap)
        for line in f1 :
            if line[0] == '#' :
                continue
            molid = int(line.split()[0]) #use int bc compare to xmap.contigQry
            if molid == previd : #get data for this mol
                molstr += line
            else : #write previous mol to output qcmap
                if molstr : #not for first mol
                    for cid in molDict[previd] :
                        outQmapFile = os.path.join(mergepath, filepref+'_contig'+str(cid)+'_q.cmap')
                        f2 = open(outQmapFile, "a")
                        f2.write(molstr)
                        f2.close()
                #prepare for next mol
                molstr = line
                previd = molid
                nmol += 1
        f1.close()
        if molDict.has_key(molid) :
            #get last molecule
            for cid in molDict[molid] :
                outQmapFile = os.path.join(mergepath, filepref+'_contig'+str(cid)+'_q.cmap')
                f2 = open(outQmapFile, "a")
                f2.write(molstr)
                f2.close()
    logOrPrintError("split_XMapQcmap_byContig: wrote %i _q.cmaps with %i molecules" % (len(xmapDict), nmol), varsP, warn=False) 
#end split_Qcmap_byContig_new


def logOrPrintError(err_msg, varsP=None, warn=True) :
    """If varsP, call updatePipeReport, otherwise print. """
    if varsP :
        varsP.updatePipeReport(err_msg+"\n")
        if warn :
            util.LogError("warning", err_msg) #assume this makes sense only if varsP
    else :
        print err_msg
#end logOrPrintError


def mergeRcmaps(outFileList, outdir, varsP=None, splitByContig=None, stageName="") :
    """Given a list of file prefixes (outFileList), append "_r.cmap" to them, and merge them
    to outdir. Report to varsP if supplied, stdout if not.
    Also support outFileList is full paths (including "_r.cmap").
    If splitByContig < 1, output each contig separately, if == 1, only output single merged cmap,
    and if > 1, do both.
    Always use stagename if supplied; if not, must supply varsP otherwise prefix is empty.
    """
    
    if not util.checkDir(outdir) :
        err_msg = "Warning in AlignModule.mergeRcmaps: could not make outdir %s, skipping copy number" % outdir
        logOrPrintError(err_msg, varsP)
        return

    if not outFileList : #just an argument check--check for presence on disk is below
        err_msg = "Warning in AlignModule.mergeRcmaps: no maps supplied"
        logOrPrintError(err_msg, varsP)
        return

    outFileList.sort() #for reproducibility with runAlignMerge.py (different order when listing dir)
    rsuf = "_r.cmap"
    #mappref = os.path.split(outFileList[0])[1] #this is just prefix, but with integer suffix--get it before -- no longer used
    #mappref = mappref[:mappref.rfind("_")+1] #remove integer suffix
    #even though outFileList should all be there, a job may have failed--check all, just existence
    present = []
    for outf in outFileList :
        target = (outf+rsuf if not outf.endswith(rsuf) else outf) #now support either
        if not util.checkFile(target) :
            err_msg = "Warning in AlignModule.mergeRcmaps: missing _r.cmap %s" % target
            logOrPrintError(err_msg, varsP)
        else :
            present.append(target)
    if not present : #no _r.cmaps found (this will also happen for empty outFileList)
        err_msg = "Warning in AlignModule.mergeRcmaps: no _r.cmaps found, skipping copy number"
        logOrPrintError(err_msg, varsP)
        return
    outFileList = present #yes, it's redundant, but now have rsuf appended

    mrgstr = (varsP.alignMolvrefMergeName if varsP else "merge")
    #mergedmappath = os.path.join(outdir, mappref+mrgstr+rsuf) #this is output merged _r.cmap -- unify with filepref

    mergedmap = mc.multiCmap(outFileList[0]) #open original, edit in memory
    #now add other maps
    for rmap in outFileList[1:] : #don't add map 0 to itself
        if mergedmap.addCovOcc( mc.multiCmap(rmap) ) : #when calling addCovOcc, check return, warn if True
            err_msg = "Warning in AlignModule.mergeRcmaps: addCovOcc call failed for map %s" % rmap
            logOrPrintError(err_msg, varsP)
    #now it's merged, but the resulting map need to be written back to disk
    filepref = (varsP.outputContigPrefix if varsP and stageName == "" else stageName) #see split_XMapQcmap_byContig
    if splitByContig < 1 or splitByContig > 1 :
        #print "\nself.varsP.outputContigPrefix", self.varsP.outputContigPrefix, "\n" #debug
        #filepref = (varsP.outputContigPrefix if varsP else stageName) #same as line in split_XMapQcmap_byContig
        mergedmap.writeAllMapsToDisk( os.path.join(outdir, filepref+'_contig'), outsuf="_r" )
        report = "mergeRcmaps: wrote %i cmaps" % len(mergedmap.cmapdict)
    if splitByContig > 0 :
        mergedmap.writeToFile( os.path.join(outdir, filepref+"_"+mrgstr+rsuf) ) #was mergedmappath
        report = "mergeRcmaps: wrote merged cmap with %i contigs" % len(mergedmap.cmapdict)
    #report result
    logOrPrintError(report, varsP, warn=False)

#end mergeRcmaps
