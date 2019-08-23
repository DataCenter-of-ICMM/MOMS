
import os

import mapClasses
#import MapClassesRev
#import RefinementModule as rm
import Multithreading as mthread
import utilities as util

"""
@package CharacterizeModule 
Get general stats and mapping stats (if reference) for contigs

"""


util.setVersion("$Id: CharacterizeModule.py 6703 2017-07-20 18:28:03Z wandrews $")


class dummyCharacterize() :
    """For getting noise parameters in case of bypassing characterize."""
    def __init__(self, varsP) :
        self.curCharacterizeFileRoots = []
        self.varsP = varsP #bc Characterize uses this for totAssemblyLenMb

        outdir = os.path.join(varsP.outputContigFolder, self.varsP.characterizeDirName) #'alignref'
        if not util.checkDir(outdir, makeIfNotExist=False) : #if this doesn't exist, we can't get what we need
            err = "characterize dir not found (%s)" % outdir
            self.varsP.updatePipeReport( "warning : "+err+"\n" )
            util.LogError("warning", err)
            return
        outfile = ""
        for qfile in os.listdir(outdir) :
            if qfile.endswith(".err") : #just take first .err file
                outfile = qfile
                break
        if not outfile : #if no .err files found, give up
            err = "characterize err file not found in dir (%s)" % outdir
            self.varsP.updatePipeReport( "warning : "+err+"\n" )
            util.LogError("warning", err)
            return
        outfile = os.path.join(outdir, outfile.replace(".err",""))
        self.curCharacterizeFileRoots.append(outfile)
        #print "dummyCharacterize: curCharacterizeFileRoots =", self.curCharacterizeFileRoots #debug
        #also want to get varsP.totAssemblyLenMb
        self.varsP.totAssemblyLenMb = mapClasses.multiCmap(varsP.latestMergedCmap, lengthonly=True).totalLength / 1e6
#end class dummyCharacterize


class Characterize(mthread.jobWrapper):
    """Align contigs to reference to characterize quality of assembly.
    If no reference, report basic contig stats only.
    
    """
    def __init__(self, varsP, argset=-1):
        '''argset is toggle between CharacterizeDefault and CharacterizeFinal argumets:
        -1 is default, 1 is final
        '''
        self.varsP = varsP
        #self.argStr = ("Final" if argset == 1 else "Default") #!=1 and !=-1 is error in generateJobList, but not here
        self.argStr = ("F" if argset == 1 else "D") #shortened for easier reading
        self.stageName = 'Characterize'+self.argStr+'_'+self.varsP.stageComplete
        util.LogStatus("progress", "stage_start", self.stageName)
        mthread.jobWrapper.__init__(self, varsP, self.stageName,clusterArgs=varsP.getClusterArgs('characterizeDefault'))
        self.xmapTarget = None
        self.curCharacterizeFileRoots = []
        outdir = self.varsP.characterizeDirName # = 'alignref'
        if argset == 1 : #this is final
            outdir += '_final'
        varsP.contigAlignTarget = os.path.join(varsP.outputContigFolder, outdir)
        if not(os.path.exists(varsP.contigAlignTarget)):
            os.mkdir(varsP.contigAlignTarget)
        self.generateJobList(argset)
        
    def generateJobList(self,argset=-1):
        if not self.varsP.ref : #no jobs if no ref
            return
        jobargs = [self.varsP.RefAlignerBin, '-ref', self.varsP.ref]
        if argset == -1 and self.varsP.argData.has_key('characterizeDefault') : # don't use nominal default
            opta = self.varsP.argsListed('characterizeDefault')
        elif argset == 1 and self.varsP.argData.has_key('characterizeFinal') : #extend (on default) -- make this default
            opta = self.varsP.argsListed('characterizeFinal')
        else : #this is an error
            self.varsP.updatePipeReport("ERROR in CharacterizeModule.generateJobList: invalid argset %s\n" % str(argset))
            return
        
        for i, curCharacterizeCmap in enumerate(self.varsP.curCharacterizeCmaps):
            if self.varsP.numCharacterizeJobs == 1:
                jobName = 'Char'+self.argStr+'_%s' % self.varsP.stageComplete
            else:
                jobName = 'Char'+self.argStr+'_%s_%d' % (self.varsP.stageComplete, i+1)
            outFileName = os.path.split(curCharacterizeCmap)[-1].replace(".cmap", "")
            outfile = os.path.join(self.varsP.contigAlignTarget,outFileName)
            self.curCharacterizeFileRoots.append(outfile)
            expectedResultFile = outfile+".xmap"
            self.xmapTarget = expectedResultFile
            currentArgs = jobargs + ["-i", curCharacterizeCmap, "-o", outfile]
            stdoutf = None
            if self.varsP.stdoutlog :
                currentArgs.extend( ['-stdout', '-stderr'] )
                stdoutf = outfile+".stdout"
            currentArgs += ['-maxthreads', str(self.varsP.nThreads)]
            currentArgs += ['-output-veto-filter', '_intervals.txt$']
            currentArgs += opta
            s1Job = mthread.singleJob(currentArgs, jobName, expectedResultFile, jobName.replace(' ',''),maxThreads=self.varsP.nThreads,clusterLogDir=self.varsP.clusterLogDir, expectedStdoutFile=stdoutf)
            self.addJob(s1Job)
            if i==0:
                self.logArguments()
        
    def runJobs(self):
        if not self.varsP.ref :
            return
        self.multiThreadRunJobs(self.varsP.nThreads, sleepTime = 0.2)
    
    def checkResults(self):
        #old heading says complete here and then summary after contig list; new says summary here
        outstr = 'Stage Summary: %s\n' % self.stageName
        if not self.varsP.ref : #still want contig stats
            infoReport = "Skipping Characterize because no reference (-r)\n"
            self.varsP.updatePipeReport(infoReport, printalso=False) #put this in pipereport just as an fyi
            infoReport += outstr
            #infoReport += 'Stage Complete: %s\n' % self.groupName #set in jobWrapper constructor
            #infoReport += MapClassesRev.ContigCharacterizationNoRef(self.varsP,self.groupName)
            infoReport += characterizeContigs(self.varsP)
            #there is no CharacterizeFinal for no ref, so this only happens once
            simple = self.varsP.stageComplete.startswith("refineFinal")
            self.varsP.updateInfoReport(infoReport + '\n', simple=simple)
            return
        self.doAllPipeReport()
        #infoReport = 'Stage Complete: %s\n' % self.groupName #set in jobWrapper constructor
        #infoReport += MapClassesRev.TopLevelCharacterization(self.varsP,self.curCharacterizeFileRoots,self.groupName)
        #infoReport += 'OLD characterize\n' #debug
        infoReport = characterizeContigs(self.varsP, self.xmapTarget)
        if infoReport : #None on failure
            simple = (self.argStr == "F" and self.varsP.stageComplete.startswith("refineFinal"))
            self.varsP.updateInfoReport(outstr + infoReport + '\n', simple=simple)
        util.LogStatus("progress", "stage_complete", self.stageName)

#end class Characterize


class referenceProcess(mthread.jobWrapper) :
    """Pre-process the reference for SV jobs, applying -mres.
    """

    def __init__(self, varsP) :
        jobName = "reference_process"
        opta_section = "referenceSvdetect"
        default_mres = "2.9"
        mres = "-mres"
        self.varsP = varsP
        usedefault = False
        if self.varsP.argData.has_key(opta_section) : #check if in optargs
            opta = self.varsP.argsListed(opta_section)
            if not mres in opta : #must have mres
                self.varsP.updatePipeReport("Warning in referenceProcess: "+mres+" missing in optArguments section "+opta_section+"\n")
                usedefault = True
        else :
            self.varsP.updatePipeReport("Warning in referenceProcess: optArguments section "+opta_section+" missing\n")
            usedefault = True
        if usedefault :
            opta = [mres, default_mres]

        mresstr = opta[opta.index(mres)+1] #get string for mres value for output name
        mresstr = mresstr.replace(".","")

        if not util.checkDir(self.varsP.refFolder) :
            self.varsP.updatePipeReport( "ERROR in referenceProcess: could not make output dir %s\n" % self.varsP.refFolder )
            return None
        refpref = os.path.basename(self.varsP.ref[:self.varsP.ref.rfind(".")]) + "_res" + mresstr
        outarg = os.path.join(self.varsP.refFolder, refpref) #refFolder is new output folder for this job
        expectedResultFile = outarg+".cmap" #if ref is spots, is this spots?
        args = [self.varsP.RefAlignerBin, '-f', '-o', outarg, '-i', self.varsP.ref, '-merge'] + opta
        stdoutf = None
        if self.varsP.stdoutlog :
            args.extend( ['-stdout', '-stderr'] )
            stdoutf = outarg+".stdout"
        args += ['-maxthreads', str(self.varsP.nThreads)]

        super(referenceProcess, self).__init__(self.varsP, jobName, clusterArgs=self.varsP.getClusterArgs("autoNoise0"))

        job = mthread.singleJob(args, jobName, expectedResultFile, jobName, maxThreads=self.varsP.nThreads, clusterLogDir=self.varsP.clusterLogDir, expectedStdoutFile=stdoutf)
        self.addJob(job)

        util.LogStatus("progress", "stage_start", jobName)
        self.varsP.runJobs(self, "referenceProcess")
        self.doAllPipeReport()
        if not self.allResultsFound() : #this is an error, but we'll continue processing without SV detect
            err = "ERROR in referenceProcess: job failed, disabling SV detect"
            self.varsP.updatePipeReport( err+"\n" )
            util.LogError("error", err)
            #self.varsP.runSV = False #no need since this class is used in SVModule
        else :
            self.varsP.refDeresed = expectedResultFile #store good result for SV detect
            self.varsP.updatePipeReport( "referenceProcess: using reference %s for svdetect\n" % self.varsP.refDeresed )
        util.LogStatus("progress", "stage_complete", jobName)            
    #end __init__
            
    def runJobs(self):
        self.multiThreadRunJobs(self.varsP.nThreads, sleepTime = 0.2)

#end class referenceProcess



#get the mapped length and the alignment parameters from the align to ref call
#mapped, params = getMappedStats(aligndir, cpath)
#return the alignParams object so that it can be used in characterizeContigs
def getMappedErrStats(aligndir, cpath):
    """Only used in characterizeContigs if listcontigs.
    """
    errfile = aligndir+os.path.split(cpath)[-1].replace(".cmap", ".err") #path checking exists in alignparams constructor
    return alignParams(errfile)


    
def characterizeContigs(varsP, xmappath=None) :
    """Log simple contigs stats, and optionally align stats from xmappath.
    """
    if not varsP.latestMergedCmap : #if this is missing, then characterize probably failed, so don't bother trying to read xmap
        err = "characterizeContigs: missing merged cmap for stage %s" % self.varsP.stageComplete
        self.varsP.updatePipeReport( "warning: "+err+"\n" )
        util.LogError("warning", err)
        return None
    #print "xmappath:", xmappath
    unitscale = 1e-6
    dorefalign = bool(xmappath) #i'm never actually calling refaligner here--this is just using xmappath
    haveref = bool(varsP.ref)

    #refcmap = mapClasses.multiCmap() #not used
    aligndir = varsP.contigAlignTarget

    try :
        #refcmap = mapClasses.multiCmap(varsP.ref)
        #reflen = refcmap.totalLength #note: total length of _all_ contigs
        reflen = mapClasses.multiCmap(varsP.ref, lengthonly=True).totalLength
        #in summary table, this is a denominator--make sure it's non-zero, don't bail (still get contig summary)
        if reflen <= 0 :
            #print "Warning in CharacterizeModule.characterizeContigs: bad reflen", reflen, "defaulting to 1" #not necessary
            reflen = 1.
    except:
        reflen = 1.

    outstr = "" #Contig Characterization:\n"

    #check for .hmaps in same dir as latestMergedCmap: if any, add a line for haploid genome size
    hmaps = util.getListOfFilesFromDir(os.path.dirname(varsP.latestMergedCmap), ".hmap")
    haplotype = (len(hmaps) > 0)
    haplotypelen = 0
    hapcontiglens = []

    totcontiglen = 0; totalignlen = 0; nmapcontigs = 0; totalignqlen = 0 #defalignlen = 0; 
    contiglens = [] #lens of all contigs in bases
    uniqueseg = {} #the argument to util.uniqueRange--stores all the map lengths which go into totalignlen--now each of these is a value, and the keys are the reference contig id or chromosome depending on dorefidchr
    for citr, cpath in enumerate([varsP.latestMergedCmap]) : #always use contigpaths
        mapi = mapClasses.multiCmap(cpath) 
        totcontiglen += mapi.totalLength
        contiglens += mapi.getAllMapLengths() #getAllMapLengths is list of all map lengths
        if haplotype :
            haplotypelen += mapi.getHaplotypeTotalMapLength()
            hapcontiglens.extend( mapi.getHaplotypeMapLengths() )

        #store a list of the contig ids in this multiCmap, then remove them if they're in the xmap
        # if they're not, print at the end
        mapids = mapi.getAllMapIds() #this is once per cmap, not once per characterizeModule call--becuase it's keys, it's already a copy, no need to copy explicitly
        ncontigs = len(mapids) #this is ncontigs in this file, ie, in mapi (see below)

        xmapobj = mapClasses.xmap() #empty map to fix xmapobj scope
        if dorefalign : #get xmap object
            if util.checkFile(xmappath, ".xmap") :
                xmapobj = mapClasses.xmap(xmappath)

        for xitr, xmapentry in enumerate(xmapobj.xmapLookup.values()) :

            #get map length from multicmap.getMapLength--returns 0 for any exception
            contiglen = mapi.getMapLength(xmapentry.contigQry)
            if contiglen <= 0 : #this strikes me as clumsy...but I don't want to return non-zero from multiCmap.getMapLength
                contiglen = 1.
            contigcov = mapi.getMapAvgCoverage(xmapentry.contigQry)

            #don't print lenr for each contig--just total them
            lenr = xmapentry.getMappedRefLen()
            lenq = xmapentry.getMappedQryLen()
            refid = xmapentry.contigRef #int

            totalignlen  += lenr
            totalignqlen += lenq

            #uniqueseg is now a dict to take into account which chromosome the query contig is on
            #note need refid bc need to separate different contigs on the _same_ chromosome
            if not uniqueseg.has_key(refid) : #if first contig on chromosome, need to init new list
                uniqueseg[refid] = []
            uniqueseg[refid].append( [xmapentry.RefStart, xmapentry.RefStop] )

            #process mapids--remove contig id (contigQry) from mapids if they're in the xmap so non-aligning contigs can be printed
            if xmapentry.contigQry in mapids :
                mapids.remove(xmapentry.contigQry)
            
        #end loop on xmap entries

        #now that all xmap entries are processed, all contigs with an alignment are removed from mapids,
        # so we can get n contigs align using this and ncontigs
        nmapcontigs += ncontigs - len(mapids) #sum multiple cmaps

    #end loop on contigs

    varsP.totAssemblyLenMb = totcontiglen*unitscale
    ncontigs = len(contiglens) #contigpaths is just files--contiglens is all contigs
    avgcontiglen = (float(totcontiglen)/ncontigs if ncontigs > 0 else 0)

    if unitscale > 1e-6 : #if not megabases
        fstr = "%9.0f"
    else : #megabases
        fstr = "%8.3f" 

    if haplotype : #new format for haplotype
        #if haplotypelen != sum(hapcontiglens) : #simply print warning in this case (do not log): ignore this bc of floating point rounding
            #print "Warning in characterizeContigs: haplotype lengths are inconsistent:", haplotypelen, sum(hapcontiglens)
        #diploid is same as else below, but names change
        outstr += "Diploid number Genome Maps: %i\n" % ncontigs
        outstr += ("Diploid Genome Map Length        (Mbp): "+fstr+"\n") % (totcontiglen*unitscale)
        outstr += ("Diploid Mean Genome Map Length   (Mbp): "+fstr+"\n") % (avgcontiglen*unitscale)
        outstr += ("Diploid Median Genome Map Length (Mbp): "+fstr+"\n") % (util.getMedian(contiglens)*unitscale)
        outstr += ("Diploid Genome Map N50           (Mbp): "+fstr+"\n") % (util.getn50(contiglens)*unitscale)
        #haploid : ignore haplotypelen, just use the list hapcontiglens
        outstr += "Haploid number Genome Maps: %i\n" % len(hapcontiglens)
        tot = sum(hapcontiglens)
        avg = ( tot/len(hapcontiglens) if len(hapcontiglens) else 0 )
        outstr += ("Haploid Genome Map Length        (Mbp): "+fstr+"\n") % (tot*unitscale)
        outstr += ("Haploid Mean Genome Map Length   (Mbp): "+fstr+"\n") % (avg*unitscale)
        outstr += ("Haploid Median Genome Map Length (Mbp): "+fstr+"\n") % (util.getMedian(hapcontiglens)*unitscale)
        outstr += ("Haploid Genome Map N50           (Mbp): "+fstr+"\n") % (util.getn50(hapcontiglens)*unitscale)
    else : #default to old format
        outstr += "Number Genome Maps: %i\n" % ncontigs
        outstr += ("Total Genome Map Length  (Mbp): "+fstr+"\n") % (totcontiglen*unitscale)
        outstr += ("Mean Genome Map Length   (Mbp): "+fstr+"\n") % (avgcontiglen*unitscale)
        outstr += ("Median Genome Map Length (Mbp): "+fstr+"\n") % (util.getMedian(contiglens)*unitscale)
        outstr += ("Genome Map N50           (Mbp): "+fstr+"\n") % (util.getn50(contiglens)*unitscale)

    if haveref :
        outstr += ("Total Reference Length   (Mbp): "+fstr+"\n") % (reflen*unitscale)
        outstr += ("Total Genome Map Length / Reference Length : "+fstr+"\n") % (totcontiglen/reflen)
    if dorefalign :
        ratio = (float(nmapcontigs)/ncontigs if ncontigs > 0 else 0)
        outstr += ("Total number of aligned Genome Maps    : %i (%.2f)\n") % (nmapcontigs, ratio)
        outstr += ("Total Aligned Length (Mbp)             : "+fstr+"\n") % (totalignlen*unitscale)
        outstr += ("Total Aligned Length / Reference Length: "+fstr+"\n") % (totalignlen/reflen)
        uniquelen = 0
        for segs in uniqueseg.values() : # need to sum on dict entries
            util.uniqueRange(segs) #this modifies list in place
            uniquelen += util.totalLengthFromRanges( segs )
        outstr += ("Total Unique Aligned Length (Mbp)      : "+fstr+"\n") % (uniquelen*unitscale)
        outstr += ("Total Unique Aligned Length / Reference Length: "+fstr+"\n") % (uniquelen/reflen)

    return outstr

#end characterizecontigs
