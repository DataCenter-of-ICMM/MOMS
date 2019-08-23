import os

import Multithreading as mthread

"""@package PairwiseModule Defines jobs for execution of pairwise comparison

Major operation modes:
Basic distributed pairwise jobs
Distributed pairwise jobs following hash result
Distributed pairwise jobs following distributed hash resuls
Distributed pairwise jobs using triangle input partitioning
"""

import math
import utilities as util
util.setVersion("$Id: PairwiseModule.py 6695 2017-07-19 21:52:42Z wandrews $")



class Pairwise(mthread.jobWrapper):
    """Populates Multithreading package for distributed pairwise jobs
    """
    def __init__(self, varsP):
        self.varsP = varsP
        stageName = 'Pairwise'
        mthread.jobWrapper.__init__(self, varsP, stageName,clusterArgs=varsP.getClusterArgs('pairwise'))
        self.generateJobList()
    
    def runJobs(self):
        self.multiThreadRunJobs(self.varsP.nThreads, sleepTime = 0.05)
    
    #__This will not work with hashing__
    # becuase hash is not supported with partial
    #triangle mode is now always on, so no need to check this
    def generateJobListLinear(self):
        """Pairwise.generateJobListLinear: This method is the old way of doing pairwise
        comparison of all molecules. It uses the -partial option to RefAligner. This
        option is _incompatible_ with the various hashing options to RefAligner.
        """
        baseArgs = self.varsP.argsListed('noise0') + self.varsP.argsListed('pairwise') 
        ct = 0
        outputTarget = os.path.join(self.varsP.alignFolder, 'exp')

        cArgs = [self.varsP.RefAlignerBin, '-i', self.varsP.bnxFile]
        for jobNum in range(1,self.varsP.nPairwiseJobs + 1):
            jobName = 'Pairwise %d of %d' % (jobNum, self.varsP.nPairwiseJobs)
            outputString = 'pairwise%dof%d' % (jobNum,self.varsP.nPairwiseJobs)
            expectedResultFile = outputTarget + outputString + '.align'
            partialArgs = ['-partial', str(jobNum), str(self.varsP.nPairwiseJobs)]
            currentArgs = cArgs + baseArgs + ['-o' , outputTarget + outputString]
            if self.varsP.stdoutlog :
                currentArgs.extend( ['-stdout', '-stderr'] )
            if self.varsP.nPairwiseJobs > 1:
                currentArgs += partialArgs
            currentArgs += ['-maxthreads', str(self.varsP.maxthreadsPW)] # NOTE : Use -jp (defaults to -T) value for pairwise jobs
	    if self.varsP.bnxStatsFile!=None:
		currentArgs += ['-XmapStatRead', self.varsP.bnxStatsFile]
            sJob = mthread.singleJob(currentArgs, 
                                     jobName, 
                                     expectedResultFile, 
                                     outputString,
                                     maxThreads=self.varsP.nThreads,
                                     clusterLogDir=self.varsP.clusterLogDir, expectedStdoutFile=outputTarget + outputString+".stdout")
            ct += 1
            self.addJob(sJob)
        self.logArguments()
        

    def generateJobListTriangle(self):
		baseArgs = self.varsP.argsListed('noise0') + self.varsP.argsListed('pairwise')

		cArgs = [self.varsP.RefAlignerBin, '-i', self.varsP.bnxFile]
		ct = 0
		outputTarget = os.path.join(self.varsP.alignFolder, 'exp')
		njobs=self.varsP.nPairwiseJobs*(self.varsP.nPairwiseJobs+1)/2
		BNX_list=[]
		for i in range(1,self.varsP.nPairwiseJobs + 1):
			file1=self.varsP.bnxFile.replace(".bnx", "_%s_of_%s.bnx" %(i, self.varsP.nPairwiseJobs))
			BNX_list.append(file1+"\n")
			for j in range(i,self.varsP.nPairwiseJobs + 1):
				file2=self.varsP.bnxFile.replace(".bnx", "_%s_of_%s.bnx" %(j, self.varsP.nPairwiseJobs))
				jobName = 'Pairwise %d of %d' % (ct+1, njobs)
				outputString = 'pairwise%dof%d' % (ct+1, njobs)
				expectedResultFile = outputTarget + outputString + '.align'
				if i==j :
					currentArgs = [self.varsP.RefAlignerBin, '-i', file1] + ['-o' , outputTarget + outputString] + baseArgs
				else :
					currentArgs = [self.varsP.RefAlignerBin, "-first", "-1", "-i", file1, "-i", file2] + ['-o' , outputTarget + outputString] + baseArgs
                                if self.varsP.stdoutlog :
                                    currentArgs.extend( ['-stdout', '-stderr'] )
				#if self.varsP.nPairwiseJobs > 1:
					#currentArgs += partialArgs
				currentArgs += ['-maxthreads', str(self.varsP.maxthreadsPW)] # NOTE : always use all threads for pairwise
				if self.varsP.bnxStatsFile!=None:
					currentArgs += ['-XmapStatRead', self.varsP.bnxStatsFile]
				#if ct == 0: #redundant with logArguments below
				#	self.pipeReport += " ".join(currentArgs) + 2 * '\n'
				sJob = mthread.singleJob(currentArgs, 
							jobName, 
							expectedResultFile, 
							outputString,
							maxThreads=self.varsP.nThreads,
							clusterLogDir=self.varsP.clusterLogDir,
							expectedStdoutFile=outputTarget + outputString+".stdout",
							)#, shell=True)
				ct += 1
				self.addJob(sJob)
		self.varsP.bnxFileList=self.varsP.bnxFile.replace(".bnx", ".list")
		f=open(self.varsP.bnxFileList, "w")
		f.writelines(BNX_list)
		f.close()
		self.logArguments()

    def generateJobList(self):
		if self.varsP.pairwiseTriangleMode :
			self.generateJobListTriangle()
		else :
			self.generateJobListLinear()

	
    def checkResults(self):
        if self.varsP.ngsBypass : #this means that pairwise is skipped completely, so do not check anything
            return #return None means success
        self.doAllPipeReport() #loops over self.jobList and calls CheckIfFileFound
        #check for align files
        if not util.checkDir(self.varsP.alignFolder, makeIfNotExist=False) :
            self.varsP.updatePipeReport("ERROR: bad alignFolder:%s\n\n" % self.varsP.alignFolder)
            return 1

        alignFiles = []
        #for sJob in self.jobList:
        #    sJob.CheckIfFileFound()
        #    alignFile = sJob.expectedResultFile
        #    if sJob.resultFound:
        #        alignFiles.append(alignFile)
        #    else:
        #        self.warning += 1
        #        self.messages += '  PW Warning Missing Expected File: %s\n' % alignFile
        #if alignFiles.__len__() == 0:
        #    self.error += 1
        #    self.messages += '  Error: PW  Missing All Align Files\n' 
        #else:

        #Above uses results in singleJob instances, below reads from disk. Either way should work
        for ifile in os.listdir(self.varsP.alignFolder) :
            if ifile.endswith(".align") :
                alignFiles.append( os.path.join(self.varsP.alignFolder, ifile) )
        if len(alignFiles) == 0 :
            self.varsP.updatePipeReport("ERROR: no align files in alignFolder %s\n\n" % self.varsP.alignFolder)
            return 1

        alignFiles.sort()
        self.varsP.writeListToFile(alignFiles, self.varsP.alignTarget)
        self.varsP.stageComplete = 'Pairwise'
       


class sortBNX(mthread.jobWrapper) :

    def __init__(self, varsP) :
        """sortBNX.__init__: this class is for sorting the input bnx
        for subsequent splitting by the splitBNX class, and eventually
        easier processing with the Pairwise class. The constructor
        (this) will call varsP.runJobs and doAllPipeReport."""
        self.stageName="SortBNX"
        self.varsP = varsP #fewer code modifications below
        self.varsP.sorted_file = self.varsP.bnxFile.replace(".bnx", "_sorted")
        #replace this with checkMinMol; this needs to use sorted file which isn't yet made
        #calculateNPairwise(self.varsP, self.varsP.bnxFile.replace(".bnx","")) #run this here bc it contains check on N mol required to start pipeline
        checkMinMol(self.varsP, self.varsP.bnxFile)
        if self.generateJobList() : #return 0 for success, 1 for skip
            if not util.checkFile(self.varsP.sorted_file+".bnx") : #this happens when accidentally using bypass but no sorted bnx exists--log error
                err = "ERROR: no sorted bnx file found (%s) (check bypass (-B) argument to Pipeline)" % (self.varsP.sorted_file+".bnx")
                self.varsP.updatePipeReport(err+"\n")
                util.LogError("critical", err)
                util.LogStatus("progress", "pipeline", "failure")
                raise RuntimeError
            #calculateNPairwise(self.varsP, self.varsP.sorted_file) #correct varsP.nPairwiseJobs -- already above
            return
        util.LogStatus("progress", "stage_start", self.stageName) #after above bc check if bypass (executeCurrentStage)
        self.varsP.runJobs(self, "SortBNX")
        self.doAllPipeReport()
        if not self.allResultsFound() :
            err = "ERROR: sortBNX failed. Check: "+self.varsP.bnxFile
            self.varsP.updatePipeReport(err+"\n")
            util.LogError("critical", err)
            util.LogStatus("progress", "pipeline", "failure")
            raise RuntimeError
        util.LogStatus("progress", "stage_complete", self.stageName)


    def runJobs(self) :
        self.multiThreadRunJobs(1)

    def generateJobList(self) :

        if not self.varsP.executeCurrentStage:
            return 1 #tell self.__init__ not to continue processing
	    
        sorted_file = self.varsP.sorted_file

        self.varsP.updatePipeReport('Sorting %s into %s\n' % (self.varsP.bnxFile, sorted_file))
	    
        expectedResultFile=sorted_file+".bnx"
        # We use assembly section here because the memory usage is higher than pairwise, while the jobs are quite short.
        # assembly stage is full host, autoNoise0 is ~1/2 host: should not need that much, just use splitting
        super(sortBNX, self).__init__(self.varsP, self.stageName, clusterArgs=self.varsP.getClusterArgs("splitting"))

        sargs = self.varsP.argsListed('bnx_sort')
        #if scan scaling, decrease minlen so that bpp>500 does not throw away any (or at least many) molecules: 10% should be sufficient
        if self.varsP.doScanScale and '-minlen' in sargs :
            i = sargs.index('-minlen')
            minlen = util.getIntFromString(sargs[i+1])
            #this means getIntFromString succeeded: if it failed, this means RefAligner job probably will fail, too
            if minlen != None : 
                sargs[i+1] = str(int(minlen*0.9))
                err = "minlen adjusted to %s for bnx_sort because scan scaling enabled (original minlen is %i)" % (sargs[i+1],minlen)
                util.LogStatus("parameter", "minlen", err) #not an error, not really a warning either, so call it a parameter
                self.varsP.updatePipeReport(err+"\n")

        cargs=[self.varsP.RefAlignerBin, '-f', '-i', self.varsP.bnxFile, "-maxthreads", str(self.varsP.nThreads),  "-merge", "-sort-idinc", "-bnx", "-o", sorted_file] + sargs
	if self.varsP.bnxStatsFile!=None:
		cargs += ['-XmapStatWrite', self.varsP.bnxStatsFile]
        if self.varsP.stdoutlog :
            cargs.extend( ['-stdout', '-stderr'] )
        self.addJob(mthread.singleJob(cargs, self.stageName, expectedResultFile, self.stageName, clusterLogDir=self.varsP.clusterLogDir, expectedStdoutFile=sorted_file+".stdout"))

        return 0 #success


class splitBNX(mthread.jobWrapper) :

    def __init__(self, varsP, splitname="SplitBNX") :
        """splitBNX.__init__: this class is for splitting the sorted bnx file into
        smaller chunks for easier processing with the Pairwise class. Like the
        sortBNX class, the constructor also calls varsP.runJobs and doAllPipeReport.
        """
        util.LogStatus("progress", "stage_start", splitname)
        self.varsP = varsP #fewer code modifications below
        self.stageName = splitname
        if not self.generateJobList() : #check return value, and runJobs only if False
            self.varsP.runJobs(self, splitname)
        self.doAllPipeReport()
        if not self.allResultsFound() :
            err = "ERROR: splitBNX failed. Check: "+self.varsP.sorted_file+".bnx"
            self.varsP.updatePipeReport(err+"\n")
            util.LogError("critical", err)
            util.LogStatus("progress", "pipeline", "failure")
            raise RuntimeError
        util.LogStatus("progress", "stage_complete", splitname)


    def runJobs(self):
        self.multiThreadRunJobs(self.varsP.nThreads)


    def generateJobList(self) :
        """splitBNX.generateJobList: submit varsP.nPairwiseJobs number of split bnx jobs. """

        sorted_file = self.varsP.sorted_file
        if not util.checkFile(sorted_file+".bnx") :
            err = "ERROR: splitBNX input file (%s) not found; exiting" % self.varsP.sorted_file
            self.varsP.updatePipeReport(err+"\n")
            util.LogError("critical", err)
            util.LogStatus("progress", "pipeline", "failure")
            raise RuntimeError

        N = calculateNPairwise(self.varsP, sorted_file) #move back here (not sortBNX) bc needs to use sorted bnx
        #N = self.varsP.nPairwiseJobs

        self.varsP.updatePipeReport('Splitting BNX\n')
        #splitJobSet=mthread.jobWrapper(self.varsP,jobName,clusterArgs=self.varsP.getClusterArgs('splitting'))
        super(splitBNX, self).__init__(self.varsP, self.stageName, clusterArgs=self.varsP.getClusterArgs('splitting'))

        #should skip the rest and return 1, like in sortBNX, here:
        if not self.varsP.executeCurrentStage:
            return 1 #tell self.__init__ not to continue processing

        self.varsP.updatePipeReport("Splitting"+(" scan-scaled" if self.varsP.doScanScale else "")+" bnx file: %s.bnx\n\n" % self.varsP.sorted_file)

        #calculate threads per job: used to be fixed at 1, now file size / 1.5 GB rounded up. This was too low, add 1.
        threads = max(1, int(math.ceil( os.path.getsize(sorted_file+".bnx")/1.5e9 ))) + 1
        if threads > 1 :
            self.varsP.updatePipeReport("Using %i threads per job\n" % threads)

        sargs = self.varsP.argsListed('bnx_sort')
        #add minlen arg from bnx_sort in case bnx 0.1 was used and ScanScaling was NOT specified or failed
        if not self.varsP.doScanScale and '-minlen' in sargs :
            i = sargs.index('-minlen')
            minlen = util.getIntFromString(sargs[i+1])
            sargs = sargs[i:i+2] if minlen != None else [] #this means getIntFromString succeeded
        else :
            sargs = []

        #the change in job partitioning breaks backward compatibility and was causing too many problems; make it conditional on refaligner version
        #this is now obsolete: assume binaries are up-to-date
        if False : #self.varsP.refaligner_version < 3995 :
            for partial in range(1,N + 1):
                output_file=self.varsP.bnxFile.replace(".bnx", "_%s_of_%s" %(partial, self.varsP.nPairwiseJobs))
                cargs=[self.varsP.RefAlignerBin, '-f', '-i', sorted_file+".bnx", "-maxthreads", str(threads), "-merge", "-subsetbin", str(partial), str(N), "-bnx", "-o",  output_file]
                if self.varsP.stdoutlog :
                    cargs.extend( ['-stdout', '-stderr'] )
                #print('%d/%d' % (partial, N), cargs)
                expectedResultFile=output_file+".bnx"
                self.addJob(mthread.singleJob(cargs, self.stageName + str(partial), expectedResultFile, self.stageName + str(partial), maxThreads=threads, clusterLogDir=self.varsP.clusterLogDir, expectedStdoutFile=output_file+".stdout"))

        else :
            #change above to single command with -subsetbin 0 N
            output_file=self.varsP.bnxFile.replace(".bnx", "")
            cargs=[self.varsP.RefAlignerBin, '-f', '-i', sorted_file+".bnx", "-maxthreads", str(threads), "-merge", "-subsetbin", "0", str(N), "-bnx", "-o",  output_file]
            if self.varsP.stdoutlog :
                cargs.extend( ['-stdout', '-stderr'] )
            self.addJob(mthread.singleJob(cargs, self.stageName, output_file+".bnx", self.stageName, maxThreads=threads, clusterLogDir=self.varsP.clusterLogDir, expectedStdoutFile=output_file+".stdout"))
    #end generateJobList

#end class splitBNX


def checkMinMol(varsP, input_file, minmol=2) :
    '''Simplified version of calculateNPairwise which just checks that there are at least minmol molecules.'''
    f=open(input_file, "r")
    count=0
    #length=0
    #site_count=0.0
    for line in f:
        if line[0] == "0":
            #x=line.split()
            count+=1
            #length+=float(x[2])
        #if line[0] == "1":
        #site_count+=len(line.split())-1
        if count > minmol : #this is all we need to check
            break
    f.close()

    #check that we have more than 1 molecule; if not, there's nothing to assemble, so exit
    if count < minmol :
        err = "ERROR in calculateNPairwise: number of molecules (%i) is too few for assembly; check bnx: %s" % (count, input_file)
        varsP.updatePipeReport(err+"\n")
        util.LogError("critical", err)
        util.LogStatus("progress", "pipeline", "failure")
        raise RuntimeError #will be caught in DNPipeline.constructData
#end checkMinMol


def calculateNPairwise(varsP, sorted_file) :
    '''Given varsPipeline object and prefix of sorted bnx file, check if N pairwise jobs should be increased based on heuristic of 80k molecules per bnx file.'''
    f=open(sorted_file+".bnx", "r")
    count=0
    length=0
    site_count=0.0
    for line in f:
        if line[0] == "0":
            x=line.split()
            count+=1
            length+=float(x[2])
        if line[0] == "1":
	    site_count+=len(line.split())-1
    f.close()

    nAuto=int(math.ceil(count/100000.0))
    nAutoSites=int(math.ceil(site_count/3e6))
    if nAutoSites>nAuto:
	nAuto=nAutoSites

    # round up to nearest multiple of 6 (since there are 6 MIC cards per host), for better load balancing during stage refineB0 etc
    if nAuto > 6 : #but if fewer than 7, this isn't necessary
        nAuto = (((nAuto + 5) / 6) * 6)

    varsP.bnxNmol     = count #store these for future use
    varsP.bnxTotLenbp = length 
    varsP.updatePipeReport("calculateNPairwise: Sorted and filtered file has %d molecules, %f sites %f length total. Auto blocks=%d\n\n" %(count, site_count, length, nAuto))
        
    if nAuto > varsP.nPairwiseJobs:
        varsP.updatePipeReport("Increasing number of blocks from %d to %d\n" % (varsP.nPairwiseJobs, nAuto))
        varsP.nPairwiseJobs=nAuto

    return varsP.nPairwiseJobs
