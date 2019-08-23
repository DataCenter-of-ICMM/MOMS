import os

import Multithreading as mthread

"""@package RefinementModule Defines jobs for refinement, extension and merging 
operations


"""


import utilities
utilities.setVersion("$Id: RefinementModule.py 5905 2017-01-23 19:54:17Z tanantharaman $")


#this class replaces the old RefineA and RefineB classes by merging them
#see comment above init
class Refine(mthread.jobWrapper):
    """refineStage is a string specifying what stage of refinement you're in
    must be 'refineA', 'refineB', 'refineNGS', or 'refineFinal' (later)
    """
    def __init__(self, refineStage, varsP):
        validstages = ['refineA', 'refineB', 'refineNGS', 'refineFinal']
        if not refineStage in validstages :
            varsP.error += 1
            varsP.message += '  Error: Refine stage name invalid: '+str(refineStage)+'\n'
            return
        self.refineStage = refineStage
        self.varsP = varsP
        utilities.LogStatus("progress", "stage_start", self.refineStage)
        #super is more pythonic than referring to the base class explicitly (only matters for multiple inheritance)
        super(Refine, self).__init__(varsP, refineStage, clusterArgs=varsP.getClusterArgs(refineStage))
        intermediateContigPrefix = self.varsP.expID + self.refineStage.replace("refine", "_r")
        self.varsP.prepareContigIO(intermediateContigPrefix, refineStage)
        #modify results of varsP.prepareContigIO for special case of refineNGS
        if self.refineStage == 'refineNGS' :
            self.varsP.inputContigPrefix = self.varsP.ngsContigPrefix
            self.varsP.inputContigFolder = self.varsP.ngsInDir
        self.generateJobList()
        
    def runJobs(self):
        self.multiThreadRunJobs(self.varsP.nThreads, sleepTime = 0.05)
        
    def writeIDFile(self, nJobs):
        f1 = open(self.varsP.idFile, 'w')
        f1.write(str(nJobs))
        f1.close()

    def generateJobList(self):
        baseArgs1 = self.varsP.argsListed(self.refineStage)
        if self.refineStage != 'refineNGS' : #noise args are in refineNGS
            baseArgs1 += self.varsP.argsListed('noise0')
        contigFiles, contigIDs = self.varsP.findContigs(self.varsP.inputContigFolder, self.varsP.inputContigPrefix)
        #nJobs = len(contigFiles)
        bnx = self.varsP.sorted_file+".bnx" #was self.varsP.bnxFile, but need sorted bc ids are different after sorting
        if self.refineStage == 'refineA' : #refineA uses assembler, all others use refaligner
            r1args = [self.varsP.AssemblerBin, '-i', bnx] #need this before -contigs
            r1args += ['-contigs', os.path.join(self.varsP.inputContigFolder, self.varsP.inputContigPrefix) + '.contigs']
        else : #should be same for refineB/NGS/Final
            r1args = [self.varsP.RefAlignerBin, '-i', bnx]
            self.writeIDFile(len(contigFiles)) #nJobs) 
        output1String = os.path.join(self.varsP.outputContigFolder, self.varsP.outputContigPrefix)
        for contigID in contigIDs :
            expectedOutputString = self.varsP.outputContigPrefix + '_contig' + contigID
            expectedResultFile = os.path.join(self.varsP.outputContigFolder, expectedOutputString + '.cmap') #refineB
            jobName = self.refineStage + ' %5s' % contigID
            if self.refineStage == 'refineA' : 
                currentArgs = 2*[str(contigID)] #this must come after r1args because it's actually an argument to -contigs
            else : #should be same for refineB/NGS/Final
                r1_cmapFile = self.varsP.inputContigPrefix + '_contig' + str(contigID) + '.cmap'
                r1_cmapFile = os.path.join(self.varsP.inputContigFolder, r1_cmapFile)
                currentArgs = ['-maxthreads', str(self.varsP.maxthreads), '-ref', r1_cmapFile, '-id', contigID]
            currentArgs = r1args + ['-o', output1String] + currentArgs + baseArgs1
            if self.varsP.stdoutlog :
                currentArgs.extend( ['-stdout', '-stderr'] )
            s1Job = mthread.singleJob(currentArgs, 
                                    jobName, 
                                    expectedResultFile, 
                                    expectedOutputString,
                                    maxThreads=self.varsP.maxthreads,
                                    clusterLogDir=self.varsP.clusterLogDir)
            self.addJob(s1Job)
        self.logArguments()
    
    def checkResults(self):
        self.varsP.stageComplete = self.refineStage
        self.varsP.mergeIntoSingleCmap()
        self.doAllPipeReport() #see Multithreading.jobWrapper
        utilities.LogStatus("progress", "stage_complete", self.refineStage)

    def endStage(self): #same as GroupedRefinementModule.Refine.endStage
        utilities.LogStatus("progress", "stage_complete", self.refineStage)

#end class Refine


class Extension(mthread.jobWrapper):
    """Generate jobs for extension phase of assembly
    """
    def __init__(self, varsP):
        self.varsP = varsP
        self.varsP.extensionCount += 1
        self.stageName = 'Extension_'+str(self.varsP.extensionCount)
        utilities.LogStatus("progress", "stage_start", self.stageName)
        mthread.jobWrapper.__init__(self,varsP, self.stageName,clusterArgs=varsP.getClusterArgs('extension'))
        extContigPrefix = self.varsP.expID + '_ext%s' % self.varsP.extensionCount
        varsP.prepareContigIO(extContigPrefix, self.stageName)
        self.generateJobList()
        
    def runJobs(self):
        self.multiThreadRunJobs(self.varsP.nThreads, sleepTime = 0.05)
        
    def generateJobList(self):
        contigFiles, contigIDs = self.varsP.findContigs(self.varsP.inputContigFolder, self.varsP.inputContigPrefix)
        curargs = [self.varsP.RefAlignerBin, '-i', self.varsP.sorted_file+".bnx"] #was bnxFile
        baseArgs = self.varsP.argsListed('noise0') + self.varsP.argsListed('extension')
        nJobs = contigFiles.__len__()    
        ct = 0
        logArguments = "" #just in case the following loop isn't entered
        for jobNum in range(1,nJobs + 1):
            contigID = contigIDs[jobNum - 1]
            #jobName = 'Extend ' + contigID + ', Job ' + str(jobNum) + ' of ' + str(nJobs)
            expContigString = self.varsP.outputContigPrefix + '_contig' + contigID
            outputString = os.path.join(self.varsP.outputContigFolder, self.varsP.outputContigPrefix)
            expectedResultFile = os.path.join(self.varsP.outputContigFolder, expContigString + '.cmap')# '_refined.cmap')
            jobName = 'Ext %s' % expContigString# + ', Job ' + str(jobNum) + ' of ' + str(nJobs)
            currentContig = contigFiles[jobNum - 1]
            currentArgs = curargs + ['-o', outputString] + baseArgs 
            currentArgs += ['-maxthreads', str(self.varsP.maxthreads), '-id', contigID, '-ref', currentContig]
            if self.varsP.stdoutlog :
                currentArgs.extend( ['-stdout', '-stderr'] )
            s1Job = mthread.singleJob(currentArgs, 
                                        jobName, 
                                        expectedResultFile, 
                                        expContigString,
                                        maxThreads=self.varsP.maxthreads, 
                                        forceForward = currentContig, 
                                        clusterLogDir=self.varsP.clusterLogDir)
            self.addJob(s1Job)
            ct += 1
        self.logArguments()
    
    def checkResults(self):
        self.varsP.stageComplete = 'Extension% 2d' % self.varsP.extensionCount
        self.varsP.mergeIntoSingleCmap()
        self.doAllPipeReport() #see Multithreading.jobWrapper
        utilities.LogStatus("progress", "stage_complete", self.stageName)


        
class Merge(mthread.jobWrapper):
    """Generate jobs for merge step. 
    
    Contigs can only merge with one other contig at a time.
    External function will call execution of this class iteratively to exhaust
    all merging possibilities
    """
    def __init__(self, varsP, stagename):
        self.varsP = varsP
        self.stageName = stagename
        self.varsP.stageName = stagename
        mthread.jobWrapper.__init__(self,varsP, self.stageName,clusterArgs=varsP.getClusterArgs('merge'))
        self.alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
        #move all of this into resetStage/generateJobList
        self.resetStage(self.stageName)
        #self.iterCount = 0
        #contigTag = '_mrg%d' % self.varsP.extensionCount
        #initialPrefix = self.varsP.expID + contigTag
        #varsP.prepareContigIO(initialPrefix)
        #self.prevPrefix = None
        #self.curPrefix = None

    def resetStage(self, stagename) :
        self.stagePrefix = self.varsP.expID + ('_mrg%d' % self.varsP.extensionCount)
        self.varsP.prepareContigIO(self.stagePrefix, stagename)
        self.iterCount = 0
        self.stageName = stagename 
        self.varsP.stageName = stagename
        self.groupName = self.stageName + self.alphabet[self.iterCount] #jobWrapper data member

    def runJobs(self):
        self.multiThreadRunJobs(self.varsP.nThreads, sleepTime = 0.05)
    
    def countContigs(self,targetFolder,prefix): 
        return len(self.varsP.findContigs(targetFolder, prefix)[0])
        
    def generateJobList(self):
        """Defines job parameters for merge. Updates variables for subsequent
        completion test in mergeComplete()
        """
        self.clearJobs()
        self.prevPrefix = self.varsP.inputContigPrefix
        self.curPrefix = self.stagePrefix + self.alphabet[self.iterCount]
        self.groupName = self.stageName + self.alphabet[self.iterCount] #jobWrapper data member
        utilities.LogStatus("progress", "stage_start", self.groupName)
        self.varsP.updatePipeReport('   PREV PREFIX %s, CUR PREFIX %s' % (self.prevPrefix, self.curPrefix))
        self.iterCount += 1
        outputString = os.path.join(self.varsP.outputContigFolder, self.curPrefix)
        currentArgs = [self.varsP.RefAlignerBin, '-o', outputString]
        #if self.varsP.stdoutlog : #always use this here bc it's the only output which should always be there
        currentArgs.extend( ['-f', '-stdout', '-stderr'] )
        currentArgs += self.varsP.argsListed('merge') 
        if self.varsP.nExtensionIter == self.varsP.extensionCount and self.varsP.argData.has_key('lastMerge') :
            currentArgs += self.varsP.argsListed('lastMerge')             
        currentArgs += ['-maxthreads', str(self.varsP.nThreads)]
        contigsTextFile = os.path.join(self.varsP.inputContigFolder, 'mergeContigs.txt')
        contigFiles, contigIDs = self.varsP.findContigs(self.varsP.inputContigFolder, self.prevPrefix, txtOutput=contigsTextFile) #this method creates the mergeContigs.txt file which is necessary for this job
        self.varsP.prefixUsed.append(self.curPrefix)
        fileArgs = ['-if', contigsTextFile]
        #expoutput = outputString+".align" #don't know which contigs will disappear, but should always get an align file -- with new arg 'pairmergeRepeat', there's no .align; use stdout
        expoutput = outputString+".stdout"
        s1Job = mthread.singleJob(currentArgs + fileArgs, self.groupName, expoutput, self.groupName, maxThreads=self.varsP.nThreads, clusterLogDir=self.varsP.clusterLogDir, expectedStdoutFile = outputString + ".stdout")
        self.addJob(s1Job)
        self.logArguments() 
            
    def mergeComplete(self):
        """Test if merge possibilities are exhaused and increment names and counters.
        If RefAligner argument -pairmergeRepeat is used, always terminate.
        """
        prevCount = self.countContigs(self.varsP.inputContigFolder, self.prevPrefix)
        curCount = self.countContigs(self.varsP.outputContigFolder, self.curPrefix)
        #self.varsP.stageComplete = 'Merge% 2d' % self.varsP.extensionCount
        self.varsP.stageComplete = self.stageName
        self.checkResults()
        contigCount = '  %s %d to %s %d' % (self.prevPrefix, prevCount, self.curPrefix, curCount)
        self.varsP.inputContigPrefix = self.curPrefix
        self.varsP.inputContigFolder = self.varsP.outputContigFolder
        self.varsP.outputContigPrefix = self.curPrefix
        utilities.LogStatus("progress", "stage_complete", self.groupName) #a stage is each merge (A, B, etc), not all
        term = "-pairmergeRepeat" in self.varsP.argsListed('merge') 
        if term or curCount <= 1 or curCount >= prevCount or self.iterCount >= len(self.alphabet)-1:
            # Terminate Merging
            contigCount += '  .. Terminate Merge ..'
            self.varsP.updatePipeReport(contigCount + '\n')
            if curCount == 0:
                self.varsP.outputContigPrefix = self.prevPrefix
            self.varsP.mergeIntoSingleCmap()
            return 1
        else:
            contigCount += '  .. Continue Merge ..'
            self.varsP.updatePipeReport(contigCount + '\n')
            return 0
    
    def checkResults(self):
        self.doAllPipeReport() #see Multithreading.jobWrapper
