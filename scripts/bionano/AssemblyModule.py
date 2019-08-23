#import subprocess
#import time
import os
#import pdb

import Multithreading as mthread

"""@package AssemblyModule For defining assembly worker (singleJob)

"""


import utilities
utilities.setVersion("$Id: AssemblyModule.py 3928 2015-07-08 01:17:25Z vdergachev $")


class Assemble(mthread.jobWrapper):
    """ Class to run the assembly phase of de novo assembly
    
    """
    
    def __init__(self, varsP):
        self.varsP = varsP
        stageName = 'Assembly'
        mthread.jobWrapper.__init__(self, varsP, stageName, clusterArgs=varsP.getClusterArgs('assembly'))
        self.contigString = varsP.expID + '_unrefined'
        varsP.prepareContigIO(self.contigString, stageName)
        self.varsP.stageName=stageName
        self.contigsFile = None
        self.generateJobList()
    
    def runJobs(self):
        self.multiThreadRunJobs(self.varsP.nThreads, sleepTime = 0.2)
    
    def generateJobList(self):
        """ Instantiate job wrapper class with queue of single jobs for assembly
        
        """
        
        if self.varsP.pairwiseTriangleMode :
            AssemblerInputFlag="-if"
            AssemblerInputFile=self.varsP.bnxFileList
        else:
            AssemblerInputFlag="-i"
            AssemblerInputFile=self.varsP.bnxFile
        cargs = [self.varsP.AssemblerBin, AssemblerInputFlag, AssemblerInputFile, '-af', self.varsP.alignTarget]
	if self.varsP.bnxStatsFile!=None:
		cargs += ['-XmapStatRead', self.varsP.bnxStatsFile]
        baseArgs = self.varsP.argsListed('noise0') + self.varsP.argsListed('assembly')
        
        logFile = os.path.join(self.varsP.localRoot, 'AssemblyLog.txt')
        errFile = os.path.join(self.varsP.localRoot, 'AssemblyLog_stderr.txt')
        
        outFile = os.path.join(self.varsP.outputContigFolder, self.varsP.outputContigPrefix) #no suffix for -o arg of Assembler
        self.contigsFile = os.path.join(self.varsP.outputContigFolder, self.varsP.outputContigPrefix+".contigs") #Assembler will append this suffix
        currentArgs = cargs + baseArgs + ['-o', outFile]
        if self.varsP.stdoutlog :
            currentArgs.extend( ['-stdout', '-stderr'] )
        logArguments = "   ".join(currentArgs) + 2 * '\n'
        jobName = 'Assembly'
        #sJob = mthread.singleJob(currentArgs, jobName, self.contigsFile, jobName, stdOutFile=logFile, stdErrOutFile=errFile,clusterLogDir=self.varsP.clusterLogDir, expectedStdoutFile=outFile+".stdout")
        sJob = mthread.singleJob(currentArgs, jobName, self.contigsFile, jobName, clusterLogDir=self.varsP.clusterLogDir, expectedStdoutFile=outFile+".stdout")
        self.addJob(sJob)
        self.logArguments()
        #self.varsP.curContigPrefix = self.varsP.expID + '_unrefined' #this is no longer needed (superceded by prepareContigIO)
        #assemblyJobSet.addJob(sJob)
    
    def checkResults(self):
        self.doAllPipeReport()
        #for sJob in self.jobList:
        #    sJob.CheckIfFileFound()
        #    if sJob.resultFound:
        #        pass
        #    else:
        #        self.error += 1
        #        self.messages += '  Error: ASSMBL  Missing Contigs File\n'
        #self.varsP.error += self.error
        #self.varsP.message += self.messages
        self.varsP.stageComplete = 'Assembly'
        self.varsP.mergeIntoSingleCmap()
        #self.pipeReport += self.makeRunReport()
        #self.pipeReport += self.makeParseReport()
        #self.varsP.updatePipeReport(self.pipeReport)

     

