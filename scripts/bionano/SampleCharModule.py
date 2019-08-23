import os
#import sys
import shutil
#from collections import OrderedDict

import utilities as util
import Multithreading as mthread

"""@package SampleCharModule Defines jobs to map bnx to reference and parses
mapping results

Optional pre-processor for De Novo assembly to characterize the single molecule
noise of the input. Reference required. Can help to inform run time noise 
arguments.
"""


util.setVersion("$Id: SampleCharModule.py 5697 2016-12-09 23:42:42Z wandrews $")


class SampleChar(mthread.jobWrapper):
    """Generates jobs for distributed single molecule mapping
    """
    def __init__(self, varsP):
        self.varsP = varsP
        stageName = '  Sample Characterization, Noise Levels'
        super(SampleChar, self).__init__(self.varsP, stageName, clusterArgs=varsP.getClusterArgs('sampleChar'))
        self.generateJobList()
    
    def runJobs(self):
        self.multiThreadRunJobs(self.varsP.nThreads, sleepTime = 0.2)
    
    def generateJobList(self):
        curArgs = self.varsP.argsListed('noise0') + self.varsP.argsListed('sampleChar')
        if util.checkFile(self.varsP.bnxTarget) : #file exists only if image processing was run
            bnxFiles = parseExperimentFile(self.varsP.bnxTarget)
            if not bnxFiles : #check that you got at least one
                errstr = "ERROR in SampleChar.generateJobList: no bnx files found in: "+self.varsP.bnxTarget
                print errstr
                self.varsP.updatePipeReport(errstr+"\n\n")
                return
            basepath = "" #os.path.split(bnxFiles[0])[0] #don't use basepath for this case
        else : #otherwise, assume this is the only bnx file
            bnxFiles = [self.varsP.bnxFile]
            #here, make a dir for the results--should really check results of checkEmptyDir for errors
            basepath = os.path.join(self.varsP.localRoot, "sampleChar")
            if self.varsP.wipe and os.path.isdir(basepath) :
                shutil.rmtree(basepath)
                #util.checkEmptyDir(basepath) #will make if not exist, but if it does, will remove and re-make -- this fn doesn't exist...
            #else :
            util.checkDir(basepath) #will make if not exist, but won't remove anything
        nJobs = len(bnxFiles)
        #for i, bnxFile in enumerate(bnxFiles):
        for bnxFile in bnxFiles :
            #bnxGroupName = '%02d' % (i+1) #get this from the path, ie, bnxFiles
            cargs = [self.varsP.RefAlignerBin, '-i', bnxFile]
            bnxname = os.path.split(bnxFile)[1].replace(".bnx","")
            jobname = 'Sample_Char_' + bnxname
            #outputTarget = os.path.join(basepath, bnxGroupName)
            if basepath : #bnx input
                outputTarget = os.path.join(basepath, bnxname)
            else : #image processing
                outputTarget = bnxFile.replace(".bnx","") + "_sampleChar"
            expectedResultFile = outputTarget + '.err' #this is used in checkResults
            currentArgs = cargs + ['-ref', self.varsP.ref, '-o' , outputTarget, '-f']
            if self.varsP.stdoutlog :
                currentArgs.extend( ['-stdout', '-stderr'] )
            currentArgs += ['-maxthreads', str(self.varsP.maxthreads)] + curArgs
            sJob = mthread.singleJob(currentArgs, jobname, expectedResultFile, jobname, clusterLogDir=self.varsP.clusterLogDir) # peStr is deprecated in favor of clusterargs
            #sJob.expTag = bnxGroupName #removed from checkResults
            self.addJob(sJob)
        self.logArguments()

    
    def checkResults(self):
        self.doAllPipeReport()
        self.infoReport = '  Sample Characterization:\n'
        for i,sJob in enumerate(self.jobList):
            if not sJob.resultFound :
                continue
            errSum = parseErrFile(sJob.expectedResultFile)
            if i == 0:
                self.infoReport += '    ExpID  ' + errSum.makeHeader() + '\n'
            #self.infoReport += '    Exp ' + sJob.expTag + ' ' + errSum.makeReport() + '\n' #no longer use expTag
            self.infoReport += sJob.jobName + ' ' + errSum.makeReport() + '\n'
        self.varsP.updateInfoReport(self.infoReport)
        

    
def parseErrFile(targetFile):
    f1 = open(targetFile)
    errResults = []
    while(True):
        line = f1.readline()
        if line == '': 
            f1.close()
            break
        if line[0] == 'I' or line[0] == '#':
            continue
        else:
            errResult = errFileResult(line)
            errResults.append(errResult)
    errSum = errFileSum(errResults)
    return errSum

    
class errFileSum():
    
    def __init__(self, errFileResults):
        self.FP = errFileResults[-1].FP
        self.FN = errFileResults[-1].FN
        self.ssd = errFileResults[-1].ssd
        self.sd = errFileResults[-1].sd
        self.bpp = errFileResults[-1].bpp
        self.nMaps = errFileResults[0].nMaps
        self.mapRate = float(errFileResults[-2].nMapsMapped) / errFileResults[-2].nMaps
    
    def makeHeader(self):
        output = '% 8s% 8s% 8s% 8s% 8s% 8s% 8s' %('NumMaps', 'MapRate', 'FP', 'FN', 'bpp', 'sd', 'ssd')
        return output
        
    def makeReport(self):
        #output = '% 8d ' % self.nMaps 
        #output += '% 8.3f   ' % self.mapRate
        #output += '% 8.3f   ' % self.FP
        #output += '% 8.4f   ' % self.FN
        #output += '% 8.1f   ' % self.bpp
        #output += '% 8.4f   ' % self.sd
        #output += '% 8.4f   ' % self.ssd
        output = '% 8d% 8.3f% 8.3f% 8.4f% 8.1f% 8.4f% 8.4f' % (self.nMaps, self.mapRate, self.FP, self.FN, self.bpp, self.sd, self.ssd)
        
        return output
        
class errFileResult():

    def __init__(self, line):
        tokens = [float(x) for x in line.strip().split('\t')]
        self.FP = tokens[1]
        self.FN = tokens[2]
        self.ssd = tokens[3]
        self.sd = tokens[4]
        self.bpp = tokens[5]
        self.nMaps = int(tokens[7])
        self.nMapsMapped = int(tokens[9])    


def groupBnxByExperimentMod(bnxFiles):
    prevExp = ''
    allGroups = []
    allGroupNames = []
    curGroup = []
    for bnxPath in bnxFiles:
        dump, bnxFile = os.path.split(bnxPath)
        allGroupNames.append(bnxFile[:4])
        allGroups.append([bnxPath])
    return allGroups, allGroupNames
    
def groupBnxByExperiment(bnxFiles):
    prevExp = ''
    allGroups = []
    allGroupNames = []
    curGroup = []
    for bnxPath in bnxFiles:
        dump, bnxFile = os.path.split(bnxPath)
        curExp = bnxFile[:2]
        if curExp == prevExp:
            curGroup.append(bnxPath)
        else:
            if curGroup.__len__() > 0:
                allGroups.append(curGroup)
            curGroup = [bnxPath]
            allGroupNames.append(curExp)
            prevExp = curExp
    if curGroup.__len__() > 0:
        allGroups.append(curGroup)            
    return allGroups, allGroupNames    


#open file targetFile, read bnx paths from it
#targetFile is output of Pipeline.writeIntermediate
def parseExperimentFile(targetFile):
    targetLocations = []
    f1 = open(targetFile)
    while(True):
        line = f1.readline()
        if line == '':
            f1.close()
            break
        if line[0] == '#':
            continue
        targetLocations.append(line.strip())
    return targetLocations

class autoNoise(mthread.jobWrapper) :

    def __init__(self, varsP) :
        """splitBNX.__init__: this class is for sorting the input bnx
        for subsequent splitting by the splitBNX class, and eventually
        easier processing with the Pairwise class. The constructor
        (this) will call varsP.runJobs and doAllPipeReport, then
        instantiate splitBNX, which will do all the splitting required
        for the Pairwise class.
        """
        self.stageName = "Autonoise0"
        self.varsP = varsP #fewer code modifications below
        
        util.LogStatus("progress", "stage_start", self.stageName) #after above bc check if bypass (executeCurrentStage)

        self.output_folder = os.path.join(self.varsP.contigFolder, "auto_noise")
        if not util.checkDir(self.output_folder) : #will make if not exist, only returns False if already exists or can't make
            print "ERROR in autoNoise: bad dir:", self.output_folder
            raise RuntimeError
	    
        # We use assembly section here because the memory usage is higher than pairwise, while the jobs are quite short.
        super(autoNoise, self).__init__(self.varsP, self.stageName, clusterArgs=self.varsP.getClusterArgs("autoNoise0"))

        bnxfile = self.varsP.bnxFile if varsP.noiseOnly else self.varsP.sorted_file+".bnx"
        #was return if generateJobListChar, but need to get readparameters if bypass
        if not self.generateJobListChar({}, bnxfile, "autoNoise0") : #return 0 for success, 1 for skip
            self.varsP.runJobs(self, "AutoNoise0")
            self.doAllPipeReport()
        if not self.allResultsFound() :
            errstr = "AutoNoise0 failed. Check: "+self.output_file+".stdout"
            self.varsP.updatePipeReport("ERROR: "+errstr+"\n")
            util.LogError("critical", errstr)
            raise RuntimeError
        util.LogStatus("progress", "stage_complete", self.stageName)

        self.varsP.noise0 = readNoiseParameters(self.output_file)
	#self.isBadErrorParams(self.varsP.noise0, 0) #DISABLE check for autoNoise0

        self.stageName = "Autonoise1"
        self.groupName = self.stageName #fix so that LogStatus call in MultiThreading.multiThreadRunJobs
        util.LogStatus("progress", "stage_start", self.stageName)

        self.clearJobs()
        
	self.varsP.replaceParam("noise0", "-readparameters", self.output_file+".errbin")

        #need to call again to set self.output_file
        if not self.generateJobListChar(self.varsP.noise0, bnxfile, "autoNoise1") : #return 0 for success, 1 for skip
            self.varsP.runJobs(self, "AutoNoise1")
            self.doAllPipeReport()
        if not self.allResultsFound() :
            errstr = "AutoNoise1 failed. Check: "+self.output_file+".stdout"
            self.varsP.updatePipeReport("ERROR: "+errstr+"\n")
            util.LogError("critical", errstr)
            raise RuntimeError
            
        self.varsP.noise1 = readNoiseParameters(self.output_file)
        
	infoReport="Automatically determined noise parameters:\n"
        #klist = ["FP", "FN", "sf", "sd", "sr", "bpp", "readparameters"] #hardcoding parameters is kind of bad, but it fixes the order without using OrderedDict.
        #for v in klist :
        for v in sorted(self.varsP.noise1.keys()) :
            #if not self.varsP.noise1.has_key(v) :
            #    continue
            param=str(self.varsP.noise1[v])
            util.LogStatus("parameter", "auto_"+v, param)
            infoReport+=v+":"+param+"\n"
            self.varsP.replaceParam("noise0", "-"+v, param)
        self.varsP.updateInfoReport(infoReport + '\n', simple=True)
        self.isBadErrorParams(self.varsP.noise1, 1)

        if self.varsP.doScanScale : #change the sorted_file to the rescaled bnx file
            rescaledbnx = self.output_file + self.varsP.rescaleSuffix #no ".bnx" in suffix
            if not util.checkFile(rescaledbnx+".bnx") : #not found--not an error if bnx 0.1 is used
                err = "Warning: scan scaled bnx not found after autoNoise1; not performing scan scaling--check that bnx 1.0 or later used in input"
                self.varsP.updatePipeReport( err+"\n\n" )
                util.LogError("warning", err)
                self.varsP.doScanScale = False
            else : #log that scan scaling is used
                self.varsP.updatePipeReport( "Using scan scaled bnx: "+rescaledbnx+".bnx\n\n" )
                util.LogStatus("parameter", "scanscaled_bnx", rescaledbnx+".bnx")
                self.varsP.sorted_file = rescaledbnx #this variable is used in splitBNX (PairwiseModule.py)
            
        util.LogStatus("progress", "stage_complete", self.stageName)


    def isBadErrorParams(self, noise, stage):
        #BAD means this:
        # for both stages: sr > 0.1 or sd > 0.1 or sf > 0.5
        # also for stage 0 : sd > 0.1 and sf > 0.35
        # also for stage 1 : sd > 0 and sf > 0.25 (this used to be for both stages)
        assert stage == 0 or stage == 1, "Error: invalid arg to autoNoise.isBadErrorParams"
        badparam = False
        if not noise :
            badparam = True
        elif stage == 0 and (noise["sd"] > 0.1 and noise["sf"] > 0.35) :
            badparam = True
        elif stage == 1 and (noise["sd"] > 0   and noise["sf"] > 0.25) :
            badparam = True
        #add not noise for case of empty dict, which readNoiseParameters will return if it can't read the .err file
        if badparam or noise["sr"] > 0.1 or noise["sd"] > 0.1 or noise["sf"] > 0.5 :
            errstr = "Failed to find usable noise parameters. Try decreasing maprate parameter and/or find a better reference. You can also try disabling auto noise (no -y, or 'Rough assembly' profile) with nominal noise parameters;"
            if noise.has_key("sf") :
                errstr += " sf=%f" % noise["sf"]
            if noise.has_key("sd") :
                errstr += " sd=%f" % noise["sd"]
            if noise.has_key("sr") :
                errstr += " sr=%f" % noise["sr"]
            self.varsP.updatePipeReport(errstr+"\n")
            util.LogError("critical", errstr)
            #util.LogStatus("progress", "pipeline", "failure") #redundant with DNPipeline.finalizePipeline
            raise RuntimeError

    def runJobs(self) :
        self.multiThreadRunJobs(1)

    def generateJobListChar(self, noise_in, input_file, optSection) :
                    
        self.output_file=os.path.join(self.output_folder, optSection) #must assign before return bc used in constructor

        if not self.varsP.executeCurrentStage:
            return 1 #tell self.__init__ not to continue processing
	    
        self.varsP.updatePipeReport('%s\n' % (optSection))
        
        expectedResultFile=self.output_file+".err"

        #cargs=[self.varsP.RefAlignerBin, '-f', '-i', input_file, "-ref", self.varsP.ref, "-maxthreads", str(self.varsP.maxthreads), "-o", self.output_file] 
        cargs=[self.varsP.RefAlignerBin, '-f', '-i', input_file, "-ref", self.varsP.ref, "-o", self.output_file] #remove maxthreads bc this is always running on its own
        if self.varsP.stdoutlog :
            cargs.extend( ['-stdout', '-stderr'] )
        cargs.extend( ['-output-veto-filter', 'intervals.txt$'] )
        for v in noise_in.keys():
		cargs.extend(["-"+v, str(noise_in[v])])
		
        cargs.extend(self.varsP.argsListed(optSection))
	if self.varsP.bnxStatsFile!=None:
		cargs += ['-XmapStatWrite', self.varsP.bnxStatsFile]
        self.addJob(mthread.singleJob(cargs, self.stageName, expectedResultFile, self.stageName, clusterLogDir=self.varsP.clusterLogDir, expectedStdoutFile=self.output_file+".stdout"))

        return 0 #success

#end class autoNoise

def readNoiseParameters(input_file):
    """Argument input_file must be path plus prefix to .err file. Return noise parameters from its last iteration."""
    #print "readNoiseParameters:", input_file #debug
    raw_line=""
    try:
        with open(input_file+".err") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if line.strip()=="":
                    continue
                raw_line=line
    except:
        err = "readNoiseParameters: error reading noise parameters from %s" % (str(input_file)+".err")
        print "error:", err
        util.LogError("error", err)
        return({})
		
    data=[]
    for x in raw_line.split("\t"):
        y=x.strip()
        if y!="":
            data.append(float(y))
    #print data
    if len(data) > 13 : #be sure that we don't get IndexError
        return(dict([("FP",data[1]), ("FN",data[2]), ("sf",data[3]), ("sd",data[4]), ("bpp",data[5]), ("res",data[6]), ("sr",data[13]), ("readparameters",input_file+".errbin")]))
    elif len(data) > 5 : #assume sr is missing in this case, but rest are correct
        return(dict([("FP",data[1]), ("FN",data[2]), ("sf",data[3]), ("sd",data[4]), ("bpp",data[5]), ("res",data[6]), ("readparameters",input_file+".errbin")]))
    else : #bad input_file, data is missing and cannot be trusted
        err = "readNoiseParameters: error parsing noise parameters from %s" % (input_file+".err")
        print "error:", err
        util.LogError("error", err)
        return({})
#end readNoiseParameters
