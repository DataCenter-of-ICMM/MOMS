#import subprocess
#import time
import os
#import shutil
#import pdb
import utilities as util
from math import *
from string import Template

import Multithreading as mthread

"""@package GroupedRefinementModule
Package runs refinement phases: RefineA, RefineB, refineNGS, and refineFinal.
Grouped refinement splits refinement into two parts in order to increase CPU
efficiency. First part does alignment of molecule maps to contigs, and second
part does refinement.

"""


util.setVersion("$Id: GroupedRefinementModule.py 6583 2017-06-28 20:17:17Z tanantharaman $")

           
#see comment above init
class Refine(mthread.jobWrapper):
    """Replaces the old RefineA and RefineB classes by merging them
    
    Combines RefineA, RefineB, refineNGS, and refineFinal
    """
    #refineStage is a string specifying what stage of refinement you're in
    #must be 'refineA', 'refineB', 'refineNGS', or 'refineFinal' (later)
    def __init__(self, StageName, varsP):
        self.refineStage = StageName
        self.multigroup = True #if False, force single group (not normally good)
        self.groupalgo = 1 #0 is old algo; 1 is new
        self.varsP = varsP
        ContigPrefix = self.varsP.expID + "_" + StageName

        if StageName=="extension0":
		self.varsP.extensionCount += 1
		
	for case in util.switch(StageName):
		if case("refine(B0|B1|Final0|Final1)", regexp=True):
			self.bunching=12
			self.ref_arg="-reff"
			break
		if case("refineA"):
			self.bunching=12
			self.ref_arg="-ref"
			break
		if case("refineNGS"):
			self.bunching=1
			self.ref_arg="-ref"
			self.varsP.inputContigPrefix = self.varsP.ngsContigPrefix
			self.varsP.inputContigFolder = self.varsP.ngsInDir
			break
		if case("extension[01]", regexp=True):
			self.bunching=12
			self.ref_arg="-reff"
			ContigPrefix = self.varsP.expID + "_"+ StageName+'_%s' % self.varsP.extensionCount
			break;
		if case():
			#varsP.error += 1 #these don't do anything
			#varsP.message += '  Error: Refine stage name invalid: '+str(StageName)+'\n'
			self.varsP.updatePipeReport("Internal error: unknown stage %s" % StageName)
			return

        clusargs = varsP.getClusterArgs(StageName) #get arguments before changing StageName, then add suffix
        
        #this block previously in generateJobList
        #if this is 0 stage, need minthreads for the 1 stage. If 1 stage, minthreads isn't used. So always get minthreads from the 1 stage.
        threadstage= (self.refineStage.replace("0","1") if self.groupalgo else self.refineStage)
        minthreads= self.varsP.getClusterArgs(threadstage, category="MinThreads")
        if minthreads:
            minthreads= Template(minthreads).substitute(maxthreads=self.varsP.maxthreads)
        else:
            minthreads=self.varsP.maxthreads
        #if minthreads:
        #    print "minthreads =", minthreads, type(minthreads), "self.varsP.maxthreads =", self.varsP.maxthreads, type(self.varsP.maxthreads)
        #    minthreads= Template(minthreads).substitute(maxthreads=self.varsP.maxthreads)
        #else:
        #    minthreads=self.varsP.maxthreads
        self.minthreads = float(minthreads)

        print self.refineStage, "maxthreads= ", self.varsP.maxthreads, " minthreads= ", minthreads, " totalThreads= ", self.varsP.nThreads #DEBUG

        StageName += (("_%i" % self.varsP.extensionCount) if StageName.startswith("extension") else "") #for status.xml only
        self.varsP.stageName=StageName
        util.LogStatus("progress", "stage_start", StageName)
        #super is more pythonic than referring to the base class explicitly (only matters for multiple inheritance)
        super(Refine, self).__init__(varsP, StageName, clusterArgs=clusargs)

        #set jobs to run on host: for 1-stages only, if groupRefineHostFr (see Multithreading.jobWrapper)
        if self.varsP.groupRefineHostFr > 0 and self.refineStage in ['refineA', 'refineB1', 'refineFinal1', 'extension1'] :
            if self.varsP.clusterArgData.has_key(self.varsP.groupRefineHostCA) : #mediumHostJob
                self.clusterArgsAlt = self.varsP.getClusterArgs(self.varsP.groupRefineHostCA) 
                # self.varsP.NumHostjobs will be set to desired number of hostJobs in generateJobList() (stage 0) OR to 0 (if we want to used fixed fraction of jobs)
            else : #fall-back for old clusterArgs
                self.clusterArgsAlt = self.varsP.getClusterArgs("autoNoise0")
        #intermediateContigPrefix = self.varsP.expID + self.StageName.replace("refine", "_r")
	self.varsP.prepareContigIO(ContigPrefix, StageName)	
        #modify results of varsP.prepareContigIO for special case of refineNGS        
        self.generateJobList()

    def runJobs(self):
        self.multiThreadRunJobs(self.varsP.nThreads, sleepTime = 0.05)
        
    def writeIDFile(self, nJobs):
        f1 = open(self.varsP.idFile, 'wb')
        f1.write(str(nJobs))
        f1.close()
        
    def groupContigName(self, i) :
        return os.path.join(self.varsP.outputContigFolder, self.varsP.outputContigPrefix+"_group"+str(i))
        
    def findGroupedContigs(self):
	print self.varsP.inputContigFolder, self.varsP.inputContigPrefix
	L=getattr(self.varsP, "count_"+self.varsP.inputContigPrefix)

	#L=[]
	#for i in range(1, count+1):
	#	L.append((str(i), os.path.join(self.varsP.inputContigFolder, self.varsP.inputContigPrefix+"_group"+str(i))))
	
	return L

    def threadWeight(self, itr) :
        return int( min( self.varsP.maxthreads, max( self.minthreads, self.varsP.nThreads / 2**itr )))
        
    def groupContigs(self):
        contigFiles, contigIDs = self.varsP.findContigs(self.varsP.inputContigFolder, self.varsP.inputContigPrefix)
        #self.writeIDFile(len(contigFiles))
        M=max(map(int, contigIDs))
        #print(M, len(contigFiles))
        if M<len(contigFiles):
		print("Failure in computing max contig ID: have %d for %d files" %(M, len(contigFiles)))
		raise Exception,("Failure in computing max contig ID: have %d for %d files" %(M, len(contigFiles)))
		
        self.writeIDFile(M) 
        
        if self.bunching==1:
		return zip(contigIDs, contigIDs, contigFiles, [1]*len(contigIDs))

	if self.refineStage == "refineA":
		# For grouping refine A it is important that contigs are in order
		contigIDs=range(1, len(contigFiles)+1)
		base=os.path.join(self.varsP.inputContigFolder, self.varsP.inputContigPrefix)
		contigFiles=[base+"_contig"+str(m)+".cmap" for m in contigIDs]
	
	# Compute contig weights
	wSum, wList = util.contigComputeWeights(contigFiles)
	wLimit=wSum*self.bunching*1.0/len(contigFiles)
	
	# Order contigs to be descending in weight
	if self.refineStage == "refineA":
		# RefineA cannot be sorted this way, however, the output from Assembler is already sorted
		contigOrder=range(0, len(contigFiles))
	else:
		contigOrder=sorted(range(0, len(contigFiles)), key=lambda i: -wList[i])

        k=0
        GroupWeight=0.0
        L=[] #vector of vectors: zeroth ele is group id, first ele is last contig id in group, second ele is groupContigName
        # third ele is used for thread boost: means logic needs to change in generateJobList
        manifest_file=open(self.groupContigName("_manifest"), "wb")
        manifest_file.write("# Manifest file for " + str(self.refineStage) + "\n")
        cnt=0
        group_file=None
        #add new file to log weights
        weightlog = open(self.groupContigName("_weights"), "wb")
        weightlog.write("# Weight log for " + str(self.refineStage) + "\n")
        weightlog.write("# contigID groupID MICthreads Weight GroupCnt GroupCumWeight GroupMaxWeight TotalCumWeight\n")

        if self.groupalgo : #new grouping algo
            #Thomas heuristic : Try to approximately equalize for all groups the estimated run time = MaxWeight / NumThreads
            # weight returned from contigComputWeights is heuristic estimate of total CPUtime on MIC ~ "run time" x NumThreads, provided NumThreads is not too large.
            # For 12Mb contig, weight 12e-6 bp^-2, so could use MaxWeight = 12e-6, or easier to read (in log file) is to multiply weights by 1e6

            # use larger values of origMaxWeight to generate fewer (longer running) jobs
            if self.refineStage == "refineA":
                MaxWeight = origMaxWeight = 240.
            elif self.refineStage == "refineB0":
                MaxWeight = origMaxWeight = 90.
            elif self.refineStage == "refineB1":
                MaxWeight = origMaxWeight = 90.
            elif self.refineStage == "refineFinal0":
                MaxWeight = origMaxWeight = 12.
            elif self.refineStage == "refineFinal1":
                MaxWeight = origMaxWeight = 12.
            else: # extension 
                MaxWeight = origMaxWeight = 90.


            MaxWeightNonHost = MaxWeight * 2 # avoid having very large contigs run on Xeon Phi, to avoid long running jobs
            MaxWeight = origMaxWeight = MaxWeight * self.varsP.GroupSize

            TotalWeight = 0.0

            scale = 1e6
            witr = 0
            gid = 1
            fname = self.groupContigName(gid)
            group_file = open(fname, "wb")

            for m in contigOrder:
                contigFile = contigFiles[m]
                contigID   = contigIDs[m]
                contigWeight = wList[m]*scale

                # start new group file, if contigID is first contig in current group gid
                if k <= 0:
                    fname = self.groupContigName(gid)
                    group_file = open(fname, "wb")

                # add contigID to current group
		group_file.write(contigFile+"\n")
                manifest_file.write(" "+str(contigID))

		k += 1
                GroupWeight += contigWeight
                TotalWeight += contigWeight

                while GroupWeight < MaxWeight and self.threadWeight(witr) > self.minthreads :
                    witr += 1
                    MaxWeight /= 2

#                print "cid= %3s, weight= %5.2f, gid= %2i, witr= %2i, k= %2i, GroupWeight= %5.2f, MaxWeight= %5.2f, nthreads= %2i, CumWeight= %5.2f" % (str(contigID), contigWeight, gid, witr, k, GroupWeight, MaxWeight, self.threadWeight(witr), TotalWeight ) #DEBUG
                weightlog.write("%s %i %i %.2f %i %0.2f %0.2f %0.2f\n" % (str(contigID), gid, self.threadWeight(witr), contigWeight, k, GroupWeight, MaxWeight, TotalWeight))

                if GroupWeight >= MaxWeight : # terminate current group

                    group_file.flush()
                    os.fsync(group_file.fileno())
                    group_file.close()
                    L.append([str(gid), contigID, fname, self.threadWeight(witr), GroupWeight, k]) #fifth ele to sort, 6th element is number of contigs in group
                    manifest_file.write("\n")

                    # start new group
                    gid += 1
                    k = 0
                    GroupWeight = 0


            #end for contigs
                    
            # terminate last group
            if k > 0:
                group_file.flush()
                os.fsync(group_file.fileno())
                group_file.close()
                L.append([str(gid), contigID, fname, self.threadWeight(witr), GroupWeight, k]) #fifth ele to sort, 6th element is number of contigs in this group
                manifest_file.write("\n")

            # flush manifest file
            manifest_file.write("#END\n")
            manifest_file.flush()
            os.fsync(manifest_file.fileno())
            manifest_file.close()

            # flush weightlog file
            weightlog.flush()
            os.fsync(weightlog.fileno())
            weightlog.close()

            # submit jobs in descending order of GroupWeight because for the smallest MaxWeight, the GroupWeight (which predicts runtime in MICs) is no longer monotonically decreasing
            if self.refineStage != "refineA": #CANNOT sort refineA because -contigs arg requires sort by contigID (which is already done)
                L = sorted(L, key = lambda x: -x[4])

            if self.varsP.groupRefineHostFr > 0 : # compute self.varsP.NumHostjobs by accumulating GroupWeight
                WeightSum = 0.0
                TargetWeight = TotalWeight * self.varsP.groupRefineHostFr
                MaxTargetWeight = 2 * TargetWeight # limit host jobs to no more than 2x the normal amount
                self.varsP.NumHostjobs = 0
            
                for job in L:
                    self.varsP.NumHostjobs += 1
                    WeightSum += job[4]
                    if WeightSum > TargetWeight and ( job[4] <= MaxWeightNonHost or WeightSum > MaxTargetWeight ) :
                        break;
                
                print "NumHostjobs= %i, TotalJobs= %i, HostWeight= %5.2f, TargetWeight= %5.2f, TotalWeight= %5.2f" % ( self.varsP.NumHostjobs, len(L), WeightSum, TargetWeight, TotalWeight )

#            print "groupID, LastContigID, name, MICthreads, GroupWeight, NumContigs:\n"+"\n".join([str(x) for x in L]) # DEBUG
            return L

        #old grouping algo
        
        self.varsP.NumHostjobs = 0 # default to fixed fraction of Host jobs if self.varsP.groupRefineHostFr > 0

        gid = 1
	# Decrease bunching if we have too few contigs to analyze
	if len(contigFiles)< self.bunching*300:
		self.bunching=len(contigFiles)/300
		if self.bunching<1:
			self.bunching=1

	# Increase bunching if we have too many contigs to avoid going past open files limit (4096)
	if self.bunching>1 and len(contigFiles)>self.bunching*1000:
		self.bunching=len(contigFiles)/1000
		if self.bunching<1:
			self.bunching=1

	# Do not group the last tail_len contigs so that the jobs finish faster
	tail_len = 0 #doesn't work bc scheduler will find slots for them before large jobs finished
        for m in contigOrder:
		contigFile=contigFiles[m]
		contigID=contigIDs[m]
		contigWeight=wList[m]

		# Make sure that extra large contigs are in one group
		if group_file and ((not self.multigroup) or (k>=self.bunching) or (GroupWeight>wLimit) or (cnt+tail_len>len(contigFiles))) :
			group_file.close()
			w=L[len(L)-1]
			L[len(L)-1]=(w[0], w[1], w[2], (GroupWeight*1.0/wLimit)**2)
			k=0
			GroupWeight=0.0
			
		if k==0:
			fname=self.groupContigName(gid)
			L.append((str(gid), contigID, fname, (contigWeight*1.0/wLimit)**2))
			gid+=1
			#print "Creating", fname
			manifest_file.write("\n"+str(contigID))
			group_file=open(fname, "wb")
		else:
			manifest_file.write(" "+str(contigID))			
		group_file.write(contigFile+"\n")
		k+=1
		GroupWeight+=contigWeight
		cnt+=1
			
	if k>0:
		group_file.close()
		w=L[len(L)-1]
		L[len(L)-1]=(w[0], w[1], w[2], (GroupWeight*1.0/wLimit)**2)
        manifest_file.write("\n")
        manifest_file.flush();
        os.fsync(manifest_file.fileno());
	manifest_file.close()
	return L
    #end groupContigs
	

    def generateJobList(self):
        baseArgs1 = self.varsP.argsListed(self.refineStage)
        
	for case in util.switch(self.refineStage):
		if case("refine(B1|Final1)", regexp=True):
			baseArgs1 += self.varsP.argsListed('noise0')
			ContigGroupList=self.findGroupedContigs()
			r1args = [self.varsP.RefAlignerBin]
			break
		if case("refine(B0|Final0)", regexp=True):
			baseArgs1 += self.varsP.argsListed('noise0')
			ContigGroupListFull=self.groupContigs()
			setattr(self.varsP, "count_"+self.varsP.outputContigPrefix, (ContigGroupListFull))
			#print self.varsP.outputContigPrefix, getattr(self.varsP, "count_"+self.varsP.outputContigPrefix)
			#r1args = [self.varsP.RefAlignerBin, '-i', self.varsP.bnxFile]
			#InputFileList=[self.varsP.bnxFile]
			r1args = [self.varsP.RefAlignerBin]
			ContigGroupList = zip(range(1,self.varsP.nPairwiseJobs + 1), range(1,self.varsP.nPairwiseJobs + 1), [self.varsP.bnxFile.replace(".bnx", "_%s_of_%s.bnx" %(x, self.varsP.nPairwiseJobs)) for x in range(1,self.varsP.nPairwiseJobs + 1)], [1]*self.varsP.nPairwiseJobs)
			break
		if case("refineA"):
			baseArgs1 += self.varsP.argsListed('noise0')
			ContigGroupList = self.groupContigs()
			print("Found %d groups for refineA" % (len(ContigGroupList)))
			#r1args = [self.varsP.AssemblerBin, '-i', self.varsP.bnxFile.replace(".bnx", "_sorted.bnx")] #need this before -contigs -- can no longer use all_sorted.bnx due to scan scaling: must refer to varsP.sorted_file
			#r1args = [self.varsP.AssemblerBin, '-i', self.varsP.sorted_file+".bnx"] #need this before -contigs
			r1args = [self.varsP.AssemblerBin, '-if', self.varsP.bnxFileList] #need this before -contigs; use split files in case splitting changed (eg due to scan scaling producing labels at < 20 bp)
			r1args += ['-contigs', os.path.join(self.varsP.inputContigFolder, self.varsP.inputContigPrefix) + '.contigs']
			break
		if case("refineNGS"):
			r1args = [self.varsP.RefAlignerBin, '-i', self.varsP.bnxFile]
			ContigGroupList=self.groupContigs()
			break
		if case("extension0"):
			baseArgs1 += self.varsP.argsListed('noise0')
			ContigGroupList=self.groupContigs()
			setattr(self.varsP, "count_"+self.varsP.outputContigPrefix, (ContigGroupList))
			#print self.varsP.outputContigPrefix, getattr(self.varsP, "count_"+self.varsP.outputContigPrefix), self.varsP.inputContigFolder, self.varsP.inputContigPrefix
			#r1args = [self.varsP.RefAlignerBin, '-i', self.varsP.bnxFile]
			#InputFileList=[self.varsP.bnxFile]
			r1args = [self.varsP.RefAlignerBin]
			ContigGroupList = zip(range(1,self.varsP.nPairwiseJobs + 1), range(1,self.varsP.nPairwiseJobs + 1), [self.varsP.bnxFile.replace(".bnx", "_%s_of_%s.bnx" %(x, self.varsP.nPairwiseJobs)) for x in range(1,self.varsP.nPairwiseJobs + 1)], [1]*self.varsP.nPairwiseJobs)
			break;
		if case("extension1"):
			baseArgs1 += self.varsP.argsListed('noise0')
			ContigGroupList=self.findGroupedContigs()
			r1args = [self.varsP.RefAlignerBin]
			break;
		if case():
			varsP.error += 1
			varsP.message += '  Error: Refine stage name invalid: '+str(StageName)+'\n'
			return


        stdarg = []
        if self.varsP.stdoutlog : #this is the same for all cases below
            stdarg = ['-f', '-stdout', '-stderr'] 

        output1String = os.path.join(self.varsP.outputContigFolder, self.varsP.outputContigPrefix)

        for m in range(0, len(ContigGroupList)):
	    contigID=ContigGroupList[m][0]
	    rawContigID=ContigGroupList[m][1]
	    contig=ContigGroupList[m][2]
	    
            if self.groupalgo : #here, ContigGroupList[m][3] is simply nthreads
                nthreads = int(ContigGroupList[m][3])
            else : #original algo: here it is the boost
                # Figure out desired number of threads to use
                threadBoost = max(ceil(ContigGroupList[m][3]), 1)
                nthreads=int(round(self.minthreads*threadBoost))
                if nthreads>self.varsP.maxthreads:
		    nthreads=self.varsP.maxthreads

            jobName = self.refineStage + ' %5s' % contigID
	    for case in util.switch(self.refineStage):
		    if case("refineA"):
                        if self.groupalgo: # here, ContigGroupList[m][1] is the last contig in the group and ContigGroupList[m][5] is the Number of Contigs in the group
                            endId = int(rawContigID)
                            rawContigID = endId - int(ContigGroupList[m][5]) + 1
                        else:
                            endId=int(rawContigID)+self.bunching-1
                            if m+1<len(ContigGroupList) :
				endId=int(ContigGroupList[m+1][1])-1
			currentArgs = [str(rawContigID), str(endId)] #this must come after r1args because it's actually an argument to -contigs
			currentArgs = r1args + currentArgs + ['-o', output1String] + ['-maxthreads', str(nthreads)] + stdarg + baseArgs1
			expectedOutputString = self.varsP.outputContigPrefix + '_contig' + str(rawContigID)
			expectedResultFile = os.path.join(self.varsP.outputContigFolder, expectedOutputString + '.cmap') #refineB
			expectedStdoutFile = output1String + "_id"+str(rawContigID)+".stdout"
			break
			
		    if case("refineB1|refineFinal1|extension1", regexp=True):
			#Inputs=zip(["-i"]*self.varsP.nPairwiseJobs, [contig.replace("_group", "_group"+str(i)+"_mapped_group")+".bnx" for i in range(1,self.varsP.nPairwiseJobs + 1)]) #this line has a rare bug if "_group" appears in the path to the group file, ie, '/home/user/abc_group1_assembly/contigs/...'
			Inputs=zip(["-i"]*self.varsP.nPairwiseJobs, [os.path.join(os.path.dirname(contig),os.path.basename(contig).replace("_group", "_group"+str(i)+"_mapped_group"))+".bnx" for i in range(1,self.varsP.nPairwiseJobs + 1)])
			Inputs=[x for t in Inputs for x in t]
                        #-id must come before -o, otherwise expectedStdoutFile is wrong
			currentArgs = ['-maxthreads', str(nthreads), '-id', str(contigID), '-o', output1String, self.ref_arg, contig]
			currentArgs = r1args + currentArgs + stdarg + baseArgs1 + Inputs 
			expectedOutputString = self.varsP.outputContigPrefix + '_contig' + str(rawContigID) + '.cmap'
                        if '-Haplotype' in baseArgs1 :
                            #update for Haplotype: check for file like exp_refineFinal1_contig96460.xmap
                            # do not check for a cmap in this case because it's not known whether it will be 0.cmap or {1,2}.cmap
                            expectedOutputString = self.varsP.outputContigPrefix + '_contig' + str(rawContigID) + '.xmap'
			expectedResultFile = os.path.join(self.varsP.outputContigFolder, expectedOutputString)
			expectedStdoutFile = output1String + "_id"+str(contigID)+".stdout"
			break
			
		    if case("refineB0|refineFinal0|extension0", regexp=True):
                        nthreads = self.varsP.maxthreadsPW # NOTE : for 0 stage, always use -jp (defaults to -T)
			currentArgs = [ '-maxthreads', str(nthreads), "-ref", os.path.join(self.varsP.inputContigFolder, util.uniquifyContigName(self.varsP.inputContigPrefix)+".cmap")]
                        outputfile = os.path.join(self.varsP.outputContigFolder, self.varsP.outputContigPrefix+'_group'+str(contigID))
                        #-id must come before -o, otherwise expectedStdoutFile is wrong
			currentArgs = r1args + ['-i', contig, '-id', str(contigID), '-o', outputfile] + stdarg + currentArgs + baseArgs1
                        currentArgs += ['-refine', '0', '-grouped', os.path.join(self.varsP.outputContigFolder, self.varsP.outputContigPrefix+'_group_manifest'), '-mapped', os.path.join(self.varsP.outputContigFolder, self.varsP.outputContigPrefix+'_group'+str(contigID)+"_mapped"), '-output-filter', ".*.bnx"]
			expectedOutputString = self.varsP.outputContigPrefix+'_group'+str(contigID)+"_mapped.bnx"
			expectedResultFile = outputfile + "_mapped_group1.bnx" 
			expectedStdoutFile = outputfile + "_id"+str(contigID)+ ".stdout"
			break
			
		    if case():
			self.varsP.updatePipeReport("Internal error: cannot handle stage %s" % (self.refineStage))
			raise ValueError
                
            if self.varsP.bnxStatsFile!=None:
		currentArgs.extend(['-XmapStatRead', self.varsP.bnxStatsFile])

            #print " ".join(currentArgs) #DEBUG
            s1Job = mthread.singleJob(currentArgs, 
                                    jobName, 
                                    expectedResultFile, 
                                    expectedOutputString,
                                    maxThreads=nthreads,
                                    clusterLogDir=self.varsP.clusterLogDir,
                                    expectedStdoutFile=expectedStdoutFile,
                                    )
            self.addJob(s1Job)
        self.logArguments()
            
    def checkResults(self, stageSuffix=""):
        '''Call jobWrapper (self) .doAllPipeReport, and varsP.mergeIntoSingleCmap.
        stageSuffix, if supplied, is appended to varsP.stageComplete in order to
        fix the stage name reported by the CharacterizeModule in the informaticsReport.
        '''
        self.doAllPipeReport()
        self.varsP.stageComplete = self.refineStage + stageSuffix
        if self.refineStage not in ['refineB0', 'refineFinal0', 'extension0'] :
		self.varsP.mergeIntoSingleCmap()
        StageName = self.refineStage + ("_%i" % self.varsP.extensionCount if self.refineStage.startswith("extension") else "") #for status.xml only
        util.LogStatus("progress", "stage_complete", StageName)
        #after stage is complete, exit Pipeline IFF a job fails _in a 0 stage_
        if (self.refineStage in ['refineB0', 'refineFinal0', 'extension0'] and
            not self.allResultsFound()) :
            err = "0-stage job failed, aborting Pipeline: stage=%s" % self.refineStage
            self.varsP.updatePipeReport( err+"\n" )
            util.LogError("critical", err)
            raise RuntimeError
    #end checkResults

    def endStage(self) :
        """Call this in place of checkResults when this stage is bypassed."""
        if self.refineStage not in ['refineB0', 'refineFinal0', 'extension0'] :
            self.varsP.mergeIntoSingleCmap()
        StageName = self.refineStage + ("_%i" % self.varsP.extensionCount if self.refineStage.startswith("extension") else "") #for status.xml only
        self.varsP.stageComplete = StageName
        util.LogStatus("progress", "stage_complete", StageName)

