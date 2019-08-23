import os, sys
import subprocess
import time
import shutil
#import pdb
import threading
from string import Template


"""@package Multithreading Controls queuing, process exectution, and process 
logging for all de novo phases

Schedules jobs on single node based on user defined thread limit
Submits jobs to drmaa scheduler for distributed execution
Checks result for expected outputs and logs computation and performance details
"""


import utilities
utilities.setVersion("$Id: Multithreading.py 6416 2017-05-25 21:55:38Z tanantharaman $")


try:
    import drmaa
except:
    pass


global my_wait

if os.name == "nt":
	import ctypes
	global nt_handles_list
	global SYNCHRONIZE
	global INFINITE
	SYNCHRONIZE=0x00100000
	INFINITE = -1

	nt_handles_list = {}
	
	
	def nt_wait():
		at=ctypes.c_long *len(nt_handles_list)
		ha=at(*nt_handles_list.keys())
		ret=ctypes.windll.kernel32.WaitForMultipleObjects(len(ha), ha, False, INFINITE)
		h=ha[ret]
		ctypes.windll.kernel32.CloseHandle(h)
		pid=nt_handles_list[h]
		print "Process %d done" % nt_handles_list[h]
		del nt_handles_list[h]
		# Fake return code as I don't know how to obtain it
		return (pid, 0)
	
	my_wait=nt_wait
	
else:
	my_wait=os.wait


global thread_list
thread_list=[]

#call start_new_thread if number of open threads is < argument, otherwise wait
#don't forget to call wait_all_threads after calling this
def start_thread_limit(threadlimit, method, args=(), kwargs={}):
    if threadlimit <= 0 :
        print "Warning in start_thread_limit: invalid threadlimit:", threadlimit, "no jobs started"
        return
    global thread_list
    while len(thread_list) >= threadlimit :
        time.sleep(1)
        #this is similar to wait_all_threads (below), but just wait on any one to finish, not all
        #thread_list.reverse() #reverse doesn't return anything; reverses in place
	for t in thread_list : 
            #if t.isAlive() :
            #    t.join()
            #else :
            if not t.isAlive() :
                thread_list.remove(t)
                break
    start_new_thread(method, args, kwargs)

	
def start_new_thread(method, args=(), kwargs={}):
	t=threading.Thread(target=method, args=args, kwargs=kwargs)
	t.start()
	global thread_list
	thread_list.append(t)
	return(t)
	
def wait_all_threads():
	global thread_list
	for t in thread_list:
		while t.isAlive():
			t.join()
	thread_list=[]
	
	
global cSession
cSession = None


class jobWrapper(object):
    """Controls queueing for single node and cluster execution of jobs
        
    """
    def __init__(self, varsP, groupName, clusterArgs=None, throttle=0):
        self.varsP = varsP
        self.jobList = []
        self.allJobsComplete = False
        self.groupName = groupName
        self.elapsedTime = 0.
        self.cpuTime = 0.
        self.nThreads = 0
        
        self.onCluster = False
        self.clusterArgs = None
        self.clusterArgsAlt = None #used only in GroupedRefinementModule
        if self.varsP.clusterLogDir:
            self.onCluster = True
            self.clusterArgs = clusterArgs
        
        self.throttle = False
        if throttle:
            self.throttle = True
            self.throttleMax = throttle

    
    def clearJobs(self):
        self.jobList = []
    
    def addJob(self,sJob):
        sJob.setParams(self.varsP) #get time/perf from varsP
        self.jobList.append(sJob)
        sJob.jobNum = len(self.jobList)
    
    def addJobList(self,jobList):
        for sJob in jobList:
            self.addJob(sJob)
    
    def logArguments(self):
        """Log the exact executable call for the first job in the queue
        
        """
        #print "logArguments:", self.groupName #debug
        pipeReport = self.groupName + "\n" #record that this stage has constructed its jobs
        if len(self.jobList) >= 1:
            sJob = self.jobList[0]
            pargs = [] #can't add to sJob.args because that changes self.jobList[0].args
            if sJob.time :
                pargs += sJob.timeArgs
            if sJob.perf :
                pargs += sJob.perfArgs
            pipeReport += ' '.join(pargs + sJob.args) + 2 * '\n'
            #print "pipereport:", self.varsP.pipeReportFile #debug
            #print "pipereport:", pipeReport #debug
        self.varsP.updatePipeReport(pipeReport, printalso=False)
    
    def multiThreadRunJobs(self, nActiveThreads, sleepTime = 0.01, threadControl=False, background=False, callLogStatus=True):
        """Main Queue script, start jobs, check for completion
        
        """

        #this is useful as a generic way to skip running a module--no jobs are submitted
        if len(self.jobList) == 0 : 
            self.varsP.updatePipeReport(" Warning in multiThreadRunJobs: number of jobs is 0, skipping stage: "+self.groupName+"\n")
            return            

        if nActiveThreads == 0 :
            self.varsP.updatePipeReport(" Error in multiThreadRunJobs: nActiveThreads must be > 0, skipping stage: "+self.groupName+"\n")
            return
        
        if background:
		start_new_thread(self.multiThreadRunJobs, (nActiveThreads, sleepTime, threadControl, False))
		return
        
        utilities.logMemory(self.varsP.memoryLogpath, self.varsP.startTime, self.groupName) #call at start and end of this method
        jobw = 30 #width of job name in printout
        print ' Starting Multi-Threaded Process:'
        print '  ' + self.groupName
        self.nThreads = nActiveThreads
        availableThreads = nActiveThreads
        startTime = time.time() 
        activeJobList = []
        nActiveJobs = 0
        nFinishedJobs = 0
        nActiveThrottle = 0
        nJobs = len(self.jobList)
        nRemainingJobs = nJobs

        global cSession
        if self.onCluster and cSession == None:
            cSession = drmaa.Session()
            cSession.initialize()
        print '  Running ' + str(nJobs) + ' jobs with ' + str(nActiveThreads) + ' threads'
        if callLogStatus :
            utilities.LogStatus("progress", "jobs_outstanding", str(nJobs), self.groupName) 
            utilities.LogStatus("progress", "stage_pct_done", "0.0", self.groupName)
	job_status=(0,nJobs)
        while True :
            if nRemainingJobs > 0:
                for i,sJob in enumerate(self.jobList):
                    if sJob.jobStarted or sJob.isRunning or sJob.isComplete:
                        continue
                    if sJob.hasContingentJob:
                        if not sJob.contingentJob.isComplete :
                            continue
                    if not(sJob.onCluster):
                        if nActiveJobs >= nActiveThreads:
                            #print "skipping1", nActiveJobs, nActiveThreads #debug
                            continue
                    if not(sJob.onCluster):
                        if availableThreads < sJob.maxThreads:
                            #print "skipping2", availableThreads, sJob.maxThreads #debug
                            continue
                    if self.throttle and sJob.throttleClass:
                            if nActiveThrottle >= self.throttleMax:
                                continue
                            nActiveThrottle += 1
                    activeJobList.append(sJob)
                    nActiveJobs += 1
                    clustarg = self.clusterArgs
                    #note: will not move any jobs to host if total number < 5
                    if self.varsP.groupRefineHostFr > 0 and self.clusterArgsAlt and self.refineStage in ['refineA', 'refineB1', 'refineFinal1', 'extension1'] :
                        if ( self.varsP.NumHostjobs <= 0 ) : 
                            self.varsP.NumHostjobs = self.varsP.groupRefineHostFr * len(self.jobList)
                        if ( i < self.varsP.NumHostjobs ) :
                            clustarg = self.clusterArgsAlt
                            self.editArgsForHost(sJob)
#                            print "Moving job to host:", sJob.jobName
                        else : # For regular jobs adjust -TotalThreads value to match -T
                            if "-TotalThreads" in sJob.args :
                                a = utilities.getIntFromString(self.varsP.nThreads)
                                if a : #getIntFromString returns None if arg not int
                                    #TotalThreads takes an optional arg: must check if it's present: if so, change, if not, add
                                    t = sJob.args.index("-TotalThreads")
                                    if len(sJob.args) > t+1 and not sJob.args[t+1].startswith("-") : #replace
                                        sJob.args[ t+1 ] = str(a)
                                    else :
                                        sJob.args.insert(t+1, str(a))
                    sJob.startJob(cSession=cSession, clusterArgs=clustarg)
                    availableThreads -= sJob.maxThreads
                    nRemainingJobs -= 1
                    statusString = ('   START %4d: % '+str(jobw)+'s, %3d Thr, %4d R, %4d T, %4d F, %4d Q') % (sJob.jobNum,sJob.jobName[:jobw], sJob.maxThreads, nActiveJobs, nJobs, nFinishedJobs, nRemainingJobs)
                    print statusString
		    sys.stdout.flush()
                    time.sleep(sleepTime) #sleep between job submission, but wait to check status

            #The block below is error prone in the case of multiple jobWrapper objects running simultaneously,
            # which we have implemented for CharacterizeModule using threading. The problem is that the
            # characterize os.wait call can steal the pid of another job, say, refinement, and then the
            # refinement job will never be marked completed. Simplest is to just wait, and inside
            # CheckIfRunning, the poll will take care of each job individually.
            time.sleep(sleepTime)
            '''
	    (pid, rc)=(-1, -1) # Defaults so the statement works for cluster jobs
            if self.onCluster : #if on cluster, no os.wait call is needed; sleep instead, then check all jobs
                time.sleep(sleepTime)
            else :
                try : #if not on cluster, use os.wait to wait for child process to finish
                    global my_wait
                    (pid, rc)=my_wait() #any child?
                except OSError, e :
                    time.sleep(sleepTime)
                    #print e

            #Set the return code of the job which was stolen by the wait call above (see comment below).
            for sJob in activeJobList: 
                if sJob.markCompleted(pid, rc) :
                    break #skip rest once correct one found
                    '''

            #The old version of this loop was dangerous because it popped from the list being iterated over.
            #So, if you skip a job due to this and that job's return code was stolen by the wait above,
            # then the job is never marked complete.
            #Though the below fix is probably sufficient, do the above also just to be safe.
            #If you iterate backwards, using reversed, removing an element will not affect the loop on the remaining elements
            for sJob in reversed(activeJobList) : 
                #sJob.markCompleted(pid, rc) #this call moved into loop above (see above comments)
                if sJob.CheckIfRunning(cSession=cSession):
                    continue
                else:
                    #activeJobList.pop(i)
                    activeJobList.remove(sJob)
                    nActiveJobs -= 1
                    nFinishedJobs += 1
                    availableThreads += sJob.maxThreads
                    availableThreads = min(nActiveThreads, availableThreads)
                    if self.throttle and sJob.throttleClass:
                        nActiveThrottle -= 1
                    statusString = ('   STOP  %4d: % '+str(jobw)+'s, %3d Thr, %4d R, %4d T, %4d F, %4d Q') % (sJob.jobNum,
                                                                                                         sJob.jobName[:jobw],
                                                                                                         sJob.maxThreads,
                                                                                                         nActiveJobs,
                                                                                                         nJobs,
                                                                                                         nFinishedJobs,
                                                                                                         nRemainingJobs)
                    
                    if sJob.cpuPercent <= 0.0 :
                        ElapsedTime = 0.0
                        print "WARNING: Job %s has cpuTime= %0.4f, cpuPercent= %0.4f\n" % ( str(sJob.jobName), sJob.cpuTime, sJob.cpuPercent ) # DEBUG
                    else :
                        ElapsedTime = sJob.cpuTime * 100.0 / sJob.cpuPercent
                    statusString += '  TotalTime=' + timeFormat1(sJob.runTime) + '  RunTime=' + timeFormat1(ElapsedTime) + '  CPUload=' + ( "%i" % int(sJob.cpuPercent) ) + '%'
                    if ( sJob.HOSTname != '' ) :
                        statusString += ' ' + sJob.HOSTname
                    else :
                        statusString += ' host=NA'
                    if ( sJob.errString != '' ) :
                        statusString += ' ' + sJob.errString
                    print statusString
                    
            #log status after above loop to calculate nFinishedJobs
            pct_done = (nFinishedJobs*100.0/nJobs if nJobs > 0 else 0)
            njr = nJobs - nFinishedJobs #num jobs remaining
	    new_status=(pct_done, njr)
	    if job_status!=new_status and callLogStatus :	
		utilities.LogStatus("progress", "jobs_outstanding", "%d" % njr, self.groupName) 
		utilities.LogStatus("progress", "stage_pct_done", "%.01f" % pct_done, self.groupName) 
		job_status=new_status
                    
            if nActiveJobs == 0 and nRemainingJobs == 0:
                break
            elif nActiveJobs < 0 or nRemainingJobs < 0:
                print "ERROR in multithreading: invalid: nActiveJobs:", nActiveJobs, "nRemainingJobs:", nRemainingJobs
                break

            #Note: you cannot check len(activeJobList) here becuase if one job takes all the threads,
            # it can finish, causing the list to be empty, but there are still more jobs to submit.
            # This is not an error.

            sys.stdout.flush()
        #end job submission - check loop

        #if self.onCluster:
        #   cSession.exit()
        self.elapsedTime = time.time() - startTime
        self.cpuTime = 0.
        for sJob in self.jobList:
            self.cpuTime += sJob.runTime
        print ' Finished Multi-Threaded Process:'
        print '  ' + self.groupName
        print #extra newline for readability
        sys.stdout.flush()
        utilities.logMemory(self.varsP.memoryLogpath, self.varsP.startTime, self.groupName) #call at start and end of this method
    #end multiThreadRunJobs


    #all modules have their own version of checkResults
    #this will partially unify, and simplify, them all
    def doAllPipeReport(self):
        """Call make run report and makeParseReport, write results to pipelineReport.

        """
        #synchronize this method with change to multiThreadRunJobs: if no jobs, report nothing
        if len(self.jobList) > 0 :
#            print "Calling makeRunReport()" # DEBUG
            report =  self.makeRunReport()
#            print "Calling makeParseReport()" # DEBUG
            report += self.makeParseReport()
            self.varsP.updatePipeReport(report, printalso=False)

    
    def makeRunReport(self):
        """Log process/time details for job set
        
        """
        if not self.varsP.stdoutRecheck :
            print "\t Not Rechecking Stdout completeness for " + self.groupName + "\n" # DEBUG

        reportString = ''
        for sJob in self.jobList:
#            print "Calling CheckIfFileFound() for job:", sJob.jobName # DEBUG
            sJob.CheckIfFileFound() #here, not in CheckIfRunning, for big bear compatibility
#            print "Calling CheckStdout() for job:", sJob.jobName # DEBUG
            if self.varsP.stdoutRecheck :
                sJob.CheckStdout(tries=0) #all jobs are done so if fail, no need to keep trying
            else :
                sJob.stdoutComplete = True # just assume job completed to save time
            reportString += sJob.makeReportString() + '\n'
        reportString += '\n'
        reportString += '  Completed %d jobs on %d threads\n' % (len(self.jobList), self.nThreads)
        reportString += '  Clock time: ' + timeFormat1(self.elapsedTime)
        reportString += '  CPU time  : ' + timeFormat1(self.cpuTime)
        return reportString
        
    def makeParseReport(self):
        """Write per-job parsable text for process/time details
        
        """
        reportString = '\n  Begin machine parsable segment for: ' + self.groupName + '\n'
        for i,sJob in enumerate(self.jobList):
            if i == 0:
                reportString += sJob.GetPerformanceString(headerOnly=True) + '\n'
            #reportString += sJob.GetPerformanceString() #not working on big bear
            reportString += sJob.parseString #set by GetPerformanceString, which is called from CheckIfRunning
        reportString += '\n'
        return reportString

    def makeSimpleReport(self):
        """Like makeRunReport, but does not report any time info, and only reports failed jobs.
        """
        if not self.varsP.stdoutRecheck :
            print "\t Not Rechecking Stdout ccompleteness for " + self.groupName + "\n" # DEBUG

        reportString = ''
        for sJob in self.jobList:
            sJob.CheckIfFileFound() #here, not in CheckIfRunning, for big bear compatibility
            if self.varsP.stdoutRecheck :
                sJob.CheckStdout(tries=1, delay=0.1) #this is for local runs, assume no nfs issues
            else :
                sJob.stdoutComplete = True # just assume job completed to save time
            reportString += sJob.makeReportStringSimple() #no newline (otherwise get newline for every job)
        return reportString

    def allResultsFound(self):
        """Return whether the result was found--ie, success--for all jobs.
        """
        #return all([job.resultFound for job in self.jobList])
        #resultFound is not as reliable as stdout--switch to this
        return all([job.stdoutComplete for job in self.jobList])

    #note this should really be a singleJob method, but singleJobs don't have varsP
    def editArgsForHost(self, job):
        """When running on host, need to change threads and memory arguments to twice the -J value (to allow for 2x host overloading)
        """
        job.maxThreads = self.varsP.maxthreadsHost
        if "-maxthreads" in job.args :
            job.args[ job.args.index("-maxthreads")+1 ] = str(self.varsP.maxthreadsHost)
        if "-TotalThreads" in job.args :
            a = utilities.getIntFromString(self.varsP.maxthreadsHost)
            if a : #getIntFromString returns None if arg not int
                #TotalThreads takes an optional arg: must check if it's present: if so, change, if not, add
                t = job.args.index("-TotalThreads")
                if len(job.args) > t+1 and not job.args[t+1].startswith("-") : #replace
                    job.args[ t+1 ] = str(a*2)
                else :
                    job.args.insert(t+1, str(a*2))
        if "-maxmem" in job.args :
            job.args[ job.args.index("-maxmem")+1 ] = self.varsP.argsListed(self.varsP.largemem)[1] #should be just ['-maxmem', '128']
        if "-maxvirtmem" in job.args :            
            job.args[ job.args.index("-maxvirtmem")+1 ] = "0"
        else :
            job.args.extend( ["-maxvirtmem", "0"] )
    #end editArgsForHost

#end class jobWrapper
        
        
class singleJob(subprocess.Popen):
    """Primary worker class: call to executable
    
    """
    
    def __init__(self, args, jobName, expectedResultFile, uniqueString, maxThreads=1, forceForward=None, throttleClass=False,stdOutFile=None, stdErrOutFile=None, clusterLogDir=None, expectedStdoutFile=None, maxRestarts=3):
        self.maxThreads = maxThreads
        self.isComplete = False
        self.jobStarted = False
        self.jobName = jobName
        self.expectedResultFile = expectedResultFile
        self.expectedStdoutFile = expectedStdoutFile
        self.resultFound = False
        self.stdoutComplete = False
        self.args = args
        self.isRunning = False
        self.hasContingentJob = False
        self.jobNum = 0
        self.runTime = -1 # total time spent on SGE queue, including elapsed wait and run time
        self.cpuTime = 0.0 # total CPU time (sum over all threads)
        self.cpuPercent = 0.0 # Average CPU percentage (average number of threads * 100%)
        self.HOSTname = '' # name of host (if known, from "Connection to phiNN closed" ssh message)
        self.errString = '' # Error message in ClusterLogs (eg "Command terminated by signal 6")
        self.uniqueString = uniqueString
        self.throttleClass = throttleClass
        self.forceForward = None
        self.stdOutFile = stdOutFile
        self.StdOutFileOpen = False
        self.stdOutObject = ''
        self.stdErrOutFile = stdErrOutFile
        self.StdErrOutFileOpen = False
        self.stdErrOutObject = ''
        self.maxRestarts=maxRestarts
        self.remainingRestarts=maxRestarts
        
        self.jobID = -1
        self.onCluster = False
        if clusterLogDir:
            self.onCluster = True
            self.clusterTarget = os.path.join(clusterLogDir, uniqueString)
            self.jTemplate = None
        self.startTime = 0
        
        self.time = False
        self.perf = False
                
        self.timeArgs = ['/usr/bin/time', '-f', '%U\t%S\t%E\t%P\t%M']
        if self.onCluster:
            self.timeArgs = ['${TIME_BINARY:=/usr/bin/time}', '-f', '"%U\t%S\t%E\t%P\t%M"']
        self.timeHeader = '%s\t%s\t%s\t%s\t%s\t%s' % ('Usr','Sys','Elpsd', 'CPU','MaxMem','EndTime')
        
        self.perfArgs = ['/usr/bin/perf', 'stat', '-x', '\\t', '--log-fd', '2']
        if self.onCluster:
            self.perfArgs = ['${PERF_BINARY:=/usr/bin/perf}', 'stat', '-x', '\\\\t', '--log-fd', '2']
        self.perfHeader = '%s' % 'perfHeader'
        self.perfString = ''
        self.parseString = ''
        
    
    def setParams(self, varsP) :
        """singleJob.setParams: some parameters make more sense to set on
        adding job to the list jobWrapper.jobList, rather than in the
        creation or starting of the job. This is mainly due to varsP
        handling.
        """
        self.time = varsP.time
        self.perf = varsP.perf


    def startJobOnCluster(self,cSession,clusterArgs=None):         
        
        clusterProcFile = self.clusterTarget + '.sh'
        self.stdErrOutFile = self.clusterTarget + '.log'
        self.returncode = -1000
        self.StdErrOutFileOpen = False
        clusterExecArgs = []
        if self.time:
            clusterExecArgs += self.timeArgs
        if self.perf:
            clusterExecArgs += self.perfArgs
        clusterExecArgs += self.args
        clusterExecArgs += ["&&", "date"] #see below
        #print clusterExecArgs #debug
        clusterExecLine = ' '.join(clusterExecArgs)
        f1 = open(clusterProcFile, 'wb')
        f1.write('#!/bin/bash\n#$ -N %s\n#$ -o %s\n#$ -j y\nhostname\nulimit -a\n' % (self.uniqueString, self.stdErrOutFile))
        f1.write('date\n')
        bashEnv = 'export O64_OMP_SET_AFFINITY=false'
        f1.write('%s\n%s\n' % (bashEnv, clusterExecLine))
        #f1.write('date\n') #this causes the exit code of this command to be exit code of shell script, which is bad
        f1.flush()
        os.chmod(clusterProcFile, 0775) 
        os.fsync(f1.fileno())
        f1.close()
        
        self.jTemplate = cSession.createJobTemplate()
        self.jTemplate.remoteCommand = clusterProcFile
        self.jTemplate.joinFiles = True
        self.jTemplate.outputPath = ':' + self.stdErrOutFile
        self.jTemplate.jobName = self.uniqueString
	self.cSession = cSession
        #jTemplate.nativeSpecification = "-pe openmp 8"  # For parallel environment
        self.clusterArgs=None
        if clusterArgs:
	    global status_log_filename
	    self.clusterArgs=Template(clusterArgs).substitute(numthreads=self.maxThreads, status_log_filename=utilities.status_log_filename, restart_count=self.maxRestarts-self.remainingRestarts)
            self.jTemplate.nativeSpecification = self.clusterArgs
            
        self.submitClusterJob()
        
    def restartJobOnCluster(self):
	self.remainingRestarts-=1
        self.jTemplate = self.cSession.createJobTemplate()
        clusterProcFile = self.clusterTarget + '.sh'
        self.jTemplate.remoteCommand = clusterProcFile
        self.jTemplate.joinFiles = True
        self.jTemplate.outputPath = ':' + self.stdErrOutFile
        self.jTemplate.jobName = self.uniqueString
        #jTemplate.nativeSpecification = "-pe openmp 8"  # For parallel environment
        if self.clusterArgs:
            self.jTemplate.nativeSpecification = self.clusterArgs
            
        self.submitClusterJob()
	
    def submitClusterJob(self):
	self.jobStarted=False
        while not self.jobStarted:
		try:
                    self.jobID = self.cSession.runJob(self.jTemplate)
		#except drmaa.errors.DeniedByDrmException:
                #in addition to above, also got this exception on cluster (so just catch 'em all):
                # drmaa.errors.DrmCommunicationException: code 2: failed receiving gdi request response for mid=65534 (can't send response for this message id - protocol error).
                except Exception, e:
                    err = "Exception encountered while submitting job: %s, nativeSpecification=\"%s\"" % (str(e),self.clusterArgs)
                    utilities.LogError("warning", err)
                    print err, "\n"
                    sys.stdout.flush()
                    time.sleep(30) #100) #reduce from 100 (this is seconds)
                else : #only execute if cSession.runJob succeeds
                    self.pid = self.jobID
                    self.startTime = time.time()
                    self.isRunning = True
                    self.jobStarted = True
    
    def startJob(self,cSession=None,clusterArgs=None):
        '''singleJob.startJob: Method to start jobs. For local runs,
        start process using base class and this object's args.
        For cluster runs, call startJobOnCluster and return.
        '''
        if self.onCluster:
            self.startJobOnCluster(cSession,clusterArgs=clusterArgs)
            return
        
        jobArgs = []
        if self.time:
            jobArgs += self.timeArgs
        if self.perf: 
            jobArgs += self.perfArgs
        jobArgs += self.args
        if self.stdOutFile:
            try:
                f1 = open(self.stdOutFile, 'w')
                stdout = f1
                self.StdOutFileOpen = True
                self.stdOutObject = f1
            except:
                stdout = open(os.devnull, 'w')
        else:
            stdout = open(os.devnull, 'w')
        
        if self.stdErrOutFile:
            try:
                e1 = open(self.stdErrOutFile, 'w')
                stderr=e1
                self.StdErrOutFileOpen = True
                self.stdErrOutObject = e1
            except:
                stderr=subprocess.PIPE
        else:
            stderr=subprocess.PIPE
        super(singleJob, self).__init__(jobArgs, stdout=stdout, stderr=stderr)
        self.startTime = time.time()
        self.isRunning = True
        self.jobStarted = True
        if os.name=="nt":
		self.nt_handle=ctypes.windll.kernel32.OpenProcess(SYNCHRONIZE, False, self.pid)
		global nt_handles_list
		nt_handles_list[self.nt_handle]=self.pid
        
    def addContingentJob(self, contingentJob):
        self.hasContingentJob = True
        self.contingentJob = contingentJob
    
    def killJob(self):
        self.kill()
        time.sleep(1)
        self.isRunning = False
        if isinstance(self.poll(),int):
            return 0
        print( ' JOB KILLED: ' + self.jobName)
        return 1

    #self.returncode is set by subprocess. However, poll or os.wait() sets it to None.
    #Reset it for use in CheckIfRunning, if the pid (also set by subprocess) matches the argument.
    #Return whether it matches.
    def markCompleted(self, pid, returnCode=None):
        if self.pid == pid :
            self.returncode=returnCode
            return True #this job's return code was set
        else :
            return False #not this job
        
    def CheckIfRunning(self,cSession=None,sleepTime=0.05):
        if self.isRunning:
            a=self.poll()
            if self.onCluster:
                hasExited = False
                exitStatus = -1
                while True:
                    try:
                        if not (cSession.jobStatus(self.jobID) in [drmaa.JobState.DONE, drmaa.JobState.FAILED]):
			    time.sleep(sleepTime)
                            return self.isRunning
                        retval = cSession.wait(self.jobID, drmaa.Session.TIMEOUT_NO_WAIT);
                        hasExited  = retval.hasExited
                        exitStatus = retval.exitStatus
                        break
                    except (drmaa.DrmCommunicationException), e:
                        # Catch drmaa communication exception, log a warning message and
                        # continue to monitor job.
                        print("*** CommERROR: %s Could not communicate with the cluster scheduler to check job status." %(str(e) ))
                        time.sleep(5)
                    except Exception, e:
                        #We found a case where an 'Exception' is raised: 
                        #   File "/usr/lib/pymodules/python2.7/drmaa/errors.py", line 90, in error_check
                        #     raise _ERRORS[code-1]("code %s: %s" % (code, error_buffer.value))
                        # Exception: code 24: no usage information was returned for the completed job
                        #This must be handled differently because must break (SGE did respond, but can't find job)
                        hasExited = True
                        exitStatus = 1000 #arbitrary
                        print("*** CommERROR: %s cluster scheduler lost job." % str(e))
                        break
                #end while True

                #if exit code != 0, restart without delay in checkStdOut bc that can cause Pipeline to wait for hours
                if self.remainingRestarts > 0 and hasExited and exitStatus != 0 :
                    if exitStatus == 255 :
                        utilities.LogError("warning", "job ssh failed (exit code %i), restarting, stdout=%s" % (exitStatus, self.expectedStdoutFile))
                        self.remainingRestarts += 1 #never fail a job because of ssh failing: may result in infinite loop, but this would only occur if the cluster cannot run the job at all
                        time.sleep(sleepTime)
                    else :
                        utilities.LogError("warning", "job has non-zero exit code (%i), restarting, stdout=%s" % (exitStatus, self.expectedStdoutFile))
                    self.restartJobOnCluster()
                    return self.isRunning
                        
                # Check whether the job is successful, if not restart it
                if self.remainingRestarts > 0 and self.expectedStdoutFile:
                    if not utilities.checkStdOut(self.expectedStdoutFile):
                        utilities.LogError("warning", "job was restarted, see stdout=\"%s\"" % self.expectedStdoutFile)
                        self.restartJobOnCluster()
                        return self.isRunning
                        
                cSession.deleteJobTemplate(self.jTemplate)
                self.StdErrOutFileOpen = True

            #elif not(isinstance(self.poll(),int)):
            elif not(isinstance(a,int)):
                #print a
                return self.isRunning
        self.endTime = time.time()
        self.runTime = self.endTime - self.startTime
        self.isRunning = False
        self.isComplete = True
        #These lines have been removed becuase it isn't necessary to call them inside
        # jobWrapper.multiThreadRunJobs (which is where this method is called). Instead, 
        # move these calls to jobWrapper.makeRunReport.
        #Revert back to original behavior because big bear doesn't like this.
        #self.CheckIfFileFound() #on third thought, big bear doesn't like this one either
        #self.GetProcBlock()
        #parseString = self.makeParseString()
        #above two are condensed into this:
        self.GetPerformanceString()
        return self.isRunning
    

    def GetPerformanceString(self, headerOnly=False):
    #def GetProcBlock(self, headerOnly=False): #rename this method
        """ Get performance stats generated by either time or perf
        from StdErrOutObject file or from self.communicate() 
        
        """

        #the header definition used to be in self.makeParseString()
        if headerOnly:
            headerString = '\tJob       \tQtime   \t'
            if self.time:
                headerString = 'time_header' + headerString + self.timeHeader
            if self.perf:
                headerString = 'time_perf_header' + headerString + self.perfHeader
            return headerString
        #end header definition

        #if this method is called twice, the second communicate call will fail
        # --don't want to overwrite self.parseString in this case
        elif self.parseString :
            return self.parseString

        nLines = 3          # take 3 lines to be sure to get time output
        if self.perf:
            nLines += 12    # take 12 more to be sure to get perf output
        
        textFound = True
        
        if self.StdErrOutFileOpen:
	    if self.stdErrOutObject!='':
		self.stdErrOutObject.close()
            self.StdErrOutFileOpen = False
            try:
                f1 = open(self.stdErrOutFile, 'r')
                try:
#                    procBlockList = [x.strip() for x in f1.readlines()[-nLines:] if x.strip()]
                    procBlockList = [x.strip() for x in f1.readlines() if x.strip()]
                except:
                    textFound = False
                f1.close()
            except:
                textFound = False
        else:
            jobCommunicate = None
            try :
                jobCommunicate = self.communicate()[1]
            except :
                textFound = False
            if not jobCommunicate :
                 textFound = False
            else:
                procBlockString = jobCommunicate
                procBlockList = [x.strip() for x in procBlockString.split('\n') if x.strip()]
        
        timeString = ""
        self.errString = ""
        if textFound:
            if self.time:
#                print 'jobName=' + self.uniqueString # DEBUG
                dateString = procBlockList.pop()
#                print 'dateString=' + dateString # DEBUG

                # check if dateString is really "Connection to phiNN closed"
                dateList = [x.strip() for x in dateString.split() if x.strip()]
#                print 'dateList[0]=%s\n' % dateList[0] # DEBUG
                if ( dateList[0] == 'Connection' ) :
                    self.HOSTname = dateList[2]
                    dateString = procBlockList.pop()
#                    print 'dateString=' + dateString # DEBUG
                else: # scan back for line starting with 'core file size'
                    lines = len(procBlockList)
                    lines += -1
                    while lines >= 0 :
                      hostLine = procBlockList[lines]
                      if "core file size" in hostLine :
                          break
                      lines += -1
                    
                    if ( lines > 0 ) :
                        self.HOSTname = procBlockList[lines-1]
                    else :
                        self.HOSTname = 'host=NA'
#                    print 'HOSTname=' + self.HOSTname # DEBUG

                cpuString = procBlockList.pop()
#                print 'cpuString=' + cpuString # DEBUG
                cpuList = [x.strip() for x in cpuString.split() if x.strip()]
                if ( cpuList[0] == 'Appending' ) : # assume there is no dateString since job failed
                    cpuString = dateString
                    cpuList = [x.strip() for x in cpuString.split() if x.strip()]
                    dateString = ""
                else : # look for message "Command terminated by signal N" 
                    if ( cpuList[0] == 'Command' ) : # assume there is no dateString since job failed
                        self.errString = cpuString
                        cpuString = dateString
                        cpuList = [x.strip() for x in cpuString.split() if x.strip()]
                        dateString = ""
                        
                if ( utilities.getFloatFromString(cpuList[0]) == None or utilities.getFloatFromString(cpuList[1]) == None ) : # give up
                    print 'Invalid cpuString=' + cpuString # DEBUG
                else :
#                print 'cpuList=' + str(cpuList) # DEBUG
#                print 'cpuList[0]=%s\n' % cpuList[0] # DEBUG
#                print 'cpuList[1]=%s\n' % cpuList[1] # DEBUG
#                print 'cpuList[3]=%s\n' % cpuList[3] # DEBUG
                    self.cpuTime = float(cpuList[0]) + float(cpuList[1])
                    cpuPercent = cpuList[3].split('%')
#                print 'cpuPercent=' + str(cpuPercent) # DEBUG
#                print 'cpuPercent[0]=%s\n' % cpuPercent[0] # DEBUG
                    if ( utilities.getFloatFromString(cpuPercent[0]) == None ) :
                        print 'Invalid cpuPercent=' + cpuPercent[0] # DEBUG
                        self.cpuPercent = 0
                    else :
                        self.cpuPercent = float(cpuPercent[0])
                    timeString = cpuString + "\t" + dateString
                    timeString = timeString.replace("\"", "") #for some reason, these now have opening and closing quotes--remove
            if self.perf:
                self.GetPerfStringFromBlockList(procBlockList[-10:])

        #code below used to be in self.makeParseString()
        startTimeInsert = '\t%3.5f' % self.runTime #note: Includes time waiting on the SGE queue, while timeString only includes elapsed time while running
        parseString = 'time_'+('perf_' if self.perf else "")+'line\t'+ ( '%-10s' % self.uniqueString ) + startTimeInsert
        if self.time:
            parseString += '\t' + timeString
        if self.perf:
            parseString += '\t' + self.perfString
        self.parseString = parseString + '\n'
        return self.parseString
        
        
    def GetPerfStringFromBlockList(self, perfBlockList):
        """ Parse the performance details from the stderr output file or pipe
        
        """
        
        self.perfString = ''
        self.perfStringHeader = ''
        for line in perfBlockList:
            tokens = line.split('\t')
            try:
                val = int(tokens[0])
                self.perfString += '%d\t' % val
            except:
                try:
                    val = float(tokens[0])
                    self.perfString += '%10.5f\t' % val
                except:
                    self.perfString += 'N/A\t'
            if len(tokens) > 1 :
                self.perfStringHeader += '%s\t' % tokens[1]
        self.perfString = self.perfString.rstrip()
        self.perfHeader = self.perfStringHeader.rstrip()        
    
    
    def CheckIfFileFound(self):
        """Was the expected output file created?
        
        """
#        print "\t Checking expectedResultFile:", self.expectedResultFile # DEBUG
        if self.expectedResultFile and os.path.exists(self.expectedResultFile):
            self.resultFound = True
        #if it doesn't exist, you can't copy it
        #elif self.forceForward:
        #    shutil.copy(self.expectedResultFile, self.forceForward)
        

    #copy the defaults from utilities.checkStdOut
    def CheckStdout(self, tries=20, delay=5):
        if self.expectedStdoutFile :
            self.stdoutComplete = utilities.checkStdOut(self.expectedStdoutFile, tries=tries, delay=delay)
            if not self.stdoutComplete:
                utilities.LogError("error", "job has not completed, see stdout=\"%s\"" % self.expectedStdoutFile)


    def makeReportString(self):
        """Report on the success or failure of this job.
        
        """
        #if not(self.isComplete):
        #    reportString = '   ' + self.jobName + ' Job Not Completed'
        timerString = '%6.1f' % self.runTime
        foundString = '  Result Found'
        if self.expectedResultFile and not self.resultFound:
            foundString = '  Result NOT Found: %s ' % self.expectedResultFile
        if self.expectedStdoutFile and not self.stdoutComplete:
            foundString += '  Stdout INCOMPLETE: %s ' % self.expectedStdoutFile
        reportString = '   ' + self.jobName + ' completed in: ' + timerString + 's'
        if self.expectedResultFile :
            reportString += foundString
        if self.expectedResultFile and not self.resultFound :
            reportString += '\n   OFFENDING ARGUMENTS:\n'
            reportString += ' '.join(self.args)
        return reportString


    def makeReportStringSimple(self):
        """Like makeReportString, but more brief: empty for success, and simple error message for fail.
        Must call checkStdout otherwise may show fail on success.
        """
        if self.expectedStdoutFile and not self.stdoutComplete : #favor stdout
            return "Job failed: %s\n" % self.expectedStdoutFile
        elif self.expectedResultFile and not self.resultFound :
            return "Job failed: %s\n" % self.expectedResultFile
        else :
            return "" #empty string            

#end class singleJob
        
def timeFormat1(inputSeconds):
    clockTimeH = int((1/3600.) * inputSeconds)
    clockTimeM = (1/60.)*(inputSeconds%3600)
    return '%2dh %0.2fm' % (clockTimeH, clockTimeM)
        
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
