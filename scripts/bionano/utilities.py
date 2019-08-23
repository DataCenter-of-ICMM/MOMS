
#utility functions for pipeline (mostly just in CharacterizeModule)

import os, sys
import re
import math
import string
#import atexit
import socket
import time
import datetime
from collections import defaultdict

"""@package utilities General methods for analysis

"""


#SVN utilties

global version
version=[]

#store argument in global version
def setVersion(v):
    global version
    version.append(v)

#parse version strings in global version; return ele which corresponds to max version
#this is only called once, in Pipeline.py
def getMaxVersion():
    global version
    maxi = maxv = 0
    for i,vstr in enumerate(version) :
        try :
            ver = int(vstr.split()[2])
        except :
            continue
        if ver > maxv :
            maxv = ver
            maxi = i
    return version[maxi]


setVersion("$Id: utilities.py 6328 2017-05-03 22:59:20Z wandrews $") #use above for this file


#general utilities

def uniquifyContigName(oname):
    iname = oname.upper()
    if iname != oname:
        return iname
    iname = oname.lower()
    if iname != oname:
        return iname
    return 'tempMerge'

#because the even len case returns a float, do the same for odd len case
def getMedian(inlist):
    if not len(inlist) :
        return 0.
    insort = sorted(inlist) #copies list
    inlen = len(insort)
    if not inlen%2 : #even number of eles
        return (insort[inlen/2-1] + insort[inlen/2]) / 2.
    return float(insort[inlen/2]) #odd number of eles


def getMAD(inlist, median=None) :
    '''MAD is median absolute deviation: get absolute difference of each element from median, MAD is median of that list.'''
    if median == None :
        median = getMedian(inlist)
    return getMedian([abs(x - median) for x in inlist])


#These two moved here from mapClasses.py
def getRMSsquared(inlist):
    if len(inlist) == 0 : #otherwise it's an exception
        return 0
    sumsquared = sum( map(lambda x: x**2, inlist) )
    return sumsquared/len(inlist)

def getRMS(inlist):
    return math.sqrt( getRMSsquared(inlist) )


#n50 algo from here:
#http://code.google.com/p/biopyscripts/source/browse/trunk/seq/N50.py
#Stuglik is apparently the guy's name
def getn50stuglik(inlist) :
    teoN50 = sum(inlist) / 2.0    
    #inlist.sort() #this is in getn50
    inlist.reverse() 
   
    testSum = 0; n50 = 0
    for leni in inlist:
        testSum += leni
        if teoN50 < testSum:
            n50 = leni
            break
    return n50


#Sorts list, calls getn50stuglik (was getn50sorted)
# Input list must be positive ints or floats only
def getn50(inlist) :
    if type(inlist) != list :
        print "Error in getn50: argument must be a list:", type(inlist)
        return 0
    inlist.sort()
    return getn50stuglik( inlist ) #much faster...
    #return getn50sorted( inlist ) #my algo sucks, performance wise


#getn50sorted :
# get the n50 for the input list
#  n50 is the value such that the sum of the elements < that value is equal to the sum >.
# for this, input list must be sorted--use getn50 above for unsorted

#Note on n50 algo for case of floats:
# It is very rare that sumL == sumR for floats. Therefore, for the rare case in which casting to an int would result in sumL == sumR, we'll end up choosing the lastL or lastR instead of the average of the two.
# This is a rare enough case that I'm going to ignore it. It only matters when there's a large difference between neighboring elements near the median of the expanded list.

#misc notes
# generalize to nX? ie, n50 is x = 50--another arg???

def getn50sorted(inlist) :
    if type(inlist) != list :
        print "Error in getn50sorted: argument must be a list:", type(inlist)
        return 0
    if len(inlist) == 0 :
        return 0 #nothing to do
    if len(inlist) == 1 :
        return inlist[0] #only one ele means return it
    printdebug = False
    sumL = inlist[0]
    sumR = inlist[-1]
    lastL = 1 #index on L--_next_ ele to use
    lastR = -1 #index on R (ele already used)--this will always count backwards
    try :
        sumall = sum(inlist) #I think these may save multiple calls to these guys
    except TypeError :
        print "Error in getn50sorted: Bad types in input list", inlist
        return 0
    lenall = len(inlist) #I think this can't fail for any list
    while True :
    #for ele in inlist[1:-1] : #exclude first, last since already added
        #try the outter loop a while true, and the inner a for, break
        #while sumL <= sumR :
        for idx, ele in enumerate(inlist[lastL:lastR]) :
            sumL += ele
            newlastL = idx #don't want to overwrite bc used in for loop
            if printdebug :
                print "sumL =", sumL, "lastL =", newlastL+lastL+1 #debug
            #require > bc makes R easier to handle--also need to update lastL if end of loop
            # len of this loop is lenall - lastL - abs(lastR)
            #  And -abs(lastR) = lastR since lastR < 0, and -1 bc index vs len
            if sumL > sumR or idx == lenall - lastL + lastR - 1 : 
                lastL += newlastL + 1 #n old + n new (need additional +1?)
                break
        #end loop to increase sumL, lastL
        neleused = lastL + abs(lastR) #for L, index of next ele is n ele used on L
        if printdebug :
            print "neleused =", neleused
        #error check--n eles
        if neleused > lenall or neleused < 0 :
            print "Error in getn50sorted: bad indicies", neleused
            return 0
        if neleused == lenall : #no eles left--recall lastR < 0
            if sumL == sumR :
                #avg of last two (recall, need to take away 1 from lastL)
                return (inlist[lastL-1] + inlist[lastR])/2 
            elif sumL > sumR :
                return inlist[lastL-1] #if exceeded, last ele is n50 (lastL is next, so -1)
            #sumL < sumR if middle happens to be low--nothing wrong with this
            elif sumL < sumR :
                return inlist[lastR] #if R exceeded, need lastR
        #more eles left--increase R
        lastR -= 1
        sumR += inlist[lastR]
        if printdebug :
            print "sumR =", sumR, "lastR =", lastR #debug
        #error check--bail if summed too much
        if sumR > sumall or sumL > sumall :
            print "Error in getn50: sum excess"
            return 0
        #error check--bail if sum < 0
        if sumR < 0 or sumL < 0 :
            print "Error in getn50: negative sum"
            return 0

#end def getn50sorted



#Given a list of pairs of ints/floats, return a new list of pairs which have no overlap.
# Pairs here means a list of two elements (if this isn't satisfied, skip it, but print an error).
# No overlap means that each pair specifies a range. So two ranges which overlap should be merged to a single range.
#
#Other notes:
# Allow modification of the input list--this will make implementation easier, and will also consume less memory. Of course, the downside is that after calling this function, you can't assume your list is the same.
#  Therefore, I will remove any items (pairs) which are invalid--see error handling at beginning of first loop
# I'm not doing any checking for negative numbers, so best not to input them
# For implementation notes, see inline comments

def uniqueRange(inlist) :
    printerr = True #this amounts to a verbose flag

    if type(inlist) != list :
        if printerr :
            print "Error in uniqueRange: input must be a list:", inlist
        return #return none is fine--list is merged in place

    #Need to check that all the pairs are sorted, otherwise, can't really trust the results of below--this is a bit time consuming, but I think it's necessary
    #store idxs to delete bc it's very very bad to delete from a list while you're looping over it
    badidxs = [] 
    for idx, ele in enumerate(inlist) :
        if type(ele) != list :
            if printerr :
                print "Warning in uniqueRange: skipping non-list ele:", ele
            badidxs.append( idx )
            continue
        if len(ele) != 2 :
            if printerr :
                print "Warning in uniqueRange: skipping list of improper size:", ele
            badidxs.append( idx )
            continue
        if not isIntorFloat(ele[0]) or not isIntorFloat(ele[1]) :
            if printerr :
                print "Warning in uniqueRange: skipping list ele of improper type:", ele[0], ele[1]
            badidxs.append( idx )
            continue
        #now that we're sure that all the pairs are good, sort them
        ele.sort()

    #delete bad eles--see first 'Other notes' above
    #I think if you reverse it, it's safe bc will change from end to beginning--assumes badidxs is sorted, but it must be because of how it's constructed, right?
    badidxs.reverse() 
    for idx in badidxs : del inlist[idx]

    #need to now sort the sublists within inlist (pretty cool, eh?)
    #and it doesn't matter if the second ele is not sorted
    inlist.sort(key=lambda x: x[0])

    #print inlist #start with this (sorted)--debug

    #now merge overlapping entries
    badidxs = [] #clear this list--this should be fine, I think
    for idx, ele in enumerate(inlist[:-1]) : #leave off last ele
        nextele = inlist[idx+1]
        #if next element starts before (or at) current one ends, merge the pairs
        #careful--need to cast to floats, otherwise comp won't work
        if float(nextele[0]) <= float(ele[1]) : 
            #print "merging", ele, "and", nextele, #debug/dev
            nextele[1] = max(float(nextele[1]), float(ele[1])) #next ends as far as it can
            #need to check if next ele starts after current--if so, it needs to start at current start
            if float(nextele[0]) >= float(ele[0]) :
                nextele[0] = ele[0]
            badidxs.append( idx ) #discard current in favor of next
            #print "to", nextele #debug/dev (new is next)

    #now remove the elements which were merged
    badidxs.reverse() 
    for idx in badidxs : del inlist[idx]

    #print inlist #this is here for development purposes--remove later
    #the input list was modified in place--no return

#end def uniqueRange(inlist)



#simple helper fn for above (not so simple...see below)
def isIntorFloat(ele) :
    #they say that checking type is a sign of bad polymorphism, but I never envision using that here
    #the new version doesn't check type, but it still has flaws: 12L fails both, but is a 'number'
    # I don't know how to do better...perhaps a cast?
    # for the case of "L" numbers, this works:
    #  import numbers
    #  isinstance(12L, numbers.Integral)
    return isinstance(ele, int) or isinstance(ele, float)
    #return type(ele) == int or type(ele) == float #old way--replace with above



#leave above, but make new version (int only, not float):
# isinstance is polymorphism safe and avoids these problems:
#  try: return ele == int(ele) #this doesn't work bc '1. == 1'; returns True
#  try: return newele = int(ele) #this doesn't throw an exception for ele is float; also returns True
def isInt( ele ) :
    return isinstance( ele, int )



#String Utilities


#check if object (string) has substring
def hasSubStr(instring, substr) :
    try :
        #find returns pos in string if found
        # ie, if found, get 1, which is not -1
        # if not found, get -1, and (-1 != -1) is False
        return ( instring.find(substr) != -1 )
    except :
        return False


#check if string represents an int using cast in try
# return the int if so; return None if not
# note that empty string will return None, not 0 (because it's not a valid int)
# also note that strings for floats will also return None because of how casting works in python
#  ie, int("1.0") raises ValueError
#  Do not put a "." in instring
#  however, you could do int(float("1.0")), but then there's also int(round(float("1.5"))) ...
def getIntFromString(instring) :
    try :
        inint = int(instring)
    except ValueError : #only catch invalid cast, which is a ValueError
        return None

    return inint #no exception
#end def getIntFromString

def getFloatFromString(instring) :
    try :
        val = float(instring)
    except ValueError : #only catch invalid cast, which is a ValueError
        return None

    return val #no exception



#get index of a file whose format is:
# 'path/filenameX.suff'
#  Note that 'path' is irrelevant--just the end of the string is used, beginning is ignored
#  and X (the index) is an integer up to (and including) 4 digits
#   If you want more than 4 digits, change the 'ndigits' local
# return the index as an int if valid, and return None if not
# depends on getIntFromString, which should return None for non-integers
# Note that multiple "." will use last, ignoring others
def getPathIndex(instring) :
    ndigits = 4 #if you want more digits, increase this

    #first, strip suff--if it doesn't exist, do nothing
    #usestr = os.path.split(instring)[1] #unnecessary
    dpos = instring.rfind(".")
    if dpos != -1 : #has a period
        usestr = instring[:dpos] #slice does not include second arg
    else :
        usestr = instring #copy so that above does not modify input (not sure if strings are passed by ref or value)

    #then, try each of 3, 2, 1 digits
    for idx in range(1,ndigits+1)[::-1] : #if ndigits is 3, do 3,2,1 (the -1 of slice goes backwards)
        useint = getIntFromString(usestr[-idx:])
        if useint != None : #this means getIntFromString succeeded
            break #since go from largest to smallest, break at first non-None

    #if getIntFromString never found an int, useint is None, but return it anyway
    return useint
#def getPathIndex(instring) :



#File Utilities



#check a file--filepath is full path, not just prefix
# if filesuff, check that filepath ends with this string
# if checkReadable, require user has read access
def checkFile(filepath, filesuff="", checkReadable=True) :
    try :
        valid = os.path.isfile(filepath)
        if filesuff != "" :
            valid = (valid and filepath.endswith(filesuff))
        if checkReadable :
            valid = (valid and os.access(filepath, os.R_OK))
    except :
        valid = False
    #print "checkFile:", filepath, valid #debug
    return valid
#end def checkFile



#check a cmap--use checkFile with filesuff == ".cmap"
#returns true if arg is a file that exists and ends in ".cmap"
def checkCmap(cmappath) :
    return checkFile(cmappath, ".cmap")
#end def checkCmap


#true if arg is file with suffix "cmap" OR "spots"
def checkCmapSpots(inpath) :
    return checkFile(inpath, ".cmap") or checkFile(inpath, ".spots")
#end def checkCmapSpots


#
# Convert string to C like representation
# to avoid printing control characters
#
def normalizeString(s):
	out=""
	for ch in s:
		if ch=="\"":
			out+="\\\""
			continue
		if ch=="\\":
			out+="\\\\"
			continue
		if ch=="\n":
			out+="\\n"
			continue
		if ch=="\r":
			out+="\\r"
			continue
		if ch=="\t":
			out+="\\t"
			continue
		if string.find(string.printable, ch)>=0:
			out+=ch
			continue
		out+="\\%03o" % tuple([ord(ch)])
	return out

# Note - the default value for NFS is to cache file attributes for a maximum of 1 minute
def checkStdOut(filename, end_pattern="END of output\n", tries=20, delay=5):
    line = ""
    try:
        f=open(filename, "rb")
        try : #a-la http://www.diveintopython.net/file_handling/file_objects.html
            f.seek(-len(end_pattern), 2) # seek relative to end of file
            line=f.read()
        finally : #this will ensure the file handle is closed even if seek or read raises an IOError
            f.close()
    except Exception,e:
		# If we could not open file, try again as there might be delay due to NFS
		if tries<2:
			LogError("warning", "Could not open file \"%s\" tries=%d exception=%s\n" %(filename, tries, e))
		if tries<1:
			return 0
		time.sleep(delay)
		return checkStdOut(filename, end_pattern, tries=tries-1, delay=delay)
    if line!=end_pattern:
        if tries>0:
            time.sleep(delay)
            return checkStdOut(filename, end_pattern, tries=tries-1, delay=delay)
        LogError("warning", "Missing end marker in \"%s\" (found \"%s\" while expecting \"%s\")\n" % (normalizeString(filename), normalizeString(line), normalizeString(end_pattern)))
        # Wait that long before restarting job
        return 0
    return 1
#end checkStdOut


global status_log
global status_log_filename 
status_log = None
status_log_filename = None
global start_time

def InitStatus(filename):
    global start_time
    start_time = time.time()
    global status_log
    global status_log_filename
    status_log_filename=filename
    status_log=open(filename, "wb")
    
    LogStatus("pipeline", "pid", "%d" % os.getpid())
    LogStatus("pipeline", "hostname", socket.gethostname())
    LogStatus("pipeline", "status_file", filename)
    LogStatus("pipeline", "user_home", os.path.expanduser('~'))
    LogStatus("pipeline", "version", getMaxVersion())
    try:
        LogStatus("pipeline", "start_iso", datetime.datetime.now().isoformat())
    except:
        pass
    try: #paranoia
        LogStatus("pipeline", "start_epoch", time.time())
    except:
        pass
    try:
        LogStatus("pipeline", "uid", os.getuid())
    except:
        pass
    #try to get ulimit -a stats in logs: this feature is not done...
    #print "ulimit"
    #retc, ulim, err = runJob(["ulimit -a"], returnstdout=True, shell=True) #bc returnstdout is True, return tuple
    #LogStatus("pipeline", "ulim", ulim) #multi-line is ok?
    #try:
    #    pass
    #except Exception,e: #can raise OSError, but catch 'em all
    #    print "Exception\n",e
    #    #pass

    #how to get the proc location: _proc_status = '/proc/%d/status' % os.getpid()

#end InitStatus


def LogStatus(category, field, value, aux1=None, aux2=None, aux3=None):
    global status_log
    global start_time
    aux=""
    if aux1:
	    aux+=" val1=\"%s\"" % aux1
    if aux2:
	    aux+=" val2=\"%s\"" % aux2
    if aux3:
	    aux+=" val3=\"%s\"" % aux3
    
    if status_log:
        status_log.write("<%s attr=\"%s\" val0=\"%s\"%s time=\"%.1f\"/>\n" % (category, field, value,aux, time.time()-start_time))
        status_log.flush()
    else:
        print "Could not log status: %s %s %s%s" % (category, field, value, aux)
	
	
global ErrorLog
ErrorLog=[]

def LogError(category, message):
	global ErrorLog
	ErrorLog.append((normalizeString(category), normalizeString(message)))

	if category=="warning":
		LogStatus("warning", "message", normalizeString(message))
	else:
		LogStatus("error", normalizeString(category), normalizeString(message))
		
def SummarizeErrors(varsP=None):
	global ErrorLog
	if len(ErrorLog)==0:
		print "No errors detected\n"
		if varsP!=None:
			varsP.updatePipeReport("No errors detected\n", printalso=False)
		return {} #RETURN DICT!
		
	categories=defaultdict(int) 
	error_count=warning_count=0
        logmsg=""
        if len(ErrorLog) :
            logmsg += "\nWarning/Error messages:\n"
	for i in range(0, len(ErrorLog)):
            logmsg += "%s : %s\n" % ErrorLog[i]
            cat=ErrorLog[i][0]
            if cat=="warning":
                warning_count+=1
            else : #if it's not a warning, it's an error
                error_count+=1
            categories[cat]+=1
        if warning_count and not error_count :
            logmsg += "\nWarning summary:\n"
        elif warning_count and error_count :
            logmsg += "\nWarning/Error summary:\n"
        elif error_count :
            logmsg += "\nError summary:\n"
	for cat in categories:
            logmsg += "\t%d %s(s)\n" % (categories[cat], cat) 

        print logmsg
	if varsP!=None:
            varsP.updatePipeReport(logmsg, printalso=False) #already printed
	#return error_count #can't distinguish error/critical here
        return categories
		
		
#return a list of full paths to files in a dir
#if suffix supplied, only files with that suffix
#if not, will return all files, including sub-dirs
#if first arg is not a dir, return empty list
def getListOfFilesFromDir(indir, suffix="") :
    if not checkDir(indir, checkWritable=False, makeIfNotExist=False) :
        print "Error in getListOfFilesFromDir: argument is not a dir:", indir
        return []
    outlist = []
    for qfile in os.listdir(indir) :
        if suffix and not qfile.endswith(suffix) :
            continue
        outlist.append( os.path.join(indir, qfile) )
    return outlist
#end def getListOfFilesFromDir


#check if a file exists (checkFile) and is executable
def checkExecutable(filepath, filesuff=""):
    if not checkFile(filepath, filesuff):
        return False
    return os.access(filepath, os.X_OK)
#end def checkExecutable


#check if directory exists; make it if not
#optionally check writability--True by default
#make creation an option--True by default
#should probably make a vebose option
def checkDir(dirpath, checkWritable=True, makeIfNotExist=True) :

    #if it exists and it's not a directory, return False
    if os.path.exists(dirpath) and not os.path.isdir(dirpath) : 
        return False

    #this check was isdir, but actually, it should be existence--if it doesn't exist, make it
    elif not os.path.exists(dirpath) : #os.path.isdir(dirpath) : 
        if makeIfNotExist :
            try :
                os.makedirs(dirpath)
            except : 
                return False
            if not os.path.isdir(dirpath) : #too bad to handle
                #print "Error in utilities.checkDir: cannot create dir:", dirpath #silent
                return False
        else : #doesn't exist, and flag to make is false--return False
            return False

    #above is check of existence. Now check writable 
    if checkWritable and not os.access(dirpath, os.W_OK) : # W_OK is for writing, R_OK for reading, etc.
        return False

    return True #all checking passed
#def checkDir


#get all sub-directories (not files) in input dir
#returns empty list if no sub-dirs
def getSubDirs(indir) :
    if not checkDir(indir, makeIfNotExist=False) :
        print "Error in getSubDirs: input is not a dir"
        return None

    targetlist = []
    for qfile in os.listdir(indir) :
        targetdir = os.path.join(indir, qfile)
        if os.path.isdir(targetdir) :
            targetlist.append( targetdir )

    return targetlist
#end def getSubDirs(indir)


#strip the first dir and filename from a path
# ie, given 'first_dir/second_dir/.../n_th_dir/filename.suf'
#  you get 'second_dir/.../nth_dir'
# if no filename, get same
def stripFirstDir( inpath ) :

    #again, the 'isinstance' vs 'try' question
    if not isinstance(inpath, str) :
        print "Error in stripFirstSubDir: input not string:", inpath
        return None

    folders = getPathList( inpath ) 

    #only strip the last ele if it's a file
    #But I don't want to check file existence, just in case this is not a real path
    #so just discard if it has '.'
    if '.' in folders[-1] : #'in' operator for strings works just like it does for lists
        del folders[-1]

    #ok, now strip the first ele in the loop statement (always)
    outstr = ""
    for ele in folders[1:] : 
        outstr = os.path.join(outstr, ele)

    #print outstr #debug
    return outstr
    
#end def stripFirstSubDir( inpath ) :


#given a path, return a list of its elements
# ie, given 'first_dir/second_dir/.../n_th_dir/filename.suf',
#  you get ['first_dir', 'second_dir', ..., 'n_th_dir', 'filename.suf']
#from http://stackoverflow.com/questions/3167154/
# but slightly modified to handle case of trailing '/'
def getPathList( inpath ) :

    #again, the 'isinstance' vs 'try' question
    if not isinstance(inpath, str) :
        print "Error in stripFirstSubDir: input not string:", inpath
        return None
        
    folders=[]
    path = inpath #copy (?)
    while True :
        path, filename = os.path.split(path)

        if filename != "" :
            folders.append(filename)
        else:
            #only append if this is the last thing to append (ie, not first ending in /)
            if path != "" and os.path.split(path)[0] == '' : 
                folders.append(path)
            if os.path.split(path)[0] == '' :
                break
        #print path, filename, folders #debug

    folders.reverse()
    #print folders #debug
    return folders

#end def getPathList()


#get the number of files in a directory with a given suffix (ie, filename ends with substr)
def getNfilesDir(indir, suffix="") :
    #was assert, change to return 0 if indir doesn't exist
    if not checkDir(indir, checkWritable=False, makeIfNotExist=False) :
        return 0

    nfile = 0
    #note on listdir: you get everything you get with "-a" *except* "." and ".."
    # this means you do get ".xyz" files and "blah.suf~" files
    # And that is what you get if suffix is not supplied, because endswith("") is apparently always True
    for filename in os.listdir(indir) :
        if filename.endswith(suffix) :
            #print filename #debug
            nfile += 1

    return nfile
#end def getNfilesDir



#if a dir exists, _remove_ all its contents. Be careful.
#dependency: shutil (this is python std lib, so should be fine)
def removeDirContents(indir) :

    #check if it exists, return if not--do not create
    if not checkDir(indir, makeIfNotExist=False) :
        return None

    import shutil

    for files in os.listdir(indir): 
        shutil.rmtree( os.path.join(indir, files) ) #don't call on indir bc only want to remove its contents

#end removeDirContents


#get number of newlines in file
def wc(infile) :
    if not checkFile(infile) :
        print "Error in wc: bad file:", infile
        return
    nlines = 0
    for line in open(infile) :
        nlines += 1
    return nlines
#end wc



#list/dict utilities



#Get the total length which corresponds to the sum of all the ranges in a list of ranges
#**No checks for any failures are done.**
# Garbage in, garbage out.
#Note: requires math
def totalLengthFromRanges(inlist) :
    printerr = True #this amounts to a verbose flag

    totlen = 0
    if type(inlist) != list : #this should probably be a try...
        if printerr :
            print "Error in totalLengthFromRanges: input must be a list:", inlist
        return 0 #return int/float

    for ele in inlist :
        elelen = 0
        try :
            elelen = math.fabs(ele[1] - ele[0])
        except :
            pass
        #print "adding", elelen #debug
        totlen += elelen

    return totlen
#end totalLengthFromRanges



#wrap the fn below to print
# assume keys can be casted to strs and values are ints
#Note: sorting ints by their string values results in errors like 10 before 2.
# but chromosomes are mixed ints and strings...
def printSortedDict(adict) :
    keys, values = sortedDict(adict)
    for i in range(len(keys)) :
        print "%3s  %3i" % (keys[i], values[i])



#return tuple of keys, values such that keys are sorted
#from http://code.activestate.com/recipes/52306-to-sort-a-dictionary/, slightly modified
#note that there is an OrderedDict class in the built-in collections module
def sortedDict(adict):
    keys = adict.keys()
    keys.sort()
    #map returns a list, applying the first arg (a fn) to the list of second arg
    #dict.get is same as [], but this is apparently faster than doing [adict[key] for key in keys]
    return keys, map(adict.get, keys)


#dict.update will clobber existing entries if the argument contains the same keys
#first check for the keys, and if found, raise an exception
# if you want a bool, wrap this in try
#if no exception, do the update (dicts are passed by ref (sort of), so modify existing)
def safeDictUpdate(indict, updict) :
    if any([ indict.has_key(upkey) for upkey in updict.keys() ]) :
        #print indict, updict #debug
        #KeyError is raised if key doesn't exist, and here's it's raised if it does
        raise KeyError, "key in update dict already in original dict" 

    #I think else isn't necessary becuase if raise, exit this fn, but it seems safer this way
    else : 
        indict.update( updict )


#take two dicts; add all values in second (adict) to first (sumdict) when keys match
#In order to prevent raising a KeyError exception, check key
# if addkeys, then put the key in sumdict
# if not, ignore any values whose keys are not in sumdict
#Also, no return because after calling, the sumdict will be changed (like a pass by reference)
def addDictToDict(sumdict, adict, addkeys=True) :
    if not adict : #do nothing for empty/None adict
        return
    for key,value in adict.iteritems() :
        if sumdict.has_key(key) :
            sumdict[key] += value
        elif addkeys :
            sumdict[key] = value

# HERE : modify weight for extension jobs : count only the first 600kb of length, then add 250kb on each end for extensions
def contigWeight(fname):
	f=open(fname, "rb")
        try :
            f.seek(-512, 2)
        except :
            pass
	ctg=f.read().split("\n")
	f.close()
	
	i=len(ctg)
        if i<1:
            return 1e-16
	while i>=1:
		i-=1
		if ctg[i]!="":
			break
	line=ctg[i].split()
	
        try :
            length1=float(line[1])
            nsites=float(line[2])
        except ValueError :
            print "ERROR in contigWeight, file=%s line=%s" % (fname, line)
            return 1e-16
	density=nsites*1.0/length1
	if density< (1.0/10e3):
		density=1.0/10e3
	return density*density*nsites

def contigComputeWeights(file_list):
	L=[]
	sizesum=0
	for fname in file_list:
		#a=os.stat(fname).st_size
		a=contigWeight(fname)
		sizesum+=a;
		L.append(a)
	return sizesum, L



#bnxfile and molecule classes

class bnxfile:

    def __init__(self, sourcefile='', thresh=None, doqual=False): #sourcefile is .bnx
        #init data members
        self.moleculelist    = [] #list of molecule objects
        self.moleculeIdlist  = []
        self.moleculeLenlist = []
        self.molstats        = {}
        self.genomesizemb    = -1. #guess genome size from sourcefile
        self.header          = ""
        #self.bnxVer          = "" #empty will default to 0.1 -- not used
        self.doQuality       = doqual
        self.sourceFile      = ""
        self.lengthThresholds = (thresh if thresh else [100, 150, 180, 200])
        if checkFile(sourcefile, ".bnx") : #must end in ".bnx" -- if migrate this out of utilities, need util.checkFile
            self.makeFromFile(sourcefile)
            self.calculateMolStats()
        elif sourcefile : #if argument supplied but the file is not found, print a warning
            print "Warning in bnxfile init: sourcefile not found:", sourcefile

    def makeFromFile(self, sourcefile) :
        self.moleculelist = []
        #try to guess genome size based on sourcefile
        if hasSubStr(sourcefile.lower(), "human") :
            self.genomesizemb = 3200.
        elif hasSubStr(sourcefile.lower(), "horse") :
            self.genomesizemb = 2800.
        elif hasSubStr(sourcefile.lower(), "zebrafish") :
            self.genomesizemb = 1400.
        elif hasSubStr(sourcefile.lower(), "chicken") :
            self.genomesizemb = 1000. #rounding down from 1064 total in the fasta
        elif hasSubStr(sourcefile.lower(), "tribolium") :
            self.genomesizemb = 200
        elif hasSubStr(sourcefile.lower(), "fly") or hasSubStr(sourcefile.lower(), "dros") :
            self.genomesizemb = 140
        elif hasSubStr(sourcefile.lower(), "coli") : #use ecoli for tests
            self.genomesizemb = 4.6
        with open(sourcefile) as f1 :
            mollines = ""
            for line in f1 :
                if line[0] == "#" : #header lines
                    self.header += line
                    #if line.lower().find("version") != -1 : #not used
                    #    self.bnxVer = line.split()[-1]
                    continue

                elif line[0] == "0" : #0h Label Channel	MapID	Length
                    if mollines : #at first molecule, this is empty
                        self.moleculelist.append( molecule(mollines) ) #store previous molecule
                    mollines = line #replace with current line

                else : #for any other line type, append to mollines
                    mollines += line

                #if len(self.moleculelist) > 3 : #test
                #    break
            #end for line in f1
            self.moleculelist.append( molecule(mollines) ) #need last molecule
        #end with open sourcefile
        self.sourceFile = os.path.split(sourcefile)[1] #store sourceFile only after loading all data
    #end makeFromFile


    #make a list of just the molecule ids from the molecule list
    def makeMolIdList(self, minlenkb=0, minnlab=0) :
        assert minlenkb >= 0 and minnlab >= 0, ("Invalid arguments: %f %f" % (minlenkb, minnlab))
        if not self.moleculelist :
            self.moleculeIdlist = []
            return
        if minlenkb == 0 and minnlab == 0 :
            self.moleculeIdlist = [x.id for x in self.moleculelist]
        elif minnlab == 0 : #minlenkb > 0
            self.moleculeIdlist = [x.id for x in self.moleculelist if x.lenkb >= minlenkb]
        elif minlenkb == 0 : #minnlab > 0
            self.moleculeIdlist = [x.id for x in self.moleculelist if len(x.labposeskb) >= minnlab]
        else : #both > 0
            self.moleculeIdlist = [x.id for x in self.moleculelist if x.lenkb >= minlenkb and len(x.labposeskb) >= minnlab]
    #end makeMolIdList(self) 


    #make a list of just the molecule ids from the molecule list;
    # filter length > minlenkb, n labels > minnlab
    def makeMolLenList(self, minlenkb=0, minnlab=0) :
        assert minlenkb >= 0 and minnlab >= 0, ("Invalid arguments: %f %f" % (minlenkb, minnlab))
        if not self.moleculelist :
            self.moleculeLenlist = []
            return
        if minlenkb == 0 and minnlab == 0 :
            self.moleculeLenlist = [x.lenkb for x in self.moleculelist]
        elif minnlab == 0 : #minlenkb > 0
            self.moleculeLenlist = [x.lenkb for x in self.moleculelist if x.lenkb >= minlenkb]
        elif minlenkb == 0 : #minnlab > 0
            self.moleculeLenlist = [x.lenkb for x in self.moleculelist if len(x.labposeskb) >= minnlab]
        else : #both > 0
            self.moleculeLenlist = [x.lenkb for x in self.moleculelist if x.lenkb >= minlenkb and len(x.labposeskb) >= minnlab]
    #end makeMolLenList(self) 


    #inits data members for statistics of molecule data
    #for length thresholds specified in thresh local
    #use dictionary whose keys are the thresholds and values are the bnxstats class
    def calculateMolStats(self) :
        self.molstats = {}

        #these are no longer used
        #if not self.moleculeIdlist :
        #    self.makeMolIdList()
        #if not self.moleculeLenlist :
        #    self.makeMolLenList()

        for thresh in self.lengthThresholds :
            #default for printing len is Mb, change here for other (see bnxstats.__str__)
            stat = bnxstats(genomesizemb=self.genomesizemb) 
            #stat.makeFromLenList( [x for x in self.moleculeLenlist if x > thresh] ) #moleculeLenlist is assumed to be in kb
            stat.makeFromMolList( [x for x in self.moleculelist if x.lenkb >= thresh] )
            self.molstats[thresh] = stat
    #end calculateMolStats
    

    #print results of calculateMolStats
    def printMolStats(self) :
        if not self.molstats :
            self.calculateMolStats()
        #for thresh, stats in self.molstats.iteritems() : #not sorted
        for thresh in sorted(self.molstats.keys()) : #loop over sorted keys
            print "Len >= %i (kb)" % thresh
            print self.molstats[thresh] #stats
    #end printMolStats

    def getTotLenMb(self) :
        if not self.molstats :
            self.calculateMolStats()
        key = self.molstats.keys()[0]
        return self.molstats[key].totlen/1e3 #totlen is in kb, convert to Mb
    #end getTotLenMb

    #check for duplicate molecule: if all lab poses are the same, it's a duplicate; print summary
    #warning: this is very slow for large bnxs
    def checkDuplicateMolecule(self, minlenkb=0., verbose=True) :
        ndupe = 0
        if minlenkb > 0 :
            filteredlist = [x for x in self.moleculelist if x.lenkb > minlenkb]
        else :
            filteredlist = self.moleculelist
        nmols = len(filteredlist)
        if verbose :
            print self.sourceFile, "N mols =", nmols
        sys.stdout.flush()
        for idx1,mol1 in enumerate(filteredlist) :
            if idx1 % 10000 == 0 :
                if verbose :
                    print self.sourceFile, ":", idx1
                sys.stdout.flush()
            for idx2 in range(idx1+1,nmols) :
                mol2 = filteredlist[idx2] #is indexing faster than slice? Slice creates copy.
                if math.fabs(mol1.lenkb - mol2.lenkb) > 10. : #len is within 10kb
                    continue
                if len(mol1.labposeskb) != len(mol2.labposeskb) : #same n labels
                    continue
                #if mol1.labposeskb == mol2.labposeskb : #instead of equality, use getMaxAbsDev
                #need a threshold: 1kb
                if mol1.getMaxAbsDev(mol2) < 1. : 
                    ndupe += 1
                    print self.sourceFile, "Duplicate molecules:", mol1.id, mol2.id
                    sys.stdout.flush()
        #end outer mol loop
        if not ndupe :
            print self.sourceFile, "No duplicates found"
        else :
            print self.sourceFile, "N duplicates found:", ndupe
    #end def checkDuplicateMolecule

    #like above, but compare two bnxfile objects -- they must be pre-filtered
    #here, all molecules in first must be compared to all in second because
    # there is no overlap between the two
    def checkDuplicateBnx(self, inbnxfile, verbose=True) :
        ndupe = 0
        if verbose :
            print self.sourceFile, "N mols1 =", len(self.moleculelist)
            print inbnxfile.sourceFile, "N mols2 =", len(inbnxfile.moleculelist)
        sys.stdout.flush()
        for idx1,mol1 in enumerate(self.moleculelist) :
            if idx1 % 10000 == 0 :
                if verbose :
                    print self.sourceFile, inbnxfile.sourceFile, ":", idx1
                sys.stdout.flush()
            for mol2 in inbnxfile.moleculelist :
                if math.fabs(mol1.lenkb - mol2.lenkb) > 10. : #len is within 10kb
                    continue
                if len(mol1.labposeskb) != len(mol2.labposeskb) : #same n labels
                    continue
                #if mol1.labposeskb == mol2.labposeskb : #instead of equality, use getMaxAbsDev
                #need a threshold: 1kb
                if mol1.getMaxAbsDev(mol2) < 1. : 
                    ndupe += 1
                    print self.sourceFile, inbnxfile.sourceFile, "Duplicate molecules:", mol1.id, mol2.id
                    sys.stdout.flush()
        #end outer mol loop
        if not ndupe :
            print self.sourceFile, inbnxfile.sourceFile, "No duplicates found"
        else :
            print self.sourceFile, inbnxfile.sourceFile, "N duplicates found:", ndupe
    #end def checkDuplicateBnx

    #compare molecule IDs in self to those in argument, print any matches
    def compareMoleculeIds(self, inbnx) :
        ndupe = 0
        inidlist = [x.id for x in inbnx.moleculelist]
        for mol in self.moleculelist :
            if mol.id in inidlist :
                ndupe += 1
                print "Molecule in both bnx:", mold.id
        if not ndupe :
            print "No duplicates found"
        else :
            print "N duplicates:", ndupe
    #end compareMoleculeIds


    #for old ImageProcessing id enconding _only_:
    # <chip><flow cell: 2 digits><mol id: 6 digits>
    #make dict of chip number to number of molecules
    # so strip the last eight digits, and that's the key
    #if verbose, warn if mol id doesn't conform to above format, ie, too short
    def getNmolProfile(self, minlenkb=0, minnlab=0, verbose=True) :
        self.makeMolIdList() #populate self.moleculeIdlist
        mprof = defaultdict(int)
        for mid in self.moleculeIdlist :
            if mid < 1e8 :
                if verbose :
                    print "Warning in getNmolProfile: molecule id %i is too small (< 1e8)" % mid
                continue
            mprof[ int(mid/1e8) ] += 1
        return mprof
    #end getNmolProfile


    #I'm currently only supporting bnx version 0.1
    #which means the header must reflect this--fix the header
    def setHeaderVersion(self, newversion="0.1"):
        verstr = "BNX File Version"
        headlines = self.header.split("\n")
        headlen = len(headlines)
        #to do this properly, you should loop over every line
        #for line in headlines :
        #    if line.find(verstr) != -1 :
        #to do it easily, assume the first line
        topline = headlines.pop(0) #removes 0th element
        if topline.find(verstr) == -1 : #if this isn't in first line, do nothing
            return
        topline = "# "+verstr+":\t0.1"
        headlines[0:0] = [topline] #insert new topline back into header
        self.header = "\n".join(headlines)
        #assert headlen == len(self.header.split("\n")) #check new vs orig len
    #end setHeaderVersion

    def writeToFile(self, filename) :
        self.setHeaderVersion() #this is only necessary for reading a version 1 bnx
        if not filename.endswith(".bnx"):
            print "Error in bnxfile.writeToFile: file to write must end in .bnx"
            return
        outfile = open(filename, "w+")
        outfile.write(self.header)
        for mol in self.moleculelist :
            outfile.write(mol.printBnx())
        outfile.close()
    #end writeToFile

#end class bnxfile


#class to store stats for bnx files

class bnxstats :

    #see __str__ for explanation of unitscale
    def __init__(self, nmol=0, totlen=0, n50=0, genomesizemb=-1, unitscale=1e3, totlab=0) :
        self.nmol         = nmol
        self.totlen       = totlen
        self.n50          = n50
        self.unitscale    = unitscale
        self.genomesizemb = genomesizemb
        self.totnlab      = totlab
        self.labdensity   = (totlab*1e2/totlen if totlen else 0) #in units of /100kb, assuming totlen is in units of kb
        self.coverage     = (self.totlen/1e3/self.genomesizemb if self.genomesizemb > 0 else 0) #totlen is in kb
        self.labdensity2  = 0 #second color channel

    #inlist is list of mol lens _in kb_
    #because this is length only, no nlab or labdensity
    def makeFromLenList(self, inlist) :
        self.nmol   = len(inlist)
        self.totlen = sum(inlist)
        self.n50    = getn50(inlist)
        self.coverage = (self.totlen/1e3/self.genomesizemb if self.genomesizemb > 0 else 0) #totlen is in kb

    #inlist is list of molecule objects
    # redundant with makeFromLenList except for totnlab, which you can't get from len alone
    def makeFromMolList(self, inlist) :
        lenlist = [x.lenkb for x in inlist]
        self.nmol       = len(lenlist)
        self.totlen     = sum(lenlist) #assume kb
        self.n50        = getn50(lenlist)
        self.totnlab    = sum( [len(x.labposeskb) for x in inlist] )
        self.labdensity = (self.totnlab*1e2/self.totlen if self.totlen else 0) #in units of /100kb
        self.coverage   = (self.totlen/1e3/self.genomesizemb if self.genomesizemb > 0 else 0) #totlen is in kb
        self.totnlab2   = (sum( [len(x.labposeskb2) for x in inlist] ) if inlist[0].labposeskb2 else 0)
        self.labdensity2= (self.totnlab2*1e2/self.totlen if self.totlen else 0) #in units of /100kb

    def __str__(self) :
        ustr = "kb" #default is kb, so unitscale is relative to this
        if self.unitscale == 1e3 : #ie, 1e3 kb is Mb
            ustr = "Mb"
        out =  "N mols: %i\n" % self.nmol
        uselen = self.totlen/(self.unitscale if self.unitscale==1e3 else 1)
        #usen50 = self.n50   /(self.unitscale if self.unitscale==1e3 else 1) #keep n50 in kb
        out += ("Total len ("+ustr+"): %10.3f\n") % uselen
        out += ("Avg len ("+"kb"+")  : %10.3f\n") % (self.totlen/self.nmol if self.nmol > 0 else 0) #always kb
        out += ("Mol N50 ("+"kb"+")  : %10.3f\n") % self.n50 #usen50
        if self.labdensity > 0 : 
            out += "Lab (/100kb)  : %10.3f \n" % self.labdensity
        if self.labdensity2 > 0 : 
            out += "Lab2(/100kb)  : %10.3f \n" % self.labdensity2
        if self.genomesizemb > 0 :
            out += "Genome Cov (x): %10.3f\n" % self.coverage
        #print "max, min len (kb):", max(lenlist), min(lenlist) #forget this
        return out

#end class bnxstats


#molecule class is for entries in a bnx file

class molecule:
    
    def __init__(self, linein='', doqual=False): #linein is all lines from bnx which make molecule
        self.doQuality = doqual #not implemented
        #self.snrQ = "QX11" #this depends on the bnx version, so this should really be an argument, but assume this (version 1.0) -- not used
        if linein :
            self.labposeskb2 = None #if one-color bnx
            self.makeFromLine(linein)
        else : #init data members
            self.id    = 0
            self.lenkb = 0.
            self.labposeskb = []
            #self.labsnr = []
            self.labposeskb2 = [] #second color channel if present


    #the argument is the four (or two, or more) lines in a bnx which contain the molecule data
    #for now, doQuality will do snr only -- not implemented
    def makeFromLine(self, linein) :
        for line in linein.split("\n") :
            if not line :
                continue
            #0 is just the id and the length
            if line[0] == "0" :
                items = line.split()
                #try here for cast?
                self.id    = int(items[1])
                self.lenkb = float(items[2])/1e3
                if len(items) > 3 : #bnx 1.0 and later
                    self.avgIntensity = float(items[3])
            #1 is the label positions
            elif line[0] == "1" :
                poses = line.split()[1:-1] #first is just '1', last is length: discard both
                self.labposeskb = [float(x)/1e3 for x in poses]
            #elif self.doQuality and line.startswith(self.snrQ) :
            #    self.labsnr = 
            #2 is labels in second color channel
            elif line[0] == "2" :
                poses = line.split()[1:-1] #first is just '1', last is length: discard both
                self.labposeskb2 = [float(x)/1e3 for x in poses]


    #lab per 100kb is number of labels over length in 100kb
    def getLabPer100kb(self) :
        return len(self.labposeskb)*100/self.lenkb

    #scale length and all label positions by argument
    # ie, for length adjustment, give bpp_in/500 as argument
    def stretchAll(self, stretch) :
        self.lenkb *= stretch
        self.labposeskb = [x*stretch for x in self.labposeskb]

    #take another molecule and calculate the 'maximum absolute deviation'
    #this is the maximum of the difference in label positions
    #only for identical label indices, so this works best for when the
    # two molecules have the same number of labels
    def getMaxAbsDev(self, inmol) :
        if len(self.labposeskb) < len(inmol.labposeskb) :
            minm = self
            maxm = inmol
        else :
            minm = inmol
            maxm = self
        return max([math.fabs(pos - maxm.labposeskb[i]) for i,pos in enumerate(minm.labposeskb)])            

    #return bnx format 0.1 string: just id, length, and positions in bp, with trailing newline
    def printBnx(self) :
        retstr =  "0\t%i\t%.1f\n" % (self.id, self.lenkb*1e3)
        retstr += "1\t"+"\t".join(["%.1f"%(x*1e3) for x in self.labposeskb])+"\n"
        return retstr

    #str method is for quick printing, not bnx format
    def __str__(self) :
        strposes = ["%.3f"%x for x in self.labposeskb]
        return "%7i  %.3f" % (self.id, self.lenkb) + "\n" + " ".join(strposes)

#end class molecule


def simpleBnxStats(bnxpath,minlen=0) :
    """The above classes are very inefficient, especially with memory. 
    Replace them with a simple loop on lines in bnx file to get the stats in the bnxstats class.
    Include only molecules whose length is >= minlen (kb).
    Return total length in Mb, other lengths in kb, density in /100kb
    """
    if not checkFile(bnxpath, ".bnx") :
        return
    bnxv = getBnxVersion(bnxpath) #bnx version 0.1 vs >= 1.0; return 0 on error
    #print "version =", bnxv #debug
    if bnxv == 0 : #warn if unable to get bnx version
        print "Warning in simpleBnxStats: could not read bnx version for file: %s" % bnxpath
    if bnxv >= 1.0 :
        nlab, mollens = simpleBnxStats10(bnxpath,minlen)
    else :
        nlab, mollens = simpleBnxStats01(bnxpath,minlen)
    sumlen = sum(mollens)
    nmol = len(mollens)
    #print nlab #debug
    labden = (nlab*1e5/sumlen if sumlen > 0 else 0) # /100kb
    return({"nmol" : nmol,
            "totlen" : sumlen/1e6, #Mb
            "avglen" : (sumlen/1e3/nmol if nmol > 0 else 0), #kb
            "n50" : getn50(mollens)/1e3, #kb
            "labdensity" : labden,
            "medianlen": getMedian(mollens)})
#end simpleBnxStats


def simpleBnxStats10(bnxpath,minlen=0) :
    """Utility for simpleBnxStats: this is the loop over the bnx file for bnx version 1.0 and newer."""
    nlab = 0
    mollens = []
    f = open(bnxpath)
    for line in f :
        #the backbone is 0, and this is all we need for this fn
        if line[0] != "0" :
            continue
        tokens = line.split()
        #the NumberofLabels field is the sixth
        if len(tokens) < 6 :
            continue
        try :
            ilab = int(tokens[5])
            ilen = float(tokens[2])
        except :
            continue
        if ilen/1e3 < minlen :
            continue
        nlab += ilab
        mollens.append(ilen)
    #end for line in f
    f.close()
    return( nlab, mollens )
#end simpleBnxStats10

def simpleBnxStats01(bnxpath,minlen=0) :
    """Utility for simpleBnxStats: this is the loop over the bnx file for bnx version 0.1."""
    nlab = 0
    mollens = []
    f = open(bnxpath)
    islab = False
    for line in f :
        #the backbone is 0 and labels are 1: we need both here
        if line[0] == "1" :
            islab = True
        elif line[0] != "0" :
            continue
        tokens = line.split()
        if islab : #1, or labels
            ilab = len(tokens)-2 #first entry in tokens is "1" for the line id, and last is length of mol
            nlab += ilab
        else : #0, or backbone
            try :
                ilen = float(tokens[2])
            except :
                continue
            if ilen/1e3 < minlen :
                continue
            mollens.append(ilen)
        #end if label or backbone
        islab = False
    #end for line in f
    f.close()
    return( nlab, mollens )
#end simpleBnxStats10

def getBnxVersion(bnxpath):
    """Utility for simpleBnxStats: get bnx file version then return.
    Return 0 if version not found."""
    ver = 0.
    f = open(bnxpath)
    for line in f :
        if line.startswith("# BNX File Version") :
            try:
                ver = float(line.split(":")[1])
            except :
                ver = 0.
            break
        elif line[0] != "#" : #header lines should be first, if done with header, failed to find version
            break
    f.close()
    return(ver)
#end getBnxVersion


def findRepeat( inlist, tolerance=0.1, minele=2 ) :
    '''Simple repeat detection: if neighboring intervals are similar, they are repeat.
    Similar means that abs(interval1-interval2)/min(interval1-interval2) <= tolerance.
    For reference, recommend tolerance = 0.01 because some are off by only a few bp.
    For single molecules, 0.1 or so is better due to noise.
    Shortest possible number of repeat intervals is 2 (minele), use more to find longer repeats.
    inlist must be sorted.
    Return list of tuples: first ele is start position, second is end position, third is number
    of repeat elements, fourth is average size of elements.'''

    #interval i starts at inlist[i] and ends at inlist[i+1] (obviously one fewer element)
    intervals = [inlist[i+1]-inlist[i] for i in range(len(inlist)-1)]
    #repeats i are in interval[i] and interval[i+1]
    #which means inlist[i] to inlist[i+2] (no tolerance applied yet)
    repeats = [abs(intervals[i+1]-intervals[i])/min(intervals[i+1], intervals[i]) for i in range(len(intervals)-1)]
    #apply tolerance as filter: eles of repeats are -1 if > tolerance, or the repeat quantity if <=
    isrepeat = lambda x: x if x <= tolerance else -1
    if not any(filter(lambda x: isrepeat(x) >= 0, repeats)) : #nothing left to do
        return
    repeats = map(isrepeat, repeats)
    
    #in order to get number of eles, need number of sequential repeats--if it's unbroken, it's the same element
    #there is a choice here: keep this simple by returning a list describing repeats
    #or make a class to describe repeats--keep it simple for now, and add the more complicated part later
    minele = max(1,minele-1) #take one away from argument (minimum 1) because it's applied to number of eles in repeat, which is number of intervals less one
    reppos = []
    isrepeat = False
    startpos = startidx = nele = 0
    #maxrep = 1500 #DEBUG -- stop after this n ele, 0 to disable
    for ir,rep in enumerate(repeats) :
        #if ir > 1400 : #DEBUG
        #    print inlist[ir], ":", intervals[ir], intervals[ir+1], ":", rep
        #rep is -1 for non-repeat intervals; >= is very important, since exact is 0
        if rep >= 0 : 
            nele += 1
            if not isrepeat : #start of repeat region
                startidx = ir
                startpos = inlist[ir]
            isrepeat = True
        elif isrepeat : #rep < 0 : store current repeat, reset isrepeat, startpos, nele
            assert startpos #must be > 0
            if nele >= minele : 
                useint = intervals[startidx:ir+1]
                avgsize = sum(useint)/len(useint) #average interval size
                #nele is n elements in repeat, but n intervals is actually one more
                #and since this is the interval after the repeat ended, we already added 1--only add one more for inlist
                reppos.append( (startpos, inlist[ir+1], nele+1, avgsize) ) 
            isrepeat = False
            startpos = startidx = nele = 0
        #if maxrep and ir > maxrep : #DEBUG
        #    break
    #if the end of the list is a repeat element, the elif isrepeat isn't hit--got to catch this
    if isrepeat and nele >= minele :
        assert startpos #dummy check
        useint = intervals[startidx:ir+2] #here need + 2 bc end of repeats list
        avgsize = sum(useint)/len(useint) #average interval size
        reppos.append( (startpos, inlist[ir+2], nele+1, avgsize) ) 
    return reppos
#def findRepeat



# See http://code.activestate.com/recipes/410692/
#
# C or Tcl-like switch statement to compensate for python shortcomings.
#
# This class provides the functionality we want. You only need to look at
# this if you want to know how this works. It only needs to be defined
# once, no need to muck around with its internals.
class switch(object):
    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration
    
    def match(self, *args, **kwargs):
	regexp=kwargs.pop("regexp", False)
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True            
        elif not regexp and self.value in args: # changed for v1.5, see below
            self.fall = True
            return True
        elif regexp : # changed for v1.5, see below
            self.fall = re.match(args[0], self.value)!=None
            return self.fall
        else:
            return False

# Replace (if present) sequence of list values
# Currently only supports key-value pair
def argumentReplaceList(origList, updatedVals):
    if updatedVals.__len__() != 2:
        return origList
    targetVal = updatedVals[0]
    try:
        i = origList.index(targetVal)
    except:
        return origList + updatedVals

    newList = list(origList) #copy of origList
    for updatedVal in updatedVals:
        if i >= len(newList):
            newList.append(updatedVal)
        else:
            newList[i] = updatedVal
        i += 1
    return newList



def getPlatform():
    '''Check for platform.system() == "Linux" or "Windows".
    If neither (or error), return None. If Linux, return 1, if Windows, return 2.
    '''
    try :
        import platform
        system = platform.system()
        if system == 'Linux' :
            return 1
        elif system == 'Windows' :
            return 2
        else :
            return None #"utilities.getPlatform: unknown platform:" + system
    except :
        return None


def checkExternalDependencies(dep, islinux=True):
    '''Check external dependencies.
    Currently supported: R, perl, perf.
    islinux is bool for linux or windows. No other platform is supported.
    Return tuple whose first element is path to executable if dependency is found,
    and None if argument is not supported or not found. The second element is
    an error message if an error occurs, or None if no errors.
    '''
    if not dep in ["R", "Rscript", "perl", "perf", "time"]:
        return None, ("ERROR in utilities.checkExternalDependencies: dependency %s is not supported" % (dep))

    if islinux :
        whichcommand = 'which'
    else :
        whichcommand = 'where' #this apparently works with server 2003 and newer, including Vista, 7

    which = runJob([whichcommand, dep], returnstdout=True)
    which = [which[0], which[1].rstrip()] #must remove whitespace before calling checkExecutable
    if which[0] : #this is return code--non-zero is bad
        return None, ("utilities.checkExternalDependencies: %s failed for dependency %s" % (whichcommand, dep))
    if checkExecutable(which[1]) : #as long as the result is executable, it should be good
        return which[1], None
    else :
        return None, ("utilities.checkExternalDependencies: %s result (%s) for dependency %s is not executable" % (whichcommand, which[1], dep))
    

def runJob(args, returnstdout=False, printstdout=False, printfail=False) :
#def runJob(args, returnstdout=False, printstdout=False, printfail=False, shell=False) :
    '''Call external process defined by args (list). Return job's return code, unless returnstdout, in which case return tuple of return code, stdout, and stderr. If printstdout, also print stdout. If printfail, print on return code != 0. Will wait for job to finish before returning.
    '''
    import subprocess
    #I'm ignoring stderr, hopefully it won't be too long--RefAligner prints to stderr when using -stderr
    # stderr is now captured also, but it's possible that it will exceed the buffer size allowed by the shell (or python) and hang
    proc = subprocess.Popen( args, stdout=subprocess.PIPE, stderr=subprocess.PIPE ) 
    #proc = subprocess.Popen( args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=shell ) 
    allout = ""
    for line in iter(proc.stdout.readline, b''): #iter will keep looping until nothing else to loop on ('')
        allout += line
    errout = ""
    for line in iter(proc.stderr.readline, b''): #apparently, the args to Popen do not move stderr from here
        errout += line
    proc.wait() #wait is necessary otherwise returncode is None

    if proc.returncode != 0 and printfail :
        print "Job return non-zero:", proc.returncode, "\n", allout, "\n", errout, "\n"
    if printstdout :
        print allout
    if returnstdout :
        return (proc.returncode, allout, errout)
    else :
        return proc.returncode
#end runJob



def getRefAlignerVersion(rabin) :
    """Call rabin (must be path to RefAligner) with -version argument, parse result,
    return version string. Return None (silently) on failure.
    """
    if not checkExecutable(rabin) :
        return None
    try :
        retc, version, err = runJob([rabin, "-version"], returnstdout=True) #bc returnstdout is True, return tuple
    except OSError : #if RefAligner.mic is used, get "OSError: [Errno 8] Exec format error"
        return None
    if retc : #0 is success
        return None
    version = version.split("\n")
    if len(version) < 3 : #need third line
        return None
    #version = version[2].split() #do not assume line number: loop to find appropriate line
    useversion = None
    for line in version :
        if line.rfind("SVNversion") != -1 or line.rfind("SVNRevision") != -1 :
            useversion = line.split()
            break
    if useversion == None :
        #print version #debug
        return None

    ret = ""
    for wi,word in enumerate(useversion) :
        #old syntax: "SVNversion=X"
        if word.startswith("SVNversion=") :
            ret = word.replace("SVNversion=","")
            break
        #new syntax: "SVNRevision: X"
        if word.startswith("SVNRevision") :
            ret = version[wi+1]
            break

    if ret == "" :
        #print version #debug
        return None

    if ret.find(":") != -1 : #remove after :
        ret = ret[:ret.find(":")]
    if len(ret) and getIntFromString(ret[-1]) == None : #last char is string; remove
        ret = ret[:-1]
    ret = (getIntFromString(ret) if getIntFromString(ret) != None else 0)

    return ret
#end getRefAlignerVersion


#these two imports are for the Memory fns below
#resource should be builtin on linux, but on windows it may be missing
#if do in local scope, need to make an argument
global have_resource
try :
    import resource
    have_resource = True
except :
    have_resource = False
global have_psutil #this module is not installed by default
try :
    import psutil
    have_psutil = True
except :
    have_psutil = False


def initMemoryLog(logpath, zerotime=None) :
    """Init log file at logpath, overwriting existing file if it exists.
    zerotime is argument to logMemory. Assume default unit args to logMemory
    and getMemoryUsage."""
    f1 = open(logpath, 'w') #will overwrite if exists; get a new file
    global have_resource
    global have_psutil
    if have_resource and have_psutil : #6 columns
        f1.write("stage_string\ttime(s)\tmemory_used(MB)\tavailable_memory(GB)\tpercent_memory_used\n") #see logMemory below
    elif have_resource : #only 3 columns
        f1.write("stage_string\ttime(s)\tmemory_used(MB)\n")
    else :
        f1.write("python package resource is missing; not logging memory\n")
    f1.close()
    logMemory(logpath, zerotime, "Pipeline_init")


def logMemory(logpath, zerotime=None, stage="") :
    """Call fn getMemoryUsage and write its returned data, plus time info, to logpath;
    also include arbitrary string stage if non-empty. If no resource package, do nothing."""
    #format, tab separated:
    #stage_string time(minus zerotime) memory_used total_memory available_memory precent_memory
    global have_resource
    if not have_resource :
        return
    f1 = open(logpath, 'a')
    m = getMemoryUsage() #see below
    t = time.time() - (zerotime if zerotime != None else 0)
    global have_psutil
    if have_psutil : #here getMemoryUsage returns a tuple
        content = "%32s\t%.1f\t%.2f\t%.2f\t%.3f\n" % (stage, t, m[0], m[1], m[2])
    else : #here getMemoryUsage returns a float
        content = "%32s\t%.1f\t%.2f\n" % (stage, t, m)
    f1.write(content) 
    f1.close()


def getMemoryUsage(units=1e3) :
    """Use python module resource to get resident memory, and psutil to get system memory stats.
    If no resource package, do nothing.
    Use units as denominator to resource.getrusage, eg, units=1 is kB, units=1e3 is MB. Units
    for the psutil data will be 1e-3 times this, eg, if getrusage is MB, this is GB."""
    #to quote the docs for resource: "may also raise error exception in unusual circumstances." (whatever that means)
    #to not quote the docs, because it's not there, the units appear to be kB, based on /proc//status comparison
    global have_resource
    if not have_resource :
        return
    try :
        maxrss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / units
    except :
        maxrss = 0
    global have_psutil #if don't have it, nothing left to do
    if not have_psutil : 
        return maxrss 
    try :
        m = psutil.virtual_memory()
        #note: m.percent is 100*(total - available)/total, ie, percent used on the system, not percent used by process
        #1e3 bc this data is bytes, and another 1e3 bc of my convention of changing units as described above
        return maxrss, m.available/(units*1e6), m.percent
    except :
        return maxrss, 0, 0, 0
