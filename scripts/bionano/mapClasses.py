#import pdb
import os
import math
import copy #for deepcopy

import utilities as util


"""@package mapClasses: Analysis for Xmap, Cmap Results

"""


util.setVersion("$Id: mapClasses.py 4889 2016-05-05 22:25:00Z wandrews $")


class xmap:

    def __init__(self, sourceFile=""):
        if util.checkFile(sourceFile, ".xmap") :
            self.makeFromFile(sourceFile)
        else :
            if sourceFile : #not empty, but doesn't exist--print warning
                print "Warning in xmap init: sourceFile is not a file:", sourceFile
            self.xmapLookup = {}
            self.header = ""

    def makeFromFile(self, sourceFile):
        self.xmapLookup = {}
        self.header = ""
        #sourceFile was already checked in constructor
        f1 = open(sourceFile)
        for line in f1 :
            if line[0] == '#':
                self.header += line
                continue
            tokens = line.split()
            ID          = int(tokens[0])
            contigRef   = int(tokens[2])
            contigQry   = int(tokens[1])
            RefStart    = float(tokens[5])
            RefStop     = float(tokens[6])
            QryStart    = float(tokens[3])
            QryStop     = float(tokens[4])
            Orientation = tokens[7]
            Confidence  = float(tokens[8])
            CigarString = tokens[9]
            qrylen      = (float(tokens[10]) if len(tokens) > 10 else None)
            reflen      = (float(tokens[11]) if len(tokens) > 11 else None)
            labch       = (util.getIntFromString(tokens[12]) if len(tokens) > 12 else None)
            align       = (      tokens[13]  if len(tokens) > 13 else None)
            dxmap = xmapEntry(ID, contigRef, contigQry, RefStart, RefStop, QryStart, QryStop, Orientation, Confidence, CigarString, qrylen, reflen, labch, align)
            self.xmapLookup[ID] = dxmap
        f1.close()


    #print all the xmapLookup entries, newlines in between
    def __str__(self) :
        retstr = self.header #remember, strings are immutable
        for idx in range(1,len(self.xmapLookup)+1) :
            retstr += str(idx) + "\t" + str(self.xmapLookup[idx]) + "\n"
        return retstr

    #like above, but add each xmapEntry's ConfByNMatch
    def strWithConfByNMatch(self, header=False) :
        retstr = ""
        if header :
            retstr += self.header
        for idx in range(1,len(self.xmapLookup)+1) :
            xentry = self.xmapLookup[idx]
            retstr += str(idx) + "\t" + str(xentry) + ("\t%.3f\n" % xentry.getConfByNMatch())
        return retstr

    #return list of aligned lengths on query
    #if qcontigid == None, do all queries; if not, do only that query id
    #unitscale is to convert to, eg, kb, if 1e3
    def getAllQryLens(self, qcontigid=None, unitscale=1e3) :
        return [x.getMappedQryLen()/unitscale for x in self.xmapLookup.values() if (not qcontigid or x.contigQry == qcontigid)]


    #like above, but ref
    def getAllRefLens(self, qcontigid=None, unitscale=1e3) :
        return [x.getMappedRefLen()/unitscale for x in self.xmapLookup.values() if (not qcontigid or x.contigQry == qcontigid)]


    #for sum, call above and use sum
    def getSumMappedQryLen(self, qcontigid=None, unitscale=1e3) :
        return sum(self.getAllQryLens(qcontigid, unitscale))


    #same as above but reference length
    def getSumMappedRefLen(self, qcontigid=None, unitscale=1e3) :
        return sum(self.getAllRefLens(qcontigid, unitscale))


    #invert an xmap: call invert method of each entry
    def invert(self) :
        for idx,xentry in self.xmapLookup.iteritems() :
            xentry.invert()


    #check if a query contig shows up in at least one xmapEntry
    def hasQueryContig(self, qcontigid) :
        for xentry in self.xmapLookup.itervalues() :
            if xentry.contigQry == qcontigid :
                return True
        return False


    #like above, but count the number of entries
    def countQueryContig(self, qcontigid) :
        count = 0
        for xentry in self.xmapLookup.itervalues() :
            if xentry.contigQry == qcontigid :
                count += 1
        return count


    #get a list of the confidences; if qcontigid supplied, just for those maps
    def getConfList(self, qcontigid=None) :
        confl = []
        for xentry in self.xmapLookup.values() :
            if qcontigid == None or xentry.contigQry == qcontigid :
                confl.append( xentry.Confidence )
        return confl


    #like above, but for N matches in cigar string
    def getNMatchList(self, qcontigid=None) :
        matchl = []
        for xentry in self.xmapLookup.values() :
            if qcontigid == None or xentry.contigQry == qcontigid :
                matchl.append( xentry.getNCigarStringMatches() )
        return matchl


    #side-length refers to the length on the sides of an SV
    # align-side-length means the aligned (as opposed to total, or 'true')
    #If one or more than two xmapEntries, where the SV is (or if there is one) is not well defined
    # return zero in this case
    #If two xmapEntries, return min of ref, qry align len in each
    # the reasoning for this is that the most conservative is the minimum
    #verbose behavior:
    # 0 : nothing
    # 1 : nothing (compatibility with various svdetect methods)
    # 2 : xmapLookup size 1 warning; left and right align lengths
    def getMinAlignSideLength(self, verbose=0) :

        if len(self.xmapLookup) <= 1 :
            if verbose > 1 : #should be high verbose
                print "Warning in xmap.getMinAlignSideLength: 0 or 1 xmapLookup entries" 
            return 0 #not None so that this can be added

        #for the purposes of the in-silico studies (svdetect), I want to know if there are more than two
        # if not for this, just combine this if with above
        #I don't store the sourcefile arg to makeFromFile, so, just do generically
        if len(self.xmapLookup) > 2 :
            #could be verbose > x, but for now, want to see it
            print "Warning in xmap.getMinAlignSideLength: more than two xmap entries" 
            return 0 #can do -1 if want to check, but that's a bit dangerous if, eg, you try to add it

        #since len == 2 after above, indicies ought to be 1 and 2 bc these are the xmap entry ids
        # note: having two entries does not guarantee the key values, though they should always be this
        if not self.xmapLookup.has_key(1) or not self.xmapLookup.has_key(2) :
            #this should never happen, so always print
            print "Error in xmap.getMinAlignSideLength: invalid xmapLookup keys"
            return 0 #can do -1 if want to check, but that's a bit dangerous if, eg, you try to add it

        left  = min( self.xmapLookup[1].getMappedQryLen(), self.xmapLookup[1].getMappedRefLen() )
        right = min( self.xmapLookup[2].getMappedQryLen(), self.xmapLookup[2].getMappedRefLen() )
        if verbose > 1 : #same level as svdetect.getTrueSideLength
            print "xmap.getMinAlignSideLength: left =", left, "right =", right
        return min( left, right )
    #end def getMinAlignSideLength

    #find chimeric contigs
    #this is designed for one-contig-per-chromosome to simplify
    # don't consider multiple alignments to the same chromosome as chimeric
    #qryid is required--only compare same qryid
    #Multiple alignments must overlap by _less_ than 75% of the total to be chimeric (maxchimoverlap)
    #consider each pair of alignments
    #return number of chimeric alignments and total chimeric length
    def findChimericContig(self, qryid, unitscale=1e3, maxchimoverlap = 0.75, verbose=0) :
        assert unitscale > 0, "Invalid argument unitscale:"+str(unitscale)
        if len(self.xmapLookup) <= 1 : #one alignment can't be chimeric
            if verbose : print "single alignment"
            return []
        keys = self.xmapLookup.keys() #use keys to loop to be sure to not double-count pairs
        nchim = 0
        chimlen = 0
        for i1 in range(len(keys)) :
            #key1 = keys[i1]
            xe1 = self.xmapLookup[keys[i1]]
            if xe1.contigQry != qryid : #only consider this qry contig
                continue
            for i2 in range(i1+1,len(keys)) :
                #key2 = keys[i2]
                xe2 = self.xmapLookup[keys[i2]]
                if xe1.contigRef == xe2.contigRef or xe2.contigQry != qryid : #same ref contig
                    continue
                #find overlap
                p1 = sorted([xe1.QryStart, xe1.QryStop])
                p2 = sorted([xe2.QryStart, xe2.QryStop])
                #no overlap happens when end of first is before start of second--this is chimeric
                if ( (p1[0] <= p2[0] and p1[1] <= p2[0]) or
                     (p1[0] >= p2[0] and p1[0] >= p2[1]) ) :
                    nchim += 1
                    chimi = min(xe1.getMappedQryLen(), xe2.getMappedQryLen())
                    chimlen += chimi
                    if verbose : print i1+1, i2+1, "no overlap: nchim:", nchim, "chimlen:", chimlen, chimi 
                #the overlap non-chimeric case is when one is a sub-alignment of the other
                elif ( (p1[0] <= p2[0] and p1[1] >= p2[1]) or
                       (p2[0] <= p1[0] and p2[1] >= p1[1]) ) :
                    if verbose : print i1+1, i2+1, "Non-chimeric sub-alignment" #debug
                #last case is partial overlap -- chimerism depends on maxchimoverlap
                elif ( (p1[0] <= p2[0] and p1[1] <= p2[1]) or #p1 < p2
                       (p1[0] >= p2[0] and p1[1] >= p2[1]) ) : #p2 < p1
                    overlap  = (p1[1] - p2[0] if p1[0] <= p2[0] else p2[1] - p1[0])
                    totalign = (p2[1] - p1[0] if p1[0] <= p2[0] else p1[1] - p2[0])
                    if overlap/totalign < maxchimoverlap :
                        nchim += 1
                        chimlen += overlap
                        if verbose : print i1+1, i2+1, "chimeric overlap: nchim:", nchim, "chimlen:", chimlen, overlap 
                    else :
                        if verbose : print i1+1, i2+1, "non-chimeric overlap:", overlap, totalign, overlap/totalign 
                else : #I think this should never happen
                    errstr = str(i1+1) + " " + str(i2+1) + " overlap case not covered: " + str(p1) + " " + str(p2)
                    print "ERROR:", errstr
                    assert False, errstr
        return [nchim, chimlen/unitscale]
    #end findChimericContig


    #header for statistics table
    def statisticsTableHeader(self, molxmap=None, refxmap=None) :
        headl = []
        if molxmap :
            headl += [
                "MCount",    #Molecule number of alignments
                "MqLen",     #Molecule sum all molecules aligned length (kb)
                "MqLen2",    #Molecule sum all molecules aligned length squared (kb^2)
                "McLen",     #Molecule sum all aligned contig length (kb)
                "McLen2",    #Molecule sum all aligned contig length squared (kb^2)
                "MavgLen",   #Molecule average all molecules aligned length (kb)
                "MavgCLen",  #Molecule average all aligned contig length (kb)
                "MmaxConf",  #Molecule maximum alignment confidence
                "MmaxCbN",   #Molecule (confidence / N matches) for maximum confidence alignment
                "MavgConf",  #Molecule average alignment confidence
                "MsumConf",  #Molecule sum alignment confidence
                "MsumConf2", #Molecule sum alignment confidence squared
                "MavgCbN",   #Molecule average (confidence / N matches)
                "MsumCbN",   #Molecule sum (confidence / N matches)
                "MmaxMtch",  #Molecule maximum N matches (ie, aligned sites)
                "MavgMtch"   #Molecule average N matches
                ]
            #these are self-expanatory--just straight out of .err file
            headl += ["M_FP", "M_FN", "M_sf", "M_sd", "M_bpp"] 
        if refxmap :
            headl += [
                "RCount",   #Reference number of alignments
                "RsumLen",  #Reference sum of aligned length (on reference) (kb)
                "RsumQLen", #Reference sum of aligned length (on query) (kb)
                "RmaxQLen", #Reference maximum of aligned length (on query) (kb)
                "RmxQLIdx", #Reference index of maximum query length alignment
                "RsumConf", #Reference sum of confidences of all alignments
                "RmaxConf", #Reference maximum confidence of all alignments 
                "RmxCfIdx", #Reference index of maximum confidence alignment
                "RsumMtch", #Reference sum of number of matches of all alignments
                "RmxMatch", #Reference maximum number of matches of all alignments  
                "RmxMtIdx", #Reference index of maximum number of matches alignment
                "RNChim",   #Reference number of chimeric pairs of alignments
                "RChimLen"  #Reference sum of lengths of chimeric portions of all chimeric pairs (kb)
                ]
            headl += ["R_FP", "R_FN", "R_sf", "R_sd", "R_bpp"]
        return "  ".join(headl) + "\n"


    #analogous to the cmap method of the similar name
    #qryid is in case the xmap has multiple queries
    # if not supplied, assert they're all the same
    #for the header, put it in the cmap because that's easier, but put as comment here:
    def statisticsTableRefRow(self, qryid=None, unitscale=1e3) :
        if not qryid :
            qryid = self.xmapLookup[1].contigQry
            assert all(qryid == x.contigQry for x in self.xmapLookup.values()), "if qryid not supplied, must be single qry contig xmap"
        aligncount = self.countQueryContig(qryid)
        rlenl = self.getAllRefLens(qryid, unitscale) #in kb
        qlenl = self.getAllQryLens(qryid, unitscale) #in kb
        #avgqlen = (sum(qlenl)/len(qlenl) if qlenl else 0) #this is redundant bc I have N and sum
        #list comprehension should return single ele list; if multi entries with max, just take first
        maxqlenidx = [i for i,xe in self.xmapLookup.iteritems() if xe.getMappedQryLen()/unitscale == max(qlenl)][0] if qlenl else 0
        confl = self.getConfList(qryid)
        #confr = (sum(confl)/len(confl) if confl else 0) #redundant
        maxconfidx = [i for i,xe in self.xmapLookup.iteritems() if xe.Confidence == max(confl)][0] if confl else 0 
        matchl = self.getNMatchList(qryid)
        #print matchl
        maxmatchidx = [i for i,xe in self.xmapLookup.iteritems() if xe.getNCigarStringMatches() == max(matchl)][0] if matchl else 0
        chiml = self.findChimericContig(qryid, unitscale) #, maxchimoverlap = 0.75) :
        datalist = []
        datalist.append( "%3i"   % aligncount ) #ACnt
        datalist.append( "%7.1f" % sum(rlenl) ) #SRefLen (in kb)
        datalist.append( "%7.1f" % sum(qlenl) ) #SQryLen (in kb)
        datalist.append( "%7.1f" % (max(qlenl) if qlenl else 0) ) #MQryLen (in kb)
        datalist.append( "%3i"   % maxqlenidx ) #MxQLIdx
        #datalist.append( "%7.1f" % avgqlen ) #AQryLen (in kb)
        datalist.append( "%7.2f" % sum(confl) ) #SumConf
        datalist.append( "%6.2f" % (max(confl) if confl else 0) ) #MaxConf
        datalist.append( "%3i"   % maxconfidx ) #MxCfIdx
        #datalist.append( "%6.2f" % confr ) #AvConf
        datalist.append( "%5i"   % sum(matchl) ) #SumMtch
        datalist.append( "%4i"   % (max(matchl) if matchl else 0) ) #MxMtch
        datalist.append( "%4i"   % maxmatchidx ) #MxMtIdx
        #datalist.append( "%6.1f" % (sum(matchl)/len(matchl) if matchl else 0) ) #AvMtch
        datalist.append( "%4i"   % (chiml[0] if chiml else 0) ) #NChim
        datalist.append( "%7.1f" % (chiml[1] if chiml else 0) ) #ChimLen (in kb)
        return "  ".join(datalist)


    #like above, but for xmaps comparing molecules to contigs
    #each contig has its own xmap, so no need for refid (qryids are the mol ids)
    def statisticsTableMolRow(self) :
        #just be sure all ref contig ids are the same, otherwise logic below won't work
        assert all(self.xmapLookup[1].contigRef == x.contigRef for x in self.xmapLookup.values()), "Not all ref contig ids are equal"
        qlenl = self.getAllQryLens() #default is kb
        rlenl = self.getAllRefLens() #default is kb
        confl = self.getConfList()
        matchl = self.getNMatchList()
        cbnl = [x/matchl[i] for i,x in enumerate(confl)] #this will only crash if confl is not empty but matchl is
        maxconf = max(confl) if confl and matchl else 0 #confl and matchl bc use matchl later (if maxconf)
        cbn_maxconf = 0
        if maxconf :
            #if there are multiple occurrences of this confidence, pick the lowest n matches to maximize conf/Nmatch
            if confl.count(maxconf) > 1 :
                maxindl = [i for i,x in enumerate(confl) if x == maxconf] #indices of these
                cbn_maxconf = maxconf/min([matchl[i] for i in maxindl])
            else : #only one max conf
                cbn_maxconf = maxconf/matchl[confl.index(maxconf)]

        datalist = []
        datalist.append( "%3i"   % len(self.xmapLookup) ) #Nmol
        datalist.append( "%9.1f" % sum(qlenl) ) #molLen
        datalist.append( "%9.1f" % sum(x**2 for x in qlenl) ) #molLen2
        datalist.append( "%9.1f" % sum(rlenl) ) #mlCLen
        datalist.append( "%9.1f" % sum(x**2 for x in rlenl) ) #mlCLen2
        datalist.append( "%6.1f" % (sum(qlenl)/len(qlenl) if qlenl else 0)) #avgMLen
        datalist.append( "%6.1f" % (sum(rlenl)/len(rlenl) if rlenl else 0)) #avgCLen
        datalist.append( "%6.1f" % maxconf) #mlMConf
        datalist.append( "%6.3f" % cbn_maxconf) #mlMCbN
        datalist.append( "%6.1f" % (sum(confl)/len(confl) if confl else 0) ) #avgConf
        datalist.append( "%6.1f" % sum(confl) ) #sumConf
        datalist.append( "%9.1f" % sum(x**2 for x in confl) ) #sumCnf2
        datalist.append( "%6.3f" % (sum(cbnl)/len(cbnl) if cbnl else 0) ) #avgCbN
        datalist.append( "%6.3f" % sum(cbnl) ) #sumCbN
        datalist.append( "%4i"   % (max(matchl) if matchl else 0) ) #mlMMtch
        datalist.append( "%6.1f" % (sum(float(x) for x in matchl)/len(matchl) if matchl else 0) ) #mlAMtch
        return "  ".join(datalist)
    #end statisticsTableMolRow


    def editHeaderMaps(self, newpath, query=True, smap=False) :
        '''Replace the path in the line 'Query Maps', for query=True,
        or 'Reference Maps', for query=False, with argument.'''
        newheader = ""
        if query and smap :
            print "Error in xmap.editHeaderMaps: cannot edit both query and smap"
            return
        for line in self.header.split('\n') :
            if not line :
                continue
            if ( (line.find("Query Maps") != -1 and query and not smap) or
                 (line.find("Reference Maps") != -1 and not query and not smap) or
                 (line.find("Smap Entries") != -1 and smap) ) :
                newheader += line.split(":")[0] + ":\t" + newpath + "\n"
            else :
                newheader += line+"\n"
        self.header = newheader


    def editHeaderQueryMaps(self, newpath) :
        self.editHeaderMaps(newpath, query=True)


    def writeToFile(self, filename):
        if not filename.endswith(".xmap"):
            print "Error in xmap.writeToFile: file to write must end in .xmap"
            return
        outfile = open(filename, "w")
        outfile.write(str(self)) #includes header and all entries with newlines
        outfile.close()

#end class xmap


            
class xmapEntry:

    def __init__(self, ID, contigRef, contigQry, RefStart, RefStop, QryStart, QryStop, Orientation, Confidence, CigarString, RefLen=None, QryLen=None, LabCh=None, Align=None):
        self.contigRef = contigRef
        self.contigQry = contigQry
        self.RefStart = RefStart
        self.RefStop = RefStop
        self.QryStart = QryStart
        self.QryStop = QryStop
        if Orientation == '+':
            self.OrientForward = True
        else:
            self.OrientForward = False
        #self.Orientation = Orientation
        self.Confidence  = Confidence
        self.CigarString = CigarString
        self.RefLen = RefLen
        self.QryLen = QryLen
        self.LabCh = LabCh
        self.Align = Align

    #this is what you get when you do print on this object
    #the xmapEntry doesn't know its ID, so ignore this field
    #no newline
    def __str__(self) :
        orientstr = ("+" if self.OrientForward else "-")
        sep = "\t"
        ret = ("%i"+sep+"%i"+sep+"%.1f"+sep+"%.1f"+sep+"%.1f"+sep+"%.1f"+sep+"%s"+sep+"%.2f"+sep+"%s") % (self.contigQry, self.contigRef, self.QryStart, self.QryStop, self.RefStart, self.RefStop, orientstr, self.Confidence, self.CigarString)
        if self.RefLen != None and self.QryLen != None and self.LabCh != None and self.Align != None :
            ret += (sep+"%.1f"+sep+"%.1f"+sep+"%i"+sep+"%s") % (self.RefLen, self.QryLen, self.LabCh, self.Align)
        return ret

    def getMappedQryLen(self) :
        return math.fabs( self.QryStart - self.QryStop )

    def getMappedRefLen(self) :
        return math.fabs( self.RefStart - self.RefStop )

    def getMinMappedLen(self) :
        return min( self.getMappedQryLen(), self.getMappedRefLen() )


    #return confidence / min_len, where min_len means min of QryLen and RefLen
    #divide by unitconv if it's > 0
    def getConfByMinLen(self, unitconv=0) :
        #qrylen = ( self.getMappedQryLen() if unitconv <= 0 else self.getMappedQryLen()/unitconv )
        #reflen = ( self.getMappedRefLen() if unitconv <= 0 else self.getMappedRefLen()/unitconv )
        #min_len = min( qrylen, reflen )
        min_len = self.getMinMappedLen() / (unitconv if unitconv > 0 else 1.)
        return self.Confidence / min_len


    #This calls the global fn (this file) CigarStringToList,
    # and then calls count on the output.
    # Note, this is very inefficient for large cigar strings (because you make the intermediate list).
    #  Then again, garbage collection should free that immediately.
    # Assume match char is "M".
    def getNCigarStringMatches(self) :
        return CigarStringToList(self.CigarString).count("M")


    #combine above two: confidence / N cigar string matches
    def getConfByNMatch(self) :
        count = CigarStringToList(self.CigarString).count("M")
        return self.Confidence / count if count else 0 #this must be > 0 if valid cigar, but just to be safe


    #invert an xmapEntry:
    #*only implemented for + orientation entries* (I think you need to force +, at least for viewing)
    # -Flip qry and ref contig ids.
    # -Flip qry and ref start and end pos
    # -Flip I and D in the cigar string
    def invert(self) :
        tmp = self.contigQry
        self.contigQry = self.contigRef
        self.contigRef = tmp

        tmp = self.QryStart
        self.QryStart = self.RefStart
        self.RefStart = tmp
        tmp = self.QryStop
        self.QryStop = self.RefStop
        self.RefStop = tmp

        newcigar = ""
        for char in self.CigarString :
            if char == "I" :
                newcigar += "D"
            elif char == "D" :
                newcigar += "I"
            else :
                newcigar += char
        self.CigarString = newcigar
    #end invert

#end class xmapEntry



#alignParams holds data from a .err file
class alignParams:
    def __init__(self, sourceFile=''): #sourceFile is .err
        if util.checkFile(sourceFile, ".err") : #must end in ".err"
            self.makeFromFile(sourceFile)
            return
        elif sourceFile : #argument given, but not a .err file
            print "Warning in alignParams.__init__, invalid .err file:", sourceFile
        
        self.header = ""
        self.filename = ""
        self.fp  = float(0) #false positive
        self.fn  = float(0) #false negative
        self.sf  = float(0) #site/fixed error/std dev (kb)
        self.sd  = float(0) #scaling error/std dev (kb)
        self.bpp = float(0) #bases per pixel--default should not matter (?)
        self.res = float(0) #resolution in pixels
        self.bppsd = 0.
        self.sr    = 0. #noise model quadratic term
        self.nmaps = 0
        self.llrm  = 0. #log10LR/Maps
        self.goodmaps = 0
        self.llrgm = 0. #log10LR/GoodMaps
        self.fprate = 0.
        self.ressd  = 0.
        self.maprate = 0. #goodmaps/nmaps
    
    def makeFromFile(self, sourceFile): #sourceFile is .err
        self.filename = os.path.split(sourceFile)[1]
        self.header = ""
        last = '' #just take last line of the .err
        nextlast = ''
        comchar = ["#", "I"] #the column header starts with "I" for Iteration
        with open(sourceFile) as f1 :
            for line in f1 :
                if line[0] in comchar :
                    self.header += line
                nextlast = last #need for goodmaps
                last = line.strip()
        
        tokens = last.split()
        #tokens[0] is iteration--not necessary to store
        self.fp = float(tokens[1]) #false positive /100kb
        self.fn = float(tokens[2]) #false negative rate
        self.sf = float(tokens[3]) #site/fixed error/std dev (kb)
        self.sd = float(tokens[4]) #scaling error/std dev (kb)
        self.bpp = float(tokens[5]) #bases per pixel
        self.res = float(tokens[6]) #resolution in pixels
        self.nmaps = int(tokens[7]) #n maps
        self.bppsd  = (float(tokens[11]) if len(tokens) > 11 else 0.)
        self.fprate = (float(tokens[12]) if len(tokens) > 12 else 0.) #if index 12 exists (length is 13), it should be the FPrate
        self.sr     = (float(tokens[13]) if len(tokens) > 13 else 0.)
        self.ressd  = (float(tokens[14]) if len(tokens) > 14 else 0.)
        tokens = nextlast.split() #nextlast line
        self.llrm     = float(tokens[8]) #log10LR/Maps
        self.goodmaps = int(tokens[9]) #good maps
        self.llrgm    = float(tokens[10])
        self.maprate = (float(self.goodmaps)/self.nmaps if self.nmaps else 0)

    
    def printParams(self):
        print "false positive                  = ", self.fp
        print "false negative                  = ", self.fn
        print "site/fixed error/std dev (kb)   = ", self.sf
        print "scaling error/std dev (root kb) = ", self.sd
        print "bases per pixel                 = ", self.bpp
        print "resolution (pixels)             = ", self.res
    
    def getParamList(self):
        #since this is meant for subprocess.Popen, need strings as args
        ret      = ["-FP" , str(self.fp ) ]
        ret.extend(["-FN" , str(self.fn ) ])
        ret.extend(["-sf" , str(self.sf ) ])
        ret.extend(["-sd" , str(self.sd ) ])
        ret.extend(["-bpp", str(self.bpp) ])
        ret.extend(["-res", str(self.res) ])
        return ret

    def getParamHeader(self, spacer = "  "):
        return ("%5s"+spacer+"%5s"+spacer+"%5s"+spacer+"%5s"+spacer+"%4s"+spacer+"%3s"+spacer+"%4s"+spacer+"%5s"+spacer+"%5s") % ("fp","fn","sf","sd","sr","bpp","res","nmaps","mapr")

    def getParamString(self, spacer = "  ", nores=True, nonmaps=True, printfilename=False):
        ret =           ("%5.3f" % self.fp )
        ret += spacer + ("%5.3f" % self.fn )
        ret += spacer + ("%5.3f" % self.sf )
        ret += spacer + ("%5.3f" % self.sd )
        ret += spacer + ("%5.3f" % self.sr )
        ret += spacer + ("%3.0f" % self.bpp)
        if not nores :
            ret += spacer + ("%4.2f" % self.res)
        if not nonmaps : #do for nmaps and maprate
            ret += spacer + ("%5i"   % self.nmaps)
            ret += spacer + ("%5.3f" % self.maprate)
        if printfilename :
            ret += spacer + ("%s" % self.filename) 
        return ret

    def __str__(self) :
        return self.getParamString(nonmaps=False) #print nmaps and maprate

    #def writeToFile(self, outfile, fpd=1., fnd=0.15, sfd=0.25, sdd=0.15, bppd=500., resd=3.3, srd=0) :
    #    """Only write two lines to output: a 0 line with defaults as specified in arguments,
    #    and a 1 line with the data in this object."""
    #on second thought, it makes more sense for 0 and 1 to be the same params so that the llrs correspond to the correct params
    def writeToFile(self, outfile) :
        sep = "\t"
        f1 = "%.6f" #FP, FN, sf, sd, llr, and sr are all these
        outf = open(outfile, "w+")
        outf.write(self.header) #should include trailing newline
        formatstr = (sep+f1+sep+f1+sep+f1+sep+f1+sep+"%.2f"+sep+"%.3f"+sep+"%i"+sep+f1+sep+"%i"+sep+f1+sep+"%.2f"+sep+f1+sep+f1+sep+"%.3f"+sep+"%.3f"+sep+"%.3f"+sep+"\n")
        #outf.write("0"+formatstr % (fpd, fnd, sfd, sdd, bppd, resd, self.nmaps, self.llrm, self.goodmaps, self.llrgm, 0, 0, srd, 0, 0, 0))
        outf.write("0"+formatstr % (self.fp, self.fn, self.sf, self.sd, self.bpp, self.res, self.nmaps, self.llrm, self.goodmaps, self.llrgm, self.bppsd, self.fprate, self.sr, self.ressd, 0, 0))
        outf.write("1"+formatstr % (self.fp, self.fn, self.sf, self.sd, self.bpp, self.res, self.nmaps, 0, 0, 0, self.bppsd, self.fprate, self.sr, self.ressd, 0, 0))
        outf.close()
    #end writeToFile

#end class alignParams



class cmap:

    #if lengthonly == True, do not load entire map--just length
    def __init__(self, sourceFile='', lengthonly=False):
        self.featurePositions = []
        self.featureLookup = {}
        self.header = ""
        self.length = 0
        self.mapid  = 0
        self.channelid = 0
        self.nsites = 0
        self.quality = 0
        exists = os.path.isfile(sourceFile)
        if exists and sourceFile.endswith('.cmap'):
            self.makeFromFile(sourceFile, lengthonly)
        elif exists and sourceFile.endswith('.spots'):
            self.makeFromSpotsFile(sourceFile, lengthonly)
        elif len(sourceFile) > 0 :
            print "Error in cmap constructor: bad sourceFile:", sourceFile
        
        
    def makeFromFile(self, sourceFile, lengthonly=False):
        self.featurePositions = []
        self.featureLookup = {}
        self.header = ""
        FirstEntry = True
        f1 = open(sourceFile)
        for line in f1 : #no readline--that will iterate over each char in line instead of line itself
            if line[0] == '#':
                self.header += line
                continue

            mapid, length, nsites, channelID, cfeat = self.parseLine(line)
            self.populateData(mapid, length, nsites, channelID, cfeat, storeFirstData=FirstEntry)
            
            if FirstEntry: #data storing now done in populateLine
                FirstEntry = False
                if lengthonly :
                    break
        #end loop on sourceFile
        f1.close()
        assert len(self.featurePositions) == len(self.featureLookup), "cmap.makeFromFile: lens of featurePositions and featureLookup disagree in "+sourceFile #this amounts to a unique label position check
        assert lengthonly or self.nsites == len(self.featurePositions), "cmap.makeFromFile: lens of featurePositions is not equal to nsites in "+sourceFile #for lengthonly, featurePositions is first entry only
    #end def makeFromFile(...)


    #helper method for makeFromFile, also for use in multiCmap class
    #given a line from a cmap, add its contents to self.featurePositions and self.featureLookup
    #return all data that needs to be stored
    def parseLine(self, line) :

        #put everything in here bc line.strip/split can fail also
        try :
            tokens = line.strip().split()
            if len(tokens) < 9 : #need all of the columns
                print "Warning in mapClasses.cmap.parseLine: bad map line", line
                return None
            mapid     = int(tokens[0])
            length    = float(tokens[1])
            nsites    = int(tokens[2]) #used for firstEntry only
            #column (index) 3 is siteID--just index in self.featurePositions + 1 (start at 1 not 0)
            channelID = int(tokens[4])
            position  = float(tokens[5])
            dev       = float(tokens[6])
            cov       = float(tokens[7])
            occur     = float(tokens[8])
            quality   = (float(tokens[9]) if len(tokens) >= 10 else None) #optional 10th column: quality
        except (AttributeError,ValueError), e : #valueerror for cast, AttributeError for line.strip/split
            print "Warning in mapClasses.cmap.parseLine: bad line:\n"+line+"\n"+e
            return None
        #return in same format populateData takes:
        return mapid, length, nsites, channelID, cmapFeature(position, dev, cov, occur, quality)
    #end def parseLine

    
    #another helper method for makeFromFile, also for use in multiCmap class
    #take output of cmap.parseLine and fill appropriate 
    #if arg storeFirstData is true, store the length, mapid, and channelid in the appropriate data members
    def populateData(self, mapid, length, nsites, channelID, cfeat, storeFirstData=False) :

        if channelID == 0:
            if cfeat.quality :
                self.quality = cfeat.quality
            return None #this is just a dummy label (ie, end of map); don't store it

        #self.featurePositions.append(position)
        #self.featureLookup[position] = cmapFeature(position, dev, cov, occur)
        self.featurePositions.append(cfeat.position)
        self.featureLookup[cfeat.position] = cfeat
        if storeFirstData :
            self.mapid = mapid #int(tokens[0])
            self.length = length #float(tokens[1])
            self.nsites = nsites #int(tokens[2])
            self.channelid = channelID #int(tokens[4])
    #end def populateData

    
    def makeFromSpotsFile(self, spotsFile, lengthonly):
        self.featurePositions = []
        self.featureLookup = {}
        self.mapid = 1 #assume this is true (it's not in the spots file)
        self.channelid = 1 #assume this also (also not in the spots file)
        for line in open(spotsFile) :
            tokens = line.strip().split('\t')
            if line.startswith("#") or line.startswith("Nick") : #NickID starts the header row
                if line.find('Reference Size') != -1:
                    self.length = int(tokens[2])
                    if lengthonly :
                        break
                continue
            position = float(tokens[2])
            self.featurePositions.append(position)
            self.featureLookup[position] = cmapFeature(position, 0,0,0)
        #end loop on sourceFile
        assert len(self.featurePositions) == len(self.featureLookup), "lens of featurePositions and featureLookup disagree in "+spotsFile #this amounts to a unique label position check
    #end def makeFromSpotsFile(...)


    #if you have an empty cmap (ie, constructor with no args), make a fresh one from a list and a length
    #if no length, use last label pos + 20.
    #do not allow labels with pos < 20 -- actually, forget this...it's not worth it, just use list as-is
    #do not allow negative label pos--may throw exception if can't cast ele of labellist to float, or can't sort--always call sort on labellist
    def makeFromList(self, labellist, length=0, header="") :
        minlengthbuffer = 20. #for start/end of featurePositions
        labellist.sort() #remember, this returns nothing...
        self.featurePositions = list( labellist ) #list will make shallow copy
        self.mapid = 1 #assume this is true
        self.channelid = 1 #assume this also
        if header :
            self.header = header #now an arg
        #this will (intentionally) throw an exception if length is not castable to float
        if length > 0 and float(length) : 
            self.length = length
        else :
            self.length = self.featurePositions[-1] + minlengthbuffer
        #lastly, fill in featureLookup same as is done in makeFromSpotsFile
        self.featureLookup.clear() #make sure it's empty first
        #note--should not check position here because need every ele in featurePositions to be in featureLookup
        for position in self.featurePositions :
            self.featureLookup[position] = cmapFeature(position, 0,0,0)


    #Coverage is defined for each site. Just average all.
    def getAvgCoverage(self) :
        if len(self.featureLookup) == 0 : #no features, so bail
            return 0
        #the list comprehension returns a list of coverages. Just sum and divide by len.
        return sum(feat.cov for feat in self.featureLookup.values())/float(len(self.featureLookup))


    def getMedianCoverage(self) :
        if len(self.featureLookup) == 0 : #no features, so bail
            return 0
        return util.getMedian([feat.cov for feat in self.featureLookup.values()])


    #Occurrence is also defined for each site. Just average all.
    def getAvgOccurrence(self) :
        if len(self.featureLookup) == 0 : #no features, so bail
            return 0
        return sum(feat.occur for feat in self.featureLookup.values())/float(len(self.featureLookup))


    def getMedianOccurrence(self) :
        if len(self.featureLookup) == 0 : #no features, so bail
            return 0
        return util.getMedian([feat.occur for feat in self.featureLookup.values()])


    #Coverage and Occurrence are defined for each site. Get per-site ratio, then average.
    def getAvgOccByCov(self) :
        if len(self.featureLookup) == 0 : #no features, so bail
            return 0
        return sum(float(feat.occur)/feat.cov for feat in self.featureLookup.values())/float(len(self.featureLookup))


    def getMedianOccByCov(self) :
        if len(self.featureLookup) == 0 : #no features, so bail
            return 0
        return util.getMedian([float(feat.occur)/feat.cov for feat in self.featureLookup.values()])


    #get distance (in bp) from first label to last label--ie, length without space at end
    # this is useful when comparing to aligned length
    def getMaxLabelLength(self) :
        if len(self.featurePositions) == 0 : #no labels, so bail
            return 0
        return self.featurePositions[-1] - self.featurePositions[0]


    #average nick density is number of sites / length
    #Technically, this is label density, not nick density, but name this like existing multiCmap method
    #unitden = 1e5 means labels per 100kb
    def getAverageNickDensity(self, unitden=1e5) :
        if not self.length : #don't divide by zero
            return 0
        return len(self.featurePositions)*unitden/self.length


    #this method name is now completely misleading since it no longer even has the capability to draw
    #the doPlot and color arguments now do nothing--what they did do is commented out--also yOffset
    def draw(self, start, stop, orientForward = True, xOffset = 0, yOffset = 0, doPlot = True, scalef = 1., color="" ):
        if orientForward:
            xDat = [x + xOffset - start for x in self.featurePositions if ((x >= start) and (x<= stop))]
        else:
            xDat = [self.length - (x + xOffset - stop) for x in self.featurePositions if ((x >= stop) and (x<= start))] #org
            #xDat = [x for x in self.featurePositions if ((x >= stop) and (x<= start))] #try... (doxoffset True in drawXmapHit)
            xDat.reverse()
        if scalef != 1. : #easier to apply after making list above
            totlen = max(xDat) #math.fabs(stop - start) #stop/start are not label poses--use pos itself
            print "ref maxpos:", totlen
            for idx, ele in enumerate(xDat) : 
                xDat[idx] = ele * scalef #* ele/totlen #apply as fraction of length
        #if doPlot:
        #    yDat = [yOffset] * len(xDat)
        #    #pl.plot([0 + xOffset, self.length + xOffset], [yOffset, yOffset]) 
        #    pl.plot(xDat, yDat, '-'+color+'o')
        #    print '  Plotting Feature Count:', len(xDat)
        return xDat
    

    #this used to call pl.fill_between (commented out below)--now return xDat and yDat instead
    def DrawCoverage(self, start, stop, orientForward = True, xOffset = 0, yOffset = 0, scaleCov = 1.):
        if orientForward:
            xDat = [x + xOffset - start for x in self.featurePositions if ((x >= start) and (x<= stop))]
            yDat = [scaleCov * self.featureLookup[x].cov + yOffset for x in self.featurePositions if ((x >= start) and (x<= stop))]
            y0Dat = [yOffset for x in self.featurePositions if ((x >= start) and (x<= stop))]
        else:
            xDat = [self.length - (x + xOffset - stop) for x in self.featurePositions if ((x >= stop) and (x<= start))]
            yDat = [scaleCov * self.featureLookup[x].cov + yOffset for x in self.featurePositions if ((x >= stop) and (x<= start))]
            y0Dat = [yOffset for x in self.featurePositions if ((x >= stop) and (x<= start))]
        #pl.fill_between(xDat, yDat, y2 = y0Dat)
        return xDat, yDat


    #utility for use in cmap.writeToFile as well as multiCmap.writeToFile
    #do not include header because the multiCmap use will not use that
    # instead, that goes in the writeToFile method; all else goes here
    def getWriteString(self) :
        outstr = ""
        #this check is redundant for cmap.writeToFile, but it is necessary for multiCmap.writeToFile
        if len(self.featurePositions) != len(self.featureLookup):
            print "Error in cmap.getWriteString: lens of featurePositions and featureLookup disagree, mapid", self.mapid
            return "" #return a string
        numsites = len(self.featurePositions)
        idx = 1 #this is for SiteID--just count
        for pos in self.featurePositions:
            #I'm not sure about an opening tab. It seems to be three spaces.
            #pos = self.featureLookup[pos].position #already have this
            dev = self.featureLookup[pos].dev
            cov = self.featureLookup[pos].cov
            occ = self.featureLookup[pos].occur
            outstr += "   %i\t%.1f\t%i\t%i\t%i\t%.1f\t%.1f\t%i\t%i\n" % (self.mapid, self.length, numsites, idx, self.channelid, pos, dev, cov, occ)
            idx += 1
        #write last line--mapid, length, numsites are same; index is numsites+1; channelID is zero, pos is len, default dev, cov, occur
        outstr += "   %i\t%.1f\t%i\t%i\t%i\t%.1f\t%.1f\t%i\t%i\n" % (self.mapid, self.length, numsites, numsites+1, 0, self.length, 0, 1, 1)
        return outstr


    def writeToFile(self, filename):
        if not filename.endswith(".cmap"):
            print "Error in cmap.writeToFile: file to write must end in .cmap"
            return
        #check that len of featurePositions and featureLookup are same, otherwise this won't really work
        elif len(self.featurePositions) != len(self.featureLookup):
            print "Error in cmap.writeToFile: lens of featurePositions and featureLookup disagree, mapid", self.mapid
            return
        #done error checking
        outfile = open(filename, "w+")
        outfile.write(self.header)
        outfile.write(self.getWriteString())
        outfile.close()


    #return list (array) of pos of labels in range, where range is specified in fraction of length
    def getLabelPosFromRange(self, start, stop):
        poslist = []
        for pos in self.featurePositions:
            if pos > self.length*stop:
                break 
            elif pos > self.length*start:
                poslist.append(pos)
        return poslist


    #similar to above, but two key differences:
    # *arguments are in units of _bases_, not fractions of length* (this means startpos and stoppos)
    # just return number of labels, not a list of positions
    # Also important that self.featurePositions is sorted--this is assumed (it should always be)
    # there is no intermediate list made--just slicing of self.featurePositions
    #  Note: will not use featureLookup at all, just featurePositions
    # start of zero means beginning, stop of zero means end
    # if a label is at exactly startpos or stoppos, it is also included
    def getNLabelsInRange(self, startpos=0, stoppos=0):

        if not util.isIntorFloat(startpos) or not util.isIntorFloat(stoppos) or startpos < 0 or stoppos < 0 :
            print "Error in cmap.getNLabelsInRange: bad arguments:", startpos, stoppos
            return None

        #this is not well-defined, so bail with error message--if stoppos is 0, use end
        if stoppos != 0 and stoppos < startpos : 
            print "Error in cmap.getNLabelsInRange: stop arg must be > start arg:", startpos, stoppos
            return None

        #this is not an error--just full length
        if startpos == 0 and stoppos == 0 :
            return len(self.featurePositions)

        #this is zero labels (unless these are a label pos)
        if startpos == stoppos :
            return (1 if startpos in self.featurePositions else 0)

        #must loop because need the index of the label adjacent to startpos and stoppos
        startidx = stopidx = 0
        for idx, pos in enumerate(self.featurePositions) :
            #print idx, pos #debug
            #do this only for _first_ label past startpos
            if startpos > 0 and startidx == 0 and pos >= startpos :
                startidx = idx
                if stoppos == 0 : #no need to find stop
                    break
            if stoppos > 0 and pos >= stoppos : #unless pos == stoppos, you go past stoppos
                stopidx = idx
                if pos == stoppos :
                    stopidx += 1 #because otherwise, this pos will be excluded in the slice
                break #always break after stoppos

        #if startpos > 0, startidx also > 0. Actually, start pos could be 1k and the first label could be 5k, so no.

        #if stoppos == 0 but startpos is non-zero, stopidx is len(self.featurePositions)
        # that's equivalent to [startidx:]
        #also, if stoppos is past end of last label, include up to end
        if stoppos == 0 or stoppos > self.featurePositions[-1] :
            stopidx = len(self.featurePositions)

        #slicing keeps the first index of the slice, but discards the second
        nlab = len(self.featurePositions[startidx:stopidx])
        #print "getNLabelsInRange: nlab =", nlab #debug

        return nlab
    #end def getNLabelsInRange


    #Adjust the positions of all the labels after the insertion or deletion (deletion arg)
    #inputs: length to change by, pos to start change, whether to add (insertion) or subtract (deletion--last arg)
    def adjustLabelPos(self, startpos, deltalength, deletion=True):
        #need to go from end to beginning because in case of multiple insertions of the same segment, there's no way to know that the next label already exists
        featcopy = copy.deepcopy( self.featurePositions ) #so that when i change featurePositions, I don't change this too
        if not deletion :
            featcopy.reverse()

        for idx, pos in enumerate( featcopy ):
        #for idx, pos in enumerate(self.featurePositions): #old way
            #always do pos > startpos: do nothing if <=
            if pos <= startpos : #same forward or backward
                continue

            if deletion :
                orgidx = idx
            else :
                orgidx = len(self.featurePositions) - 1 - idx #-1 bc max idx == len -1

            if deletion: #if a segment is deleted, pos decreases by length of deleted segment
                self.featurePositions[orgidx] -= deltalength
            else: #if a segment is inserted, pos increases by length of inserted segment
                self.featurePositions[orgidx] += deltalength

            #Editing dictionaries in python is a pain in the ass. You need to add a new entry with the new key (old value) then delete the old entry.
            value = self.featureLookup[pos]
            if deletion:                                    #if a segment is deleted, pos decreases by length of deleted segment            
                self.featureLookup[pos-deltalength] = value #insert new entry which is old value with new key
            else:                                           #if a segment is inserted, pos increases by length of inserted segment
                self.featureLookup[pos+deltalength] = value #same, but different place
            del self.featureLookup[pos]                     #remove old entry

            #print "adjustLabelPos:", idx, self.featurePositions[orgidx], featcopy[idx], len(self.featurePositions), len(self.featureLookup) #debug

        #This should perhaps be an assert
        if len(self.featurePositions) != len(self.featureLookup): #just a dummy check
            print "Error in adjustLabelPos: lens of featurePositions and featureLookup disagree"
            #raise Exception("Length Exception") #debug (nice, eh?)

    #end def adjustLabelPos(...)


    #Invert the map. This means flip it around.
    # Change is pos_i = len - oldpos_i. 
    #  All old labels are removed and replaced with new ones with same attributes but different poses.
    def invertMap(self): #could take a range, but that's more complicated. Assume all.
        #if there are no labels, there's nothing to invert
        if len(self.featurePositions) == 0 :
            return

        newfeatlook = {} #fix rare case where the newpos of label N equals pos of label M
        for idx, pos in enumerate(self.featurePositions):
            newpos = self.length - pos

            #I think that it's safe to change the value of the current index so long as you don't delete it
            self.featurePositions[idx] = newpos

            #Add new label to new dict
            newfeatlook[newpos] = self.featureLookup[pos]

        del self.featureLookup #remove old
        self.featureLookup = newfeatlook #and replace with new
        
        #Need to call reverse because otherwise the beginning is at the end
        self.featurePositions.reverse()
        
        #print len(self.featurePositions), len(self.featureLookup) #debug
        #This should perhaps be an assert--same line as in adjustLabelPos
        if len(self.featurePositions) != len(self.featureLookup): #just a dummy check
            print "Error in cmap.invertMap: lens of featurePositions and featureLookup disagree"


    #just print the map in the range beg, end, where these are specified in fractions of the length
    def printRange(self, beg=0.0, end=1.0):
        for idx, pos in enumerate(self.featurePositions):
            if pos < self.length*beg :
                continue
            elif pos > self.length*end :
                break
            print "   %i\t%.1f\t%i\t%i\t%i\t%.1f\t%.1f\t%i\t%i\n" % (self.mapid, self.length, numsites, idx, self.channelid, pos, dev, cov, occ)
    #end def printRange


    #take another cmap as an argument, and if all sites are the same,
    #add the coverage and occurrence of the argument to the current map
    #map length must also be the same, as a sanity check
    #return True on failure (to distinguish from success, which is None)
    def addCovOcc(self, incmap) :
        if incmap.length != self.length :
            print "Error in cmap.addCovOcc: map lengths are not equal: %s %s" % (incmap.length, self.length)
            return True
        if len(incmap.featureLookup) != len(self.featureLookup) :
            print "Error in cmap.addCovOcc: n sites are not equal: %s %s" % (len(incmap.featureLookup), len(self.featureLookup))
            return True
        newlookup = {}
        for pos,site in self.featureLookup.iteritems() :
            #if sites have different positions, this operation is not reliable, so must bail
            if not incmap.featureLookup.has_key(pos) :
                print "Error in cmap.addCovOcc: site missing in arg: %i" % pos
                return True
            #dev is std dev--just copy from original
            newlookup[pos] = cmapFeature(pos, site.dev, site.cov + incmap.featureLookup[pos].cov,
                                         site.occur + incmap.featureLookup[pos].occur)
        assert len(newlookup) == len(self.featureLookup), ("Invalid length of new featureLookup: %s" % len(newlookup)) #not necessary...
        self.featureLookup = newlookup #replace old w new, and that's it
    #end addCovOcc


    #check that all sites in featurePositions are >= 20.0 and <= lengh-20.0
    #check that values of featureLookup are all cmapFeatures and keys are 1:1 with eles of featurePositions
    #nothing more, nothing less
    #will not modify self
    def isWellFormed(self) :
        minoffset = 20. #this is the minimum site position, and max is length minus this

        #make list of keys of featureLookup--remove from list when found in featurePositions
        # at the end, it should be an empty list
        flu_keys = self.featureLookup.keys()

        for pos in self.featurePositions :
            #check pos in range stated above
            if pos < minoffset or pos > self.length - minoffset :
                print "mapClasses.cmap.isWellFormed: invalid featurePosition:", pos
                return False

            #check that pos is a key in featureLookup
            if not pos in flu_keys :
                print "mapClasses.cmap.isWellFormed: featurePosition absent from featureLookup:", pos
                return False

            del flu_keys[flu_keys.index(pos)] #remove pos from flu_keys--will raise ValueError if not found

            if not isinstance(self.featureLookup[pos], cmapFeature) :
                print "mapClasses.cmap.isWellFormed: invalid featureLookup value at:", pos
                return False

        #end loop on featurePositions

        #all keys of featureLookup should have been checked and therefore removed from flu_keys
        if len(flu_keys) :
            print "mapClasses.cmap.isWellFormed: key of featureLookup absent from featurePositions:", flu_keys
            return False            

        return True
    #end def isWellFormed


    #header of statistics table
    def statisticsTableHeader(self) :
        headl = [
            "idx ", #row number in table
            "CmapID", #Cmap ID
            "Clength", #Cmap length
            "CNsites", #Cmap number sites
            "CavgDen", #Cmap average density
            "CavgCov", #Cmap average coverage
            "CmedCov", #Cmap median coverage
            "CavgOcc", #Cmap average occurrence
            "CmedOcc", #Cmap mediant occurrence
            "CaObC", #Cmap average (occurrence / coverage)
            "CmObC" #Cmap median (occurrence / coverage)
            ]
        return "  ".join(headl)


    #each row of statistics table is a cmap--return the relevant information as a string
    #lenunit is denominator for lengths--1e3 means kb
    def statisticsTableRow(self, lenunit=1e3) :
        datalist = []
        datalist.append( "%6s"   % self.mapid )
        datalist.append( "%9.1f" % (self.length/lenunit) )
        datalist.append( "%4i"   % self.nsites )
        datalist.append( "%4.1f" % self.getAverageNickDensity() ) #use default arg unit
        datalist.append( "%5.1f" % self.getAvgCoverage() )
        datalist.append( "%5.1f" % self.getMedianCoverage() )
        datalist.append( "%5.1f" % self.getAvgOccurrence() )
        datalist.append( "%5.1f" % self.getMedianOccurrence() )
        datalist.append( "%5.3f" % self.getAvgOccByCov() )
        datalist.append( "%5.3f" % self.getMedianOccByCov() )
        return "  ".join(datalist)

    #end def statisticsTabelRow

#end class cmap



#based on cmap.getWriteString (above), the number of columns is fixed:
# correct the header based on the number of fields, currently this is always 9
# this must be the value of the nfield argument
def fixHeaderColumns(header, nfield = 9) :
    if not header :
        return header
    newh = ""
    for line in header.split("\n")[:-1] : #split with "\n" will produce '' as last element--discard it
        if line.startswith("#h") :
            tokens = line.split()
            if len(tokens) > nfield : #otherwise nothing to fix
                #want #h<space><tab-separated fields>
                line = tokens[0]+" "+("\t".join(tokens[1:(nfield+1)])) # +1 for the #h
        newh += line+"\n" #restore newline removed by split
    return newh
#end fixHeaderColumns


class cmapFeature:
    
    def __init__(self, position, dev, cov, occur, quality=None): #quality is only optional arg
        self.position = position
        self.dev = dev
        self.cov = cov
        self.occur = occur
        self.quality = quality

#end class cmapFeature

        
def CigarStringToList(CigarString):
    CigarList = []
    curCount = ''
    for val in CigarString:
        if ord(val) <= 57:
            curCount += val
        else:
            CigarList += [val] * int(curCount)
            curCount = ''
    return CigarList



#class for handling multi-contig cmaps

class multiCmap :

    def __init__(self, infile='', lengthonly=False, warnNoInfile=True) :
        self.cmapdict = {} #keys are contig ids of cmap, values are cmap objects
        self.header = "" #reset in makeFromFile
        self.totalLength = 0 #add all cmap lengths; reset in makeFromFile
        if util.checkCmap(infile) : #infile is a cmap
            self.makeFromFile(infile, lengthonly)
        elif util.checkFile(infile, filesuff=".spots") : #infile is a spots--this means a single contig reference
            self.cmapdict[1] = cmap(infile) #assign spots file contig id 1, cmap constructor handles the rest
            self.header      = self.cmapdict[1].header #copy header
            self.totalLength = self.cmapdict[1].length
        elif infile and warnNoInfile :
            print "Error in cmap constructor: bad sourceFile:", infile

    def makeFromFile(self, sourceFile, lengthonly=False):
        self.cmapdict = {} #must clear cmapdict
        self.header = ""
        self.totalLength = 0 #add all cmap lengths
        FirstEntry = True
        currid = -1 #need to store the current map id, then make a new map each time it changes
        cdum = cmap() #just to call parseLine--sort of dumb, but whatever
        infile = open(sourceFile)
        for line in infile : 
            if line[0] == '#':
                self.header += line
                continue

            ret = cdum.parseLine(line) #read current line
            if ret == None : #if exception is raised in parseLine, None is returned
                continue
            mapid, length, nsites, channelID, cfeat = ret
            if mapid != currid : #first time this mapid has been encountered
                FirstEntry = True
                currid = mapid
                self.cmapdict[mapid] = cmap()
                self.totalLength += length
            else :
                FirstEntry = False

            if lengthonly :
                continue #here, don't break, because I want the length of all the contigs
            
            #store the data
            self.cmapdict[mapid].populateData(mapid, length, nsites, channelID, cfeat, storeFirstData=FirstEntry)
            
        #end loop over sourceFile
        #i should probably call 'isWellFormed' on each cmap...
        infile.close()
    #end def makeFromFile


    #return a list of the lengths of all the cmaps in self.cmaplist
    def getAllMapLengths(self) :
        return map(lambda x: x.length, self.cmapdict.values()) #simple, eh?

    
    #call getn50 on list of all map lengths as obtained from self.getAllMapLengths
    def getAllMapN50(self) :
        return util.getn50( self.getAllMapLengths() )


    #return a list of the contig ids of all the cmaps in self.cmaplist--they must be unique
    def getAllMapIds(self) :
        return self.cmapdict.keys() #though trivial, better to wrap for compatibility


    #total average nick density is total number of labels over total length
    #unitden = 1e5 means labels per 100kb
    def getTotalAverageNickDensity(self, unitden=1e5) :
        nlab   = sum( map(lambda x: len(x.featurePositions), self.cmapdict.values()) )
        totlen = sum( self.getAllMapLengths() ) / unitden
        return nlab/totlen


    #return the length of the cmap with index contigid
    #the index is just the key, so this would be a one liner, except...
    #if the key doesn't exist, you get a KeyError exception which I can't let out of CharacterizeModule
    #so have to avoid this
    def getMapLength(self, contigid) :
        try :
            return self.cmapdict[contigid].length
        except : #most conservative is to catch all exceptions (?)
            print "Error in multiCmap.getMapLength for contigid =", contigid
            return 0


    #follow example of getMapLength
    def getMapAvgCoverage(self, contigid) :
        try :
            return self.cmapdict[contigid].getAvgCoverage()
        except : #most conservative is to catch all exceptions (?)
            print "Error in multiCmap.getMapAvgCoverage for contigid =", contigid
            return 0


    #haplotype contigs end in 0 (no split) or 1 and 2 (for a split pair):
    #get the size of all 0 maps plus the average of each 1 and 2 pair
    def getHaplotypeTotalMapLength(self) :
        totlen = 0
        for cid in self.cmapdict.keys() :
            if str(cid)[-1] == '0' :
                totlen += self.cmapdict[cid].length
            else : #since we only need sum, no need to average: just add half
                totlen += self.cmapdict[cid].length/2
        return totlen

    #like above, but return list
    def getHaplotypeMapLengths(self) :
        hapcmapdict = {} #keys are cmapdict keys transformed to haploid (changed to str and removed last digit)
        for cid in self.cmapdict.keys() : #cid are integers
            base,suf = str(cid)[:-1],str(cid)[-1] #base are keys of hapcmapdict
            if suf == '0' : 
                hapcmapdict[base] = self.cmapdict[cid].length
            elif hapcmapdict.has_key(base) : #must check if key present
                hapcmapdict[base].append( self.cmapdict[cid].length )
            else :
                hapcmapdict[base] = [self.cmapdict[cid].length]
        #end cmapdict.keys()
        #print "getHaplotypeMapLengths:", hapcmapdict #debug
        ret = [] #rather than just return lens after collapsing, make new list because modifying object which is looped on worries me
        for cid,lens in hapcmapdict.iteritems() :
            if isinstance(lens, list) :
                ret.append(sum(lens)/len(lens)) #len should always be 2, but just in case...
            else :
                ret.append(lens)
        #print "getHaplotypeMapLengths:", ret #debug
        return ret
    #end getHaplotypeMapLengths

    #write all maps in self.cmaplist to disk
    #outpath should be a path _plus_ a prefix like ./mydir/base,
    # eg, '../yh_1216_refguide/contigs2/chr6_split/yh_1216_chr6_contig'
    #  then the path '../yh_1216_refguide/contigs2/chr6_split/' must exist and be writable,
    #  and the file prefix will be 'yh_1216_chr6_contig'
    #  appended to the suffix will be the mapid of each map in cmaplist
    #  note this will not work for just a suffix--must include path such that os.path.split
    #   will return path, suffix
    # note that no uniqueness check is done for entries in that list--duplicates will clobber each other
    #mincontiglen will filter out contigs whose length is less than the argument if it's supplied
    # since cmap.length is in units of bases, this is too
    def writeAllMapsToDisk(self, outpath, mincontiglen=0, outsuf="") :

        if len(self.cmapdict) == 0 : #nothing to write
            return None #silent for now

        dirpath, suffix = os.path.split(outpath)
        if not util.checkDir(dirpath, checkWritable=True, makeIfNotExist=False) :
            print "Error in mapClasses.multiCmap.writeAllMapsToFile: bad path arg:", outpath
            return None

        self.header = fixHeaderColumns(self.header)
        for cmap in self.cmapdict.values() :
            if mincontiglen and cmap.length < mincontiglen :
                continue
            #copy header of original to each contig
            cmap.header = self.header
            cmap.writeToFile(outpath+str(cmap.mapid)+outsuf+".cmap")

    #end def writeAllMapsToDisk


    #unlike above, only a single file is created--another multiCmap--with all contigs in it
    def writeToFile(self, filename):

        if len(self.cmapdict) == 0 : #nothing to write
            return None #silent for now

        if not filename.endswith(".cmap"):
            print "Error in cmap.writeToFile: file to write must end in .cmap"
            return None

        outfile = open(filename, "w+")
        self.header = fixHeaderColumns(self.header)
        outfile.write(self.header)
        #multiCmaps are not necessarily sorted by mapid, but I think they often are
        #therefore, sort this also by mapids
        for cid in sorted(self.cmapdict.keys()) :
            outfile.write(self.cmapdict[cid].getWriteString())
        outfile.close()

        
    #inmcmap is another multiCmap--call method of this same name on each cmap in cmapdict
    #return True on failure, None on success
    def addCovOcc(self, inmcmap) :
        for cids, cmap in self.cmapdict.iteritems() :
            if not inmcmap.cmapdict.has_key(cids) :
                print "Error in multiCmap.addCovOcc: map id not found in argument: %i" % (cids)
                return True
            if cmap.addCovOcc( inmcmap.cmapdict[cids] ) :
                print "Error in multiCmap.addCovOcc: error in map id %i" % (cids)
                return True
    #end addCovOcc

    
    #Wrapper to the global function intoToChrPos which returns a dict as follows.
    #Keys of the dict are chromosome string as returned by intoToChrPos
    #Values are a pair (list of two eles).
    # The first item in the pair is zero becuse it's filled elsewhere (CharacterizeModule)
    # The second item is the total length of all the contigs in the same chromosome
    def makeChrSummaryDict(self) :
        retdict = {}
        for cids, cmap in self.cmapdict.iteritems() :
            chrs, pos = intToChrPos(cids)

            #make sure that the return is usable
            if not chrs :
                continue

            if retdict.has_key(chrs) :
                retdict[chrs][1] += cmap.length
            else :
                retdict[chrs] = [0, cmap.length]

        return retdict
    #end def makeChrSummaryDict(self) 


    #Call the cmap method statisticsTableRow for all cmaps
    #if xmappath supplied, look for xmaps whose filename is as follows
    #xmappath must have the xmap filename prefix, eg:
    # $PATH/exp_mrg5C_contig
    #also, if xmappath, check if .err file also exists
    def getAllMapsStatTable(self, molxmappath="", refxmappath="") :
        rets = ""
        for it,key in enumerate(self.cmapdict.keys()) : #keys are contig ids as ints
            if it == 0 :
                rets += self.cmapdict[key].statisticsTableHeader() + "  "
                rets += xmap().statisticsTableHeader(bool(molxmappath),bool(refxmappath)) #has trailing newline
            rets += ("%4i  "%(it+1)) + self.cmapdict[key].statisticsTableRow() #add index; no trailing newline
            #molecule-based alignments
            xmapfile = molxmappath+str(key)+".xmap"
            if molxmappath and util.checkFile(xmapfile) : #file exists
                rets += "  "+xmap(xmapfile).statisticsTableMolRow()
            elif molxmappath :
                print "Warning in getAllMapsStatTable: missing molecule xmap:", xmapfile
            errfile = molxmappath+str(key)+".err"
            if molxmappath and util.checkFile(errfile) : #file exists
                rets += "  "+alignParams(errfile).getParamString() #for err, no contig id is recorded
            elif molxmappath :
                print "Warning in getAllMapsStatTable: missing molecule err:", errfile
            #reference-based alignments
            xmapfile = refxmappath+str(key)+".xmap"
            if refxmappath and util.checkFile(xmapfile) : #file exists
                rets += "  "+xmap(xmapfile).statisticsTableRefRow(key) #for xmap, give contig id as arg
            elif refxmappath :
                print "Warning in getAllMapsStatTable: missing ref xmap:", xmapfile
            errfile = refxmappath+str(key)+".err"
            if refxmappath and util.checkFile(errfile) : #file exists
                rets += "  "+alignParams(errfile).getParamString() #for err, no contig id is recorded
            elif refxmappath :
                print "Warning in getAllMapsStatTable: missing ref err:", errfile
            rets += "\n"
        return rets


    #save a stat table to disk -- see getAllMapsStatTable
    def writeStatTableToDisk(self, filename, molxmappath="", refxmappath="", printtable=False) :
        table = self.getAllMapsStatTable(molxmappath, refxmappath)
        if printtable :
            print table
        outfile = open(filename, "w+")
        outfile.write(table)
        outfile.close


#end class multiCmap 



#Global fn which maps an integer to chromosome number (str) and start position (int) as follows.
#The first 'nchrdigits' digits in the integer encode the chromosome.
# They are equal to the chromosome unless 'ishuman' == True, in which case let
#  23 -> X and 24 -> Y and > 24 prints a warning and sets chromosome to zero
#The remaining digits are treated as the start position of the contig
# the units are 'posunits', ie, 1e4 = 10000 = 10kb.
# return this in units of 'retunits', default 1kb
#Note inidx must be base 10... (if you do '045', it's octal). This should never happen.
def intToChrPos(inidx, nchrdigits=2, posunits=1e4, retunits=1e3, ishuman=True, verbose=True) :
    if not util.isInt(inidx) :
        if verbose :
            print "Error in intToChrPos: invalid input (inidx):", inidx
        return "",0

    #if indix has equal or fewer digits than nchrdigits, nothing to do
    if inidx <= 10**nchrdigits - 1 :
        if verbose :
            print "Error in intToChrPos: input (inidx) too short:", inidx, nchrdigits
        return "",0

    #get chr value--first nchrdigits in inidx
    #This is probably slower, but cast to a string and then slice it--much simpler than doing arithmetic
    idxstr = str(inidx) #use for both chrn and posraw
    chrs   = idxstr[:nchrdigits] #string
    posraw = int(idxstr[nchrdigits:])

    #check chrn if ishuman
    if ishuman :
        if chrs == "23" :
            chrs = "X"
        elif chrs == "24" :
            chrs = "Y"
        elif int(chrs) > 24 : #humans (at least most of us) have only 24 chromosomes (depending on how you count)
            print "Warning in intToChrPos: chromosome number too large for human:", chrs, "using default 0"
            chrs = "0"

    #get chromosome position: remaining digits are position in units of posunits, convert to units of retunits
    retpos = posraw*posunits/retunits

    return chrs, retpos
#end intToChrPos
