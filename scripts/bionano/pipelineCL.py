#!/usr/bin/env python

import sys
import atexit
import os

"""
@package pipelineCL 
Entry point of Pipeline for command line

Usage: python pipelineCL.py -h
"""

    
#this is moved here from utilities/Pipeline.py because those files are used outside of the pipeline
@atexit.register
def on_exit():
	try:
            util.LogStatus("pipeline", "exit", "%d" % os.getpid())
	except:
            pass


if __name__ == "__main__":
    import Pipeline
    varsP = Pipeline.varsPipeline()
    varsP.prerunChecks()
    
    print('  Prerun Tests:\n\t%d ERRORS\n\t%d WARNINGS\n' % (varsP.error, varsP.warning))
    if varsP.error or varsP.warning:
        #print(varsP.message)
	varsP.printMessage()
    if varsP.error:
        print('  EXITING: See errors') 
        sys.exit()
    
    dnpipeline = Pipeline.DNPipeline()
    dnpipeline.run(varsP)
    
