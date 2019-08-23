#!/bin/bash

/usr/bin/time -f "Elapsed time: %e s\nMemory usage: %M K bytes"  Rscript /home/bionano/tools/pipeline/1.0/HybridScaffold/1.0/runTGH.R -R /home/bionano/tools/pipeline/1.0/RefAligner/1.0/RefAligner -b1 EXP_REFINEFINAL1_BSPQI.cmap -b2 EXP_REFINEFINAL1_BSSSI.cmap -N ecoli-contigs.fa -e1 GCTCTTC -e2 CACGAG -t cur_results.tar.gz -s status.txt -f hybridScaffold_two_enzymes.xml -O output 2>&1 | tee run-TGH.log
