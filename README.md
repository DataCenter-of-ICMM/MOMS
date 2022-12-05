MOMS
=====
<b>moms</b> is a Multi-channel Optical Map Scaffolder, which utilizes multiple Bionano &reg; optical maps in distinct channels simultaneously to scaffold a primary assembly generated from 2nd/3rd generation sequencing reads, with an aim to achieve chromosome-level assemblies.

### Background
Due to the existence of various repetitive sequences (such as centromeres, telomeres, etc.) and long spacer regions in eukaryotic genomes, the shortcomings of existing sequencing technologies are increasingly apparent. Although next-generation sequencing (NGS) technology can obtain a large amount of sequencing data cost-effectively, it cannot span longer repeat sequences due to read length limitations, which greatly reduces the continuity and integrity of the assembly. In addition, genome-wide structural variations cannot be easily and intuitively seen using NGS technology. Third-generation sequencing technology (TGS) greatly compensates for the problems caused by short read lengths in NGS due to the technical advantages of long read lengths, such as the assembly of repetitive sequences and the identification of larger structural variants. However, both NGS and TGS face the same problems caused by repetitive regions but to different degrees.

Single-molecule optical mapping of the genome can easily break through the limitations of sequencing technology, thereby laying the foundation for genome assembly at the chromosome level. A genome single-molecule optical map is an ordered genome-wide map of restriction enzyme cleavage sites derived from a single DNA molecule. It provides macroscopic framework, reflecting the macroscopic structural information of the entire genome. Therefore, it can assist in assembly, and improve the accuracy and completeness of results. At present, single-molecule optical mapping technology is mainly used for assisted genome assembly and/or large-scale structural variation detection.

There are two main technology providers of single-molecule optical mapping technology, *OpGen* and *BioNano*. *OpGen* launched the Argus&reg; system in 2010, which cuts single-molecule DNA immobilized on the surface area of MapCard DNA in situ through restriction enzymes, so that the sequence of the cut DNA fragments remains unchanged. The DNA fragments were stained with fluorescent dyes and placed under a fluorescence microscope to collect information on the size and sequence of each restriction endonuclease fragment, and obtain a genome-wide restriction endonuclease cleavage site map based on splicing. The *Irys&reg;/Saphyr&reg;* system launched by *BioNano&reg;* uses endonucleases to recognize and cut DNA and label fluorescence, and then use extremely fine capillary electrophoresis to straighten the DNA molecules, linearize each DNA single molecule, and carry out ultra-long single molecular high-resolution fluorescence imaging, that is, to generate a map of the distribution of enzyme cleavage sites. At present, *OpGen* has basically withdrawn from the market, and the only existing single-molecule physical optical mapping technology on the market is *BioNano&reg;*.

![saphyr](https://source.acexy.cn/view/YNogfkH)

![omaps](https://source.acexy.cn/view/YNohhCG)

Two important factors for the success of single-molecule optical mapping experiments are: the preparation of ultra-long, high-molecular-weight DNA molecules; and the choice of nucleic acid restriction enzymes. For the same sample, different restriction endonucleases will produce different optical maps. Breakpoints are likely to be generated in regions with too dense restriction sites, and insufficient information can be provided in regions with too sparse restriction sites., resulting in the inability to relocate base sequences.

In the early days of the introduction of single-molecule optical maps, due to the high cost of experiments, a typical protocol usually uses a computer program to assess the digestion density of each endonuclease on-the-shelf according to the assembly results from NGS/TGS. This selects the appropriate endonuclease. However, due to the inherent limitations of single-enzyme optical maps, when the budget allows, people began to consider the complementarity between multiple-enzyme maps to further improve the connectivity and integrity of the assembly results.

In order to solve the problem of multi-enzyme map assembly, *Bionano&reg;* company launched TGH software. In an example application, TGH performed assisted/hybrid assembly using two optical maps of BSPQI enzyme and BSSSI enzyme and NGS sequencing assembly results, achieving a 3.7-fold to 9.7-fold improvement in samples from three species including human, goat, and maize. The advantages of multi-channel optical mapping are demonstrated. However, TGH still has limitations. It cannot handle optical maps of three or more channels simultaneously for scaffolding.

As there are still plenty of large and complex genomes such as those of traditional chinese herbal medicines need to be de novo sequenced, there is an urgent need for the industry to provide more complete bioinformatics analysis tools for multi-enzyme optical mapping analysis. To this end, we have compiled an automated, high-performance, and scalable multi-channel optical map assisted scaffolding pipeline. It integrates some existing software such as alignment tool between optical maps and base sequences, as well as algorithms and scripts developed by our own, into a convenient and easy-to-use integrated software tool.

### Workflow

The final program consists of the following main steps:
1. Prepare the input file, automatically detect the endonuclease used, and use it for subsequent analysis;
2. According to the previously detected sequence of nucleotide bases recognized by the restriction enzyme, electronically digest the base sequence in the contig FASTA file to generate the NGS CMAP file;
3. For the BNG optical maps of multiple different channels from different enzyme digestions, perform a preliminary comparison with the NGS maps corresponding to the base sequences to obtain an appropriate scaling factor, and then adjust the coordinates in BNG optical maps accordingly;
4. Detect chimeric maps and interrupt them at appropriate positions;
5. Use the maps with all the chimeric sites broken to do single-enzyme map assisted hybrid scaffolding;
6. Mediated by NGS contigs, a modified Kruskal’s algorithm for maximum spanning tree (MST) is used to build the skeleton of the multi-channel map;
7. Compare the BNG cmaps with the skeleton of the multi-channel cmaps, and adjust the skeleton accordingly;
8. Compare the NGS cmaps with the skeleton of the multi-channel cmap to determine the anchor position of the contigs in the skeleton;
9. Generate the final scaffolding result and save it in FASTA format; meanwhile, save the linkage information between contigs in AGP format.

The workflow chart of MOMS is as follows:
![wkflow](https://source.acexy.cn/view/YNu7mOI)

### Core Features
1. Directed Node Graph (DNG)

MOMS uses a directed node graph to represent nucleic acid fragments and the relationships between nucleic acid fragments. A directed node graph is defined as a network structure consisting of directed nodes and edges. The node is directed, that is, a node is divided into two endpoints. Entering the node from different endpoints corresponds to the different orientations of the node in the path. In the process of traversing the graph, entering from one endpoint must come out from the other. Once the traversal order is determined, the single strand of the nucleic acid fragment corresponding to the node is also determined. The edge in the graph is undirected, and it expresses the adjacency relationship between nodes. According to the combination of connections between different endpoints, there are four adjacency relationships between two nodes, i.e. FF, RR, FR, and RF., corresponding to Forward-Forward, Reverse-Reverse, Forward-Reverse, and Reverse-Forward, respectively.

By traversing the directed node graph, under certain optimal criteria (such as the least conflict and the longest path length), we can reconstruct the positional relationship between nodes or nucleic acid fragments, thereby reconstructing the target genome.

Using this data structure, the transformation relationship between the marker coordinates of the optical maps of different channels can be limited to each local node, so as to avoid the accumulation of errors in a path traversal. Whether the coordinate position is correct or not directly affects the subsequent analysis of sequence replies, which is a key point in the entire analysis process.

2. Maximum Spanning Tree (MST) Algorithm

Genome scaffolding turned out to be an NP problem, i.e. a problem without any efficient polynomial solution. In the assembling process assisted by optical maps, as the number of consistent optical maps still reaches hundreds or thousands, it is impossible to obtain the solution in an effective time range according to the exhaustive solution method. Therefore, MOMS adopts a heuristic algorithm, a modified Kruskal’s algorithm for maximum spanning tree, to solve the problem. The actual assembly results prove that the assembly method combining the directed node graph and the maximum spanning tree is superior to the existing optical map assembly algorithms.

3. Multi-channel optical map alignment algorithm

After genome scaffolds being built, the base sequences in the NGS contigs need to be anchored to scaffolds. This is what a multi-channel optical map alignment algorithm can do. This algorithm, which is also implemented by dynamic programming, is an extension of the single-channel optical map alignment algorithm.

By utilizing multi-channel optical map alignments, we were able to put NGS contigs of smaller sizes back that would not be aligned and reposted in single-channel alignments due to insufficient information.

4. Multi-threading support

The time consuming stages of the pipeline, especially those contain optical map alignment, provide interface for multi-thread running.

### Dependencies
* RefAligner v7437 (or above) in Bionano Solve v3.2.1 (or above) (https://bionanogenomics.com/support/software-downloads/)
* Python 2.7.5 (https://www.python.org/downloads/release/python-275/)
* Perl 5.16 (https://www.perl.org/get.html)

### Operating System
* CentOS Linux release 7.6.1810

### Installation
Untar moms.tar.gz to your favourate directory. For example:
```bash
    $ mkdir -p ~/bin
    $ (cd ~/bin/; tar zxvf ~/moms.tar.gz)
    $ echo 'export PATH=~/bin/moms:$PATH' >> ~/.bashrc
    $ source ~/.bashrc
    $ (cd ~/bin/moms/scripts; ./setup.sh)
```
The file structure after installation is as follows:
```
.
├── bin
├── cpp
│   ├── common.h
│   ├── main.cpp
│   ├── Makefile
│   ├── omap.cpp
│   ├── omap.h
│   ├── worker.cpp
│   └── worker.h
├── moms.conf
├── moms.py
├── queue.conf
├── scripts
│   ├── align-final-bng.sh
│   ├── align-final-ngs.sh
│   ├── align-final-ngs-uni.sh
│   ├── assemble.sh
│   ├── bionano
│   │   ├── align_bnx_to_cmap.py
│   │   ├── AlignModule.py
│   │   ├── AssemblyModule.py
│   │   ├── CharacterizeModule.py
│   │   ├── compressPipeline.sh
│   │   ├── filter_SNR_dynamic.pl
│   │   ├── GroupedRefinementModule.py
│   │   ├── mapClasses.py
│   │   ├── Multithreading.py
│   │   ├── PairwiseModule.py
│   │   ├── pipelineCL.py
│   │   ├── Pipeline.py
│   │   ├── RefinementModule.py
│   │   ├── runAlignMol.py
│   │   ├── SampleCharModule.py
│   │   ├── SVModule.py
│   │   ├── utilities.py
│   │   └── xml
│   │       ├── alignArguments.xml
│   │       ├── assembleArguments.xml
│   │       ├── clusterArguments.xml
│   │       ├── conflictsArguments.xml
│   │       ├── globalArguments.xml
│   │       └── mergeArguments.xml
│   ├── check-files.sh
│   ├── conflictReport.sh
│   ├── fa2cmap.sh
│   ├── hybrid-scaffold.sh
│   ├── perl
│   │   ├── BNG
│   │   │   ├── refAlignerRun.pm
│   │   │   └── Utility.pm
│   │   ├── bnx_aligner.pl
│   │   ├── bnx_assembler.pl
│   │   ├── bnx_rescale.pl
│   │   ├── calc_cmap_stats.pl
│   │   ├── calc_fasta_stats.pl
│   │   ├── calc_xmap_stats.pl
│   │   ├── clusterAndLayout.pl
│   │   ├── cmap_aligner.pl
│   │   ├── cmap_aligner_two_passes.pl
│   │   ├── cmap_cutter.pl
│   │   ├── cmap_filter.pl
│   │   ├── cmap_fullset.pl
│   │   ├── cmap_gap_filler.pl
│   │   ├── cmap_joiner.pl
│   │   ├── cmap_merger_multi_color.pl
│   │   ├── cmap_merger.pl
│   │   ├── cmap_subset.pl
│   │   ├── cmap_uni.pl
│   │   ├── collect_evidences.pl
│   │   ├── collectUnusedContigs.pl
│   │   ├── conflicts_cutter.pl
│   │   ├── conflicts_identifier.pl
│   │   ├── edges_allocator.pl
│   │   ├── evidence_convertor.pl
│   │   ├── evidence_filter.pl
│   │   ├── ExportAGP.pl
│   │   ├── ExportAGP_TwoEnzyme.pl
│   │   ├── fa2cmap_multi_color.pl
│   │   ├── fitsLinearModel.pl
│   │   ├── msg.pm
│   │   ├── reviseArguments.pl
│   │   ├── xmap_associator.pl
│   │   ├── xmap_checker.pl
│   │   ├── xmap_converter.pl
│   │   ├── xmap_filter_by_cosine.pl
│   │   ├── xmap_filter.pl
│   │   ├── xmap_idmap.pl
│   │   └── xmap_to_OpGen_xml.pl
│   ├── prepare-input.sh
│   ├── python
│   ├── queue
│   │   ├── default.sh
│   │   └── main.sh
│   ├── report.sh
│   ├── rescale-cmaps.sh
│   ├── resolve-chimeras.sh
│   ├── sandwich-scaffold.sh
│   └── util.sh
└── test
    ├── ecoli-contigs.fa
    ├── EXP_REFINEFINAL1_BBVCI.cmap
    ├── EXP_REFINEFINAL1_BSPQI.cmap
    ├── EXP_REFINEFINAL1_BSSSI.cmap
    ├── run-moms.sh
    ├── test-omaligner
    │   ├── DLE1-BNG.bnx
    │   ├── DLE1-NGS.cmap
    │   └── run-omaligner.sh
    └── test-TGH
        ├── ecoli-contigs.fa -> ../ecoli-contigs.fa
        ├── EXP_REFINEFINAL1_BSPQI.cmap -> ../EXP_REFINEFINAL1_BSPQI.cmap
        ├── EXP_REFINEFINAL1_BSSSI.cmap -> ../EXP_REFINEFINAL1_BSSSI.cmap
        ├── hybridScaffold_two_enzymes.xml
        └── run-TGH.sh

```

### Usage
usage: moms.py [-h] -i FASTA -b CMAPS [CMAPS ...] -o OUTPUT [-f FORCE]
               [-m [MANCUTS [MANCUTS ...]]] [-v] [--cluster]
               [-t NUM_THREADS] [--conf CONF] [--version]
#### notes
MOMS automatically speculates internal parameters, such as those used for optical mapping alignment, as much as possible. Therefore, different from similar scaffolding tools, it only provides indispensable options in a succinct way. For those options/parameters that seldom changes between different runs, they can be specified in a configuration file.

#### required arguments
1. input FASTA file
    -i assembly.fasta : the genome assembly sequences in FASTA format assembled from NGS/TGS sequencing reads

2. input CMAP files
    -b file1.cmap file2.cmap ... : multiple Consensus optical MAPs (CMAPs) assembled from molecules nicking sequences, with each corresponding to a single nicking enzyme

3. output folder
    -o OUTPUT : the path of the output folder or directory

#### optional argument
1. Overriding flag
    -f [yes|no] : the flag of whether or not to override the output folder, with a default value of 'no'

2. Chimera-cutting files
    -m file1.txt file2.txt ... : the human investigated chrimera-cutting files for breaking chrimeral contigs/cmaps

3. Validation flag
    -v : the flag is set to run assembly validation only without OM-based scaffolding, with a default value of not set

4. SGE cluster flag
    --cluster : the flag is set to run the alignment in a SGE cluster, with a default value of not set

5. number of threads
    -t num : the number of threads for running the pipeline, whose default value is 8.

6. pipeline configuration file
    --conf file.conf : the configuration file which specifies paths of scripts as well as parameters of various stages

7. verion information
    --version : print the version of MOMS being used

### File formats
#### Input FASTA file
FASTA format is a text-based format for representing nucleotide sequences or peptide sequences, where base pairs or amino acids are represented using single [IUPAC](https://www.bioinformatics.org/sms/iupac.html) codes. Its sequence begins with a single-line description, followed by lines of sequence data. The description line, which begins with '>', gives a name and/or a unique identifier for the sequence, and may also contain additional information. An example of the file content is as follows:
```
>contig1
CTGTGGAACACTGTCTATGACAGCACTTTATGTGTAGCCTCATTTTCAAATGAGAATGCATTTACTGTATCACCATGAAGTAAGGTTGCAGTAGAGTTTTTGATAGATATTCTTCATTTGGTTAAGgaaattttcttccatctctcATGTGCTGATGATTTTTTTACAATTA
>contig2
aaaaaataaaataaaataaaaaactagccaggtgtacggtgcatgtctgtagtctgagctacttaggaggctgaggcaggaggatttcttgagcttagttcaaggctgcagtgagctatgatcatgccactgcactctagcctaggtaacCCTTTGAttattttcatcttccttttttaGTACCGAAATCCTAATTCTGATATGGCTATCCAGTGGCTA

```
#### Input CMAP files
BNG map file in cmap format, one file corresponds to one enzyme
Bionano Genomics&reg; CMAP files are raw data files that provide information on the location of marker sites in the genome map or the results of electronic digestion of reference or sequence data. It is a tab-delimited text-based file that can be opened in Microsoft Excel for easy reading, or in any text-based editor. This format file contains two parts: a CMAP information header, describing the format of the data; and a genome map information block containing the data values.

Among them, the CMAP information header contains: #CMAP file version, #tag channel, #enzyme digestion recognition motif, #number of consensus maps, #h and #f; the genome map information block contains: the first tag site in the map, the second tag site in the map, ..., [repeat until all tag sites], the last tag site is the end of the entire consensus genome map.

An example of the contents of a cmap file is as follows:
```
# hostname=bionano1b
# CMAP File Version:    0.2
# Label Channels:   1
# Nickase Recognition Site 1:   gctcttc;green_01
# Number of Consensus Maps: 4387
# Values corresponding to intervals (StdDev, HapDelta) refer to the interval between current site and next site
#h CMapId   ContigLength    NumSites    SiteID  LabelChannel    Position    StdDev  Coverage    Occurrence  ChimQuality SegDupL SegDupR FragileL    FragileR    OutlierFrac ChimNorm    Mask
#f int  float   int int int float   float   float   float   float   float   float   float   float   float   float   Hex
260 1623778.7   182 1   1   20.5    88.8    5.0 5.0 -1.00   -1.00   -1.00   23.08   0.00    0.00    -1.00   0
260 1623778.7   182 2   1   7592.8  89.9    5.0 5.0 -1.00   -1.00   -1.00   0.00    0.00    0.00    -1.00   0
260 1623778.7   182 3   1   17222.9 75.6    7.0 7.0 -1.00   -1.00   -1.00   0.00    0.00    0.00    -1.00   0
… …
```
#### Aligned XMAP files
The Bionano Genomics&reg; XMAP file is a cross-comparison between two maps. It reports the comparison derived from the alignment between a query CMAP file and a reference CMAP file. The data line displays the map start and end coordinates and the locations of the labels on the map using a tab-delimited text based file.

The XMAP file presents the information in two sections: the XMAP information header, which describes the specific format of the data; and the map alignment information block, which contains the data rows.

When imported into Bionano Access &trade;, the XMAP file is automatically filtered and ready for downstream analysis. XMAP files can also be opened in Excel for easy readability or in any tab-delimited, text-based editor.
![xmap](https://source.acexy.cn/view/YP0DWJa)

The XMAP file contains two sections: a XMAP information header, describing meta data as well as the format of the data; and an alignment information block containing specific data values.. 

Among them, the XMAP information header contains: #XMAP file version, #Reference maps, #Query maps, #h and #f; the alignment information block specifies the alignment information in accordance with the data format description in header. To be more specific, after the 3 IDs, is the first alignment of a reference map label to a query map label with orientation and confidence; then the (pseudo)-CIGAR string displays in HitEnum, followed by query and reference length and label channel; the final string shows the alignment label site in the map and is repeated for all label sites indexed per label color channel.

An example of the contents of a xmap file is as follows:
```
# hostname=bionano1b
# XMAP File Version:  0.2
# Label Channels:   1
# Reference Maps From:  /home/share2/hjiang/
# Query Maps From: 
#h XmapEntryID  QryContigID RefContigID QryStartPos QryEndPos   RefStartPos RefEndPos   Orientation Confidenence  HitEnum QryLen  RefLen  LabelChannel    Alignment
#f int          int         int         float       float       float       float       string      float       string  float   float   int	string
13  500 63  292146.8    123996.0    773333.0    941852.0    -   17.34   8M1D1M2D2M1D4M1D1M1D1M  295488.4    1996164.0   1   (110,38)(111,38)(112,37)(113,36)(114,35)(115,34)(116,33)(117,32)(118,31)(119,30)(120,30)(123,29)...
39  231 146 163243.6    371467.7    56646.0 256006.0    +   17.19   1M1D6M2D2M3I1M1D4M1D4M1D2M1D1M  512605.2    336000.0    1   (8,14)(10,15)(11,16)(12,17)(13,18)(14,19)(15,20)(18,21)(19,22)(20,26)(21,27)(22,27)(23,28)...
… …
```
#### Conflicts file
One advantage of the optical mapping technique is that it provides complementary information to sequencing technologies so that it can be used to detect chimeric contigs/scaffolds in a given assembly.

An example is illustrated in the following figure:
![conflict](https://source.acexy.cn/view/YP0YYqQ)
where a NGS contig and a BNG cmap are inconsistent before the sites noted by a red rectangle, suggesting that there may be a chimeric join in either the NGS contig or the BNG cmap.

MOMS uses the so-called conflicts file to store conflicts information. A conflicts file is a tab-delimited file, which contains two sections: a header line with a leading '#' char; and a conflicts block, each line of which presents a single conflict between a query map and a reference map.

An example of a conflicts file is as follows:
```
# xMapId    refQry  refId   leftRefBkpt rightRefBkpt    alignmentOrientation    refQry  qryId   leftQryBkpt rightQryBkpt    alignmentOrientation
13  ref 63  -1  941882  -   qry 500 -1  123966  -
39  ref 146 56616   256036  +   qry 231 163213.6    371497.7    +
40  ref 146 56616   -1  -   qry 248 272600.7    -1  -
58  ref 232 -1  495634  +   qry 19  -1  214571.1    +
59  ref 232 495574  -1  +   qry 144 94513.2 -1  +
… …
```
#### Configuration file
A MOMS configuration file is a formatted text file that includes paragraphs whose names are strings enclosed in square brackets. The paragraphs contain key-value pairs separated by equal signs '=', specifying the relevant system property values. These properties are for values that change infrequently and are used in conjunction with arguments entered on the command line.

The default configuration file is located in the same directory as moms.py, and the file name is: moms.conf. The default content is:
```
[parameters]
fa2cmap.minLabels=0
fa2cmap.minLength=0
assemble.minLen = 150
assemble.minSites = 8
assemble.pvalue = 1e-5
scaffold.single.fillgap=1
scaffold.mono.fillgap=1
anchor.unichannel=1

[paths]
scripts.dir=scripts
aligner.dir=scripts/bionano/binary
```
### Running example
````shell

```
<p style="color:blue">MOMS (Multiple-channel Optical Map Scaffolder) -- version 0.1.54</p>

```
Copyright (C) 2018-2022 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences. All Rights Reserved.

 Parameter used:
                FASTA file (-i): ecoli-contigs.fa
             Output folder (-o): .
 BNG cmap file for enzyme0 (-b): /home/bionano/MOMS/test/EXP_REFINEFINAL1_BSPQI.cmap
 BNG cmap file for enzyme1 (-b): /home/bionano/MOMS/test/EXP_REFINEFINAL1_BSSSI.cmap
            Num of threads (-t): 12
```
<p style="color:green">Step 1: Prepare input files</p>

```
[14:31:43] Run Command: mkdir -p /home/bionano/MOMS/test/Step-01_Input
[14:31:43] Run Command: /home/bionano/MOMS/scripts/prepare-input.sh /home/bionano/MOMS/test/ecoli-contigs.fa /home/bionano/MOMS/test/EXP_REFINEFINAL1_BSPQI.cmap,/home/bionano/MOMS/test/EXP_REFINEFINAL1_BSSSI.cmap /home/bionano/MOMS/test/Step-01_Input
```
<p style="color:orange">[2021-10-13 02:31:43 PM]</p>

<p style="color:green">Step 2: Encode NGS contigs to cmap file(s)</p>
```
[14:31:43] Run Command: mkdir -p /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding
[14:31:43] Run Command: /home/bionano/MOMS/scripts/fa2cmap.sh /home/bionano/MOMS/test/Step-01_Input/ecoli-contigs.fa 1:BSPQI:GCTCTTC,2:BSSSI:CACGAG /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding
Encoding with all enzyme channels ...
/home/bionano/MOMS/scripts/perl/fa2cmap_multi_color.pl -i /home/bionano/MOMS/test/Step-01_Input/ecoli-contigs.fa -e BSPQI:GCTCTTC 1 BSSSI:CACGAG 2  -o /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding
Encoding with single enzyme channels ...
/home/bionano/MOMS/scripts/perl/fa2cmap_multi_color.pl -i /home/bionano/MOMS/test/Step-01_Input/ecoli-contigs.fa -e BSPQI:GCTCTTC 1 -o /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding
/home/bionano/MOMS/scripts/perl/fa2cmap_multi_color.pl -i /home/bionano/MOMS/test/Step-01_Input/ecoli-contigs.fa -e BSSSI:CACGAG 1 -o /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding
```
<p style="color:orange">[2021-10-13 02:31:47 PM]</p>

<p style="color:green">Step 3: Rescale BNG cmap file(s)</p>
```
[14:31:47] Run Command: mkdir -p /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling
[14:31:47] Run Command: /home/bionano/MOMS/scripts/rescale-cmaps.sh /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs /home/bionano/MOMS/test/Step-01_Input/assembly /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling 12
Beginning initial NGS CMAP to BioNano CMAP alignment ...
Running command: /home/bionano/MOMS/scripts/perl/cmap_aligner.pl -s first -r /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs_BSPQI.cmap -q /home/bionano/MOMS/test/Step-01_Input/assembly_BSPQI.cmap -o /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/align0/BSPQI  -t 12 > /dev/null
Appending stderr to /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/align0/BSPQI.stdout
Running command: /home/bionano/MOMS/scripts/perl/cmap_aligner.pl -s first -r /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs_BSSSI.cmap -q /home/bionano/MOMS/test/Step-01_Input/assembly_BSSSI.cmap -o /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/align0/BSSSI  -t 12 > /dev/null
Appending stderr to /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/align0/BSSSI.stdout
Initial alignment complete in .460 s.

Rescaling BioNano CMAP ...
Running command: /home/bionano/MOMS/scripts/bionano/binary/RefAligner -merge -i /home/bionano/MOMS/test/Step-01_Input/assembly_BSPQI.cmap -readparameters /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/align0/BSPQI.errbin -o /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/align0/assembly_BSPQI_adjusted -stdout -stderr > /dev/null
Appending stderr to /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/align0/assembly_BSPQI_adjusted.stdout
Running command: /home/bionano/MOMS/scripts/bionano/binary/RefAligner -merge -i /home/bionano/MOMS/test/Step-01_Input/assembly_BSSSI.cmap -readparameters /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/align0/BSSSI.errbin -o /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/align0/assembly_BSSSI_adjusted -stdout -stderr > /dev/null
Appending stderr to /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/align0/assembly_BSSSI_adjusted.stdout
Rescaling complete in .074 s.
```
<p style="color:orange">[2021-10-13 02:31:48 PM]</p>

<p style="color:green">Step 4: Detect and resolve chimeral cmaps</p>
```
[14:31:48] Run Command: mkdir -p /home/bionano/MOMS/test/Step-04_Chimeras_resolution
[14:31:48] Run Command: /home/bionano/MOMS/scripts/resolve-chimeras.sh /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/assembly /home/bionano/MOMS/test/Step-04_Chimeras_resolution 12
Beginning initial NGS CMAP to rescaled BioNano CMAP alignment ...
Running command: /home/bionano/MOMS/scripts/perl/cmap_aligner.pl -r /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs_BSPQI.cmap -q /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/assembly_BSPQI_adjusted.cmap -o /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSPQI  -t 12 > /dev/null
Appending stderr to /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSPQI.stdout
      15 alignments are found for enzyme BSPQI between      7 BNG contigs and     12 NGS contigs.
Running command: /home/bionano/MOMS/scripts/perl/cmap_aligner.pl -r /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs_BSSSI.cmap -q /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/assembly_BSSSI_adjusted.cmap -o /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSSSI  -t 12 > /dev/null
Appending stderr to /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSSSI.stdout
      16 alignments are found for enzyme BSSSI between      5 BNG contigs and     14 NGS contigs.
Initial rescaled alignment complete in .444 s.

Beginning conflicts identification ...
Running command: /home/bionano/MOMS/scripts/perl/conflicts_identifier.pl -i /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSPQI.xmap -r /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSPQI_r.cmap -q /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSPQI_q.cmap -r0 /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs_BSPQI.cmap -q0 /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/assembly_BSPQI_adjusted.cmap -o /home/bionano/MOMS/test/Step-04_Chimeras_resolution/conflicts/BSPQI >/dev/null
    NONE of      7 BNG contigs have been flagged as conflicting for enzyme BSPQI.
    NONE of    195 NGS contigs have been flagged as conflicting for enzyme BSPQI.
Running command: /home/bionano/MOMS/scripts/perl/conflicts_identifier.pl -i /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSSSI.xmap -r /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSSSI_r.cmap -q /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSSSI_q.cmap -r0 /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs_BSSSI.cmap -q0 /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/assembly_BSSSI_adjusted.cmap -o /home/bionano/MOMS/test/Step-04_Chimeras_resolution/conflicts/BSSSI >/dev/null
    NONE of      6 BNG contigs have been flagged as conflicting for enzyme BSSSI.
    NONE of    195 NGS contigs have been flagged as conflicting for enzyme BSSSI.
Conflicts identification complete in .150 s.

Beginning conflicts cutting ...
Running command: /home/bionano/MOMS/scripts/perl/conflicts_cutter.pl -i /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSPQI.xmap -r /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSPQI_r.cmap -q /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSPQI_q.cmap -r0 /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs_BSPQI.cmap -q0 /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/assembly_BSPQI_adjusted.cmap -c /home/bionano/MOMS/test/Step-04_Chimeras_resolution/conflicts/BSPQI_conflicts.txt -s 1,2 -o /home/bionano/MOMS/test/Step-04_Chimeras_resolution/BSPQI  >/dev/null
       7 BNG contigs found for enzyme BSPQI after conflicts cutting.
     195 NGS contigs found for enzyme BSPQI after conflicts cutting.
Running command: /home/bionano/MOMS/scripts/perl/conflicts_cutter.pl -i /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSSSI.xmap -r /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSSSI_r.cmap -q /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSSSI_q.cmap -r0 /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs_BSSSI.cmap -q0 /home/bionano/MOMS/test/Step-03_BNG_CMAP_rescaling/assembly_BSSSI_adjusted.cmap -c /home/bionano/MOMS/test/Step-04_Chimeras_resolution/conflicts/BSSSI_conflicts.txt -s 2,2 -o /home/bionano/MOMS/test/Step-04_Chimeras_resolution/BSSSI  >/dev/null
       6 BNG contigs found for enzyme BSSSI after conflicts cutting.
     195 NGS contigs found for enzyme BSSSI after conflicts cutting.
Conflicts cutting complete in .550 s.
```
<p style="color:orange">[2021-10-13 02:31:51 PM]</p>

<p style="color:green">Step 5: Perform single-enzyme hybrid scaffolding using chimera-resolved cmaps</p>
```
[14:31:51] Run Command: mkdir -p /home/bionano/MOMS/test/Step-05_Pre-scaffold
[14:31:51] Run Command: /home/bionano/MOMS/scripts/hybrid-scaffold.sh /home/bionano/MOMS/test/Step-04_Chimeras_resolution/ecoli-contigs /home/bionano/MOMS/test/Step-04_Chimeras_resolution/assembly /home/bionano/MOMS/test/Step-01_Input/ecoli-contigs.fa /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs /home/bionano/MOMS/test/Step-05_Pre-scaffold 1 12
Beginning to merge NGS cmaps and BN cmaps to construct single-enzyme scaffolds ...
Running command: /home/bionano/MOMS/scripts/perl/cmap_merger.pl -f /home/bionano/MOMS/test/Step-04_Chimeras_resolution/ecoli-contigs_BSPQI_cut.cmap -s /home/bionano/MOMS/test/Step-04_Chimeras_resolution/assembly_BSPQI_adjusted_cut.cmap -e /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSPQI.errbin -o /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSPQI  -t 12 > /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSPQI/cmap_merge.stdout
Appending stderr to /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSPQI/Mrg.stdout
       4 super-contigs found for enzyme BSPQI.
Running command: /home/bionano/MOMS/scripts/perl/cmap_merger.pl -f /home/bionano/MOMS/test/Step-04_Chimeras_resolution/ecoli-contigs_BSSSI_cut.cmap -s /home/bionano/MOMS/test/Step-04_Chimeras_resolution/assembly_BSSSI_adjusted_cut.cmap -e /home/bionano/MOMS/test/Step-04_Chimeras_resolution/align1/BSSSI.errbin -o /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSSSI  -t 12 > /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSSSI/cmap_merge.stdout
Appending stderr to /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSSSI/Mrg.stdout
       3 super-contigs found for enzyme BSSSI.
Single-enzyme scaffolds construction completed in .877 s.

Beginning alignment of NGS cmap to Hybrid CMAP ...
Running command: /home/bionano/MOMS/scripts/perl/cmap_aligner_two_passes.pl -r /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSPQI.hybrid.cmap -q /home/bionano/MOMS/test/Step-04_Chimeras_resolution/ecoli-contigs_BSPQI_cut.cmap -o /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSPQI/BSPQI-NGS  -t 12 > /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSPQI/BSPQI-NGS-align.log
Appending stderr to /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSPQI/BSPQI-NGS_1st_pass.stdout
Appending stderr to /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSPQI/BSPQI-NGS_2nd_pass.stdout
      13 alignments found for enzyme BSPQI between     13 NGS contigs and      3 super-contigs.
Running command: /home/bionano/MOMS/scripts/perl/cmap_aligner_two_passes.pl -r /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSSSI.hybrid.cmap -q /home/bionano/MOMS/test/Step-04_Chimeras_resolution/ecoli-contigs_BSSSI_cut.cmap -o /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSSSI/BSSSI-NGS  -t 12 > /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSSSI/BSSSI-NGS-align.log
Appending stderr to /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSSSI/BSSSI-NGS_1st_pass.stdout
Appending stderr to /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSSSI/BSSSI-NGS_2nd_pass.stdout
      18 alignments found for enzyme BSSSI between     18 NGS contigs and      3 super-contigs.
align sequences to hybrid scaffolds completed in .693 s.

Beginning alignment of BNG cmap to Hybrid CMAP ...
Running command: /home/bionano/MOMS/scripts/perl/cmap_aligner.pl -r /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSPQI.hybrid.cmap -q /home/bionano/MOMS/test/Step-04_Chimeras_resolution/assembly_BSPQI_adjusted_cut.cmap -o /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSPQI/BSPQI-BNG  -t 12 > /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSPQI/BSPQI-BNG-align.log
Appending stderr to /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSPQI/BSPQI-BNG.stdout
       5 alignments found for enzyme BSPQI between      5 BNG contigs and      4 super-contigs.
Running command: /home/bionano/MOMS/scripts/perl/cmap_aligner.pl -r /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSSSI.hybrid.cmap -q /home/bionano/MOMS/test/Step-04_Chimeras_resolution/assembly_BSSSI_adjusted_cut.cmap -o /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSSSI/BSSSI-BNG  -t 12 > /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSSSI/BSSSI-BNG-align.log
Appending stderr to /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSSSI/BSSSI-BNG.stdout
       5 alignments found for enzyme BSSSI between      5 BNG contigs and      3 super-contigs.
align Bionano genome maps to hybrid scaffolds completed in .457 s.

Merging Hybrid CMAP with NGS not participated in the hybrid scaffold ...
Running command: /home/bionano/MOMS/scripts/perl/cmap_joiner.pl -f /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSPQI.hybrid.cmap -s /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSPQI/ecoli-contigs_BSPQI.non_used.cmap -o /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSPQI/step2.hybrid-NGS-non_used
Running command: /home/bionano/MOMS/scripts/perl/cmap_joiner.pl -f /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSSSI.hybrid.cmap -s /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSSSI/ecoli-contigs_BSSSI.non_used.cmap -o /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSSSI/step2.hybrid-NGS-non_used
Merging Hybrid CMAP with naive NGS CMAP complete in .076 s.

Merging Hybrid CMAP with BNG CMAP not participated in the hybrid scaffold ...
Running command: /home/bionano/MOMS/scripts/perl/cmap_joiner.pl -f /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSPQI.hybrid.cmap -s /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSPQI/assembly_BSPQI.non_used.cmap -o /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSPQI/step2.hybrid-BNG-non_used
Running command: /home/bionano/MOMS/scripts/perl/cmap_joiner.pl -f /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSSSI.hybrid.cmap -s /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSSSI/assembly_BSSSI.non_used.cmap -o /home/bionano/MOMS/test/Step-05_Pre-scaffold/CMAPs/BSSSI/step2.hybrid-BNG-non_used
Merging Hybrid CMAP with naive NGS CMAP complete in .056 s.

Beginning construction of AGP and FASTA file of the scaffolded and unscaffolded sequences ...
Running command: /home/bionano/MOMS/scripts/perl/ExportAGP.pl -i /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSPQI/BSPQI-NGS.xmap -c /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSPQI/BSPQI-NGS_r.cmap -s /home/bionano/MOMS/test/Step-01_Input/ecoli-contigs.fa -m /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs_BSPQI_key.txt -t /home/bionano/MOMS/test/Step-04_Chimeras_resolution/BSPQI_auto_cut_NGS_coord_translation.txt -o /home/bionano/MOMS/test/Step-05_Pre-scaffold/FASTAs/BSPQI/BSPQI > /home/bionano/MOMS/test/Step-05_Pre-scaffold/FASTAs/BSPQI/export.log
Running command: /home/bionano/MOMS/scripts/perl/ExportAGP.pl -i /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSSSI/BSSSI-NGS.xmap -c /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSSSI/BSSSI-NGS_r.cmap -s /home/bionano/MOMS/test/Step-01_Input/ecoli-contigs.fa -m /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs_BSSSI_key.txt -t /home/bionano/MOMS/test/Step-04_Chimeras_resolution/BSSSI_auto_cut_NGS_coord_translation.txt -o /home/bionano/MOMS/test/Step-05_Pre-scaffold/FASTAs/BSSSI/BSSSI > /home/bionano/MOMS/test/Step-05_Pre-scaffold/FASTAs/BSSSI/export.log
AGP and FASTA generation complete in .569 s.
```
<p style="color:orange">[2021-10-13 02:31:57 PM]</p>

<p style="color:green">Step 6: Perform multi-enzyme scaffolding mediated by NGS contigs</p>

```
[14:31:57] Run Command: mkdir -p /home/bionano/MOMS/test/Step-06_Final_Scaffold
[14:31:57] Run Command: /home/bionano/MOMS/scripts/sandwich-scaffold.sh /home/bionano/MOMS/test/Step-05_Pre-scaffold hybrid /home/bionano/MOMS/test/Step-04_Chimeras_resolution/ecoli-contigs /home/bionano/MOMS/test/Step-06_Final_Scaffold
Beginning to build relationship between the BN contigs from different enzyme assemblies ...
Running command: /home/bionano/MOMS/scripts/perl/xmap_associator.pl -f /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSPQI.xmap -s /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSSSI.xmap -t1 /home/bionano/MOMS/test/Step-04_Chimeras_resolution/BSPQI_auto_cut_NGS_coord_translation.txt -t2 /home/bionano/MOMS/test/Step-04_Chimeras_resolution/BSSSI_auto_cut_NGS_coord_translation.txt -o /home/bionano/MOMS/test/Step-06_Final_Scaffold/sandwich/BSPQI_BSSSI/glues_BSPQI_BSSSI
Running command: /home/bionano/MOMS/scripts/perl/fitsLinearModel.pl -i /home/bionano/MOMS/test/Step-06_Final_Scaffold/sandwich/BSPQI_BSSSI/glues_BSPQI_BSSSI.tsv -o /home/bionano/MOMS/test/Step-06_Final_Scaffold/sandwich/BSPQI_BSSSI/edges_BSPQI_BSSSI
Scaffolds/Contigs relationship building complete in .070 s.

Beginning to make the layout of the final assembly ...
Running command: /home/bionano/MOMS/scripts/perl/edges_allocator.pl /home/bionano/MOMS/test/Step-06_Final_Scaffold/sandwich/*/edges_*.tsv /home/bionano/MOMS/test/Step-06_Final_Scaffold/sandwich/*/glues_*.tsv /home/bionano/MOMS/test/Step-06_Final_Scaffold/sandwich
Running command: /home/bionano/MOMS/scripts/perl/clusterAndLayout.pl  -r /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSPQI/BSPQI-NGS_r.cmap -r /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSSSI/BSSSI-NGS_r.cmap -e /home/bionano/MOMS/test/Step-06_Final_Scaffold/sandwich/edges.tsv -g /home/bionano/MOMS/test/Step-06_Final_Scaffold/sandwich/glues.tsv -o /home/bionano/MOMS/test/Step-06_Final_Scaffold/sandwich/paths
Running command: /home/bionano/MOMS/scripts/perl/cmap_merger_multi_color.pl  -r /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSPQI/BSPQI-NGS_r.cmap -r /home/bionano/MOMS/test/Step-05_Pre-scaffold/anchor/BSSSI/BSSSI-NGS_r.cmap -l /home/bionano/MOMS/test/Step-06_Final_Scaffold/sandwich/paths.tsv -o /home/bionano/MOMS/test/Step-06_Final_Scaffold/multicolors
Running command: /home/bionano/MOMS/scripts/perl/cmap_subset.pl -i /home/bionano/MOMS/test/Step-06_Final_Scaffold/multicolors.cmap -c 1 -o /home/bionano/MOMS/test/Step-06_Final_Scaffold/mono/BSPQI
Running command: /home/bionano/MOMS/scripts/perl/cmap_subset.pl -i /home/bionano/MOMS/test/Step-06_Final_Scaffold/multicolors.cmap -c 2 -o /home/bionano/MOMS/test/Step-06_Final_Scaffold/mono/BSSSI
Cmap assembly complete in .083 s.
```
<p style="color:orange">[2021-10-13 02:31:58 PM]</p>

<p style="color:green">Step 7: Align BNG data to scaffolds</p>

```
[14:31:58] Run Command: mkdir -p /home/bionano/MOMS/test/Step-07_BNG_anchoring
[14:31:58] Run Command: /home/bionano/MOMS/scripts/align-final-bng.sh /home/bionano/MOMS/test/Step-06_Final_Scaffold /home/bionano/MOMS/test/Step-04_Chimeras_resolution/assembly adjusted_cut /home/bionano/MOMS/test/Step-07_BNG_anchoring 1 12
Beginning alignment of BNG contigs to the final scaffold ...
Running command: /home/bionano/MOMS/scripts/perl/cmap_aligner.pl -s BNG -r /home/bionano/MOMS/test/Step-06_Final_Scaffold/mono/BSPQI.cmap -q /home/bionano/MOMS/test/Step-04_Chimeras_resolution/assembly_BSPQI_adjusted_cut.cmap -o /home/bionano/MOMS/test/Step-07_BNG_anchoring/single/BSPQI/BSPQI  -t 12 > /home/bionano/MOMS/test/Step-07_BNG_anchoring/single/BSPQI/BSPQI-align.log
Appending stderr to /home/bionano/MOMS/test/Step-07_BNG_anchoring/single/BSPQI/BSPQI.stdout
    5 alignments found for enzyme BSPQI between      5 BNG contigs and      3 super-contigs.
Running command: /home/bionano/MOMS/scripts/perl/cmap_aligner.pl -s BNG -r /home/bionano/MOMS/test/Step-06_Final_Scaffold/mono/BSSSI.cmap -q /home/bionano/MOMS/test/Step-04_Chimeras_resolution/assembly_BSSSI_adjusted_cut.cmap -o /home/bionano/MOMS/test/Step-07_BNG_anchoring/single/BSSSI/BSSSI  -t 12 > /home/bionano/MOMS/test/Step-07_BNG_anchoring/single/BSSSI/BSSSI-align.log
Appending stderr to /home/bionano/MOMS/test/Step-07_BNG_anchoring/single/BSSSI/BSSSI.stdout
    5 alignments found for enzyme BSSSI between      5 BNG contigs and      3 super-contigs.
Alignment completed in .375 s.

Beginning to fill gaps using non-aligned BNG contigs ...
Running command: grep -v "^#" /home/bionano/MOMS/test/Step-07_BNG_anchoring/single/BSPQI.xmap | cut -f2 | sort | uniq | /home/bionano/MOMS/scripts/perl/cmap_filter.pl -c /home/bionano/MOMS/test/Step-04_Chimeras_resolution/assembly_BSPQI_adjusted_cut.cmap -i - -exclude -o /home/bionano/MOMS/test/Step-07_BNG_anchoring/gapfilled/BSPQI.non_used > /home/bionano/MOMS/test/Step-07_BNG_anchoring/gapfilled/BSPQI.nonused.log
Running command: /home/bionano/MOMS/scripts/perl/cmap_gap_filler.pl -i /home/bionano/MOMS/test/Step-06_Final_Scaffold/mono/BSPQI.cmap -a /home/bionano/MOMS/test/Step-07_BNG_anchoring/gapfilled/BSPQI.non_used.cmap --addonly -o /home/bionano/MOMS/test/Step-07_BNG_anchoring/gapfilled/BSPQI.filled  -t 12 > /home/bionano/MOMS/test/Step-07_BNG_anchoring/gapfilled/BSPQI.filled.log
Running command: grep -v "^#" /home/bionano/MOMS/test/Step-07_BNG_anchoring/single/BSSSI.xmap | cut -f2 | sort | uniq | /home/bionano/MOMS/scripts/perl/cmap_filter.pl -c /home/bionano/MOMS/test/Step-04_Chimeras_resolution/assembly_BSSSI_adjusted_cut.cmap -i - -exclude -o /home/bionano/MOMS/test/Step-07_BNG_anchoring/gapfilled/BSSSI.non_used > /home/bionano/MOMS/test/Step-07_BNG_anchoring/gapfilled/BSSSI.nonused.log
Running command: /home/bionano/MOMS/scripts/perl/cmap_gap_filler.pl -i /home/bionano/MOMS/test/Step-06_Final_Scaffold/mono/BSSSI.cmap -a /home/bionano/MOMS/test/Step-07_BNG_anchoring/gapfilled/BSSSI.non_used.cmap --addonly -o /home/bionano/MOMS/test/Step-07_BNG_anchoring/gapfilled/BSSSI.filled  -t 12 > /home/bionano/MOMS/test/Step-07_BNG_anchoring/gapfilled/BSSSI.filled.log
Gap-filling completed in .115 s.

Running command: /home/bionano/MOMS/scripts/perl/cmap_fullset.pl -i /home/bionano/MOMS/test/Step-07_BNG_anchoring/mono/BSPQI.cmap -i /home/bionano/MOMS/test/Step-07_BNG_anchoring/mono/BSSSI.cmap -o /home/bionano/MOMS/test/Step-07_BNG_anchoring/multicolors
Running command: /home/bionano/MOMS/scripts/perl/cmap_uni.pl -i /home/bionano/MOMS/test/Step-07_BNG_anchoring/multicolors.cmap -o /home/bionano/MOMS/test/Step-07_BNG_anchoring/unified
```
<p style="color:orange">[2021-10-13 02:32:00 PM]</p>

<p style="color:green">Step 8: Align NGS contigs to scaffolds</p>
```
[14:32:00] Run Command: mkdir -p /home/bionano/MOMS/test/Step-08_NGS_anchoring
[14:32:00] Run Command: /home/bionano/MOMS/scripts/align-final-ngs-uni.sh /home/bionano/MOMS/test/Step-07_BNG_anchoring /home/bionano/MOMS/test/Step-04_Chimeras_resolution/ecoli-contigs cut /home/bionano/MOMS/test/ecoli-contigs.fa 1:BSPQI:GCTCTTC,2:BSSSI:CACGAG /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding /home/bionano/MOMS/test/Step-08_NGS_anchoring 12
Beginning alignment of NGS contigs to the final scaffold using single channel ...
Running command: /home/bionano/MOMS/scripts/perl/cmap_aligner_two_passes.pl -r /home/bionano/MOMS/test/Step-07_BNG_anchoring/mono/BSPQI.cmap -q /home/bionano/MOMS/test/Step-04_Chimeras_resolution/ecoli-contigs_BSPQI_cut.cmap -o /home/bionano/MOMS/test/Step-08_NGS_anchoring/single/BSPQI/BSPQI  -t 12 > /home/bionano/MOMS/test/Step-08_NGS_anchoring/single/BSPQI/BSPQI-align.log
Appending stderr to /home/bionano/MOMS/test/Step-08_NGS_anchoring/single/BSPQI/BSPQI_1st_pass.stdout
Appending stderr to /home/bionano/MOMS/test/Step-08_NGS_anchoring/single/BSPQI/BSPQI_2nd_pass.stdout
      13 alignments are found for enzyme BSPQI between     13 NGS contigs and      3 super-contigs.
Running command: /home/bionano/MOMS/scripts/perl/cmap_aligner_two_passes.pl -r /home/bionano/MOMS/test/Step-07_BNG_anchoring/mono/BSSSI.cmap -q /home/bionano/MOMS/test/Step-04_Chimeras_resolution/ecoli-contigs_BSSSI_cut.cmap -o /home/bionano/MOMS/test/Step-08_NGS_anchoring/single/BSSSI/BSSSI  -t 12 > /home/bionano/MOMS/test/Step-08_NGS_anchoring/single/BSSSI/BSSSI-align.log
Appending stderr to /home/bionano/MOMS/test/Step-08_NGS_anchoring/single/BSSSI/BSSSI_1st_pass.stdout
Appending stderr to /home/bionano/MOMS/test/Step-08_NGS_anchoring/single/BSSSI/BSSSI_2nd_pass.stdout
      18 alignments are found for enzyme BSSSI between     18 NGS contigs and      3 super-contigs.
Single-channel aligment complete in .660 s.

Beginning alignment of unused NGS contigs to the final scaffold using unified channel ...
Running command: /home/bionano/MOMS/scripts/perl/collectUnusedContigs.pl -i /home/bionano/MOMS/test/ecoli-contigs.fa -m /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs_BSPQI_BSSSI_key.txt -t /home/bionano/MOMS/test/Step-08_NGS_anchoring/combined_NGS_coord_translation.txt -u /home/bionano/MOMS/test/Step-08_NGS_anchoring/single/used_NGS_id.txt -min 20000 -max 200000 -o /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/unused_fit_contigs > /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/contig-filtering.log
Running command: /home/bionano/MOMS/scripts/perl/fa2cmap_multi_color.pl -i /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/unused_fit_contigs.fa -e BSPQI:GCTCTTC 1 BSSSI:CACGAG 2  -o /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/encoding
Running command: /home/bionano/MOMS/scripts/perl/cmap_uni.pl -i /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/unused_fit_contigs_multicolors.cmap -o /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/unused_fit_contigs_unified
Running command: /home/bionano/MOMS/scripts/perl/cmap_aligner_two_passes.pl -r /home/bionano/MOMS/test/Step-07_BNG_anchoring/unified.cmap -q /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/unused_fit_contigs_unified.cmap -o /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/aligned --skip  -t 12 > /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/aligned.log
Appending stderr to /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/aligned_1st_pass.stdout
Running command: /home/bionano/MOMS/scripts/perl/xmap_filter_by_cosine.pl -x /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/aligned.xmap -r /home/bionano/MOMS/test/Step-07_BNG_anchoring/multicolors.cmap -q /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/unused_fit_contigs_multicolors.cmap > /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/filtered.xmap
Running command: /home/bionano/MOMS/scripts/perl/xmap_idmap.pl -i /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/filtered.xmap -k /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/encoding/unused_fit_contigs_BSPQI_BSSSI_key.txt -o /home/bionano/MOMS/test/Step-08_NGS_anchoring/unified/final.xmap
       4 alignments are found for between      4 NGS contigs and      3 super-contigs.
Unified-channel aligment complete in .074 s.

```
<p style="color:orange">[2021-10-13 02:32:03 PM]</p>

<p style="color:green">Step 9: Report the final scaffolds in AGP/FASTA format</p>
```
[14:32:03] Run Command: mkdir -p /home/bionano/MOMS/test/Step-09_Report
[14:32:03] Run Command: /home/bionano/MOMS/scripts/report.sh /home/bionano/MOMS/test/Step-08_NGS_anchoring/combined.xmap /home/bionano/MOMS/test/Step-01_Input/ecoli-contigs.fa /home/bionano/MOMS/test/Step-06_Final_Scaffold/multicolors.cmap /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs /home/bionano/MOMS/test/Step-09_Report/ecoli-contigs
Beginning to export AGP and FASTA files for multi-channel CMAPs ...
Running command: /home/bionano/MOMS/scripts/perl/ExportAGP_TwoEnzyme.pl -i /home/bionano/MOMS/test/Step-08_NGS_anchoring/combined.xmap -c /home/bionano/MOMS/test/Step-06_Final_Scaffold/multicolors.cmap -s /home/bionano/MOMS/test/Step-01_Input/ecoli-contigs.fa -o /home/bionano/MOMS/test/Step-09_Report/ecoli-contigs -m /home/bionano/MOMS/test/Step-02_NGS_CMAP_encoding/ecoli-contigs_BSPQI_BSSSI_key.txt -t /home/bionano/MOMS/test/Step-08_NGS_anchoring/combined_NGS_coord_translation.txt 2>&1 > /home/bionano/MOMS/test/Step-09_Report/export.log
Exporting AGP and FASTA files complete in .122 s.

Running command: ls /home/bionano/MOMS/test/Step*/*.cmap | xargs /home/bionano/MOMS/scripts/perl/calc_cmap_stats.pl -f -o /home/bionano/MOMS/test/Step-09_Report/cmap_file_stats.txt
Loading file: EXP_REFINEFINAL1_BSPQI.cmap
Loading file: EXP_REFINEFINAL1_BSSSI.cmap
Loading file: ecoli-contigs_BSPQI_BSSSI.cmap
Loading file: ecoli-contigs_BSPQI.cmap
Loading file: ecoli-contigs_BSSSI.cmap
Loading file: assembly_BSPQI_adjusted.cmap
Loading file: assembly_BSSSI_adjusted.cmap
Loading file: assembly_BSPQI_adjusted_cut.cmap
Loading file: assembly_BSSSI_adjusted_cut.cmap
Loading file: ecoli-contigs_BSPQI_BSSSI.cmap
Loading file: ecoli-contigs_BSPQI_cut.cmap
Loading file: ecoli-contigs_BSSSI_cut.cmap
Loading file: BSPQI.hybrid.cmap
Loading file: BSSSI.hybrid.cmap
Loading file: multicolors.cmap
Loading file: multicolors.cmap
Loading file: unified.cmap

Running command: ls -1 /home/bionano/MOMS/test/Step-09_Report/ecoli-contigs*.fa* | xargs /home/bionano/MOMS/scripts/perl/calc_fasta_stats.pl -f > /home/bionano/MOMS/test/Step-09_Report/fasta_file_stats.txt
Loading file: ecoli-contigs_NCBI.fasta
Loading file: ecoli-contigs_NOT_SCAFFOLDED.fasta

Running command: /home/bionano/MOMS/scripts/perl/calc_xmap_stats.pl -q /home/bionano/MOMS/test/Step-08_NGS_anchoring/combined.xmap -r /home/bionano/MOMS/test/Step-06_Final_Scaffold/multicolors.cmap > /home/bionano/MOMS/test/Step-09_Report/xmap_file_stats.txt
Loading file: /home/bionano/MOMS/test/Step-06_Final_Scaffold/multicolors.cmap
Loading file: /home/bionano/MOMS/test/Step-08_NGS_anchoring/combined.xmap
Calculating coverage ... done
```
<p style="color:orange">[2021-10-13 02:32:04 PM]</p>

```
Elapsed time: 21.35 s
Memory usage: 171648 K bytes
````
### Citation
