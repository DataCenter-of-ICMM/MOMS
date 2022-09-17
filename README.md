MOMS
====

<b>moms</b> is a Multi-channel Optical Map Scaffolder, which utilizes multiple Bionano (R) optical maps in distinct channels simultaneously to scaffold a primary assembly generated from 2nd/3rd generation sequencing reads, with an aim to achieve chromosome-level assemblies.

### Citation

### Features
* Implemented as an user-friendly command line pipeline
* Modularized as encoding, rescaling, hybrid-scaffolding, sandwitch-scaffolding, anchoring, and reporting modules
* A Directed Node Graph (DNG) is used for representing association between optical maps
* Multi-threading support

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
```
### Usage
usage: moms.py [-h] -i FASTA -b CMAPS [CMAPS ...] -o OUTPUT [-f FORCE]
               [-m [MANCUTS [MANCUTS ...]]] [-c CLUSTER]
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

3. SGE cluster flag
  -c [yes|no] : the flag of whether to run the alignment in a SGE cluster, with a default value of 'yes'

4. number of threads
  -t num : the number of threads for running the pipeline, whose default value is 8.

5. pipeline configuration file
  --conf file.conf : the configuration file which specifies paths of scripts as well as parameters of various stages

6. verion information
  --version : print the version of MOMS being used

#### content of an example configuration file
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

