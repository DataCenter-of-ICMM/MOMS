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

    $ mkdir -p ~/bin
    $ (cd ~/bin/; tar zxvf ~/moms.tar.gz)
    $ echo 'export PATH=~/bin/moms:$PATH' >> ~/.bashrc
    $ source ~/.bashrc
