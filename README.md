MOMS
====

<b>moms</b> is a Multi-channel Optical Map Scaffolder, which utilizes multiple Bionano (R) optical maps in distinct channels simultaneously to scaffold a primary assembly generated from 2nd/3rd generation sequencing reads, with an aim to achieve chromosome-level assemblies.

### Citation

### Features
* Implemented as an user-friendly command line pipeline
* Modularized as encoding, rescaling, hybrid-scaffolding, sandwitch-scaffolding, anchoring, and reporting modules
* A Directed Node Graph (DNG) is used for representing association between optical maps
* Multi-threading support

### Installation
Untar moms.tar.gz to your favourate directory. For example:

    $ mkdir -p ~/bin
    $ (cd ~/bin/; tar zxvf ~/moms.tar.gz)
    $ echo 'export PATH=~/bin/moms:$PATH' >> ~/.bashrc
    $ source ~/.bashrc
