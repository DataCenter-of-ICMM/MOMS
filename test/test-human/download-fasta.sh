## download ref data
wget -c --no-check-certificate https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/013/985/GCA_001013985.1_ASM101398v1/GCA_001013985.1_ASM101398v1_genomic.fna.gz
## Break scaffolds into contigs
seqtk cutN -n 1 GCA_001013985.1_ASM101398v1_genomic.fna.gz >pacbio.fasta
