In-silico PCR on complete genomes.

Installing
----------

`sudo apt-get install -y build-essential g++ cmake git-all liblzma-dev zlib1g-dev libbz2-dev liblzma-dev`

`git clone --recursive https://github.com/gear-genomics/silica.git`

`cd silica`

`make all`

Running
-------

`./src/silica -g hg19.fa.gz sequences.fasta`
