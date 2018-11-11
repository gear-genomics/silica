# Silica
In-silico PCR on complete genomes.

Installing
----------

`sudo apt-get install -y build-essential g++ cmake git-all liblzma-dev zlib1g-dev libbz2-dev liblzma-dev libboost-date-time-dev libboost-program-options-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev`

`git clone --recursive https://github.com/gear-genomics/silica.git`

`cd silica`

`make all`

Install Genome Index
--------------------

Currently done by Indigo, place in fm folder.

Running
-------

`./src/silica -g hg19.fa.gz sequences.fasta`


Setup and run the server
------------------------

The server runs in a terminal

Install the dependencies:

`sudo apt install python python-pip`

`pip install flask flask_cors`

Start the server:

`cd PATH_TO_SILICA/silica`

`python server/server.py`

Setup and run the client
------------------------

The client requires a different terminal

Install the dependencies:

`cd PATH_TO_SILICA/silica/client`

`sudo apt install npm`

`sudo npm install`

Start the client:

`cd PATH_TO_SILICA/silica/client`

`npm run dev`


