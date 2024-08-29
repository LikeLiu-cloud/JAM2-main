# JAM2: a Monte-Carlo event generator for high-energy nuclear collisions

This code simulates high energy nuclear collisions based on the transport theoretical models. 
Models include hadronic cascade model, relativistic quantum molecular dynamics (RQMD.RMF, RQMDs, RQMDv), and hydrodynamics.
Hybrid simulation of hydrodynamics + cascade model is available.

First download pythia8.307 code from
https://pythia.org/download/pythia83/pythia8307.tar.bz2
from https://pythia.org/releases/

Then Install Pythia8 to your preferred directory:

./configure --prefix=$HOME/lib/pythia8

make install

Then go to JAM source directory and in case there is no file 'configure' do

autoreconf -i

./configure PYTHIA8=$HOME/lib/pythia8 --prefix=$HOME/lib

make

make install

In case you set the environment variable 'PYTHIA8DATA', it should be identical to
the directory:'$HOME/lib/pythia8/share/Pythia8/xmldoc'.

1. Example files are in jam2/input/ 
You can find the options for cascade and meanfield and hydro modes in the file Manual.txt 
and the other options in the file JAM.cxx under jamDefaultParam() function.

2. You will need to change the random seed in the input file each time you run to have different set of events

3. There is also a file makeTree.C to copy the output to a root tree(Where you can choose particle spices you want to save to the root tree)

4. once you compile you should get an executable called jam in the folder jam2-main/jam2
Test as ./jam -f input_file

5. You should generate random seed when you submit jobs. Usually this can be just the JOBINDEX when you use star-submit. Just make sure the seed is incremented accordingly for each set of jobs you run so that same JOBINDEX will give a seed thats different from the previous set

6. submit example in /jam2-main/submitAll
