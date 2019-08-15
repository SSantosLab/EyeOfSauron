# general conda setup:
export CONDA_DIR=/cvmfs/des.opensciencegrid.org/fnal/anaconda2
source $CONDA_DIR/etc/profile.d/conda.sh
# this environment has the geckodriver installed: 
conda activate des18a

# run data preparation script
#python prep_data.py

# launch the notebook
theusername=`whoami`
thismachine=`hostname`
echo "Launching jupyter notebook with no browser."
echo "To interact with the notebook on your local machine, do:"
echo "ssh -N -L localhost:8888:localhost:8889 $theusername@$thismachine"
jupyter notebook --no-browser --port=8889

# old instructions:
#
# To run this notebook:
#   * the geckodriver executable needs to be in the environmental variable $PATH
export PATH=$PATH:/data/des30.a/data/annis/dae-haven/py-lib/lib/python2.7/site-packages/geckodriver/
#  * you'll need to run 
###conda activate des18a
#  * and set the python path:
PYTHONPATH=$PYTHONPATH:/data/des30.a/data/annis/dae-haven/py-lib/lib/python;export PYTHONPATH;
PYTHONPATH=$PYTHONPATH:/data/des30.a/data/annis/dae-haven/py-lib/lib/python2.7/site-packages;export PYTHONPATH;
       
# *(to run a noetbook remotely, see 
# http://home.fnal.gov/~kadrlica/fnalstart.html

