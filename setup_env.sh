#!/bin/bash

VENV=ps2ff
PYVER=3.5

DEPARRAY=(numpy scipy pandas pytest pytest-cov)

# turn off whatever other virtual environment user might be in
source deactivate
    
# remove any previous virtual environments called pager
CWD=`pwd`
cd $HOME;
conda remove --name $VENV --all -y
cd $CWD
    
# create a new virtual environment called $VENV with the below list of dependencies installed into it
conda create --name $VENV --yes --channel conda-forge python=$PYVER ${DEPARRAY[*]} -y

#activate the new environment
source activate $VENV

#install some items separately
#conda install -y sqlalchemy #at the time of this writing, this is v1.0, and I want v1.1
conda install -y psutil

#tell the user they have to activate this environment
echo "Type 'source activate ${VENV}' to use this new virtual environment."
