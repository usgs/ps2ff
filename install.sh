#!/bin/bash

VENV=ps2ff
PYVER=3.5

# Is conda installed?
conda=$(which conda)
if [ ! "$conda" ] ; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -O miniconda.sh;
    bash miniconda.sh -f -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
fi

conda update -q -y conda
conda config --prepend channels conda-forge

DEPARRAY=(numpy=1.11 scipy=0.19.1 pandas=0.20.3 pytest=3.2.0 pytest-cov=2.5.1 configobj=5.0.6)

# Is the Travis flag set?
travis=0
while getopts t FLAG; do
  case $FLAG in
    t)
      travis=1
      ;;
  esac
done

# Append additional deps that are not for Travis CI
if [ $travis == 0 ] ; then
    DEPARRAY+=(ipython=6.1.0 spyder=3.2.1 jupyter=1.0.0 seaborn=0.8.0 \
        sphinx=1.6.3)
fi

# Turn off whatever other virtual environment user might be in
source deactivate
    
# Remove any previous virtual environments called pager
CWD=`pwd`
cd $HOME;
conda remove --name $VENV --all -y
cd $CWD
    
# Create a new virtual environment
conda create --name $VENV -y python=$PYVER ${DEPARRAY[*]}

# Activate the new environment
echo "Activating the $VENV virtual environment"
source activate $VENV

# Install some items separately
pip install https://github.com/usgs/earthquake-impact-utils/archive/master.zip

# This package
echo "Installing ps2ff..."
pip install -e .

# Tell the user they have to activate this environment
echo "Type 'source activate ${VENV}' to use this new virtual environment."
