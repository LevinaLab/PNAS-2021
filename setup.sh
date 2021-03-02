#!/bin/bash
#Bash strartup
#installs nest,nestml,pip3
sudo apt-get update
sudo apt-get install libssl-dev
#Update Cmake
wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null |
    sudo apt-key add -
sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
sudo apt update
sudo apt install cmake
# Install python
sudo apt install python3-pip
alias python=python3
alias pip=pip3
#Make sure the default python is python3 
sudo update-alternatives --install /usr/bin/python python /usr/bin/python3.6 2
sudo python3 -m pip install ipython\
	numpy\
	scipy\
	sortedcontainers\ 
	matplotlib\
	seaborn\
	cython\
	networkx\
	wheel\
	nose\
	numba\
	h5py\
    xlrd\
    colorama==0.3.9

sudo apt-get install -y build-essential\
    cmake libltdl7-dev libreadline6-dev \
    libncurses5-dev libgsl0-dev libtool-bin\
    openmpi-bin libopenmpi-dev

#Install NEST
pip install pygsl
mkdir nest-build &&\
mkdir nest-install

git clone  https://github.com/nest/nest-simulator.git  -b v2.20.0
cd /home/$USER/nest-build
cmake  -DCMAKE_INSTALL_PREFIX:PATH=/home/$USER/.local /home/$USER/nest-simulator
make -j 18
make install
sudo echo "source /home/$USER/.local/bin/nest_vars.sh" >> ~/.bashrc
sudo echo "export LD_LIBRARY_PATH=/home/$USER/.local/lib/" >> ~/.bashrc
#Install MyModel (adLIF)
cd /home/$USER/cleanEI/nestmodels/
cmake -Dwith-nest=/home/$USER/.local/bin/nest-config .
make
make install
#echo "export LD_LIBRARY_PATH=/home/$USER/.local/lib/" >> ~/.bashrc
