# ======================================
# Installation guide for OpenFOAM6 on Ubuntu18.04
# Author: Sheng Pan (DLUT)
# ======================================

# Step 1: Souece replacement by domestic mirrors for fast download speed:(Only needed for Chinese users)
  $ cd /etc/apt
  $ cp sources.list sources.list.bak
  $ sudo vim sources.list
  # Change all http://cn.archive.ubuntu.com with http://mirrors.aliyun.com
  source sources.list

# Step 2: System update & upgrade
  $ sudo apt-get update && sudo apt-get upgrade

# Step 3: Openmpi installation 
  # Prequisites:
  $ sudo apt-get install gcc g++ gfortran make libc6-dev dpkg-dev cmake 
  $ sudo apt-get update
  $ sudo apt-get install build-essential
  
  # Installation (It's highly recommanded to install openmpi with apt-get!):
  $ sudo apt-get install openmpi-bin openmpi-doc libopenmpi-dev

# Step 4: OpenFOAM installation 
  $ sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key | apt-key add -"
  $ sudo add-apt-repository http://dl.openfoam.org/ubuntu
  $ sudo apt-add-repository universe
  $ sudo apt-get update
  $ sudo apt-get -y install openfoam6
  # Then OpenFOAM6 and Paraview5.4 will be installed in /opt directory
