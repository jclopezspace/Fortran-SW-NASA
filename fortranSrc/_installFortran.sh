#!/bin/bash

# On Linux you need to run this first to ensure the script has execute permissions:
# chmod +x _installFortran.sh

# This script gets the latest apps to install and then installs gfortran
sudo apt-get update
sudo apt-get install -y gfortran