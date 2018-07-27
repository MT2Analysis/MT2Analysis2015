### This README is work-in-progress

## Setup and Installation
### MT2Analysis2015 must work within the cmsenv, which is set through a valid release of CMSSW
cmsrel CMSSW_8_0_12 # later versions are not tested
cd CMSSW_8_0_12/src
cmsenv
### You can now do the following steps from any directory you like
### Add a "link" to the remote repository where the code is
git remote add ana-mt2 https://github.com/MT2Analysis/MT2Analysis2015.git
### Then clone it locally from the remote
git clone -o ana-mt2 -b mg-data2017 https://github.com/MT2Analysis/MT2Analysis2015.git myMT2Analysis

### After installation, every time you log in
cd CMSSW_8_0_12/src
cmsenv

## Compilation
make <name-of-file-you-want-to-compile>
### Note: everytime you change interface/mt2.h you have to:
make clean
# and then recompile

# Run - time !
## Control Region Plots
### First you need to run the estimates
./regionEventYields <cfg-file-name> <data/mc/signal>
./llepControlRegion <cfg-file-name> <data/mc/signal>
./zllControlRegion <cfg-file-name> <data/mc/signal>
### Then you can make the Data/MC plots
./drawLostLeptonControlRegionDataMC <cfg-file-name> <lumi/shape>
./drawZllControlRegion <cfg-file-name> <lumi/shape>

### Next step: run the background estimation !
