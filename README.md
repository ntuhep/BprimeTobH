#BprimeTobH
#==========
#Bprime To b H  Analysis 
#Please execute the following commans in lxplus5

setenv SCRAM_ARCH slc5_amd64_gcc462
cmsrel CMSSW_5_3_11
cd  CMSSW_5_3_11/src
cmsenv
# QGTagger
cvs co -r v1-2-3 -d QuarkGluonTagger/EightTeV UserCode/tomc/QuarkGluonTagger/EightTeV
mkdir BprimebHAnalysis
cd BprimebHAnalysis 
git clone https://github.com/ntuhep/BprimeTobH.git
cd ..
scram b -j 8




