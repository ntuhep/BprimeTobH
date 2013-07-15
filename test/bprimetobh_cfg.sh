#!/bin/bash

RUNDIR='/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_9_BpbH/src/BprimebHAnalysis/BprimeTobH/test/'
BATCHDIR=${PWD}

lean_up () {
	#try to recover log files and root files
	echo try to recover log files and root files ...
	cp -p *.root $RUNDIR
	cp -p *.log  $RUNDIR
	cp -p *.out  $RUNDIR
	exit
}
#LSF signals according to http://batch.web.cern.ch/batch/lsf-return-codes.html
trap clean_up HUP INT TERM SEGV USR2 XCPU XFSZ IO

cd ${RUNDIR}
echo Setting up ${PWD} as CMSSW environment.
eval `scramv1 runtime -sh` 
cd ${BATCHDIR}
echo The running directory is ${PWD}.
cp ${RUNDIR}/bprimetobh_cfg.py ${BATCHDIR}
time cmsRun bprimetobh_cfg.py maxEvents=-1 
cp -pu *.root *.log *.out ${RUNDIR} 



