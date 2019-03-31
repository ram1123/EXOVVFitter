#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess
import tarfile
import datetime
import commands

#os.system("sed -i 's/0.8,0.,1.0/0.5,0.0,0.8/' g1_exo_doFit_class_new.py")
#OUTDIR = 'WWTree_CleanedCode_Isolated_NaNFixed_Btag30GeV_2018_03_16_00h13_BothLSBUSB_UpDownVarWjet'
#OUTDIR = 'WWTree_CommonNtuple_For1and2Lepton_MuonPtScale_2018_07_09_18h38_LedLep50GeV'
#OUTDIR = 'WWTree_CWC_FirstComment_ZV'
OUTDIR = 'WWTree_After_CWR_ZV'
changes = raw_input("\n\nWrite change summary: ")

print "==> ",changes


currentDir = os.getcwd();
CMSSWDir =  currentDir+"/../";
print "PWD = ",currentDir

TestRun = 0

if OUTDIR == "":
	JobName = "job"
else:
	JobName = OUTDIR
# Get date and time for output directory
## ADD "test" IN OUTPUT FOLDER IF YOU ARE TESTING SO THAT LATER YOU REMEMBER TO WHICH DIRECTORY YOU HAVE TO REMOVE FROM EOS
os.system('xrdfs root://cmseos.fnal.gov/ mkdir /store/user/rasharma/BackgroundEstimation/' + OUTDIR)
if TestRun:
	outputFolder = "/store/user/rasharma/BackgroundEstimation/"+OUTDIR+'/'+datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M')+"_TEST/";
	OutputLogPath = "Logs/"+OUTDIR+'/' + datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M') + "_TEST";
else:
	outputFolder = "/store/user/rasharma/BackgroundEstimation/"+OUTDIR+'/'+datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M');
	OutputLogPath = "Logs/"+OUTDIR+'/' + datetime.datetime.now().strftime('%Y_%m_%d_%Hh%M');


print "Name of output dir: ",outputFolder
print "Name of Log dir: ",OutputLogPath
# create a directory on eos
os.system('xrdfs root://cmseos.fnal.gov/ mkdir ' + outputFolder)
os.system('xrdfs root://cmseos.fnal.gov/ mkdir ' + outputFolder + '/cards_em_HP')
os.system('xrdfs root://cmseos.fnal.gov/ mkdir ' + outputFolder + '/plots_em_HP')
# create directory in pwd for log files
os.system('mkdir -p ' + OutputLogPath)

def exclude_function(filename):
    if filename.endswith('.root'):
            return True
    else:
            return False

## Function to create a tar file
def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
            tar.add(source_dir, arcname=os.path.basename(source_dir),  exclude=exclude_function)

# Get CMSSW directory path and name
cmsswDirPath = commands.getstatusoutput('echo ${CMSSW_BASE}')
CMSSWRel = os.path.basename(cmsswDirPath[1])

print "CMSSW release used : ",CMSSWRel

# create tarball of present working CMSSW base directory
os.system('rm CMSSW*.tgz')
make_tarfile(CMSSWRel+".tgz", cmsswDirPath[1])

# send the created tarball to eos
os.system('xrdcp -f ' + CMSSWRel+".tgz" + ' root://cmseos.fnal.gov/'+outputFolder+"/" + CMSSWRel+".tgz")
os.system('git diff g1_exo_doFit_class_new.py > mypatch.patch')
os.system("sed -i '1s/^/Changes Summay : "+changes+"\\n/' mypatch.patch")
os.system('git log -1 --format="%H" >> mypatch.patch ')
os.system('xrdcp -f mypatch.patch root://cmseos.fnal.gov/'+outputFolder+'/mypatch.patch')



outJDL = open("runstep2condor_ZV.jdl","w");


outJDL.write("Executable = runstep2condor_ZV.sh\n");
outJDL.write("Universe = vanilla\n");
#outJDL.write("Requirements =FileSystemDomain==\"fnal.gov\" && Arch==\"X86_64\"");
outJDL.write("Notification = ERROR\n");
outJDL.write("Should_Transfer_Files = YES\n");
outJDL.write("WhenToTransferOutput = ON_EXIT\n");
#outJDL.write("include : list-infiles.sh |\n");
#outJDL.write("Transfer_Input_Files = "+inputlist+"\n");
outJDL.write("x509userproxy = $ENV(X509_USER_PROXY)\n");

outJDL.write("Output = "+OutputLogPath+"/"+JobName+".stdout\n");
outJDL.write("Error  = "+OutputLogPath+"/"+JobName+".stdout\n");
outJDL.write("Log = "+OutputLogPath+"/"+JobName+".log\n");
outJDL.write("Queue\n");
	    
outJDL.close();

#command = "python g1_exo_doFit_class_new.py -b -c em --mass 600 --category HPW --sample Signal_aQGC --jetalgo PuppiAK8_jet_mass_so --type vbf"
command = "python g1_exo_doFit_class_new_ZV.py -b "

outScript = open("runstep2condor_ZV.sh","w");
outScript.write('#!/bin/bash');
outScript.write("\n"+'source /uscmst1/prod/sw/cms/bashrc prod');
outScript.write("\n"+'echo "Starting job on " `date`');
outScript.write("\n"+'echo "Running on: `uname -a`"');
outScript.write("\n"+'echo "System software: `cat /etc/redhat-release`"');
outScript.write("\n"+'source /cvmfs/cms.cern.ch/cmsset_default.sh');
outScript.write("\n"+'### copy the input root files if they are needed only if you require local reading');
outScript.write("\n"+'xrdcp -s root://cmseos.fnal.gov/'+outputFolder+"/" + CMSSWRel+".tgz  .");
outScript.write("\n"+'tar -xf '+ CMSSWRel +'.tgz' );
outScript.write("\n"+'rm '+ CMSSWRel +'.tgz' );
outScript.write("\n"+'cd ' + CMSSWRel + '/src/EXOVVFitter' );
outScript.write("\n"+'rm -rf plots_*');
outScript.write("\n"+'rm -rf cards_*');
outScript.write("\n"+'echo "============================================" ');
outScript.write("\n"+'echo "====> List output files : " ');
outScript.write("\n"+'ls');
outScript.write("\n"+'echo "============================================" ');
outScript.write("\n"+'echo "pwd = $PWD"');
outScript.write("\n"+'echo "===========\t Project Rename 	=====================" ');
outScript.write("\n"+'scramv1 b ProjectRename');
outScript.write("\n"+'echo "===========\t Load CMS environment	=====================" ');
outScript.write("\n"+'eval `scram runtime -sh`');
outScript.write("\n"+'echo "===========\t scram clean and scram b ===================" ');
outScript.write("\n"+'scramv1 b clean; scramv1 b');
outScript.write("\n"+'echo "============================================" ');
outScript.write("\n"+'export PATH=${PATH}:'+currentDir);
outScript.write("\n"+'cd PDFs');
outScript.write("\n"+'rm *.pcm *.d *.so');
outScript.write("\n"+'echo "====> List output files : " ');
outScript.write("\n"+'ls');
outScript.write("\n"+'echo "============================================" ');
outScript.write("\n"+'root -l -b -q compilePdfs.C');
outScript.write("\n"+'echo "============================================" ');
outScript.write("\n"+'cd ../');
outScript.write("\n"+'echo "===========\t Run main script	=====================" ');
outScript.write("\n"+command);
outScript.write("\n"+'echo "====> List output files : " ');
outScript.write("\n"+'ls plots_em_HP/');
outScript.write("\n"+'echo "xrdcp output for condor"');
outScript.write("\n"+'xrdcp -f g1_exo_doFit_class_new.py root://cmseos.fnal.gov/'+outputFolder+'/');
outScript.write("\n"+'xrdcp -r -fs plots_em_HP root://cmseos.fnal.gov/'+outputFolder+'/plots_em_HP/');
outScript.write("\n"+'xrdcp -r -fs cards_em_HP root://cmseos.fnal.gov/'+outputFolder+'/cards_em_HP/');
outScript.write("\n"+'for FILE in *.root');
outScript.write("\n"+'do');
outScript.write("\n"+'\techo "xrdcp -r -f ${FILE} root://cmseos.fnal.gov/'+outputFolder+'/"');
outScript.write("\n"+'\txrdcp -r -f ${FILE} root://cmseos.fnal.gov/'+outputFolder+'/ 2>&1');
outScript.write("\n"+'done');
outScript.write("\n"+'cd ${_CONDOR_SCRATCH_DIR}');
outScript.write("\n"+'rm -rf ' + CMSSWRel);
outScript.write("\n");
outScript.close();
os.system("chmod 777 runstep2condor_ZV.sh");

print "===> Set Proxy Using:";
print "\tvoms-proxy-init --voms cms --valid 168:00";
print "\"condor_submit runstep2condor_ZV.jdl\" to submit";
os.system("condor_submit runstep2condor_ZV.jdl")
