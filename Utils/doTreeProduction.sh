#!/bin/bash

############################## configuration lines (TO DO: consider to move this into a separate file) #################################

#MEMO of the command necesasry to list files on T2: 
# env -i X509_USER_PROXY=~/.x509up_u`id -u` gfal-ls gsiftp://storage01.lcg.cscs.ch/pnfs/lcg.cscs.ch/cms/trivcat/store/user/mschoene/crab

# NB: Insert path starting from /store/user/... The rest will be added automatically

#inputProductionFolder="/store/user/mschoene/crab/8_0_26/reMini2017_runC_Feb28"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/reMini2017_runB_Feb28"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/reMini2017_runD_Feb28"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/reMini2017_runE_Feb28"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/reMini2017_runF_Feb28"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/reMini2017_runG_Feb28"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/reMini2017_runH_Feb28"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/MT2_gg_15Mar_CutBased15Id"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/MT2gg_22Mar_Signal"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/MT2gg_29Mar_Signal_bugged/"
#inputProductionFolder="/store/user/mschoene/crab/8_0_11/signal2017_Jan11/"
#inputProductionFolder="/store/user/mschoene/crab/8_0_11/mc2016_Nov15"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/gg_MC_Apr18/"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/gg_data2016_Sep04/"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/gg_Sig_May29/"

#inputProductionFolder="/store/user/mschoene/crab/9_2_4/data2017_reReco_jan16"
#inputProductionFolder="/store/user/mschoene/crab/9_2_4/data2017rR_jan21_newPhotonId/"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/data2016_05Feb18_newId/"
#inputProductionFolder="/store/user/mschoene/crab/9_4_1/data2017_05Feb18/"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/mc2016_05Feb18_newId/"
#inputProductionFolder="/store/user/mschoene/crab/9_4_1/mc2017_05Feb18"
#inputProductionFolder="/store/user/mschoene/crab/9_4_1/sig2017_05Feb18"

#inputProductionFolder="/store/user/mschoene/crab/8_0_26/sig2017_15Feb18/"

#inputProductionFolder="/store/user/mschoene/crab/8_0_26/data2016_27Feb18_hgg"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/sig2016_27Feb18_hgg"
#inputProductionFolder="/store/user/mschoene/crab/8_0_26/mc2016_27Feb18_hgg"
inputProductionFolder="/store/user/mschoene/crab/8_0_26/data2017_27Feb18_hgg"

era=2017    #this will change the input file and the output folder location

# In case you want to run the same production twice, adding a post-fix may help
#postFix=""
#postFix="_signal_pP_Apr19"
#postFix="_VH_pP_Jun19"
#postFix="_VH_pP_Test3Jun22"
#postFix="_pPSep08_smaller"
postFix="_pP_Mar01"

# For reading input from T2 (default):
site="lcg.cscs.ch"
se="storage01"
# or alternatively for reading from T3 (for legacy)
#site="psi.ch"
#se="t3dcachedb03"

# You should uncomment only one of the two, because data and MC production usually require different settings 
listOfSamplesFile="postProcessing2017-Data.cfg"  #for data inputs
#listOfSamplesFile="postProcessing2016-Data.cfg"  #for data inputs
#listOfSamplesFile="postProcessing2016-MC.cfg"   # for MC inputs

isCrab=1
inputPU="MyDataPileupHistogram.root"

#GoldenJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
GoldenJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt"


#"$PWD/gold_runF.txt"  #produced, for example for runE, with: filterJSON.py --min=276831 --max=277420 --output=gold_runE.txt gold_json.txt

###CHANGE
doSkimmingPruning=1 #1 as default; 0 for *_forQCD datasets (in data), which don't contain the necessary info to run the skimming and which are already pruned
applyJSON=1     #0 for MC
doFilterTxt=0   #0 for MC
doAllSF=0       #1 for MC
doPreProc=1     #0 (only 1 for large MC samples (almost all of them now!) or if you want to split MC samples, then run ./doTreeProduction pre first)




# You rarely need to edit the following cfg parameters
doSilver=0      #0 for MC
SilverJSON=$GoldenJSON
treeName="mt2"
fileExt="_post.root"
PUvar="nTrueInt"

# --- in current implementation one also needs to change this in runSkimmingPruning.sh ------
# current default values should be ok and you shouldn't need to touch anything here
useXRD="false"
gfalProtocol="gsiftp" # if useXRD disabled, use gfal via the given protocol
#gfalProtocol="srm" # alternative to gsiftp (gsiftp supposed to be more stable)
# -------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# NOTA BENE: the following files should also be (may have to) edited according 
# to the type of production (for example MT2 vs ZGamma, data vs MC):
#
# postProcessing2016-MC.cfg/postProcessing2016-Data.cfg
# skimmingPruning.cfg, skimmingPruningQCD.cfg, skimmingPruningMonoJet.cfg
#
# TO DO: add switches here in doTreeProduction.sh instead of force editing by hand
# ------------------------------------------------------------------------------


######################################### END OF CONFIGURATION PART #########################################################################





# ------------- Initialization -----------------------

host="${se}.${site}"

inputFolder="/pnfs/"$site"/cms/trivcat"$inputProductionFolder

productionName="$(basename $inputFolder)$postFix" 

outputFolder="/pnfs/psi.ch/cms/trivcat/store/user/`whoami`/MT2production/"$era"/PostProcessed/"$productionName"/"



jobsLogsFolder="./${productionName}"
workingFolder="/scratch/`whoami`/"$productionName

# warning: these 3 lines need to be replicated also in all other bash-script created on the fly (see next)
shopt -s expand_aliases
alias semkdir="env -i X509_USER_PROXY=~/.x509up_u`id -u` gfal-mkdir -p"
alias secp="env -i X509_USER_PROXY=~/.x509up_u`id -u` gfal-copy"




if [[ "$#" -eq 0 ]]; then
    echo "Relaunch the script with one of the following options: "
    echo "./doTreeProduction.sh pre       # pre-processing"
    echo "./doTreeProduction.sh post      # post-processing"
    echo "./doTreeProduction.sh mergeData # merge data and remove duplicates (not implemented yet)"
    echo "./doTreeProduction.sh postCheck # check post-processing"
    echo "./doTreeProduction.sh addAllSF  # add all scale factor weights"
    echo "./doTreeProduction.sh addISR    # add isr weights variables(not fully implemented yet for the standalone version)"
    echo "./doTreeProduction.sh addBtag   # add b-tagg weights variables"
    echo "./doTreeProduction.sh addLepSF  # add lepton weights variables"
    echo "./doTreeProduction.sh clean     # clean (not implemented yet)"
    exit;
fi;


##### beginning pre processing
if [[ "$1" = "pre" ]]; then
    if [ -d "$jobsLogsFolder" ]; then 
	echo "ERROR: the logFolder" $jobsLogsFolder " already exists."
	echo "Delete it and start from a clean area, or redirect the logs in a different place."
	echo "Exiting ..."
	exit
    else
	mkdir  $jobsLogsFolder
    fi

    if [ $CMSSW_BASE ]; then
	myCMSSW=$CMSSW_BASE  
    else
	myCMSSW=/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw-patch/CMSSW_7_4_12_patch4
    fi

    env -i X509_USER_PROXY=~/.x509up_u`id -u` gfal-ls ${gfalProtocol}://t3se01.psi.ch$outputFolder &> /tmp/checkOutputDir


    while read line; 
    do 
	case "$line" in \#*) continue ;; esac; #skip commented lines
	case "$line" in *"_ext"*) continue ;; esac; #skip extensions
	if [ -z "$line" ]; then continue; fi;  #skip empty lines
	id=`echo $line |awk '{print $1}'`
	name=`echo $line |awk '{print $2}'`

	fileList=inputChunkList_${name}.txt
	
	fileList=$jobsLogsFolder/inputChunkList_${name}.txt

	if [ -e inputChunkList_${name}.txt ]; then
	    echo "deleting the old file list"
	    rm $fileList
	fi;

	crabExt=""
	if [ ${isCrab} = 1 ]; then
	#crabExt=$(ls $inputFolder/$name/)
	    echo $(xrdfs $host ls $inputFolder/$name/)
	    crabExt=$(xrdfs $host ls $inputFolder/$name/)
	fi;
	echo "crabExt: " $crabExt;

	

	if [ ${isCrab} = 1 ]; then
	    for ((i=0; i<10; i++)); do  #BM: where this 10 is coming from ?? 
		xrdfs $host ls $crabExt/000${i}/ &> /dev/null;
		remoteDirectoryExist=`echo $?`
		if [ "$remoteDirectoryExist" -eq "0" ]; then
		    remoteFiles=`xrdfs $host ls $crabExt/000${i}/ |grep mt2`
		    for f in $remoteFiles; do
			echo $f>>$fileList
		    done;
		else
		    break;
		fi;
	    done;
	else
	    for f in $inputFolder/$name/mt2*.root; do
		echo $f>>$fileList
	    done;
	fi;


	numFiles=$(wc -l $fileList | awk '{print $1}')
	echo "number of files = " $numFiles

        ###And now for the extensions if they exist	
	crabExt=""
	if [ ${isCrab} = 1 ]; then
	#crabExt=$(ls $inputFolder/$name/)
	    echo "input to xrdfs: " xrdfs $host ls $inputFolder/${name}_ext/
	    echo $(xrdfs $host ls $inputFolder/${name}_ext/)
	    crabExt=$(xrdfs $host ls $inputFolder/${name}_ext/)
	fi;
	echo "crabExt: " $crabExt;

	
	if [ ${isCrab} = 1 ]; then
	    for ((i=0; i<10; i++)); do  #BM: where this 10 is coming from ?? 
		xrdfs $host ls $crabExt/000${i}/ &> /dev/null;
		remoteDirectoryExist=`echo $?`
		if [ "$remoteDirectoryExist" -eq "0" ]; then
		    remoteFiles=`xrdfs $host ls $crabExt/000${i}/ |grep mt2`
		    for f in $remoteFiles; do
			echo $f>>$fileList
		    done;
		else
		    break;
		fi;
	    done;
	else
	    for f in $inputFolder/$name/mt2*.root; do
		echo $f>>$fileList
	    done;
	fi;

	numFiles=$(wc -l $fileList | awk '{print $1}')
	echo "number of files = " $numFiles
    
        ###Seriously, the logic here could be improved
        ###And now for the SECOND extensions if they exist	
	crabExt=""
	if [ ${isCrab} = 1 ]; then
	#crabExt=$(ls $inputFolder/$name/)
	    echo "input to xrdfs: " xrdfs $host ls $inputFolder/${name}_ext2/
	    echo $(xrdfs $host ls $inputFolder/${name}_ext2/)
	    crabExt=$(xrdfs $host ls $inputFolder/${name}_ext2/)
	fi;
	echo "crabExt: " $crabExt;

	
	if [ ${isCrab} = 1 ]; then
	    for ((i=0; i<10; i++)); do  #BM: where this 10 is coming from ?? 
		xrdfs $host ls $crabExt/000${i}/ &> /dev/null;
		remoteDirectoryExist=`echo $?`
		if [ "$remoteDirectoryExist" -eq "0" ]; then
		    remoteFiles=`xrdfs $host ls $crabExt/000${i}/ |grep mt2`
		    for f in $remoteFiles; do
			echo $f>>$fileList
		    done;
		else
		    break;
		fi;
	    done;
	else
	    for f in $inputFolder/$name/mt2*.root; do
		echo $f>>$fileList
	    done;
	fi;

	numFiles=$(wc -l $fileList | awk '{print $1}')
	echo "number of files = " $numFiles



	if [ "$site" == "psi.ch" ]; then
	    sed -e "s#^#dcap://t3se01.psi.ch:22125/#" $jobsLogsFolder/inputChunkList_${name}.txt > $jobsLogsFolder/temp_${name}_dcap.txt
	else
    	    sed -e "s#^#root://$host/#" $jobsLogsFolder/inputChunkList_${name}.txt > $jobsLogsFolder/temp_${name}_dcap.txt
	fi
	mv $jobsLogsFolder/temp_${name}_dcap.txt $jobsLogsFolder/chunkPart_${name}.txt


	fileListTemp=$jobsLogsFolder/chunkPart_${name}.txt

	scriptName=batchScript_${name}.sh

	outputFile=${name}_pre #old with .cfg

	cat <<EOF > $scriptName

#!/bin/bash
shopt -s expand_aliases
alias semkdir="env -i X509_USER_PROXY=~/.x509up_u`id -u` gfal-mkdir -p"
alias secp="env -i X509_USER_PROXY=~/.x509up_u`id -u` gfal-copy"


#### The following configurations you should not need to change
# Job name (defines name seen in monitoring by qstat and the
#     job script's stderr/stdout names)
#$ -N preP_${name}_`whoami`

### Specify the queue on which to run
#$ -q short.q

# Change to the current working directory from which the job go
# submitted . This will also result in the job report stdout/stderr being
# written to this directory, if you do not override it (below).
#$ -cwd

# here you could change location of the job report stdout/stderr files
#  if you did not want them in the submission directory
####v$ -o $jobsLogsFolder/${name}.out
####c$ -e $jobsLogsFolder/${name}.err
#$ -o $jobsLogsFolder/${name}.out
#$ -e $jobsLogsFolder/${name}.err

source $VO_CMS_SW_DIR/cmsset_default.sh
#source /mnt/t3nfs01/data01/swshare/glite/external/etc/profile.d/grid-env.sh
export SCRAM_ARCH=slc6_amd64_gcc491
export LD_LIBRARY_PATH=/mnt/t3nfs01/data01/swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH
echo "Loading your CMSSW release"
echo "from $myCMSSW"
cd $myCMSSW
eval `scramv1 runtime -sh`
cd -
echo "preProcessing(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",$id,\"$fileListTemp\");"
echo "gSystem->Load(\"BTagCalibrationStandalone_cc.so\");gROOT->LoadMacro(\"preProcessing.C+\"); preProcessing(\"$name\",\"$inputFolder\",\"$outputFile\",\"$treeName\",$id,\"$fileListTemp\"); gSystem->Exit(0);" |root.exe -b -l ;


EOF

###qsub  -q short.q -l h_vmem=5g batchScript_${name}.sh;

#	qsub -q long.q -l h_vmem=5g $scriptName; 
	qsub -q all.q $scriptName;
	rm $scriptName;
	
    done < $listOfSamplesFile

fi;
##### done pre processing






########## post-processing ######################
if [[ "$1" = "post" ]]; then


#xrdfs t3dcachedb.psi.ch ls $outputFolder &> ./checkOutputDir
env -i X509_USER_PROXY=~/.x509up_u`id -u` gfal-ls ${gfalProtocol}://t3se01.psi.ch$outputFolder &> /tmp/checkOutputDir

# --- check the existence of outputFolder on SE ---
if [ -n "`cat /tmp/checkOutputDir|grep 'No such file or directory'`"  ]; then
    :
else
    echo "WARNING: output directory " $outputFolder "already exists."
    read -r -p "Are you sure you want to write there ? [y/N] " response
    if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
    then
	echo "Ok. Will proceed then ... "
    else
	echo "Exiting ..."
	exit;
    fi
fi
# --- 

if [ -d "$jobsLogsFolder" ]; then 
    echo "ERROR: the logFolder" $jobsLogsFolder " already exists."
    echo "Delete it and start from a clean area, or redirect the logs in a different place."
    echo "Exiting ..."
    exit
else
    mkdir  $jobsLogsFolder
fi

#gfalProtocol="gsiftp" # if useXRD disabled, use gfal via the given protocol
semkdir ${gfalProtocol}://t3se01.psi.ch/$outputFolder
python $PWD/convertGoodRunsList_JSON.py $GoldenJSON >& goodruns_golden.txt
#python $PWD/convertGoodRunsList_JSON.py $GoldenJSON >& $jobsLogsFolder/goodruns_golden.txt

echo $outputFolder
xrdfs t3dcachedb03.psi.ch cp $GoldenJSON $outputFolder/
secp file://$GoldenJSON ${gfalProtocol}://t3se01.psi.ch/$outputFolder/

if [ $doSilver -eq 1 ]; then
    secp file://$SilverJSON ${gfalProtocol}://t3se01.psi.ch/$outputFolder/
    python $PWD/convertGoodRunsList_JSON.py $SilverJSON >& $jobsLogsFolder/goodruns_silver.txt
fi

echo "Location of log files is: " $jobsLogsFolder
echo "Location of final files on SE is: " $outputFolder
echo "Working folder on working-node is: " $workingFolder


if [ $CMSSW_BASE ]; then
    myCMSSW=$CMSSW_BASE  
else
    myCMSSW=/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw-patch/CMSSW_7_4_12_patch4
fi

#get the pileUp histogram
if [ ! -f MyDataPileupHistogram.root ]; then
    pileupCalc.py -i $GoldenJSON --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69000 --maxPileupBin 50 --numPileupBins 50  MyDataPileupHistogram.root
fi


echo "gROOT->LoadMacro(\"goodrunClass.cc+\"); gSystem->Exit(0);" |root.exe -b -l ;
echo "gROOT->LoadMacro(\"leptonSF.cc+\"); gSystem->Exit(0);" |root.exe -b -l ;
echo "gROOT->LoadMacro(\"BTagCalibrationStandalone.cc+\"); gSystem->Exit(0);" |root.exe -b -l ;
echo "gSystem->Load(\"goodrunClass_cc.so\"); gSystem->Load(\"BTagCalibrationStandalone_cc.so\"); gROOT->LoadMacro(\"postProcessing.C+\"); gSystem->Exit(0);" |root.exe -b -l ;

while read line; 
do 
    case "$line" in \#*) continue ;; esac; #skip commented lines
    if [ -z "$line" ]; then continue; fi;  #skip empty lines
    id=`echo $line |awk '{print $1}'`
    name=`echo $line |awk '{print $2}'`
    xsec=`echo $line |awk '{print $3}'`
    filter=`echo $line |awk '{print $4}'`
    kfactor=`echo $line |awk '{print $5}'`

    doPruning="true"
    if [ $id -lt 10 ]; then
	doPruning="false"
	doAllSF="false"
    fi;

    if [[ $id -lt 10 && $doFilterTxt == 1 ]]; then
	# this is the merged (and duplicate removed) file including Nov15 txts and Oct 15 txts for all datasets
	# filterTxt=/shome/casal/eventlist_Nov14/filter_cscNov15_ecalscnNov15_cscOct15_sortu.txt
	# latest one from december
	filterTxt=/shome/mmasciov/CMSSW_7_4_7_MT2PostProcessing/src/analysisCode/Utils/allFilters_19Jan.txt
	outputFilteredFile=${workingFolder}/${name}_filtered$fileExt;
    fi;
    
    crabExt=""
    if [ ${isCrab} = 1 ]; then
	#crabExt=$(ls $inputFolder/$name/)
	echo "input to xrdfs: " xrdfs $host ls $inputFolder/$name/
	crabExt=$(xrdfs $host ls $inputFolder/$name/)
    fi;
    echo "crabExt: " $crabExt;

    
    #default is to not do the splitting into 10 files
    #also this should NEVER be done for MC as some scale factors need normalization
    #unless you did the preprocessing

    fileList=$jobsLogsFolder/inputChunkList.txt

    if [ -e  $fileList ]; then
	echo "deleting the old file list"
	rm $fileList
    fi;

    if [ ${isCrab} = 1 ]; then
	for ((i=0; i<10; i++)); do  #BM: where this 10 is coming from ?? laziness
	    #if [ -d $inputFolder/${name}/$crabExt/000${i}/ ]; then
	    xrdfs $host ls $crabExt/000${i}/ &> /dev/null;
	    remoteDirectoryExist=`echo $?`
	    if [ "$remoteDirectoryExist" -eq "0" ]; then
		echo "adding files in 000"${i}" folders"
		remoteFiles=`xrdfs $host ls $crabExt/000${i}/ |grep mt2`
		for f in $remoteFiles; do
		    echo $f>>$fileList
		done;
	    else
		break;
	    fi;
	done;
    else
	for f in $inputFolder/$name/mt2*.root; do
	    echo $f>>$fileList
	done;
    fi;
    
    numFiles=$(wc -l $fileList | awk '{print $1}')
    echo "number of files = " $numFiles
    
    ###And now for the extensions if they exist########################
    crabExt=""
    if [ ${isCrab} = 1 ]; then
	echo "input to xrdfs: " xrdfs $host ls $inputFolder/${name}_ext/
	crabExt=$(xrdfs $host ls $inputFolder/${name}_ext/)
    fi;
    echo "crabExt: " $crabExt;

    
    if [ ${isCrab} = 1 ]; then
	for ((i=0; i<10; i++)); do  #BM: where this 10 is coming from ?? 
	    xrdfs $host ls $crabExt/000${i}/ &> /dev/null;
	    remoteDirectoryExist=`echo $?`
	    if [ "$remoteDirectoryExist" -eq "0" ]; then
		remoteFiles=`xrdfs $host ls $crabExt/000${i}/ |grep mt2`
		for f in $remoteFiles; do
		    echo $f>>$fileList
		done;
	    else
		break;
	    fi;
	done;
    else
	for f in $inputFolder/$name/mt2*.root; do
	    echo $f>>$fileList
	done;
    fi;


    #numFiles=$(wc -l inputChunkList.txt | awk '{print $1}')
    numFiles=$(wc -l $fileList | awk '{print $1}')
    echo "number of files = " $numFiles

    ###For the second extentions if they exist, logic to be improved..
    ###And now for the SECOND extensions if they exist########################
    # crabExt=""
    # if [ ${isCrab} = 1 ]; then
    # 	echo "input to xrdfs: " xrdfs $host ls $inputFolder/${name}_ext2/
    # 	crabExt=$(xrdfs $host ls $inputFolder/${name}_ext2/)
    # fi;
    # echo "crabExt: " $crabExt;

    
    if [ ${isCrab} = 1 ]; then
	for ((i=0; i<10; i++)); do  #BM: where this 10 is coming from ?? 
	    xrdfs $host ls $crabExt/000${i}/ &> /dev/null;
	    remoteDirectoryExist=`echo $?`
	    if [ "$remoteDirectoryExist" -eq "0" ]; then
		remoteFiles=`xrdfs $host ls $crabExt/000${i}/ |grep mt2`
		for f in $remoteFiles; do
		    echo $f>>$fileList
		done;
	    else
		break;
	    fi;
	done;
    else
	for f in $inputFolder/$name/mt2*.root; do
	    echo $f>>$fileList
	done;
    fi;


    #numFiles=$(wc -l inputChunkList.txt | awk '{print $1}')
    numFiles=$(wc -l $fileList | awk '{print $1}')
    echo "number of files = " $numFiles

    maxNfiles=200
    counter=-1
    preProcFile=""
    if [[ $doPreProc ]]
    then    
	preProcFile=${name}_pre;
	echo $preProcFile
	echo "YOU SHOULD NOT BE HERE FFS" 
    fi

    if [[ (( $doPreProc -eq 1 && (( $numFiles -gt $maxNfiles )) ))  || (( $id -lt 10 )) ]]; then
   	echo "File will be split into multiple files for speed and memory limit purposes"
	counter=0;
	maxNfiles=200;
	if [[ $id > 10 ]]; then
	    preProcFile=${name}_pre;
	fi;
    fi;
    
    if [[ $doPreProc ]]
    then    
#	preProcFile=${name}_pre;
	echo $preProcFile
	echo "this better has the pre proc file name above" 
    fi

    while (( (( (( $numFiles  )) > (($counter * $maxNfiles)) )) || $(($counter < 0 )) )); do 

        counterFile=$jobsLogsFolder/chunkPart_${name}_${counter}.txt

        ##the input file list for the current range
	if [[ $counter -gt -1 ]]; then
	    echo "running on chunks $((counter*maxNfiles+1)) to $(((counter+1)*maxNfiles))"
 	    sed -n $((counter*maxNfiles+1)),$(((counter+1)*maxNfiles))p  $fileList > $jobsLogsFolder/chunkPart_${name}_${counter}.txt
	    if [ "$site" == "psi.ch" ]; then
		sed -e "s#^#dcap://t3se01.psi.ch:22125/#" $jobsLogsFolder/chunkPart_${name}_$counter.txt > $jobsLogsFolder/chunkPart_${name}_${counter}_dcap.txt
	    else
		sed -e "s#^#root://$host/#" $jobsLogsFolder/chunkPart_${name}_$counter.txt > $jobsLogsFolder/chunkPart_${name}_${counter}_dcap.txt
	    fi
	    mv $jobsLogsFolder/chunkPart_${name}_${counter}_dcap.txt $jobsLogsFolder/chunkPart_${name}_${counter}.txt
	fi;

	echo "Submitting the job for part" $counter; 

	#if [[ ${doPreProc} == 0 ]]; then
	if [[ $counter == -1 ]]; then
	    if [ "$site" == "psi.ch" ]; then
		sed -e "s#^#dcap://t3se01.psi.ch:22125/#" $jobsLogsFolder/inputChunkList.txt > $jobsLogsFolder/temp_${name}_${counter}_dcap.txt
	    else
    		sed -e "s#^#root://$host/#" $jobsLogsFolder/inputChunkList.txt > $jobsLogsFolder/temp_${name}_${counter}_dcap.txt
	    fi
	    mv $jobsLogsFolder/temp_${name}_${counter}_dcap.txt $jobsLogsFolder/chunkPart_${name}_${counter}.txt
	fi;

	counterName=${name}_${counter};

	if [[ $counter == -1 ]]; then
	    outputFile=${workingFolder}/${name}$fileExt;
	    counterName=${name}
	    scriptName=batchScript_${name}.sh
	else
	    scriptName=batchScript_${counter}_${name}.sh;
	    outputFile=${workingFolder}/${counterName}$fileExt;
	    if [[ $doFilterTxt == 1 ]]; then
		outputFilteredFile=${workingFolder}/${counterName}_filtered$fileExt;
	    fi;
	fi;

	cat <<EOF > $scriptName
#!/bin/bash
shopt -s expand_aliases
alias semkdir="env -i X509_USER_PROXY=~/.x509up_u`id -u` gfal-mkdir -p"
alias secp="env -i X509_USER_PROXY=~/.x509up_u`id -u` gfal-copy"

#### The following configurations you should not need to change
# Job name (defines name seen in monitoring by qstat and the
#     job script's stderr/stdout names)
#$ -N pP_${name}_`whoami`

### Specify the queue on which to run
#$ -q short.q

# Change to the current working directory from which the job go
# submitted . This will also result in the job report stdout/stderr being
# written to this directory, if you do not override it (below).
#$ -cwd

# here you could change location of the job report stdout/stderr files
#  if you did not want them in the submission directory
####v$ -o $jobsLogsFolder/${name}.out
####c$ -e $jobsLogsFolder/${name}.err
#$ -o $jobsLogsFolder/${counterName}.out
#$ -e $jobsLogsFolder/${counterName}.err

#source /mnt/t3nfs01/data01/swshare/psit3/etc/profile.d/cms_ui_env.sh
#source $VO_CMS_SW_DIR/cmsset_default.sh
#source //mnt/t3nfs01/data01/swshare/ROOT/thisroot.sh 
#eval \`scramv1 runtime -sh\`



source $VO_CMS_SW_DIR/cmsset_default.sh
#source /mnt/t3nfs01/data01/swshare/glite/external/etc/profile.d/grid-env.sh
export SCRAM_ARCH=slc6_amd64_gcc491
export LD_LIBRARY_PATH=/mnt/t3nfs01/data01/swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH
echo "Loading your CMSSW release"
echo "from $myCMSSW"
cd $myCMSSW
eval `scramv1 runtime -sh`
cd -


mkdir -p $workingFolder

semkdir ${gfalProtocol}://t3se01.psi.ch/$outputFolder

echo "postProcessing(\"$name\",\"$counterFile\",\"$outputFile\",\"$treeName\",$filter,$kfactor,$xsec,$id,\"$crabExt\",\"$inputPU\",\"$PUvar\",\"$jobsLogsFolder/goodruns_golden.txt\",\"$jobsLogsFolder/goodruns_silver.txt\",$applyJSON,$doAllSF,$doSilver,\"$preProcFile\"); gSystem->Exit(0);"

echo "gSystem->Load(\"goodrunClass_cc.so\"); gSystem->Load(\"leptonSF_cc.so\");  gSystem->Load(\"BTagCalibrationStandalone_cc.so\"); gROOT->LoadMacro(\"postProcessing.C\"); postProcessing(\"$name\",\"$counterFile\",\"$outputFile\",\"$treeName\",$filter,$kfactor,$xsec,$id,\"$crabExt\",\"$inputPU\",\"$PUvar\",$applyJSON,$doAllSF,$doSilver,\"$preProcFile\"); gSystem->Exit(0);" |root.exe -b -l ;

####echo "gSystem->Load(\"goodrunClass_cc.so\");  gSystem->Load(\"BTagCalibrationStandalone_cc.so\"); gROOT->LoadMacro(\"postProcessing.C\"); postProcessing(\"$name\",\"$counterFile\",\"$outputFile\",\"$treeName\",$filter,$kfactor,$xsec,$id,\"$crabExt\",\"$inputPU\",\"$PUvar\",\"$jobsLogsFolder/goodruns_golden.txt\",\"$jobsLogsFolder/goodruns_silver.txt\",$applyJSON,$doAllSF,$doSilver,\"$preProcFile\"); gSystem->Exit(0);" |root.exe -b -l ;



if [[ $id -lt 10 && $doFilterTxt == 1 ]]; then 
   echo "filterFromTxt(\"$filterTxt\",\"$outputFile\",1,\"$outputFilteredFile\",\"mt2\",\"\")"
   echo "gROOT->LoadMacro(\"filterFromTxt.C\"); filterFromTxt(\"$filterTxt\",\"$outputFile\",1,\"$outputFilteredFile\",\"mt2\",\"\"); gSystem->Exit(0);" |root.exe -b -l ;
   mv $outputFilteredFile $outputFile;
fi;

####### mv $outputFile $outputFolder
secp file://$outputFile ${gfalProtocol}://t3se01.psi.ch/$outputFolder

rm $outputFile

if [[ $doSkimmingPruning == 1 ]]; then
#Normal skim
skimmingPruningCfg="${workingFolder}/skimmingPruning_${counterName}.cfg"
     cat skimmingPruning.cfg |grep -v \# | sed  "s#INPUTDIR#${outputFolder}#" |sed "s#INPUTFILTER#${counterName}_#" \
 	| sed "s#OUTPUTDIR#${outputFolder}/skimAndPrune#" | sed "s#DOPRUNING#${doPruning}#" > \$skimmingPruningCfg
./runSkimmingPruning.sh \$skimmingPruningCfg
rm \$skimmingPruningCfg

# #qcd skim
# skimmingPruningCfgQCD="${workingFolder}/skimmingPruningQCD_${counterName}.cfg"
#     cat skimmingPruningQCD.cfg |grep -v \# | sed  "s#INPUTDIR#${outputFolder}#" |sed "s#INPUTFILTER#${counterName}_#" \
# 	| sed "s#OUTPUTDIR#${outputFolder}/QCDskimAndPrune#" | sed "s#DOPRUNING#${doPruning}#" > \$skimmingPruningCfgQCD
# ./runSkimmingPruning.sh \$skimmingPruningCfgQCD
# rm \$skimmingPruningCfgQCD

# #qcd skim for monojet
# skimmingPruningCfgMonoJet="${workingFolder}/skimmingPruningMonoJet_${counterName}.cfg"
#     cat skimmingPruningMonoJet.cfg |grep -v \# | sed  "s#INPUTDIR#${outputFolder}#" |sed "s#INPUTFILTER#${counterName}_#" \
# 	| sed "s#OUTPUTDIR#${outputFolder}/QCDMonoJetSkimAndPrune#" | sed "s#DOPRUNING#${doPruning}#" > \$skimmingPruningCfgMonoJet
# ./runSkimmingPruning.sh \$skimmingPruningCfgMonoJet
# rm \$skimmingPruningCfgMonoJet

else
 echo "skipping skimming/pruning steps"
fi;

#echo "is anything left in working folder? workingFolder: " 
#ls $workingFolder


EOF

        #if you have a big file and no time to change the code to be smoother: qsub  -q short.q -l h_vmem=5g batchScript_${name}.sh;

	qsub -q all.q $scriptName;
#	qsub -q short.q  $scriptName;
#	qsub -q short.q   -l h_vmem=5g $scriptName;
	rm $scriptName;

	if (($counter < 0)); then
	    break;
	else
	counter=$[$counter+1]; #	    counter=$[$counter+1]; 
	fi;
	
	#rm $counterFile

    done;

    ### BM: re-check logic of chunks creation. 
    #rm -f /tmp/testChunks
    #for x in $jobsLogsFolder/chunkPart_${name}_*.txt; do cat $x >> /tmp/testChunks; done;    
    #echo "number files in inputList, sum of chuncks: "
    #cat $fileList |wc -l
    #cat /tmp/testChunks |wc -l


done < $listOfSamplesFile




###rm -f postProcessing_C.d postProcessing_C.so;

fi;

if [[ "$1" = "postCheck" ]]; then
    #logsErr=`cat ${jobsLogsFolder}/*.err`
    # use the following until a solution is found to remove dictionary warning messages, 
    # which appear only when running skimming/pruning:
    logsErr=`cat ${jobsLogsFolder}/*.err | grep -v "found in libCore.so  is already in libDataFormatsStdDictionaries.so" | grep -v "BTagCalibrationReader found in libCondToolsBTau.so  is already in libCondFormatsBTauObjects.so"`    
    if [[ -z "$logsErr" ]]; then
	echo "there were no errors. Zipping all logs and copying them to the SE"	 
	cd $jobsLogsFolder
	tar -czvf logs.tgz  *
	secp file://`pwd`/logs.tgz ${gfalProtocol}://t3se01.psi.ch/$outputFolder
	cd ..
	rm $jobsLogsFolder/*
	rmdir $jobsLogsFolder
	./doTreeProduction.sh clean 
    else
	echo "ERROR: something went wrong. Check your logs in " $jobsLogsFolder
    fi

fi





if [[ "$1" = "mergeData" ]]; then

    seString="root://t3se01.psi.ch/"

    # It assumes that the 'doTreeProduction mergeData' script is run after the 'doTreeProduction post' step.
    # If this is not the case, the 'input' variable here may need to be set properly by hand
    # no automated yet to merge the three skim flavours... (un)comment out as necessary
    # input="${outputFolder}/skimAndPrune/"
    #input="${outputFolder}/"
    input="${outputFolder}/QCDskimAndPrune/"
    #input="${outputFolder}/QCDMonoJetSkimAndPrune/"


    echo "InputFolder for mergeData script: " $input

    # create tmp folder in the scratch area before copying everything to the SE
    tmpOutputDir="/scratch/$USER/$RANDOM/"
    mkdir $tmpOutputDir
    inputFilesList="${tmpOutputDir}/fileList.txt"

    # Add other relevant strings here if you want to merge more than these 3 datasets
    #datasets="MET HTMHT JetHT"
    #datasets="MET HTMHT JetHT SingleElectron SingleMuon SinglePhoton DoubleEG DoubleMuon MuonEG"
    datasets="DoubleEG DoubleMuon HTMHT JetHT MET MuonEG SingleElectron SingleMuon SinglePhoton"

    rootFileName="merged"
    for d in $datasets; do
	rootFileName=${rootFileName}_$d
    done
    rootFileName=${rootFileName}.root
    tmpOutputFile=$tmpOutputDir/$rootFileName
    outputFile=$input/$rootFileName

    for x in $datasets; do
	prefix=$x
	for x in $input/${prefix}_*.root; do echo $seString$x >> $inputFilesList ; done;
    done
    sed -i 's#pnfs/psi.ch/cms/trivcat/##g' $inputFilesList

    echo "merged file will be: " $outputFile
    echo "gROOT->LoadMacro(\"removeDuplicates.C\"); removeDuplicatesFromChain(\"$inputFilesList\",1,\"$tmpOutputFile\"); gSystem->Exit(0);" |root.exe -b -l ;

    xrdcp -vf $tmpOutputFile root://t3dcachedb.psi.ch:1094//$input/$rootFileName
    rm $tmpOutputFile
    rm $inputFilesList
    rmdir $tmpOutputDir
fi


if [[ "$1" = "addLepSF" ]]; then
#the outputfolder of postprocessing is the input file for the adding of the scale factors
    ./doLepSF.sh $outputFolder 
fi

if [[ "$1" = "addAllSF" ]]; then
#the outputfolder of postprocessing is the input file for the adding of the scale factors
    ./doAllSF.sh $outputFolder 
fi

if [[ "$1" = "addBtag" ]]; then
    ./doBTagSF.sh
fi

if [[ "$1" = "addISR" ]]; then
    ./doISRSF.sh
fi

if [[ "$1" = "clean" ]]; then
    rm -f postProcessing_C*;
    rm -f goodrun_cc*;
    rm -f goodrunClass_cc*;
    rm -f BTagCalibrationStandalone_cc*;
fi


