#!/bin/bash


usage() {
    echo "`basename $0` -f filelist -c configfile -n njobs(default=1) [-i]"
    echo "  -f filelist: filelist of input (e.g. etc/data/ExampleDataDiphoton8TeVSkimEOS.list)"
    echo "  -c configfile: config file (e.g. etc/config/mvaAnalysis_Legacy2013.config)"
    echo "  -i: interactive mode (by default job sent to batch system)"
    echo "  -o: output directory (ovewrite the one defined in config file)"
    echo "  -n: number of jobs per input file"
    echo "  -a: automatic setting of the number of jobs with 200000 evts max per job by default"
}

mainConfigFile=notSpecified
inFileList=notSpecified
outputDir=notSpecified
nJobs=1
interactive=0
autonjobs=0
nevtperjobauto=-1

if ! options=$( getopt -o hc:f:n:o:i:a -l help -- "$@" )
then
    # something went wrong, getopt will put out an error message for us
    exit 1
fi

eval set -- $options
while [ $# -gt 0 ]; do
    case "$1" in
	-h | --help) usage; exit 0;;
	-c) mainConfigFile=$2; shift;;
	-f) inFileList=$2; shift;;
	-n) nJobs=$2; shift;;
	-o) outputDir=$2; shift;;
	-i) interactive=1;;
	-a) autonjobs=1;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; usage >> /dev/stderr; exit 1;;
	(*)  break;;
    esac
    shift
done

###### print and check configuration
echo "================================="
echo " - config file: " $mainConfigFile
echo " - file list  : " $inFileList
echo " - njobs      : " $nJobs 
if [ $autonjobs -eq 1 ]; then
    echo " - automatic number of jobs set with max: " $nevtperjobauto " evts per job"
fi
echo " - outputDir  : " $outputDir
echo " NEWWWWWW: automatic setting of nEvtsPerJob such that:  nJobs x nEvtPerJobs > nTot"

if [ $mainConfigFile == "notSpecified" ]; then
    echo no good: $mainConfigFile
    usage; 
    exit 0;
elif [ ! -f $mainConfigFile ]; then
    echo " config file does not exist... bailout"
    exit 0;
fi

if [ $inFileList == "notSpecified" ]; then
    usage; 
    exit 0;
elif [ ! -f $inFileList ]; then
    echo " input filelist does not exist... bailout"
    exit 0;
fi
echo "================================="

#exit 0;

optionsub="-q 8nh"
#optionsub="-q 1nh"


dirafs=`pwd`
dirscript=${dirafs}/tmpBatchOut/
mkdir -p $dirscript

###### obtain grid certificate to run over T2
source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh
voms-proxy-init
userid=`echo $KRB5CCNAME | awk -F"_" '{print $2}'`
pathToGridCert=/tmp/x509up_u${userid}
gridCert=`basename $pathToGridCert`
cp -f $pathToGridCert ${dirscript}/$gridCert
pathToGridCert=${dirscript}/$gridCert
echo $pathToGridCert

config() {
    file=$1
    script=${2}
    ijob=$3
    config=$4
    nEvts=-1

    outputDirOpt=""
    if [ $outputDir == "notSpecified" ]; then
	outputDirOpt=""
    else
	outputDirOpt=" -o "$outputDir
    fi
    
    if [ $# -ge 5 ]; then
	nEvts=$5
    fi
    
    exe="xAnaExe ${config} --iJob $ijob --inputFile $file --nEvtsPerJob $nEvts"$outputDirOpt
    echo "$exe"

    
    cat > $script<<EOF 
#!/bin/bash
cd $dirafs
source etc/scripts/setup_xAna.sh
export X509_USER_PROXY=$pathToGridCert
echo $exe
$exe
EOF
    chmod +x $script
}



configfile=$mainConfigFile
nJobsTot=$nJobs
nEvtsPerJob=-1

let "nJobs=nJobs-1"
for f in $( cat $inFileList ); do
    #JM add automatic splitting 
    nJobsAuto=0
    if [ $autonjobs -eq 1 ]; then
	filelog=/tmp/${USER}/nEvtsFile.log
	rm -f $filelog
	getNumberOfEvents.py $f  > $filelog
	cat $filelog
	nEvtsTot=`cat $filelog | grep nEvtsInTree | awk '{print $2}'`
	nJobsAuto=$nEvtsTot/$nevtperjobauto
	let "nJobsAuto=nJobsAuto+1"
	nJobsTot=$nJobsAuto
	let "nJobs=nJobsAuto-1"
	echo  "nJobsAuto: $nJobsAuto nEvtsTot: $nEvtsTot nJobs:$nJobs"
    fi	

    if [ $nJobsTot -eq 1 ]; then
	nEvtsPerJob=-1
    else
	filelog=/tmp/${USER}/nEvtsFile.log
	rm -f $filelog
	getNumberOfEvents.py $f $nJobsTot > $filelog
	cat $filelog
	nEvtsInFile=`cat $filelog | grep nEvtsInTree | awk '{print $2}'`
	nEvtsPerJob=`cat $filelog | grep nEvtsPerJob | awk '{print $2}'`
    fi

    for ijob in  $( seq 0 1 $nJobs ); do	
	filename=`basename $f`
	script=script$$
	script=${script}_${filename/root/sh}_${ijob}
	cd $dirscript/
	config ${f} $script $ijob  $configfile $nEvtsPerJob
	echo "-> created script: "$script
	if [ $interactive -eq 1 ]; then
	    source $script
	else
	    bsub $optionsub $script 
	fi
	cd -
    done
done


