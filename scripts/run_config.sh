#!/bin/bash
#                 ______ __  __ _____ _      _____ 
#                |  ____|  \/  |_   _| |    |_   _|
#                | |__  | \  / | | | | |      | |  
#                |  __| | |\/| | | | | |      | | 
#                | |____| |  | |_| |_| |____ _| |_ 
#                |______|_|  |_|_____|______|_____|
#
#   irace run configurator for EMILI version 1.0
set -e

function check_file 
{
 if [ -e $1 ]
 then
    printf "OK\n"
 else
   printf "ERROR\n file $1 does not exists! \n"
   exit 1
 fi
 return 0
}

if [ $# -lt 6 ]
then


echo "                 ______ __  __ _____ _      _____ "
echo "                |  ____|  \/  |_   _| |    |_   _|"
echo "                | |__  | \  / | | | | |      | |  "
echo "                |  __| | |\/| | | | | |      | |  "
echo "                | |____| |  | |_| |_| |____ _| |_ "
echo "                |______|_|  |_|_____|______|_____|"
echo "                                                  "
echo " irace run configurator for EMILI version 1.0     "
echo "                                                  "
echo "usage: $0 run_path grammar_path runtime irace_budget training_set_path description [options]"
echo " "
echo "  run_path            : name that will be used to create the directory for the run"
echo "  grammar_path        : path to the grammar file"
echo "  runtime             : ro value you want to use"
echo "  irace_budget        : budget for irace"
echo "  training_set_path   : path to the dir that contains the training instances"
echo "  description         : description of the run between \" \" "
echo "options: "
echo "  --executbable path  : path to the EMILI executable (if not specified the default path will be used) "
echo "  --recursion level   : number of times a recursive rule is allowed to be expanded"
echo "  --candidates path   : the final candidates of the irace run at path " 
echo "                        will be used as intial candidates for the current run"
echo ""
else
RUNDIR=$1
IRACE_VER=`irace --version | grep -i version | cut -d':' -f2` 
GRAMMAR=$2
RUNTIME=$3
BUDGET=$4
TRAININGSET=$5
COMMENT=$6
G2CODE='/lustre/home/fpagnozzi/bin/grammar2code'
EMILI_BIN='/home/fpagnozzi/emili/build/EMILI'
TEMPLATE_DIR='/home/fpagnozzi/noccirace/templates'
DEFAULT_HOOK='hook_run'
DEFAULT_TARGET='target-runner'
DEFAULT_TUNE='/home/fpagnozzi/noccirace/templates/tune-main-cluster-mpi-em'
RECURSION=3
FIRST_LEVEL=""
shift 6
#if [ $# -eq 7 ]
#then
#  EMILI_BIN=$7
#fi
PARAMS=
while [ $# -gt 0 ]; do
    case "$1" in
        --executable) shift; EMILI_BIN="$1"; shift;;
        --recursion) shift; RECURSION="$1"; shift;;
        --candidates) shift; FIRST_LEVEL="$1"; shift;;
        *) PARAMS="$PARAMS $1"; shift;;# terminate case
  esac
done

#IRACEHOOK="target-runner"
IRACEHOOK=$DEFAULT_TARGET
IRACE_PARAMS="scenario.txt"
if [ `echo $IRACE_VER | cut -d'.' -f1` -lt 2 ]
then
#	IRACEHOOK="hook-run"
	IRACEHOOK=$DEFAULT_HOOK
	IRACE_PARAMS="tune-conf"
fi
EMILI_VER=`$EMILI_BIN | grep commit | cut -d':' -f 2`

# check that rundir exists
printf "checking destination path..."
if [ ! -d $RUNDIR ]
then
   printf "creating dir..."
   mkdir $RUNDIR
fi

check_file $RUNDIR

# copy iracehook template to rundir
printf "checking $IRACEHOOK..."
check_file "$TEMPLATE_DIR/$IRACEHOOK"

cp "$TEMPLATE_DIR/$IRACEHOOK" $RUNDIR

# copy scenario.txt if present
printf "checking $IRACE_PARAMS..."
if [ -e $IRACE_PARAMS ]
then
  cp $IRACE_PARAMS $RUNDIR
fi
printf "OK\n"

printf "checking tune-main-cluster-mpi-em..."
#check_file tune-main-cluster-mpi-em
check_file $DEFAULT_TUNE

#cp tune-main-cluster-mpi-em $RUNDIR
cp $DEFAULT_TUNE $RUNDIR

#move to rundir
pushd $RUNDIR > /dev/null

FULL_DIR=`pwd`
# create parameters file
printf "creating parameters file with $RECURSION levels of recursion..."
$G2CODE $GRAMMAR -d $RECURSION -p parameters.txt 1> /dev/null 2> /dev/null
printf "\nchecking file..."
check_file parameters.txt

# modify iracehook
#sed -i "s/^GRAMMAR_DIR=.*/GRAMMAR_DIR=$GRAMMAR/g" $IRACEHOOK
printf "modifying $IRACEHOOK..."
sed -i "s,^GRAMMAR_DIR=.*,GRAMMAR_DIR='$GRAMMAR',g" $IRACEHOOK
sed -i "s/^TIMER_SETTING=.*/TIMER_SETTING='-ro $RUNTIME'/g" $IRACEHOOK
sed -i "s,^ERR_OUT=.*,ERR_OUT='$FULL_DIR/ERROR',g" $IRACEHOOK
sed -i "s,^EMILI_BIN=.*,EMILI_BIN='$EMILI_BIN',g" $IRACEHOOK
printf "OK\n"
printf "Using the executable $EMILI_BIN\n"

# create scenario.txt
if [ ! -e $IRACE_PARAMS ]
then
	printf "creating $IRACE_PARAMS..."
	printf "trainInstancesDir <- \"$TRAININGSET\"\n maxExperiments <- $BUDGET\n" > $IRACE_PARAMS
	printf "OK\n"
else
	printf "modifying $IRACE_PARAMS..."
	sed -i "s,^trainInstancesDir.*,trainInstancesDir <- \"$TRAININGSET\",g" $IRACE_PARAMS
	sed -i "s/^maxExperiments.*/maxExperiments <- $BUDGET/g" $IRACE_PARAMS
	printf "OK\n"
fi

# setup second level run

if [ ! -z $FIRST_LEVEL ]
then
  printf "creating configurations file from previous irace execution..."
  iu_build_canditates_file irace2 $FIRST_LEVEL/irun/irace.Rdata > initconfs
  check_file initconfs
  #printf "OK\n"
  printf "updating scenario file..."
  echo 'configurationsFile="initconfs"' >> scenario.txt
  printf "OK\n"
fi

# create rundir
printf "creating execution dir..."
mkdir irun
printf "OK\n"

#print description
printf "creating description file..."
printf "Run description\n IraceVersion: $IRACE_VER\n Grammar: $GRAMMAR\n EMILI: $EMILI_VER\n Runtime: $RUNTIME\n Comment: \"$COMMENT\"\n\n" > runDescription
printf "OK\n"
popd > /dev/null
printf "\nDONE\n"
fi
