#!/bin/bash                                                                    
# Initialize variables for command line arguments
MAX_EVENTS=""
PEDESTALS=false
FINESIGMA=false
ONLINE=false
HELP=false
# Function to show usage
usage() {
	echo "Usage: $0 [--max-events=number] [--pedestals] [--pedestals-fine] [--online] [--help]"
	exit 1
}
# Parse command line arguments
for arg in "$@"; do
	case $arg in
		--max-events=*)
			MAX_EVENTS="${arg#*=}"
			shift # Remove argument from processing
			;;
		--pedestals)
			PEDESTALS=true
			shift # Remove argument from processing
			;;
		--fine-sigma)
			FINESIGMA=true
			shift # Remove argument from processing
			;;
		--online)
			ONLINE=true
			shift # Remove argument from processing
			;;
		--help)
			HELP=true
			shift # Remove argument from processing
			;;
		*)
			usage
			;;
	esac
done
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/u/land/fake_cvmfs/10/cvmfs_fairroot_v18.8.0_fs_nov22p1_extra/ucesb/hbook
#source /home/land/root/build/bin/thisroot.sh
cd /u/land/r3broot/202402_s091_s118/R3BParams_s091_s118/macros/exp/online/foot/event_display_dev/src
clear ; clear
make clean
make
# Forming the main command with the parsed options
CMD="./event_display"
if [ -n "$MAX_EVENTS" ]; then
	CMD+=" --max-events=$MAX_EVENTS"
fi

if [ "$PEDESTALS" = true ]; then
	CMD+=" --pedestals"
fi

if [ "$FINESIGMA" = true ]; then
	CMD+=" --fine-sigma"
fi

if [ "$ONLINE" = true ]; then
	CMD+=" --online"
fi

if [ "$HELP" = true ]; then
	CMD+=" --help"
fi

# Execute the command
echo
echo "Executing: $CMD"
eval $CMD
