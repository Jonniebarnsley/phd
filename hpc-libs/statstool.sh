#!/bin/bash

usage() { echo "Usage: $0 [-m <mask>] [-b <basin id>] <stats exec> <plotfile>" 1>&2; exit 1; }

# handle options
while getopts ":m:b:" option; do
    case "${option}" in
        m)  # specify mask
            MASK=${OPTARG}
            ;;
        b)  # specify basin number
            BASIN_ID=${OPTARG}
            ;;
        *)  # invalid option
            usage
            ;;
    esac
done
shift $((OPTIND-1))

# check args
if [ "$#" -ne 2 ]; then
    usage
fi

STATS="$1"
FILE="$2"

# constants
RHO_ICE=918.0 # density of ice
RHO_OCEAN=1028.0 # density of seawater
GRAVITY=9.81 # gravitational constant
SEA_LEVEL=0.0

# now apply the stats tool to the plotfile

# if no mask
if [ -z "${MASK}" ]; then
    $STATS $FILE $RHO_ICE $RHO_OCEAN $GRAVITY $SEA_LEVEL
    exit
fi

# if mask given but incorrect filepath
if [ ! -f "$MASK" ]; then
    echo "error: $MASK does not exist" 
    exit
fi

# otherwise
$STATS $FILE $RHO_ICE $RHO_OCEAN $GRAVITY $SEA_LEVEL $MASK $BASIN_ID