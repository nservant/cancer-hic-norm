#!/bin/bash

## Nicolas Servant
## This script is used to simulate Cancer HiC data
## Starting from a diploid map, it simulated CNV following the description in ${cnvbed} files
##

## PATH
R_PATH=$(which R)

function usage
{
    echo -e "usage : segment_hic_data.sh -i INPUT -b INPUT_BED -s BKP_STRINGENCY -o ODIR [-h]"
    echo -e "BKP_STRINGENCY must be normal, low, high"
    echo -e "Use option -h|--help for more information"
    exit 1
}

## Get options
while getopts ":i:b:c:o:h" opt; do
    case $opt in
	i) INPUT_MATRIX=$OPTARG ;;
	b) INPUT_BED=$OPTARG ;;
	s) BKP_STRINGENCY=$OPTARG ;;
	o) OUT_PATH=$OPTARG ;;
	h) usage ;;
	\?)
	    echo "Invalid option: -$OPTARG" >&2
	    usage
	    ;;
    esac
done


function annotate
{
    cmd="${R_PATH}/R --no-save --no-restore CMD BATCH \"--args input='$1' xgi='$2' odir='$3'\" annotate_hicdata.R annotate_hicdata.Rout"
    echo $cmd
    eval $cmd
}

function segment
{
    local prefix=$(basename(${INPUT_MATRIX}) | sed -e 's/.RData//')
    ${R_PATH}/R --no-save --no-restore CMD BATCH "--args input.annot='$1' output_prefix='$prefix' glad.stringency='$2'"  segment_simulated_data.R segment_simulated_data.Rout
}


if [ $# == 0 ]; then
    usage
fi

if [[ ! -e ${INPUT_MATRIX} || ! -e ${INPUT_BED} ]]; then
    echo "Input(s) not found. Exit."
    exit 1
fi

annotate ${INPUT_MATRIX} ${INPUT_BED} ${OUT_PATH}

inannot=$(basename(${INPUT_MATRIX}) | sed -e 's/.mat\(rix\)*/.RData/')
if [ ! -e $inannot ]; then
    echo -e "Annotated file not found. Exit."
fi

segment ${inannot} ${BKP_STRINGENCY} && echo "Segmentation run with success ..."

exit 0