#!/bin/bash

## Nicolas Servant
## This script is used to simulate Cancer HiC data
## Starting from a diploid map, it simulated CNV following the description in ${cnvbed} files
##

## PATH
R_PATH=$(which R)

function usage
{
    echo -e "usage : simulateCNV.sh -i INPUT -b INPUT_BED -c CNV_FILE -o OUTPUT [-h]"
    echo -e "Use option -h|--help for more information"
    exit 1
}

## Get options
while getopts ":i:b:c:o:h" opt; do
    case $opt in
	i) DIPLOID_MATRIX=$OPTARG ;;
	b) DIPLOID_BED=$OPTARG ;;
	c) CNV_BED=$OPTARG ;;
	o) OUT_PATH=$OPTARG ;;
	h) usage ;;
	\?)
	    echo "Invalid option: -$OPTARG" >&2
	    usage
	    ;;
    esac
done

#DIPLOID_MATRIX=./data/rao/IMR90_CCL186_${RES}.matrix
#DIPLOID_BED=./data/rao/IMR90_CCL186_${RES}_abs.bed
#OUT_PATH=./results/
#CNV_BED=cnv1.bed

function run_simulation
{
    local outfile=${4}/$(basename $1 | sed -e 's/mat\(rix\)*//')_$(basename $3 | sed -e 's/.bed//').matrix
    
    ## Run simulation
    echo "Simulate CNVs from ${3} file ..."
    cmd="${R_PATH} --no-save --no-restore CMD BATCH \"--args input='$1' xgi='$2' cnv='${3}' output_prefix='${outfile}'\" simulateCNV.R simulateCNV_${4}.Rout"
    echo $cmd
    eval $cmd
}

if [ $# == 0 ]; then
    usage
fi

if [[ ! -e ${DIPLOID_MATRIX} || ! -e ${DIPLOID_BED} || ! -e ${CNV_BED} ]]; then
    echo "Input(s) not found. Exit."
    exit 1
fi

run_simulation ${DIPLOID_MATRIX} ${DIPLOID_BED} ${CNV_BED} ${OUT_PATH}
exit 0
