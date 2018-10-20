#!/bin/bash

## CRISPR-Cas9 large deletions detector (call v1.0.0)
## Author: Rui Chen
################################################################# PARSE ARGUMENTS ##########################################################

## input
threads_default=False 
while getopts "t:f:s:b:p:o:a:g:h" OPT; do
{
    case ${OPT} in
        t)
            threads="${OPTARG}"
            threads_default=True;; 
        f) 
            reference="${OPTARG}";;
        s)
            target="${OPTARG}";;
        b)
            bam="${OPTARG}";;
        p)
            prefix="${OPTARG}";;
        o)
            output="${OPTARG}";;
        a)
            parents="${OPTARG}";;
        g)
            offspring="${OPTARG}";;
        h)
            echo "Usage: clld [options] -f <reference> -s <bed> -b <bam> -p <prefix> -o <output_dir> -a <parents> -g <offspring>" && exit 
    esac
} done

## multi input handle
bam_array=(${bam//,/ })
number=${#bam_array[@]}

## load the config file 
source ./clld_io.config
source ./clld_software.config

## check the directory
source ./clld_check_dir.sh

## option
if [ ${threads_default} = "False" ]
    then threads=24
fi

################################################################# FUNCTIONS ##########################################################

call_manta()
{
    check_file_exist "${manta_path}/runWorkflow.py"

    ## combine the multi input
    para=""
    while [[ ${number} -gt 0 ]]
    do
        number=$(expr ${number} - 1)
        single_bam=${bam_array[${number}]}
        para="${para}--bam ${single_bam} "
    done
    
    ## set the workflow
    configManta.py \
    ${para} \
    --referenceFasta "${reference}" --runDir "${manta_path}" --generateEvidenceBam
    
    ## run the workflow
    echo "run the manta to call structure variants\n"
    time "${manta_path}/runWorkflow.py" -m local -j ${threads}
    gzip -d "${manta_path}/results/variants/diploidSV.vcf.gz"
    cp "${manta_path}/results/variants/diploidSV.vcf" "${result_path}/${prefix}.diploidSV.vcf"
    cp -r "${manta_path}/results/evidence/" "${result_path}"
    time python clld.py "${result_path}/${prefix}" "${result_path}/${prefix}.diploidSV.vcf" "${target}" "${parents}" "${offspring}"
}

################################################################# MAIN #####################################################################

## call large deletion
call_manta > "${log_path}/${prefix}.manta.log" 2>&1 &
