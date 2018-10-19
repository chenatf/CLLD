#!/bin/bash

## CRISPR-Cas9 large deletions detector (mapping v1.0.0)
## Author: Rui Chen
################################################################# PARSE ARGUMENTS ##########################################################

## handle input
threads_default=False
while getopts "t:r:f:1:2:p:o:h" OPT; do
{
    case ${OPT} in
        t)
            threads="${OPTARG}"
            threads_default=True;;
        r)
            flag="${OPTARG}";; 
        f) 
            reference="${OPTARG}";;
        1)
            rawdata_r1="${OPTARG}";;
        2)
            rawdata_r2="${OPTARG}";;
        p)
            prefix="${OPTARG}";;
        o)
            output="${OPTARG}";;
        h)
            echo "Usage: csvd_mapping -r <flag> -f <reference> -1 <pair1> -2 <pair2> -p <prefix> -o <output_dir>" && exit 
    esac
} done

## multi input handle
rawdata_r1=(${rawdata_r1//,/ })
rawdata_r2=(${rawdata_r2//,/ })
number=${#rawdata_r1[@]}

## load the config file
source ./clld_io.config
source ./clld_software.config

## check the directory
source ./clld_check_dir.sh

## options
fastp_threads=8
mosdepth_threads=4
if [ ${threads_default} = "False" ]
    then threads=12
fi

################################################################# FUNCTIONS ##########################################################

mapping_ratio()
{
    all_result=$(samtools view -c "${1}")
    no_alighnment=$(samtools view -c -f 4 "${1}")
    no_primany=$(samtools view -c -f 2304 "${1}")
    all_read=$(expr ${all_result} - ${no_primany})
    mapping_reads=$(expr ${all_read} - ${no_alighnment})
    ratio=$(expr 100 \* ${mapping_reads} / ${all_read})
    echo "All reads are:${all_read}  Mapping reads:${mapping_reads}  mapping ratio:${ratio}%"
}

mapping()
{
    while [[ ${number} -gt 0 ]]
    do
        ## get the data number
        number=$(expr ${number} - 1)
	data_r1=${rawdata_r1[${number}]}
        data_r2=${rawdata_r2[${number}]}

        ## quality control
        echo "step1: run fastp to imporve the quality of ${data_r1}, ${data_r2}\n"
        time fastp -w ${fastp_threads} -c \
        -h "${result_path}/${prefix}_${number}.html" -j "${fastp_path}/${prefix}_${number}.json" \
        -i "${data_r1}" -o "${fastp_path}/${prefix}_${number}_r1.fq.gz" \
        -I "${data_r2}" -O "${fastp_path}/${prefix}_${number}_r2.fq.gz"

        ## mapping
        echo "step2: run bwa and samtools to mapping the reads to the reference\n" 
        time ${speedseq} align -o "${speedseq_path}/${prefix}_${number}_" -t ${threads} \
        -R "${flag}" "${reference}" \
        "${fastp_path}/${prefix}_${number}_r1.fq.gz" "${fastp_path}/${prefix}_${number}_r2.fq.gz"

        ## calculate mapping ratio
        echo "step3: run samtools to calculate mapping ratio\n"
        time mapping_ratio "${speedseq_path}/${prefix}_${number}_.bam"
    done

    ## merge batches
    echo "step4: run samtools to merge batches and calculate the total mapping ratio\n"
    merge_bam=($(ls "${speedseq_path}/${prefix}"_*_markdup_.bam))
    time samtools merge -f "${samtools_path}/${prefix}.bam" "${merge_bam[@]}"
    time samtools index "${samtools_path}/${prefix}.bam"
    time mapping_ratio "${samtools_path}/${prefix}.bam" && rm "${merge_bam[@]}"

    ## coverage analysis
    echo "step5: run mosdepth to calucate coverage\n"
    time mosdepth -t ${mosdepth_threads} "${mosdepth_path}/${prefix}" "${samtools_path}/${prefix}.bam"
    time python3 "${plot-dist}" -o "${result_path}/${prefix}.dist.html" "${prefix}" \
    "${mosdepth_path}/${prefix}.mosdepth.global.dist.txt"
}

################################################################# MAIN #####################################################################

## mapping 
mapping > "${log_path}/${prefix}.mapping.log" 2>&1 &
