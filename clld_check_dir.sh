#!/bin/bash

check_dir_exist()
{
    ## check the dir exist
    if [ ! -d "$1" ]
        then  mkdir "$1"
    fi
}

check_file_exist()
{
    ## check the file exist
    if [ -r "$1" ]
        then  rm "${1}"
    fi
}

## check all the path exist
check_dir_exist "${metadata_path}"
check_dir_exist "${result_path}"
check_dir_exist "${log_path}"
check_dir_exist "${fastp_path}" 
check_dir_exist "${samtools_path}"
check_dir_exist "${speedseq_path}"
check_dir_exist "${mosdepth_path}"
check_dir_exist "${manta_path}"