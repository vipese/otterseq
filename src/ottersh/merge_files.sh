#!/bin/bash 

# Default value
prefix="merged_snps"

# Parse named arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --merge-file)
        merge_file="$2"
        shift # past argument
        shift # past value
        ;;
        --outpath)
        outpath="$2"
        shift # past argument
        shift # past value
        ;;
        --prefix)
        prefix="$2"
        shift # past argument
        shift # past value
        ;;
        *)    # unknown option
        echo "Unknown option $key"
        exit 1
        ;;
    esac
done

outfile="${outpath}"/"${prefix}"
plink --merge-list "${merge_file}" --make-bed --out "${outfile}" --silent