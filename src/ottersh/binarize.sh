#!/bin/bash 

# Parse named arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --files)
        FILES="$2"
        shift # past argument
        shift # past value
        ;;
        --outpath)
        OUTPATH="$2"
        shift # past argument
        shift # past value
        ;;
        *)    # unknown option
        echo "Unknown option $key"
        exit 1
        ;;
    esac
done

# For each file, binarize
for file in "${FILES[@]}" ; do
    filename=$(basename "$file")
    outfile=${OUTPATH}/${filename}
    plink2 --pedmap "${file}" --make-bed --out "${outfile}" --silent
done 

# Print 
echo "Files binarized"


