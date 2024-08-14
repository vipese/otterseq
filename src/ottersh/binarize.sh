#!/bin/bash 

# Parse named arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --file)
        file="$2"
        shift # past argument
        shift # past value
        ;;
        --outpath)
        outpath="$2"
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
filename=$(basename "$file")
outfile=${outpath}/${filename}
plink2 --pedmap "${file}" --make-pgen  --out "${outfile}" --silent

# Print 
echo "Files binarized"


