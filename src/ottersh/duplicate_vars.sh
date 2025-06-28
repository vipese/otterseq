#!/bin/bash

# Parse named arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --bfile)
        bfile="$2"
        shift # past argument
        shift # past value
        ;;
        *)    # unknown option
        echo "Unknown option $key"
        exit 1
        ;;
    esac
done

# Duplicated by variant and coordinates
bin/plink --bfile "${bfile}" --list-duplicate-vars suppress-first \
    --allow-no-sex --out "${bfile}" --silent 
