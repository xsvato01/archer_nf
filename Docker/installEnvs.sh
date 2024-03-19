#!/bin/bash

# Loop through all .yml files in the current directory
for file in *.yml; do
    # Check if the file exists
    if [ -f "$file" ]; then
        # Execute the conda command to create environment from the YAML file
        mamba env create -f "$file"
    fi
done
