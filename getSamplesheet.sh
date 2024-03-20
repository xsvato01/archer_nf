#!/bin/bash
# Download the CSV file from a URL
wget "https://docs.google.com/spreadsheets/d/1WOktQDMH13d_zr0g8be1BqysMd7TIyy4LrzPywfvjhA/export?format=csv" -O "samplesheet.csv"
current_dir=$(pwd)

# Convert the CSV file to JSON using a Python script
python $current_dir/project/xsvato01/TP53_nf/scripts/CsvToJson.py $current_dir/samplesheet.csv $current_dir/samplesheet.json
#rm samplesheet.csv