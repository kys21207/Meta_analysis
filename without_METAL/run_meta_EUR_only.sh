#!/bin/bash

# Directory containing the tsv.gz files
INPUT_DIR="/mnt/project/project/publically_available_summary_statistics/meta/"

# Directory for output files (using the same directory)
OUTPUT_DIR="/mnt/project/project/publically_available_summary_statistics/meta/"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Make R script executable
chmod +x run_meta_analysis.R

# Find all tsv.gz files and process them
echo "Starting meta-analysis of all tsv.gz files in $INPUT_DIR"
echo "----------------------------------------------------"

# Get count of files
FILE_COUNT=$(find "$INPUT_DIR" -name "*.tsv.gz" | wc -l)
echo "Found $FILE_COUNT files to process"

# Initialize counter
COUNTER=0

# Process each file
find "$INPUT_DIR" -name "*.tsv.gz" | while read -r file; do
    COUNTER=$((COUNTER + 1))
    filename=$(basename "$file")
    echo "[$COUNTER/$FILE_COUNT] Processing $filename"
    
    # Run the R script on this file
    ./run_meta_analysis.R "$file"
    
    # Check if successful
    if [ $? -eq 0 ]; then
        echo "✓ Successfully processed $filename"
    else
        echo "✗ Error processing $filename"
    fi
    echo "----------------------------------------------------"
done

echo "All files processed. Results saved to $OUTPUT_DIR"
