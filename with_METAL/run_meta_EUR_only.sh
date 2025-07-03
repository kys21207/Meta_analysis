#!/bin/bash

# Directory containing the tsv.gz files
INPUT_DIR="/mnt/project/publically_available_summary_statistics/meta/"
DATA_DIR="/opt/notebooks/input_data"
# Directory for output files (using the same directory)
OUTPUT_DIR="/opt/notebooks/meta_eur"

chmod +x ./generic-metal/metal

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Function to create meta-analysis configuration file
create_meta_config() {
    local trait_name=$1
    local config_file="$OUTPUT_DIR/${trait_name}_meta_config.txt"
    
    cat > "$config_file" << EOF
# This file includes a series of comments. Each comment is marked 
# by a # character as the first character in the line.
#
# This is a comment!

# Meta-analysis weighted by standard error does not work well
# when different studies used very different transformations.
# In this case, some attempt was made to use similar trait
# transformation and you can request a standard error based
# analysis by uncommenting the following line:
 SCHEME   STDERR
 AVERAGEFREQ ON
 
# Describe and process the DGI input files
MARKER   SNP
WEIGHT   ukbb_n
ALLELE   ALT REF 
FREQ     ukbb_af_alt
EFFECT   ukbb_beta
STDERR   ukbb_sebeta
PVAL     ukbb_pval

PROCESS ${DATA_DIR}/${trait_name}_ukb_EUR.tsv.gz 

# Describe and process the FUSION input files
MARKER   SNP
WEIGHT   MVP_n
ALLELE   ALT REF 
FREQ     MVP_EUR_af_alt
EFFECT   MVP_EUR_beta
STDERR   MVP_EUR_sebeta
PVAL     MVP_EUR_pval

PROCESS ${DATA_DIR}/${trait_name}_mvp_EUR.tsv.gz

# Describe and process the SardiNIA input files
MARKER   SNP
WEIGHT   fg_n
ALLELE   ALT REF 
FREQ     fg_af_alt
EFFECT   fg_beta
STDERR   fg_sebeta
PVAL     fg_pval

PROCESS ${DATA_DIR}/${trait_name}_fg_EUR.tsv.gz

OUTFILE ${OUTPUT_DIR}/${trait_name}_meta_ .tsv

# Execute meta-analysis
ANALYZE HETEROGENEITY
EOF

    echo "Created meta-analysis configuration file: $config_file"
}

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
#file="/mnt/project/publically_available_summary_statistics/meta/G6_MS_meta_out.tsv.gz"
    COUNTER=$((COUNTER + 1))
    filename=$(basename "$file")
    echo "[$COUNTER/$FILE_COUNT] Processing $filename"
    
    # Extract trait name from filename (assuming format is trait_name_something.tsv.gz)
    # This extraction method may need adjustment based on your actual filename format
    trait_name=$(basename "$filename" _meta_out.tsv.gz)
    
    # Run the R script on this file
#    Rscript ./meta_analysis/pullout_each_EUR_data.R "$file" "$DATA_DIR"
    
    # Check if successful
    if [ $? -eq 0 ]; then
        echo "✓ Successfully processed $filename"
        
        # Generate meta-analysis configuration file
        create_meta_config "$trait_name"
    else
        echo "✗ Error processing $filename"
    fi
    echo "Run METAL"
    ./generic-metal/metal "$OUTPUT_DIR/${trait_name}_meta_config.txt"
done

echo "All files processed. Results saved to $OUTPUT_DIR"
