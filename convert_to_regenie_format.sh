#!/bin/bash

input_path="./meta_eur"
endpoints_file="/mnt/project/publically_available_summary_statistics/meta/endpoints.tsv"

for input_file in ${input_path}/*.tsv; do
#input_file="./meta_eur/G6_MS_meta_1.tsv"
phenocode=$(basename "$input_file" _meta_1.tsv)
output_file="${phenocode}_regenie.tsv.gz"

# Calculate N based on phenocode from endpoints.tsv
N=$(awk -F'\t' -v pc="$phenocode" '
    NR > 1 && $2 == pc {
        N = $6 + $7 + $8 + $9 + $12 + $13
        print N
        exit
    }
' "$endpoints_file")

if [ -z "$N" ]; then
    echo "Error: Phenocode $phenocode not found in $endpoints_file"
    exit 1
fi

echo "Using N = $N for phenocode $phenocode"

# Process the input file and generate the regenie format
awk -v N="$N" '
    BEGIN {
        OFS="\t"
        print "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "CHISQ", "LOG10P", "INFO"
    }
    
    # Skip header row
    NR == 1 { next }
    
    {
        # Check if HetDF starts with 2 (i.e., HetDF==2)
        if (substr($12, 1, 1) == "2") {
            # Extract CHROM and GENPOS from MarkerName
            split($1, marker_parts, ":")
            chrom = marker_parts[1]
            genpos = marker_parts[2]
            
            # Map columns to regenie format
            id = $1
            allele0 = $2  # Allele1
            allele1 = $3  # Allele2
            a1freq = $4  # Freq
            beta = $6    # Effect
            se = $7      # StdErr
            pvalue = $8  # P-value
            
            # Calculate LOG10P
            log10p = -log(pvalue) / log(10)
            
            # Calculate CHISQ (beta/se)^2
            chisq = (beta/se)^2
            
            # Set A1FREQ and INFO to NA
            info = "NA"
            
            print chrom, genpos, id, allele0, allele1, a1freq, N, beta, se, chisq, log10p, info
        }
    }
' "$input_file" | gzip > "$output_file"

done
echo "Conversion complete. Results saved to $output_file"
