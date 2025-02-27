#!/bin/bash

# Starting with the HicContactCaller outputt file, this script finds adaptive hubs. To make it work you need:
# a) HicContactCaller output file
# b) List of DE genes names (i.e evm.model.Herato1201.1", it needs to be 1 column, originating from RNA-seq DEseq)
# c) Selective sweep bed file (Sweepfinder2 output) This file is edited to include only the top 0.5% sweeps with highest CLR values, and has 2kb up and downstream (ie a 4kb range in each position).

# With these files this script will print a file based on the hicContactCaller output file that but that only contains those contacts that have a pvalue < 0.05, and whose genes match those that are differentially
# expressed (as per the DE genes file), and whos chipPeak is under strong selection (as per selective sweep file). Then the last filter makes sure only those CREs who are functionally connected to two or more DE genes are printed. 
# This last filter is satisfy the last criteria of the Adaptive Hubs Definintion of J. Lewis in Science Advances (2020)

echo "Processing adaptive hubs..."

########################  

# Ensure correct paths
cd /scratch/Counterman_Lab_Shared/Jamie/contact_caller_hi-c/by_chunks

# File paths
file1="HW_FW_Combined_0.05_sig_DE_genes_clean_sort.csv" # This is the RNAseq DEseq file, its based on differential expression between postman and rayed, done separetely for HW and FW and then combined afterwards.
file2="Output_chunks_combined.txt" #This is the hic Contact Caller output file that i combined after processing it in chunks (to bypass the everpresent memory issues in invoking R within the HPC)
file3="all.amalfreda.topSS.bed" # Selective sweep file

output_file="output_HiC_DE_p0.05_SS_all.amalfreda.txt"
sorted_output_file="sorted_${output_file}"
output_file_unique_loci="adaptive_hubs_all.amalfreda.txt"
output_file_unique_loci_bed="adaptive_hubs_all.amalfreda.bed" # Additional BED format output

#########################

# Use awk to count, then filter based on multiple conditions including gene count per loci, considering only statistically significant associations
awk '
    BEGIN { FS = OFS = "\t" }
    # Load gene names from file1 into an array for easy lookup
    NR == FNR { genes[$1]; next }
    # Load selective sweep ranges from file3
    FNR < NR && FILENAME == ARGV[2] {
        sweeps[$1] = (sweeps[$1] ? sweeps[$1] ";" : "") $2 "-" $3;
        next
    }
    # Count genes per loci in file2, considering only those that pass the p-value filter of 0.05, ofc this figure can be changed below but also must be changed in the next section.
    FNR < NR && FILENAME == ARGV[3] {
        gene_name = gensub(/[";]/, "", "g", $10);
        loci = $1 FS $2;
        if (gene_name in genes && $6 < 0.05) {
            loci_gene_count[loci]++;
            loci_genes[loci][gene_name]; # Track unique genes per loci that pass the p-value filter
        }
    }
    # Second pass through file2 to print lines that meet all criteria
    FILENAME == ARGV[3] {
        gene_name = gensub(/[";]/, "", "g", $10);
        loci = $1 FS $2;
        # Ensure gene is in file1, p-value < 0.05, and loci has >= 2 unique genes that passed the p-value filter
        if (gene_name in genes && $6 < 0.05 && loci in loci_gene_count) {
            count_unique_genes = 0;
            for (g in loci_genes[loci]) count_unique_genes++;
            if (count_unique_genes >= 2) {
                # Check if loci falls within any selective sweep range
                split(sweeps[$1], ranges, ";");
                for (i in ranges) {
                    split(ranges[i], bounds, "-");
                    if ($2 >= bounds[1] && $2 <= bounds[2]) {
                        print;
                        break;
                    }
                }
            }
        }
    }
' "$file1" "$file3" "$file2" "$file2" > "$output_file"

#########################

# Sort the output file and replace the original with the sorted one
sort -k1,1 -k2,2n "$output_file" > "$sorted_output_file"
mv "$sorted_output_file" "$output_file"

# Create a second output file with unique loci (columns 1 and 2), sorted and duplicates removed.
cut -f1,2 "$output_file" | sort -k1,1 -k2,2n | uniq > "$output_file_unique_loci"

# Additionally, create a BED format file with the specified adjustments
awk 'BEGIN {OFS="\t"} {
    start = $2 - 2000; 
    end = $2 + 2000;
    if (start < 1) start = 1; 
    print $1, start, end
}' "$output_file_unique_loci" | sort -k1,1 -k2,2n | uniq > "$output_file_unique_loci_bed"

echo "Processing complete. Output files generated:"
echo "- Main output: $output_file"
echo "- Unique loci output: $output_file_unique_loci"
echo "- BED format output: $output_file_unique_loci_bed"
