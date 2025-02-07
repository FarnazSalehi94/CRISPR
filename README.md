# CRISPR

## First step

### preprocessing 
```
#!/bin/bash

# Root directory where the individual folders are located
ROOT_DIR="/lizardfs/salehi/crispr"

# Loop through all subdirectories in the root directory
for DIR in "$ROOT_DIR"/*/; do
    # Enter each subdirectory
    cd "$DIR" || exit
    
    # Find the .fa.gz file
    FILE=$(find . -name "*.fa.gz" -print -quit)
    
    if [[ -n "$FILE" ]]; then
        # Get the base name without extension
        BASENAME=$(basename "$FILE" .fa.gz)
        
        # Gunzip the file while keeping the original and index it
        gunzip -c "$FILE" > "$BASENAME.fa"
        samtools faidx "$BASENAME.fa"
    else
        echo "No .fa.gz file found in $DIR"
    fi
done
```

## Second step

### wfmash
In here, we do map a set of query sequences against a reference genome using ```wfmash``` .

```
#!/bin/bash

# Root directory where the individual folders are located
ROOT_DIR="/lizardfs/salehi/crispr"

# Reference fasta file
REFERENCE_FA="/lizardfs/salehi/crispr/grch38.fa"

# Loop through all subdirectories in the root directory
for DIR in "$ROOT_DIR"/*/; do
    cd "$DIR" || exit
            
    # Extract the basename of the current directory
    BASENAME=$(basename "$DIR")

    # Run wfmash using sbatch
    sbatch -p workers -c 48 --wrap "wfmash ${DIR}${BASENAME}.fa $REFERENCE_FA  > ${DIR}${BASENAME}-grch38.paf"
done
```

## Third step

### Hamming-fasta

```hamming-fasta``` searches for a specific sequence in a FASTA file. 
Our Target sequence is ```CTAACAGTTGCTTTTATCACNGG``` .


```
#!/bin/bash

ROOT_DIR="/lizardfs/salehi/crispr"

for DIR in "$ROOT_DIR"/*/; do
    cd "$DIR" || continue  # Skip if cd fails

    BASENAME=$(basename "$DIR")

    # Convert PAF to chain format
    paf2chain -i "$BASENAME-grch38.paf" > "$BASENAME-grch38.chain"

    # Run hamming-fasta
    hamming-fasta --fasta "$BASENAME.fa" --sequence CTAACAGTTGCTTTTATCACNGG -d 5 > "$BASENAME-output5.bed"
    
    # Process bed file
    awk 'NR>1 { $2=""; print }' "$BASENAME-output5.bed" > "$BASENAME-modified_output5.bed"

    # Use liftOver
    liftOver "$BASENAME-modified_output5.bed" "$BASENAME-grch38.chain" "$BASENAME-liftover5.bed" "$BASENAME-unlifted.bed" || echo "liftOver failed for $BASENAME"
done
```

## Fourth step

### Matrix

The output is a matix showing the frequency of interested sequnces with 5 or less mismatches in FASTA file among all individuals. 

```
#!/bin/bash

# Define the directory containing the files
dir="/lizardfs/salehi/crispr"

# Process all relevant files and merge them
awk '
{
    key = $1;  # Chromosome identifier
    start = $2; # Start position
    end = $3;   # End position
    seq = $4;   # Sequence
    mismatches = $5; # Mismatches

    # Initialize entry if not exists
    if (!(key in data)) {
        data[key] = start "\t" end "\t" seq "\t" mismatches;
        haplotype[key] = "";
    } else {
        # Update start and end positions
        if (start < split(data[key], arr, "\t")[1]) {
            data[key] = start "\t" arr[2] "\t" arr[3] "\t" arr[4];
        }
        if (end > split(data[key], arr, "\t")[2]) {
            data[key] = arr[1] "\t" end "\t" arr[3] "\t" arr[4];
        }
    }
    # Append haplotype count (assuming 5th column is same for all)
    haplotype[key] = haplotype[key] ? haplotype[key] "\t" mismatches : mismatches;
} END {
    print "Chromosome\tStart\tEnd\tSequence\tmismatches\tHap1\tHap2...";
    for (key in data) {
        print key "\t" data[key] "\t" haplotype[key];
    }
}' "$dir"/*/*-liftover5.bed > merged_output.txt

echo "Output written to merged_output.txt"
```
