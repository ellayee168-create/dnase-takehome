#!/usr/bin/env python
"""
Convert BAM files to base-resolution BigWigs for DNase-seq.

Replicates ChromBPNet's reads_to_bigwig.py methodology exactly.
"""

import os
import subprocess
import tempfile
import pyfaidx

GENOME_FA = "reference/hg38.fa"
CHROM_SIZES = "reference/hg38.chrom.sizes"

# DNase: AlphaGenome says shifts are 0/+1
# ChromBPNet auto-detects either (0,0) or (0,1)
# We assume the standard case: detected shift = (0, 1)
# So: plus_shift_delta = -0 = 0, minus_shift_delta = 1 - 1 = 0
# 
# But AlphaGenome paper says "Shifts of 0/+1 for DNase-seq were applied"
# This means the FINAL positions should be: +0 for plus, +1 for minus
# So if detected=(0,0), delta=(0,1) -> applies +1 to minus strand end
PLUS_SHIFT_DELTA = 0
MINUS_SHIFT_DELTA = 1  # This gets added to the END coordinate for - strand

def load_valid_chroms_from_fasta(genome_fa):
    """Load valid chromosome names from reference FASTA (ChromBPNet method)."""
    with pyfaidx.Fasta(genome_fa) as g:
        return set(g.keys())

def convert_bam(bam_path, output_prefix, valid_chroms):
    """
    Convert BAM to base-resolution BigWig for DNase-seq.
    
    Matches ChromBPNet's reads_to_bigwig.py exactly:
    1. BAM -> tagAlign (BED6) via bedtools bamtobed (NO quality filtering)
    2. Filter chromosomes not in reference FASTA
    3. Apply shifts: + strand start += delta, - strand end += delta
    4. bedtools genomecov -bg -5 (count 5' ends)
    5. Sort and convert to BigWig
    """
    output_bw = f"{output_prefix}_unstranded.bw"
    
    if os.path.exists(output_bw):
        print(f"  Already exists: {output_bw}")
        return output_bw
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bedGraph', delete=False) as tmp_bg:
        bedgraph_path = tmp_bg.name
    
    try:
        # ChromBPNet's awk command (reformatted for clarity):
        # For + strand: print $1, $2+plus_shift_delta, $3, $4, $5, $6
        # For - strand: print $1, $2, $3+minus_shift_delta, $4, $5, $6
        #
        # Then bedtools genomecov -5 counts the 5' end:
        # - For + strand reads: position = start
        # - For - strand reads: position = end - 1 (BED is 0-based half-open)
        
        # Build the command matching ChromBPNet exactly
        awk_cmd = """awk -v OFS="\\t" '{{if ($6=="+"){{print $1,$2{plus:+},$3,$4,$5,$6}} else if ($6=="-") {{print $1,$2,$3{minus:+},$4,$5,$6}}}}'""".format(
            plus=PLUS_SHIFT_DELTA, 
            minus=MINUS_SHIFT_DELTA
        )
        
        cmd = f"""
        bedtools bamtobed -i {bam_path} | \
        {awk_cmd} | \
        sort -k1,1 | \
        bedtools genomecov -bg -5 -i stdin -g {CHROM_SIZES} | \
        LC_COLLATE="C" sort -k1,1 -k2,2n
        """
        
        print(f"  Running ChromBPNet-style conversion...")
        
        # Run with chromosome filtering (matching ChromBPNet's stream_filtered_tagaligns)
        p1 = subprocess.Popen(
            ["bedtools", "bamtobed", "-i", bam_path],
            stdout=subprocess.PIPE
        )
        
        # Filter and shift in one pass
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_bed:
            tmp_bed_path = tmp_bed.name
            
            for line in p1.stdout:
                fields = line.decode('utf-8').strip().split('\t')
                chrom = fields[0]
                
                # Filter: only keep chromosomes in reference FASTA
                if chrom not in valid_chroms:
                    continue
                
                start = int(fields[1])
                end = int(fields[2])
                name = fields[3]
                score = fields[4]
                strand = fields[5]
                
                # Apply shifts exactly like ChromBPNet's awk
                if strand == "+":
                    new_start = start + PLUS_SHIFT_DELTA
                    new_end = end
                else:  # strand == "-"
                    new_start = start
                    new_end = end + MINUS_SHIFT_DELTA
                
                tmp_bed.write(f"{chrom}\t{new_start}\t{new_end}\t{name}\t{score}\t{strand}\n")
        
        p1.wait()
        
        # Sort, genomecov, sort (matching ChromBPNet)
        cmd2 = f"""
        sort -k1,1 {tmp_bed_path} | \
        bedtools genomecov -bg -5 -i stdin -g {CHROM_SIZES} | \
        LC_COLLATE="C" sort -k1,1 -k2,2n > {bedgraph_path}
        """
        
        subprocess.run(cmd2, shell=True, check=True)
        os.remove(tmp_bed_path)
        
        # Convert to BigWig
        print(f"  Converting to BigWig...")
        subprocess.run([
            "bedGraphToBigWig",
            bedgraph_path,
            CHROM_SIZES,
            output_bw
        ], check=True)
        
        print(f"  Created: {output_bw}")
        return output_bw
        
    finally:
        if os.path.exists(bedgraph_path):
            os.remove(bedgraph_path)

def main():
    os.makedirs("bigwigs", exist_ok=True)
    
    print(f"Loading valid chromosomes from {GENOME_FA}...")
    valid_chroms = load_valid_chroms_from_fasta(GENOME_FA)
    print(f"  Found {len(valid_chroms)} chromosomes")
    
    with open("metadata.tsv") as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split("\t")
            file_acc = parts[2]
            
            bam_path = f"bams/{file_acc}.bam"
            output_prefix = f"bigwigs/{file_acc}"
            
            if not os.path.exists(bam_path):
                print(f"BAM not found: {bam_path}")
                continue
            
            print(f"\nConverting {file_acc}...")
            convert_bam(bam_path, output_prefix, valid_chroms)

if __name__ == "__main__":
    main()