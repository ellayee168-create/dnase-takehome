#!/usr/bin/env python
"""Download BAM files listed in metadata.tsv"""

import os
import subprocess

def main():
    # Create download directory
    os.makedirs("bams", exist_ok=True)
    
    with open("metadata.tsv") as f:
        header = next(f)  # Skip header
        print(f"Header: {header.strip()}")
        
        for line in f:
            parts = line.strip().split("\t")
            # Columns: cell_line, experiment_accession, file_accession, bam_url, spot1_score, read_length, biological_replicates
            cell_line = parts[0]
            file_acc = parts[2]
            bam_url = parts[3]
            
            output_path = f"bams/{file_acc}.bam"
            
            if os.path.exists(output_path):
                print(f"Already exists: {output_path}")
                continue
            
            print(f"Downloading {file_acc} ({cell_line})...")
            subprocess.run([
                "wget", "-O", output_path, bam_url
            ], check=True)

if __name__ == "__main__":
    main()