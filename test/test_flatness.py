import os
import sys
import subprocess
import numpy as np
from scipy.stats import pearsonr

def create_files():
    np.random.seed(42)
    os.makedirs('test_adjS_flat/in', exist_ok=True)
    
    n_regions = 200
    cpgs_per_region = 100
    coverage = 20

    ref_file = open('test_adjS_flat/in/reference.tsv', 'w')
    sim_file = open('test_adjS_flat/in/simulations.tsv', 'w')
    bed_file = open('test_adjS_flat/in/regions.bed', 'w')
    
    start_pos = 1
    
    regions_info = []

    for i in range(n_regions):
        # Sample a true probability p for this region
        p = np.random.uniform(0.1, 0.9)
        regions_info.append(p)
        
        region_start = start_pos
        
        for j in range(cpgs_per_region):
            pos = start_pos + j * 2
            ref_file.write(f"chr1\t{pos}\t{pos+1}\n")
            
            # binomial sampling
            meth = np.random.binomial(coverage, p)
            sim_file.write(f"chr1\t{pos}\t{meth}\t{coverage}\t{meth/coverage:.4f}\n")
            
        region_end = start_pos + cpgs_per_region * 2
        bed_file.write(f"chr1\t{region_start}\t{region_end}\n")
        
        start_pos += cpgs_per_region * 2 + 100
        
    ref_file.close()
    sim_file.close()
    bed_file.close()

def run_yamet():
    cmd = [
        "../method/build/yamet",
        "-c", "test_adjS_flat/in/simulations.tsv",
        "-r", "test_adjS_flat/in/reference.tsv",
        "-i", "test_adjS_flat/in/regions.bed",
        "--norm-det-out", "current_test_adjS.out",
        "--cores", "1"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print("yamet failed!")
        print(res.stderr)
        sys.exit(1)

def evaluate():
    adjS_vals = []
    avg_meth_vals = []
    
    with open("current_test_adjS.out", "r") as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 6:
                # Based on file_map.cpp exportNormDetOut
                # columns might be: chr, start, end, <file_ids>..., shannon_norm, avg_meth
                # If no cluster, there's 1 file, then shannon_norm, avg_meth
                # wait, let me look at exportNormDetOut header again
                # "chr\tstart\tend\tsimulations.tsv\tshannon_norm\tavg_meth"
                try:
                    sampen_norm = float(parts[3])
                    avg_meth = float(parts[5])
                    adjS_vals.append(sampen_norm)
                    avg_meth_vals.append(avg_meth)
                except ValueError:
                    pass

    # Calculate Pearson correlation with distance from 0.5 to detect U-shape
    # If S has an inverted U-shape, it will be highly negatively correlated with |avg_meth - 0.5|
    dist_from_half = [abs(m - 0.5) for m in avg_meth_vals]
    corr, pval = pearsonr(dist_from_half, adjS_vals)
    print(f"Sample size: {len(adjS_vals)}")
    print(f"Correlation between adjS (sampen_norm) and |avg_meth - 0.5|: {corr:.4f} (p={pval:.4g})")
    
    # Check if relationship is flat: correlation should be small
    if abs(corr) > 0.15:
        print("FAIL: The relationship is not flat, absolute correlation with distance from 0.5 is large (U-shape detected).")
        sys.exit(1)
    else:
        print("PASS: The relationship is flat.")
        sys.exit(0)

if __name__ == "__main__":
    import os
    # Move to root dir of repo to run (same directory root as github actions)
    # The script should be executed from yamet/test
    # But pathing above assumes yamet root is CWD.
    if hasattr(sys, 'argv') and len(sys.argv)>1:
        os.chdir(sys.argv[1])
    
    create_files()
    run_yamet()
    evaluate()
