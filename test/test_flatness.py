import os
import sys
import subprocess
import math
import random

def create_files(sparsity=1.0):
    random.seed(42)
    os.makedirs('test_adjS_flat/in', exist_ok=True)
    
    n_regions = 200
    cpgs_per_region = 100
    coverage = 20

    with open('test_adjS_flat/in/reference.tsv', 'w') as ref_file, \
         open('test_adjS_flat/in/simulations.tsv', 'w') as sim_file, \
         open('test_adjS_flat/in/regions.bed', 'w') as bed_file:
         
        start_pos = 1
        
        for i in range(n_regions):
            p = float(random.uniform(0.1, 0.9))
            region_start = float(start_pos)
            
            for j in range(cpgs_per_region):
                pos = start_pos + j * 2
                ref_file.write(f"chr1\t{pos}\t{pos+1}\n")
                
                # Apply sparsity (downsampling)
                if random.random() <= sparsity:
                    # binomial sampling
                    meth = sum(1 for _ in range(coverage) if random.random() < p)
                    sim_file.write(f"chr1\t{pos}\t{meth}\t{coverage}\t{meth/coverage:.4f}\n")
                
            region_end = float(start_pos + cpgs_per_region * 2)
            bed_file.write(f"chr1\t{int(region_start)}\t{int(region_end)}\n")
            
            start_pos += int(cpgs_per_region * 2 + 100)

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

def evaluate(label=""):
    adjS_vals = []
    avg_meth_vals = []
    
    with open("current_test_adjS.out", "r") as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 6:
                try:
                    sampen_norm = float(parts[3])
                    avg_meth = float(parts[5])
                    # Skip sentinel values (-1 encodes "no valid windows" in yamet)
                    if sampen_norm == -1:
                        continue
                    adjS_vals.append(sampen_norm)
                    avg_meth_vals.append(avg_meth)
                except ValueError:
                    pass

    n = len(adjS_vals)
    if n == 0:
        print(f"[{label}] No output parsed")
        return

    dist_from_half = [float(abs(m - 0.5)) for m in avg_meth_vals]
    mean_dist = sum(dist_from_half) / n
    mean_adjs = sum(adjS_vals) / n
    
    num = sum((dist_from_half[i] - mean_dist) * (adjS_vals[i] - mean_adjs) for i in range(n))
    den_dist = sum((dist_from_half[i] - mean_dist)**2 for i in range(n))
    den_adjs = sum((adjS_vals[i] - mean_adjs)**2 for i in range(n))
    
    corr = num / math.sqrt(den_dist * den_adjs) if (den_dist * den_adjs) > 0 else 0
    
    print(f"[{label}] Sample size: {n}")
    print(f"[{label}] Correlation between adjS and |avg_meth - 0.5|: {corr:.4f}")
    if abs(corr) > 0.15:
        print(f"[{label}] FAIL: High correlation -> Shape is not flat (U-shape detected).")
    else:
        print(f"[{label}] PASS: The relationship is flat.")

if __name__ == "__main__":
    if hasattr(sys, 'argv') and len(sys.argv)>1:
        os.chdir(sys.argv[1])
        
    print("Testing Dense Data (100% density)")
    create_files(sparsity=1.0)
    run_yamet()
    evaluate(label="Dense")
    
    print("\nTesting Sparse Data (20% density)")
    create_files(sparsity=0.20)
    run_yamet()
    evaluate(label="Sparse")
