Every row in the parameter file is the parameterization of one simulated datasets. See `./input/parameters.tsv` for an example.    
The following parameters can be set:    

  - `n_cells`: Number of cells.   
  - `n_cpgs`: Number of CpGs to be simulated.   
  - `mode`: Mode for simulating lengths of missing and covered stretches, either `nb` (sampling from a negative binomial) or `rand` (repeated bernoulli).    
  - `covParams`: Coverage levels, can be either `low`, `medium`, `high` or `low_real`.    
 Correspond to different parameterizations of negative binomial distributions for sampling lengths of missing and covered stretches. `low_real` leads to estimating parameters for sampling lengths of covered (and missing) stretches from E7.5 cells acquired in the study of [Argelaguet et al., 2019](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121708).   
 Only used in mode `nb`.
 - `transMat`: Transition matrix for simulating sequences of single CpG methylation states (0, 1), levels representing transition matrices of the form $\begin{array}{c|ccc}
                                    & \text{0} & \text{1}\\
                                    \hline
                                    \text{0} & p_{0,0}  & p_{0,1}   \\
                                    \text{1} & p_{1,0}  & p_{1,1} \\
                                    \end{array}$ are:     
                                    \
      - `lmr` $\begin{bmatrix}0.8 & 0.2 \\ 0.8 & 0.2 \end{bmatrix}$            
      lowly methylated
      - `hmr` $\begin{bmatrix}0.2 & 0.8 \\ 0.2 & 0.8 \end{bmatrix}$        
      highly methylated 
      - `imr_cons` $\begin{bmatrix}0.8 & 0.2 \\ 0.2 & 0.8 \end{bmatrix}$        
      intermediatley methylated, conserved methylation across neighboring CpGs
      - `imr_rand` $\begin{bmatrix}0.5 & 0.5 \\ 0.5 & 0.5 \end{bmatrix}$          
      intermediatley methylated, uniform transition of methylation states across neighboring CpGs       
\        
\      
Resulting simulated datasets are named the following:    
`sim_<n_cells>_<n_cpgs>_<mode>_<covParams>_<transMat>`  

for example: 
`sim_10_50_nb_medium_hmr.tsv` or 
`sim_10_50_nb_<low_real>_<imr_rand>.tsv` 

Further look-up tables with the CpG-positions of the simulated files are generated named:
`cpgPositions_<n_cpgs>.tsv`