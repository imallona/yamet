Every row in the parameter file is the parameterization of one simulated datasets. See `./input/parameters.tsv` for an example.    
The following parameters can be set:    

  - `n_cells`: Number of cells.   
  - `n_cpgs`: Number of CpGs to be simulated.   
  - `mode`: Mode for simulating lengths of missing and covered stretches, either `nb` (sampling from a negative binomial) or `rand` (repeated bernoulli).    
  - `covParams`: Coverage levels, can be either `low`, `medium`, `high` or `low_real`.    
 Correspond to different parameterizations of negative binomial distributions for sampling lengths of missing and covered stretches. `low_real` leads to estimating parameters for sampling lengths of covered (and missing) stretches from E7.5 cells acquired in the study of [Argelaguet et al., 2019](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121708).   
 Only used in mode `nb`.
 - `transMat`: Transition matrix for simulating sequences of single CpG methylation states (0, 1), levels representing transition matrices with entries: $p_{0,0}, p_{0,1}, p_{1,0}, p_{1,1}$ are:
      - `lmr`: $p_{0,0}=0.8,p_{0,1}=0.2, p_{1,0}=0.8, p_{1,1}=0.2$        
      lowly methylated
      - `hmr` $p_{0,0}=0.2, p_{0,1}=0.8, p_{1,0}=0.2, p_{1,1}=0.8$    
      highly methylated   
      - `imr_cons` $p_{0,0}=0.8, p_{0,1}=0.2, p_{1,0}=0.2, p_{1,1}=0.8$            
      intermediatley methylated, conserved methylation across neighboring CpGs
      - `imr_rand` $p_{0,0}=0.5, p_{0,1}=0.5, p_{1,0}=0.5, p_{1,1}=0.5$        
      intermediatley methylated, uniform transition of methylation states across neighboring CpGs       

     
Resulting simulated datasets are named the following:    
`sim_<n_cells>_<n_cpgs>_<mode>_<covParams>_<transMat>`  

for example: 
`sim_10_50_nb_medium_hmr.tsv` or 
`sim_10_50_nb_<low_real>_<imr_rand>.tsv` 

Further look-up tables with the CpG-positions of the simulated files are generated named:
`cpgPositions_<n_cpgs>.tsv`