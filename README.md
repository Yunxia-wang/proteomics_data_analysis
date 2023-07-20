# Proteomics data analysis

# Contents
- [1. Generate a PCA plot](#generate_PCA)
- [2. Calculate ratios](#calculate_ratios)
- [3. Identify proteins that change in expression over 2 fold](#identify_proteins)
- [4. perform GO-term enrichment analysis](#perform_GO_term)

# 1. Generate a PCA plot

The results of the PCA plot is "PCA for samples.png" in the directory "proteomics_data_analysis/Result_png/".

# 2. Calculate ratios

The calculation results of this calculating ratios are in the directory "proteomics_data_analysis/Result_data/".
The ratios of each experimental channel to the average of expanded channels within each donor is the "ratios of each experimental channel.csv". The mean ratio between the donors is the "mean ratio between the donors.csv".

# 3. Identify proteins that change in expression over 2 fold

Identification of changed proteins is in the directory "proteomics_data_analysis/Result_data/". All changed proteins including up-regulated and down-regulated are in "DEG_M.csv". Among them, the down-regulated ones are in "DOWN_DEG.csv". The down-regulated ones are in "UP_DEG.csv". Volcano plot of changed genes is "Volcano plot for samples.png" in "proteomics_data_analysis/Result_png".

# 4. perform GO-term enrichment analysis

The result of performing GO-term enrichment analysis on the changing proteins is "EnrichGO_result.csv" in "proteomics_data_analysis/Result_data/".

