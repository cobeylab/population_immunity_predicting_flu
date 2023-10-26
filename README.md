# population_immunity_predicting_flu

## Data
- Data/sera/data contains HA and NA titers
- Data/sequences/data contains HA and NA reference virus sequences

## Analyses

### Genealogy of HA and NA
Analyses/genealogy
- `sample_sequences_for_genealogy.py` downsamples HA and NA sequences and saves the sampled sequences as “*.fasta” files. 
- `write_meta.py` takes sampled sequences for genealogy and assign them to clades, save files as  “*_clade.csv”. 
- `tree_201718_shaded.R` makes plots of HA genealogy for season 16-17 and 17-18.
- `na_tree_201718_shaded.R` makes plots of NA genealogy for season 16-17 and 17-18.

### Forecasting susceptibilities
Analyses/forecasting/src
- `binning.R` As pre-processing before calculating susceptibility, ages are binned and the age bins are saved as age_bin_ha.rds and age_bin_na.rds using this code. 
- `call_susceptibility.R` sets the options for calculating relative susceptibility and calls `susceptibility.R,` which calculates and saves the results of relative susceptibility according to the options. 
- `plot_susceptibility_rank.R` takes the results of relative susceptibility using cutoff and makes susceptibility plots such as Figure 2. 
- `plot_susceptibility_gmt_rank.R` takes results of relative susceptibility using GMT and makes susceptibility plots. 
- `titer_distribution.R` makes plots of titer distribution.

### Frequency and age distribution of clades
Analyses/frequency/src
- `write_age_and_allele.py` takes sequences and meta data and assign sequences to clades. Uses `age_distribution_by_clade_functions.py` and `define_virus.py` and make “*_clade_assigned_season_*.csv” files.
- `make_plot_frequency.R` takes “*_clade_assigned_season_*.csv” files and makes plots of frequencies by clade.
- `difference_in_proportion_test.R` performs chi-square test to test the difference in the proportion of 3C.2A2 between children and adults.

Analyses/age_distribution_of_clades/src
- `age_distribution.R` takes “*_clade_assigned_season_*.csv” files in frequency/data_clade_assigned folders and makes age distribution plots.

### Correlation between pairs of viruses
Analyses/correlation_and_clustering/src
- `correlation_analysis_boot_rmv_all_undetectable.R` calculates correlation coefficients for each reference virus pairs in each age group. It also tests significant difference of the correlations using bootstrapping. This code uses `analysis_uitl.R`, `correlation_functions.R`, and `correlation_function_bootstrap.R`.
- `plot_correlation.R` takes the results from `correlation_analysis_boot_rmv_all_undetectable.R` and make heat maps.
- `correlation_analysis_boot_equal_bin.R` calculates correlation coefficients for each reference virus pairs in each age group which is binned with equal number of years. It also draws plot of this result.
- `correlation_analysis_detail_plot.R` makes plots showing correlation between titers to different reference viruses.

### Cosine similarity between pairs of people
Analyses/correlation_and_clustering/src
- `cosine_similarity_continuous_1000.R` calculates cosine similarity between pairs of people using HA titer vectors and makes result plots.
- `cosine_similarity_continuous_hana_1000.R` calculates cosine similarity between pairs of people using HA-NA titer vectors and makes result plots.

### Clustering analysis
Analyses/correlation_and_clustering/src
- `clustering_analysis_using_cosine_smilarity.R` performs clustering analysis and makes result plots.
- `correlation_analysis_using_clusters.R` calculate correlation coefficient for each clustering group and make result plots.
- `larger_effect_in_targeting.R` performs ANOVA to test if between-age-group titer difference is larger than within-age-group titer difference.

### HA NA analysis
Analyses/ha_na_analysis
- `hana_distribution_analysis.R` makes plot of age distribution of HA and NA detectable titers. Performs logistic regression to test if fraction of detectable titers change by age. 
- `ha_na_correlation.R`, `ha_na_correlation_by10years.R` make plots of correlation between HA and NA titers.






