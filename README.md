README
=======

<p align="center">
<img src="https://raw.githubusercontent.com/befriendabacterium/communityinvasion/main/inputs/external/images/schematic.jpg" width="100%">
</p>

## Description

* author = Matt Lloyd Jones
* web = https://github.com/befriendabacterium/
* date = July 1st, 2021
* description = This repository contains a pipeline the reproduce the analyses in the preprint 'Relationships between community composition, productivity and invasion resistance in semi-natural bacterial microcosms' (https://doi.org/10.1101/2019.12.18.881102).

## Pipeline (must be run from RStudio)

All the scripts used in the pipeline are in the directory `src`. Input data for the analysis is in the `inputs` directory - with the raw composition, growth and invasion data and in the 1_raw sub-directory; wrangled/munged but not yet ready for analysis versions of these data in the `2_wrangled` sub-directory; and ready-for-analysis data in the `3_ready` sub-directory. Additionally in `inputs`, there's a sub-directory `4_semmodelspecifications`, containing the structural equation model specifications as text files (but no data) and another called `external`, which contains some external inputs needed that were produced outside of R. Outputs of the analyses are in the `outputs` directory - with the random forest regressions used to compare (compositional) dimensionality reduction approaches, and the final forests using the selected dimensionality reduction and the growth variables (`\full`) contained in the `randomforest` sub-directory; the structural equation models outputs in the `sem` sub-directory; the figures in the `figures` directory.

The whole pipeline can be rerun by running `Z_run_wholepipeline.R`, which empties the directories to contain only the data needed to start the whole analysis, and runs each script in correct order.

### 0. Check what packages are required and download them

  * **script**: `0_acquirepackages` - Run this script to check what packages are required across all the scripts in the `src` directory, compare them to your installed packages, and install and load the ones you're missing so you can run the analysis easily.

### 0. Download data from OSF

  * **script**: `0_downloaddata.R` - Run this script to download the data from OSF.

  * **inputs**: There are three options for downloading the 'inputs' and 'outputs' data, related to where in the pipeline you want to start running the code.

   * `1_start`: These inputs and outputs folders contain only the 'rawest' data i.e. the data needed to run the analysis from the start to the end. If you want to run the whole pipeline, after downloading this, just run `0_run_wholepipeline.R`, though it'll take a while so perhaps go cook something.

   * `2_preanalysis`: These inputs and outputs folders contain the data after pipeline steps 1-6 have been run i.e. time-consuming pre-analysis has been done. This allows you to run the analysis from the `8_randomforests.R` point; the point described in the main body of the manuscripts' Results.

   * `3_end`: These inputs and outputs folders contain the end result of running all of the steps/scripts in the pipeline.

You can choose which option you want by changing the value of the 'whichpoint' R object to the relevant string (i.e. '1_start', '2_preanalysis', '3_end').

  * **outputs**:

   * Same as inputs; the script is just downloading the data from OSF.

### 1. Cleaning of explanatory variables

  * **script**: `1_cleanexplanatoryvariables.R` - Cleaning of the genotypic and phenotypic explanatory variables; removal of rare sequences for the genotypic data (OTU table) and removal of assay measurements not used for this analysis for the phenotypic data.

  * **inputs**:
    * `inputs/1_raw/composition/Functional_Redundancy_Treehole_genera.csv` - The raw OTU table for the inoculated resident communities obtained from the sequencing company MRDNA, that has been through their proprietary, but standard bioinformatic pipeline (described in Methods: Computational techniques in the manuscript [1]). This is the basis of all compositional explanatory variables used in further analyses.
    * `inputs/1_raw/growth/20151016_Functions_remainder.csv` - The raw phenotypic data for the resident communities obtained from the laboratory assays of their cell numbers, respiration, enzyme and ATP activity in the growth period prior to invasion at 14 days. This is the basis of all phenotypic explanatory variables used in further analyses. N.B There are 4 pseudoreplicates of each assay, which will be averaged in step 3.

  * **outputs**:
    * `inputs/2_wrangled/composition/composition_otu_cleaned.csv` - OTU table with rare sequences removed.
    * `inputs/2_wrangled/growth/growthdata_cleaned.csv` - Dataframe of phenotypic measurements with irrelevant measurements removed.

### 2. Aggregation/averaging (mean) of the pseudoreplicated assays

  * **scripts**: `2_aggregate_assays.R` - This script aggregates the 4 pseudoreplicates of the phenotypic measurements (explanatory variables) and the invasion assays (response variables).

  * **inputs**:
    * `inputs/1_raw/invasion/invasiondata_luminescence_4assays.csv` - Raw invasion assays data obtained from the lab.
    * `inputs/2_wrangled/growth/growthdata_cleaned.csv`- Cleaned phenotypic assays data from Step 2.

  * **outputs**:
    * `inputs/2_wrangled/growth/growthdata_cleaned_assaymeans.csv` - Phenotypic data averaged across pseudoreplicates.
    * `inputs/2_wrangled/invasion/invasiondata_luminescence_assaymeans.csv` - Invasion data averaged across psuedoreplicates.

### 3. Match all the wrangled datasets by community (end of data preparation)

  * **scripts**: `3_matchdatasets.R` - Matches all the wrangled datasets - genotypic (explanatory), phenotypic (explanatory), invasion (response) - according to the Community ID (some communities were not assayed in all assay types).

  * **inputs**:
    * `inputs/2_wrangled/composition/composition_otu_cleaned.csv` - OTU table with rare sequences removed from Step 1.
    * `inputs/external/functionInk/SamplePropsOFtaxaClus_Time0_NL_Average_StopStep-91_AbundMean.dat` - The abundances of each of the functional groups calculated with the functionInk pipeline in Alberto Pascual Garcia's repository (https://github.com/apascualgarcia/Invasion)
    * `inputs/growth/growthdata_cleaned_assaymeans.csv` - Phenotypic data averaged across pseudoreplicates, from Step 3.
    * `inputs/2_wrangled/invasion/invasion/invasiondata_luminescence_assaymeans.csv` - Phenotypic data averaged across pseudoreplicates.

  * **outputs**:
    * `inputs/3_ready/composition/composition_otu_matched.csv` - Cleaned OTU table and matched to contain only communities assayed in phenotypic and invasion assays.
    * `inputs/3_ready/growth/growthdata_assaymeans_matched.csv` - Phenotypic data averaged across pseudoreplicates, and matched to contain only communities assayed in genotypic and invasion assays.
    * `inputs/3_ready/invasion/invasion/invasiondata_luminescence_assaymeans_matched.csv` - Invasion data averaged across pseudoreplicate, and matched to contain only communities assayed in genotypic and phenotypic dataframes.

### 4. Conversion of luminescence measurements into estimated cells/ml for invasion assays

  * **script**: `4_convertinvasionlumtocfu.R` - Converts the luminescence assays for invader abundance estimation into CFU/ml units.

  * **inputs**:
    * `inputs/1_raw/invasion/invadercalibrationgrowthcurve_96h_8reps.csv` - Data from assays performed with the invaders in monoculture to calibrate luminescence against CFU/ml.
    * `inputs/3_ready/invasion/invasion/invasiondata_luminescence_assaymeans_matched.csv` - Invasion data averaged across psuedoreplicates, from Step 3.

  * **outputs**:
    * `inputs/2_wrangled/invasion/invasiondata_cfu_assaymeans_matched.csv` - Invasion data averaged across pseudoreplicates.

### 5. Calculation of various dimensionality reductions of the OTU data

  * **scripts**: `5_calculatedimensionalityreductions.R` - This script A) aggregates the abundances of the cleaned OTU table at each level to provide tables at each taxonomic rank resolution to explore as a dimensionality reduction technique; B) Reads in the functional group abundances, matches them to the Communities (they are in a different order); C) Calculates a PCoA on the OTU table

  * **inputs**:
    * `inputs/2_wrangled/composition/composition_otu_cleaned.csv` - OTU table with rare sequences removed.
    * `inputs/external/taxonomy/taxonomyNCBI.csv` - A table of the taxonomy of each OTU at the Genus, Family, Order, Class and Phylum levels, fetched from NCBI.
    * `inputs/external/functionInk/SamplePropsOFtaxaClus_Time0_NL_Average_StopStep-91_AbundMean.dat` =  The abundances of each of the functional groups calculated with the functionInk pipeline in Alberto Pascual Garcia's repository (https://github.com/apascualgarcia/Invasion).

  * **outputs**:
    * `/inputs/3_ready/composition/composition_genus_matched.csv` - Cleaned OTU table aggregated at the assigned-genus level.
    * `/inputs/3_ready/composition/composition_family_matched.csv` - Cleaned OTU table aggregated at the assigned-genus level.
    * `composition/composition_order_matched.csv` - Cleaned OTU table aggregated at the assigned-genus level.
    * `/inputs/3_ready/composition/composition_class_matched.csv` - Cleaned OTU table aggregated at the assigned-genus level.
    * `/inputs/3_ready/composition/composition_phylum_matched.csv` - Cleaned OTU table aggregated at the assigned-genus level.
    * `/inputs/3_ready/composition/composition_funcgroups_matched.csv` - The abundances of each of the functional groups calculated with the functionInk pipeline in Alberto Pascual Garcia's repository (https://github.com/apascualgarcia/Invasion), in the right Community order.
    * `/inputs/3_ready/composition/composition_pcoa.RDS` - PCoA of OTU abundance data based on Jensen-Shannon distance of communities.

### 6. Find the best dimensionality reduction technique for the genotypic data

  * **scripts**: `6_find_bestdimreduction_prerunranfs.R` - These scripts compare the performance (for explaining invasion success) of aggregation of abundances by taxonomic rank aggregation, aggregation of abundances by functional group, and PCoA as dimensionality reduction techniques for the compositional data. The former script `...prerunranfs.R` reads in pre-run random forests (for convenience as they take a long time to run), the latter script `...prerunranfs.R` reruns the random forests and overwrites the pre-run ones (though output should be identical).

  * **inputs**:
      * `/inputs/3_ready/composition/composition_TAXONOMICRANK_matched.csv` - Each of the cleaned OTU tables, aggregated at each taxonomic level, and matched to contain only communities assayed in phenotypic and invasion assays.
      * `/inputs/3_ready/composition/composition_funcgroups_matched.csv` - The abundances of each of the functional groups calculated with the functionInk pipeline in Alberto Pascual Garcia's repository (https://github.com/apascualgarcia/Invasion), in the right Community order.
      * `inputs/3_ready/composition/PCoA/composition_pcoa.RDS` - PCoA of the OTU table.

 * **outputs**:
      * `outputs/randomforest/composition/DIMREDUCTION/INVADERNAME.TIMEOFSAMPLING.ranf_DIMREDUCTION.RDS` - Each of the random forests for each dimensionality reduction technique (DIMREDUCTION) filed into a relevant folder containing the 6 models (2 invaders (INVADERNAME) x 3 times of sampling (TIMEOFSAMPLING) after invasion).
      * `outputs/randomforest/extracted_varexps_df.csv` - the variance explained by each of the models, formatted as a table containing metadata about the model.
      * `outputs/randomforest/all_randomforests_list.RDS` - A big RDS file containing all the randomforests, saved as an R list.

### 7. Calculate OTU-level diversity metrics as additional dimensionality reductions of compositional data

  * **scripts**: `7_calculatediversitymetrics.R` - Calculates Simpson's diversity, Rao's diversity and Community-Invader phylogenetic distance metrics.

  * **inputs**:
    * `inputs/3_ready/composition/composition_otu_matched.csv` - Each of the cleaned OTU tables, aggregated at each taxonomic level, and matched to contain only communities assayed in phenotypic and invasion assays.
    * `inputs/external/geneious/tree101_rooted_invdistNAMES.newick` - The phylogenetic tree produced from the 16S sequences of the OTUs alongside those of the two invaders using Geneious and the methods detailed in the manuscript [1].
    * `inputs/external/taxonomy/sp.names.final.csv` - A simple table mapping the OTU names/numbers to the taxonomic assignment names in the raw OTU table.
    * `inputs//3_ready/invasion/invasion/invasiondata_cfu_assaymeans_matched.csv` - Invasion data averaged across pseudoreplicate, and matched to contain only communities assayed in genotypic and phenotypic dataframes.

  * **outputs**: Working directory `inputs/3_ready/composition`
      * `outputs/diversitymetrics.csv` - A table containing Simpson's diversity., Rao's diversity and Community-Invader phylogenetic distance (for each invader, in 2 separate columns) of each resident Community.

### 8. Run full random forests

  * **scripts**: `8_randomforests.R` - Script to run the full random forest regressions for the main analysis, using the chosen dimensionality reduction technique for the genotypic/OTU data and the phenotypic assays as explanatory variables, and the invader success measurements as response variables.

  * **inputs**:
    * See `data_prepper.R` (this script is sourced to pull in and format all input data, so inputs are same)

  * **outputs**:
    * `outputs/randomforest/full/extracted_varexps_df.csv` - Table of variance explained (%) by each of the 6 random forests (2 invaders x 3 sampling points)
    * `outputs/randomforest/full/all_randomforests_list.RDS` - RDS file containing an R list of all 6 random forests (2 invaders x 3 sampling points)

### 9. Structural equation modelling

  * **scripts**:
    * This script (`9_launch_SEM_multiple_simpsonsdiversity.R`) runs the No mediation, Partial mediation and Complete mediation SEM models, outputting the models and a model comparison. The user must define the SEM models in lavaan format in external files located in the folder `inputs/4_semmodelspecifications`.Then, the variable `selectModel` should point to that file. The input data is located in `inputs/3_raw` and called with `data_prepper.R`. Note that the function to retrieve the files works only with *RStudio*, and it should be manually included otherwise. Other functions are coded in different files, see details in the respective scripts.

  * **inputs**:
    * `inputs/4_semmodelspecifications` - The model specifications encoded in a text file, with the ad-hoc extension .mod
    * Experimental data - See data_prepper.R (this script is sourced to pull in and format all input data, so inputs are same)

  * **output**: `outputs/sem`
    * The script will run each model and store it in a .mod file, and a summary of results of the fit and the MI indexes is stored in another file (.fit). It will then create visual representations (PDF Imag) with the script `SEM_PathPlot.R` for the plain model (excluding variances), the structural model and the full model. If several models are analysed with the "multiple" version, an analysis of all the models can be made.

### 10. Extract model coefficients for interpretation & plotting

 * **scripts**:
   * `10_coeffs_extractor.R` - This script extracts the coefficients from whichever set of 3 SEM models have been run (i.e. with or without diversity), in order to estimate the relationships between composition, productivity and invasion.

 * **inputs**:
    * See 9. Structural equation modelling - this script sources one of these scripts for all inputs (by default `src/9_launch_SEM_multiple_simpsonsdiversity.R`)

 * **output**:
    * `outputs/sem/modelcomparison.csv` - Model comparison of the No mediation, Partial mediation and Complete mediation models, in a table.
    * `outputs/sem/all_model_coeffs.csv` - Coefficients for all the relationships in the No mediation, Partial mediation and Complete mediation models, in a table.

### X. Data preperation for main analyses

  * **scripts**: `X_data_prepper.R` - This script reads in the ready-for-analysis data and does some formatting and scaling of the data for the main analyses in the paper.

  * **inputs**:
    * `inputs/3_ready/composition/diversitymetrics.csv` - A table containing Simpson's diversity., Rao's diversity and Community-Invader phylogenetic distance (for each invader, in 2 separate columns) of each resident Community.
    * `inputs/3_ready/composition/funcgroups/SamplePropsOFtaxaClus_Time0_NL_Average_StopStep-91_ZscoreMean.dat` - The abundances of each of the functional groups calculated with the functionInk pipeline in Alberto Pascual Garcia's repository (https://github.com/apascualgarcia/Invasion).
    * `inputs/3_ready/growth/growthdata_assaymeans_matched.csv` - Phenotypic data averaged across pseudoreplicates, and matched to contain only communities assayed in genotypic and invasion assays.
    * `inputs/3_ready/invasion/invasion/invasiondata_cfu_assaymeans_matched.csv` - Invasion data averaged across pseudoreplicate, and matched to contain only communities assayed in genotypic and phenotypic dataframes.

  * **outputs**:
    * `inputs/3_ready/divandgrowth_data.csv` - A table containing the diversity and phenotypic explanatory variables for use in the random forest models (composition dimensionality reductions are left out to enabling easy substitution of a different technique).

### X. Sensitivity analysis

  * **scripts**: `X_sensitivityanalysis.R` - This script runs a sensitivity analysis to check how robust the results of the SEM model comparison are to changes in the data below the strict detection limit of the luminescence assay for measuring invasion success. Basically, it just randomises the values below the detection limit for each invader within the range of possible values and re-runs the SEMs and model comparison 999 times, checking how often the result (i.e. best fitting model) changes.

  * **inputs**:
    * Experimental data - See `data_prepper.R` (this script is sourced to pull in and format all input data, so inputs are same)
    * See 9. Structural equation modelling - this script sources one of these scripts for all inputs (by default `src/9_launch_SEM_multiple_simpsonsdiversity.R`)

* **outputs**
    * `outputs/sem/aicranks_sensitivityanalysis` - A .csv file containing the model comparison tables of the 999 randomisations/iterations, bound into one dataframe.

### Figure 1

  * **scripts**: `Figure1.R` - This script makes Figure 1 in the manuscript, which contains a schematic explaining the pre-experiment processing of the communities and the experimental set-up, and 2 graphs showing the broad patterns of invader survival across the 3 post-invasion sampling points.

  * **inputs**:
    * `inputs/external/images/schematic.jpeg` - Schematic describing steps before and during experiment.
    * `inputs/1_raw/invasion/invadercalibrationgrowthcurve_96h_8reps.csv` - Data from assays performed with the invaders in monoculture to calibrate luminescence against CFU/ml.
    * `inputs/3_ready/invasion/invasion/invasiondata_luminescence_assaymeans_matched.csv` - Invasion data averaged across psuedoreplicates, from Step 3.
    * `inputs/1_raw/invasion/detectionlims.RDS` - RDS file containing the detection limits of each of the 2 invaders for the invasion assay.

  * **outputs**:
    * `Figure1.tiff` - TIFF image of Figure 1.

### Figure 2

  * **scripts**: `Figure2.R` - This script makes Figure 2 in the manuscript, which summarises the variance explained by each of the full random forests and the variable importance of each of the explanatory variables.

  * **inputs**:
    * `outputs/randomforest/full/varimps.csv` - Variable importance matrices for all 6 random forests (2 invaders x 3 sampling points)
    * `outputs/randomforest/full/extracted_varexps_df.csv` - Table of variance explained (%) by each of the 6 random forests (2 invaders x 3 sampling points)
    * `outputs/randomforest/full/all_randomforests_list.RDS` - RDS file containing an R list of all 6 random forests (2 invaders x 3 sampling points)

  * **outputs**: Working directory `outputs/figures`
    * `Figure2.tiff` - TIFF image of Figure 2.

### Figure 3

  * **scripts**: `Figure3.R` - This script makes Figure 3 in the manuscript, a panel figure summarising some of the main relationships between _P. fluorescens_ invasion success and selected explanatory variables.

  * **inputs**:
    * See `data_prepper.R` (this script is sourced to pull in and format all input data, so inputs are same)

  * **outputs**:
    * `outputs/figures/Figure3.tiff` - TIFF image of Figure 3.

### Figure 4

  * **scripts**: `Figure4.R` - This script makes Figure 4 in the manuscript, which shows the structure and fit of the 3 main structural equation models used in the MS, using the extracted fit and coefficients information.

  * **inputs**:
    * `outputs/sem/modelcomparison.csv` - Model comparison of the No mediation, Partial mediation and Complete mediation models, in a table.
    * `outputs/sem/all_model_coeffs.csv` - Coefficients for all the relationships in the No mediation, Partial mediation and Complete mediation models, in a table.

  * **outputs**:
    * `outputs/figures/Figure4.tiff` - TIFF image of Figure 4.

### Supplementary Figure 1

  * **scripts**: `SupplementaryFigure1.R` - This script makes Supplementary Figure 1 in the manuscript, which compares the performance of various dimensionality reductions of composition when permuted in random forests.

  * **inputs**:
    * `outputs/randomforest/extracted_varexps_df.csv` - Variance explained extracted from all random forests
    * `outputs/randomforest/all_randomforests_list.RDS` - All random forests stored as an R list/RDS.

  * **outputs**:
    * `outputs/figures/SupplementaryFigure1.tiff` - TIFF image of Supplementary Figure 1

### Supplementary Figure 2

  * **scripts**: `SupplementaryFigure2.R` - This script makes Supplementary Figure 2 in the manuscript, which shows the relationships between _P. fluorescens_ invader survival and the four enzyme assays.

  * **inputs**:
    * See `inputs/3_ready/data_prepper.R` (this script is sourced to pull in and format all input data, so inputs are same)

  * **outputs**:
    * `outputs/figures/SupplementaryFigure2.tiff` - TIFF image of Supplementary Figure 2

### Supplementary Figure 3

  * **scripts**: `SupplementaryFigure3.R` - This script makes Supplementary Figure 3 in the manuscript, which is a rank-abundance plot of the OTUs in the communities.

  * **inputs**:
    * See `data_prepper.R` (this script is sourced to pull in and format all input data, so inputs are same)

  * **outputs**:`
    * `outputs/figures/SupplementaryFigure3.tiff` - TIFF image of Supplementary Figure 2

### Supplementary Figure 4

  * **scripts**: `SupplementaryFigure4.R` - This script makes Supplementary Figure 4 in the manuscript, which is a plot of the first two dimensions of the PCoA.

  * **inputs**:
    * `inputs/3_ready/composition/composition_pcoa.RDS` - PCoA of the OTU abundance data from Step 5.

  * **outputs**: Working directory `outputs/figures`
    * `SupplementaryFigure4.tiff` - TIFF image of Supplementary Figure

### Supplementary Figures 5-7

  * **scripts**: `SupplementaryFigure5to7.R` - This script makes Supplementary Figure 2 in the manuscript, which shows the relationships between _P. fluorescens_ invader survival and the four enzyme assays.

  * **inputs**:
    * None, just renames and moves files

  * **outputs**:
    * `outputs/figures/SupplementaryFigure5.tiff` - TIFF image of Supplementary Figure 5
    * `outputs/figures/SupplementaryFigure6.tiff` - TIFF image of Supplementary Figure 6
    * `outputs/figures/SupplementaryFigure7.tiff` - TIFF image of Supplementary Figure 7

## REFERENCES

[1] Jones, M. L., Rivett, D. W., Pascual-García, A., & Bell, T. (2021). Relationships between community composition, productivity and invasion resistance in semi-natural bacterial microcosms. bioRxiv.

[2] Pascual‐García, A., & Bell, T. (2020). functionInk: An efficient method to detect functional groups in multidimensional networks reveals the hidden structure of ecological communities. Methods in Ecology and Evolution, 11(7), 804-817.

[3] Pascual-García, A., & Bell, T. (2020). Community-level signatures of ecological succession in natural bacterial communities. Nature communications, 11(1), 1-11.

[4] Rivett, D. W. & Bell, T. Abundance determines the functional role of bacterial phylotypes in complex communities. Nat. Microbiol. 3, 767 (2018).
