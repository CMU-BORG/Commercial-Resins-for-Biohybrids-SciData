# README for the Supplemental Files to:

“Cytotoxicity and Characterization of 3D-Printable Resins Using a Low-Cost Printer for Muscle-based Biohybrid Devices”

AS Liao, K Dai, AB Irez, A Sun, MJ Bennington, S Schaffer, B Chopra, JM Seok, R Adams, YJ Zhang, VA Webster-Wood

2024

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14014330.svg)](https://doi.org/10.5281/zenodo.14014330)

Preprint Available on bioRxiv: [![Static Badge](https://img.shields.io/badge/DOI-10.1101%2F2024.10.31.621115-blue?link=https%3A%2F%2Fdoi.org%2F10.1101%2F2024.10.31.621115)](https://doi.org/10.1101/2024.10.31.621115)

Living Machines 2024 Conference Paper: [![Static Badge](https://img.shields.io/badge/DOI-10.1007%2F978--3--031--72597--5__27-blue)](https://doi.org/10.1007/978-3-031-72597-5_27)

Liao, A.S. et al. (2025). Biocompatibility of Asiga Dental Resins Using a Low-Cost Printer for Biohybrid Actuator Applications. In: Szczecinski, N.S., Webster-Wood, V., Tresch, M., Nourse, W.R.P., Mura, A., Quinn, R.D. (eds) Biomimetic and Biohybrid Systems. Living Machines 2024. Lecture Notes in Computer Science(), vol 14930. Springer, Cham. https://doi.org/10.1007/978-3-031-72597-5_27

## Overview:

These are the supplemental files associated with the article. The article encompasses 3 major investigations of multiple 3D-printable resins:

- Cytotoxicity

- Mechanical Testing

- Print Fidelity

Each directory encompasses files related to these studies.

### Directory: [Cytotoxicity](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/tree/main/Supplemental/Cytotoxicity)
The Excel spreadsheets in these files are the raw data exported from the plate reader with some organization in the subsequent sheets. The Minitab files encompass the statistical analyses and figure generation from the plate reader data.

#### File List

- CytotoxicityPlateReader_allExceptFL_asl_20240209_resin_72h_edited.xlsx: Raw plate reader data for the initial study, which included all resin samples except for the follow-up cytotoxicity analysis for Formlabs Silicone 40A

- CytotoxicityPlateReader_FollowUp_asl_20240402_resin_72h_followUp_withReRuns_aslEdit.xlsx: Raw plate reader data for the follow-up cytotoxicity analysis for Formlabs Silicone 40A

- CytotoxicityStats_allExceptFL_MINITAB_PLATEREADER_72HVWW_ASL_V06.mpx: Minitab data analysis file for the initial study, which included all resin samples except for the follow-up cytotoxicity analysis for Formlabs Silicone 40A

- CytotoxicityStats_followUp_72h.mpx: Minitab data analysis file for the follow-up cytotoxicity analysis for Formlabs Silicone 40A



### Directory: [CAD Files for 3D Printing Samples](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/tree/main/Supplemental/CAD%20Files%20for%203D%20Printing%20Samples)
This directory includes the SolidWorks part files that house the designs for the samples used in the study.

#### File List

- 96wellspecimen.SLDPRT: Design for the cytotoxicity analysis samples - these printed parts were cultured with the cells in a 96 well plate

- Compression_ASTM_D575_40percentscale.SLDPRT: Design for the compression test sample pucks for the elastomeric resins

- Compression_ASTM_D695_20percentscale.SLDPRT: Design for the compression test sample pucks for the rigid resins

- Dogbone_ASTM_D412C_halfsize.SLDPRT: Design for the tensile test dogbones for the elastomeric resins

- Dogbone_ASTM_D638V_31percentscale.SLDPRT: Design for the tensile test dogbones for the rigid resins

- resolutionCoupon_merged.SLDPRT: Design for the samples printed to assess print fidelity



### Directory: [Print Fidelity](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/tree/main/Supplemental/Print%20Fidelity)
This directory houses the raw data and the subsequent analyses for the print fidelity assessment. Images of the print fidelity resolution coupons can be found within the BiocompatibilityResolutionCoupons. These images were assessed manually via a Google Form (Print Fidelity Scoring for supplemental to paper - Google Forms.pdf). The raw data for the resulting assessment can be found in the Excel spreadsheet. The statistical analyses can be found in the Minitab files. Personally identifiable information has been recoded into a numerical value.

#### File List

- AppendixC_PrintFidelity_BoxPlots_mixed model (1)_ASL (2) ASL_20240812_withBar.mpx: Minitab data analysis file for generating box plots for all of the print fidelity manual assessment data

- print fidelity form data.xlsx: Raw data from the manual assessment of print fidelity using the Google Form survey

- Print Fidelity mixed effects model.mpx: Minitab data analysis file for assessing the manual assessment of print fidelity and generation of heatmaps

- Print Fidelity Scoring for supplemental to paper - Google Forms.pdf: PDF copy of the Google Form used as part of the manual assessment of print fidelity

- Sub-Directory: [BiocompatibilityResolutionCoupons](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/tree/main/Supplemental/Print%20Fidelity/BiocompatibilityResolutionCoupons)

  - Original microscope images of each print fidelity resolution coupon used as part of the manual assessment of print fidelity

  - Overview of structure: BiocompatibilityResolutionCoupons/[RESIN]/[STERILIZATION]/{images}

    - Sub-Directory: <RESIN>: This sub-directory houses more directories, each labeled for the specific resin sample images

      - Sub-Directory: <STERILIZATION>: Within each resin directory, there are sub-directories for each sterilization condition

        - Within the sterilization folder are the images of the samples. Each image is of a different, distinct sample



### Directory: [Mechanical Testing Data and Modeling](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling)
The Excel document in this folder contains the mechanical data for all tensile and compressive tests. Excess data after specimens failed are removed in the Compession_Rem and Tensile_Rem tabs for the compression and tensile tests respectively. The MATLAB files contain all of the analysis code used to process the force-length data into stress-strain data, fit the mechanical models according to the methods laid out in the main manuscript, generate the mechanical analysis results figures, and calculate the numerical metrics reported in the manuscript. The .tgn file is used to generate the table of Yeoh model parameters (Table 10). Finally, the MATLAB Outputs folder contains the generated outputs of the MATLAB scripts. Within these folders, the two .mat files contain the results of the mechanical model fits, and the Excel file contains the compile scalar metrics and Yeoh model parameters. The Output Figures folder contains the rendered data plots generated by the mechanical analysis MATLAB scripts.

#### File List

- masterRawData_v08.xlsx: Contains all tensile and compression testing data for both rigid and elastomeric samples

- MaterialModelFitting_MJB.m: MATLAB script that conducts all mechanical analysis and modeling fitting reported in the manuscript. 

- Bootstrap_Model_Fit.m: MATLAB function that conducts the bootstrapping parameter estimation described in the elastomeric resin modeling section

- RoundToSigFigs.m: MATLAB helper function that rounds numerical values to an appropriate number of significant figures

- YeohModelParameters.tgn: Latex Table file to generate Table 10

- Sub-Directory: [MATLAB Outputs](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/MATLAB%20Outputs)

  - Sub-Directory: [FinalAnalysis](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/MATLAB%20Outputs/FinalAnalysis) (regenerated data reporting the same data as [20240523](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/MATLAB%20Outputs/20240523))

    - Sub-Directory: Output Figures: All data plots generated from the mechanical analysis

    - Material Properties.xlsx: Calculated scalar metrics and model parameters for all tested samples. Data are labeled for conducting statistical analysis in MiniTab

    - RigidMaterialsModels.mat: MATLAB struct containing the mechanical models and scalar metrics for the rigid resin samples

    - SoftMaterialModels.mat: MATLAB struct containing the mechanical models and scalar metrics for the elastomeric resin samples

