# README for the Supplemental Files to:

“Cytotoxicity and Characterization of 3D-Printable Resins Using a Low-Cost Printer for Muscle-based Biohybrid Devices”

AS Liao, K Dai, AB Irez, A Sun, MJ Bennington, S Schaffer, B Chopra, JM Seok, R Adams, YJ Zhang, VA Webster-Wood

2024

> [!NOTE]
> If using this dataset, please cite using the following DOI for this dataset:
> 
> [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14014330.svg)](https://doi.org/10.5281/zenodo.14014330)
>
> Liao, A. S., Dai, K., irez, A. B., Sun, A., Bennington, M. J., Schaffer, S., Chopra, B., Seok, J. M., Adams, R., Zhang, Y. J., & Webster-Wood, V. (2024). CMU-BORG/Commercial-Resins-for-Biohybrids-2024: Cytotoxicity and Characterization of 3D-Printable Resins Using a Low-Cost Printer for Muscle-based Biohybrid Devices (supplemental-data) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.14014331

## Associated Work:
Preprint Available on bioRxiv: [![Static Badge](https://img.shields.io/badge/DOI-10.1101%2F2024.10.31.621115-blue?link=https%3A%2F%2Fdoi.org%2F10.1101%2F2024.10.31.621115)](https://doi.org/10.1101/2024.10.31.621115)

> Liao, A.S. et al. (2024). Cytotoxicity and Characterization of 3D-Printable Resins Using a Low-Cost Printer for Muscle-based Biohybrid Devices. bioRxiv. https://doi.org/10.1101/2024.10.31.621115

Living Machines 2024 Conference Paper: [![Static Badge](https://img.shields.io/badge/DOI-10.1007%2F978--3--031--72597--5__27-blue)](https://doi.org/10.1007/978-3-031-72597-5_27)

> Liao, A.S. et al. (2025). Biocompatibility of Asiga Dental Resins Using a Low-Cost Printer for Biohybrid Actuator Applications. In: Szczecinski, N.S., Webster-Wood, V., Tresch, M., Nourse, W.R.P., Mura, A., Quinn, R.D. (eds) Biomimetic and Biohybrid Systems. Living Machines 2024. Lecture Notes in Computer Science(), vol 14930. Springer, Cham. https://doi.org/10.1007/978-3-031-72597-5_27

## Overview:

These are the supplemental files associated with the article. The article encompasses 3 major investigations of multiple 3D-printable resins:

- Cytotoxicity

- Mechanical Testing

- Print Fidelity

Each directory encompasses files related to these studies.

### Directory: [Cytotoxicity](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/tree/main/Supplemental/Cytotoxicity)
The Excel spreadsheets in these files are the raw data exported from the plate reader with some organization in the subsequent sheets. The Minitab files encompass the statistical analyses and figure generation from the plate reader data.

The Excel spreadsheets each have 2 sheets: 
- RawData
  - This sheet is directly exported from the Agilent Synergy H1 plate reader and the relevant metadata for the plate reader.
- TableFormat
  - This takes the data exported from the plate reader (in 'RawData') and reformats it into a table format with the metadata associated with the raw data

*<a name="Table-1"></a>Table 1: Columns and expected values for the 'TableFormat' sheet in cytotoxicity raw data spreadsheets ([CytotoxicityPlateReader_allExceptFL.xlsx](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/blob/41b6656c9c1e08af4327d9a3535cf96991ccb853/Supplemental/Cytotoxicity/CytotoxicityPlateReader_allExceptFL.xlsx) and [CytotoxicityPlateReader_FollowUp_asl_20240402_resin_72h_followUpFL.xlsx](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/blob/41b6656c9c1e08af4327d9a3535cf96991ccb853/Supplemental/Cytotoxicity/CytotoxicityPlateReader_FollowUp_asl_20240402_resin_72h_followUpFL.xlsx))*
<!--header-->
<table>
  <thead>
    <th> Column Name </th>
    <th> Possible Values </th>
    <th> Description </th>
  </thead>
  <tbody>
    <tr>
      <td> Plate ID </td>
      <td> Numeric Integer (1-3) </td>
      <td> ID used by plate reader to determine which plate was read </td>
    </tr>
    <tr>
      <td> Plate Row </td>
      <td> String (A-H) </td>
      <td> Row ID used to identify which cell the specific reading was from </td>
    </tr>
    <tr>
      <td> Plate Column </td>
      <td> Numeric Integer (1-12) </td>
      <td> Column ID of the plate used to identify which cell the specific reading was from </td>
    </tr>
    <tr>
      <td> Excel Row ID </td>
      <td> Numeric Integer </td>
      <td> Row ID of the Excel sheet in 'RawData' used to identify which cell the specific reading was from for the Calcein reading </td>
    </tr>
    <tr>
      <td> Excel Column ID </td>
      <td> String </td>
      <td> Column ID of the Excel sheet in 'RawData' used to identify which cell the specific reading was from </td>
    </tr>
    <tr>
      <td rowspan=15> Resin </td>
      <td> Live - EthD Only </td>
      <td> Negative Control Well (No treatment to cells) - Ethidium Homodimer-1 Dye Only </td>
    </tr>
    <tr>
      <td> Live - Cal Only </td>
      <td> Negative Control Well (No treatment to cells) - Calcein AM Dye Only </td>
    </tr>
    <tr>
      <td> Dead - EthD Only </td>
      <td> Positive Control Well (70% Ethanol treatment to cells, no resin) - Ethidium Homodimer-1 Dye Only </td>
    </tr>
    <tr>
      <td> Dead - Cal Only </td>
      <td> Positive Control Well (70% Ethanol treatment to cells, no resin) - Calcein AM Dye Only </td>
    </tr>
    <tr>
      <td> Blank </td>
      <td> Empty well </td>
    </tr>
    <tr>
      <td> Dye Only </td>
      <td> No cells, no resin samples, only Ethidium Homodimer-1/Calcein AM solution </td>
    </tr>
    <tr>
      <td> PDMS </td>
      <td> Cells exposed to PDMS samples </td>
    </tr>
    <tr>
      <td> FormLabs Silicone 40A - IPA/BuOAc </td>
      <td> Cells exposed to <a href="https://formlabs.com/store/materials/silicone-40a-resin/">FormLabs Silicone 40A</a> postprocessed with an 80% IPA and 20% <i>n</i>-butyl acetate (BuOAc) solution  </td>
    </tr>
    <tr>
      <td> FormLabs Silicone 40A - IPA Only </td>
      <td> Cells exposed to <a href="https://formlabs.com/store/materials/silicone-40a-resin/">FormLabs Silicone 40A</a> postprocessed with 100% IPA </td>
    </tr>
    <tr>
      <td> 3D Resyn Bioflex A10 IPA Wash </td>
      <td> Cells exposed to <a href="https://www.3dresyns.com/products/3dresyn-bioflex-a10-mb-monomer-based">3Dresyn Bioflex A10 MB Monomer Based</a> postprocessed with 100% IPA </td>
    </tr>
    <tr>
      <td> 3D Resyn Bioflex A10 Proprietary Wash </td>
      <td> Cells exposed to <a href="https://www.3dresyns.com/products/3dresyn-bioflex-a10-mb-monomer-based">3Dresyn Bioflex A10 MB Monomer Based</a> postprocessed with the manufacturer's proprietary wash solution, <a href="https://www.3dresyns.com/products/cleaning-fluid-unw2-bio-ultra-non-whitening-biocompatible-cleaner-with-medium-viscosity-for-ultra-gloss-and-transparency-of-prints">UNW2</a> </td>
    </tr>
    <tr>
      <td> Asiga DentaGuide </td>
      <td> Cells exposed to <a href="https://www.asiga.com/materials-dental/">Asiga DentaGUIDE</a> </td>
    </tr>
    <tr>
      <td> Asiga DentaGum </td>
      <td> Cells exposed to <a href="https://www.asiga.com/materials-dental/">Asiga DentaGUM</a> </td>
    </tr>
    <tr>
      <td> Liqcreate Biomed Clear </td>
      <td> Cells exposed to <a href="https://www.liqcreate.com/product/bio-med-clear-biocompatible-resin/">Liqcreate Bio-Med Clear</a> </td>
    </tr>
    <tr>
      <td> Phrozen AquaGray 8K </td>
      <td> Cells exposed to <a href="https://phrozen3d.com/products/aqua-8k-resin">Phrozen Aqua-Gray 8K</a> </td>
    </tr>
    <tr>
      <td rowspan=3> Sterilization </td>
      <td> N/A </td>
      <td> Sterilization type not applicable for this 'Resin' condition </td>
    </tr>
  <tr>
    <td> 70% Ethanol </td>
    <td> Resin samples were sterilized by submersion in 70% ethanol </td>
  </tr>
  <tr>
    <td> Autoclave </td>
    <td> Resin samples were sterilized using an autoclave </td>
  </tr>
  <tr>
    <td> Calcein Reading (485, 526) </td>
    <td> Numeric (Relative Fluorescence Units) </td>
    <td> Fluorescence reading at 485 nm excitation, 526 nm emission (for calcein AM detection) </td>
  </tr>
  <tr>
    <td> Tx Red Reading (584, 625) </td>
    <td> Numeric (Relative Fluorescence Units) </td>
    <td> Fluorescence reading at 584 nm excitation, 624 nm emission (for ethidium homodimer-1 detection) </td>
  </tr>
  </tbody>
</table>

<!--/header-->

#### File List

- [CytotoxicityPlateReader_allExceptFL.xlsx](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/blob/41b6656c9c1e08af4327d9a3535cf96991ccb853/Supplemental/Cytotoxicity/CytotoxicityPlateReader_allExceptFL.xlsx):
  - Raw plate reader data for the initial study, which included all resin samples except for the follow-up cytotoxicity analysis for Formlabs Silicone 40A (IPA/BuOAc)
  - See [Table 1](#Table-1) for descriptions on the 'TableFormat' sheet

- [CytotoxicityPlateReader_FollowUp_asl_20240402_resin_72h_followUpFL.xlsx](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/blob/41b6656c9c1e08af4327d9a3535cf96991ccb853/Supplemental/Cytotoxicity/CytotoxicityPlateReader_FollowUp_asl_20240402_resin_72h_followUpFL.xlsx):
  - Raw plate reader data for the follow-up cytotoxicity analysis for Formlabs Silicone 40A (IPA/BuOAc)
  - See [Table 1](#Table-1) for descriptions on the 'TableFormat' sheet

- CytotoxicityStats_allExceptFL_MINITAB_PLATEREADER_72HVWW_ASL_V06.mpx: Minitab data analysis file for the initial study, which included all resin samples except for the follow-up cytotoxicity analysis for Formlabs Silicone 40A

- CytotoxicityStats_followUp_72h.mpx: Minitab data analysis file for the follow-up cytotoxicity analysis for Formlabs Silicone 40A



### Directory: [CAD Files for 3D Printing Samples](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/tree/main/Supplemental/CAD%20Files%20for%203D%20Printing%20Samples)
This directory includes the SolidWorks part files that house the designs for the samples used in the study.

#### File List

- [96wellspecimen.SLDPRT](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/blob/41b6656c9c1e08af4327d9a3535cf96991ccb853/Supplemental/CAD%20Files%20for%203D%20Printing%20Samples/96wellspecimen.SLDPRT)
  - Design for the cytotoxicity analysis samples - these printed parts were cultured with the cells in a 96 well plate

- [Compression_ASTM_D575_40percentscale.SLDPRT](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/blob/41b6656c9c1e08af4327d9a3535cf96991ccb853/Supplemental/CAD%20Files%20for%203D%20Printing%20Samples/Compression_ASTM_D575_40percentscale.SLDPRT)
  - Design for the compression test sample pucks for the elastomeric resins

- [Compression_ASTM_D695_20percentscale.SLDPRT](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/blob/41b6656c9c1e08af4327d9a3535cf96991ccb853/Supplemental/CAD%20Files%20for%203D%20Printing%20Samples/Compression_ASTM_D695_20percentscale.SLDPRT)
  - Design for the compression test sample pucks for the rigid resins

- [Dogbone_ASTM_D412C_halfsize.SLDPRT](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/blob/41b6656c9c1e08af4327d9a3535cf96991ccb853/Supplemental/CAD%20Files%20for%203D%20Printing%20Samples/Dogbone_ASTM_D412C_halfsize.SLDPRT)
  - Design for the tensile test dogbones for the elastomeric resins

- [Dogbone_ASTM_D638V_31percentscale.SLDPRT](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/blob/41b6656c9c1e08af4327d9a3535cf96991ccb853/Supplemental/CAD%20Files%20for%203D%20Printing%20Samples/Dogbone_ASTM_D638V_31percentscale.SLDPRT)
  - Design for the tensile test dogbones for the rigid resins

- [resolutionCoupon_merged.SLDPRT](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/blob/41b6656c9c1e08af4327d9a3535cf96991ccb853/Supplemental/CAD%20Files%20for%203D%20Printing%20Samples/resolutionCoupon_merged.SLDPRT)
  - Design for the samples printed to assess print fidelity



### Directory: [Print Fidelity](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/tree/main/Supplemental/Print%20Fidelity)
This directory houses the raw data and the subsequent analyses for the print fidelity assessment. Images of the print fidelity resolution coupons can be found within the BiocompatibilityResolutionCoupons. These images were assessed manually via a Google Form (Print Fidelity Scoring for supplemental to paper - Google Forms.pdf). The raw data for the resulting assessment can be found in the Excel spreadsheet. The statistical analyses can be found in the Minitab files. Personally identifiable information has been recoded into a numerical value.

#### File List

- printFidelity_rawData.xlsx
  - Raw data from the manual assessment of print fidelity using the Google Form survey
  - The column names are associated with the questions in the questionnaire used for the manual scoring ([Print Fidelity Scoring for supplemental to paper - Google Forms.pdf](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/blob/41b6656c9c1e08af4327d9a3535cf96991ccb853/Supplemental/Print%20Fidelity/Print%20Fidelity%20Scoring%20for%20supplemental%20to%20paper%20-%20Google%20Forms.pdf))
  - The 'Email' collected from the questionnaire was used as a scorer ID. In the raw data spreadsheet, it has been recoded as a number to protect scorer privacy (Column B)

- AppendixC_PrintFidelity_BoxPlots_mixed model (1)_ASL (2) ASL_20240812_withBar.mpx
  - Minitab data analysis file for generating box plots for all of the print fidelity manual assessment data

- print fidelity form data.xlsx
  - Raw data from the manual assessment of print fidelity using the Google Form survey with the pairwise comparison results from the Minitab files

- Print Fidelity mixed effects model.mpx
  - Minitab data analysis file for assessing the manual assessment of print fidelity and generation of heatmaps

- [Print Fidelity Scoring for supplemental to paper - Google Forms.pdf](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/blob/41b6656c9c1e08af4327d9a3535cf96991ccb853/Supplemental/Print%20Fidelity/Print%20Fidelity%20Scoring%20for%20supplemental%20to%20paper%20-%20Google%20Forms.pdf)
  - PDF copy of the Google Form used as part of the manual assessment of print fidelity
  - For the resin type and sterilization method questions (Questions 1 and 2), these match the subdirectory names (see [Sub-Directory: BiocompatibilityResolutionCoupons](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024?tab=readme-ov-file#sub-directory-biocompatibilityresolutioncoupons) for details)
  - For the plate ID and the sample ID (Questions 4 and 5), 

#### Sub-Directory: [BiocompatibilityResolutionCoupons](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/tree/main/Supplemental/Print%20Fidelity/BiocompatibilityResolutionCoupons)

- Original microscope images of each print fidelity resolution coupon used as part of the manual assessment of print fidelity

- Overview of structure: BiocompatibilityResolutionCoupons/[RESIN]/[STERILIZATION]/{images}

  - Sub-Directory: [RESIN]: This sub-directory houses more directories, each labeled for the specific resin sample images (See [Table 2](#Table-2) below)

    - Sub-Directory: [STERILIZATION]: Within each resin directory, there are sub-directories for each sterilization condition (See [Table 3](#Table-3) below)

      - {images} Within the sterilization folder are the images (*.png) of the samples. Each image is of a different, distinct sample and labeled as follows: [resinType\]_[sterilizationType\]_Plate[pID]\_[sampleID].png.
        - The resin type [resinType] is in a shortened format. Refer to the [RESIN] directory housing the image for the longer name
        - The sterilization type [sterilizationType] is in a shortened format. Refer to the [STERILIZATION] directory housing the image for the longer name
        - The plate ID [pID] corresponds to the printing round that sample is part of since multiple prints were needed to fabricate all of the samples in the study
        - The sample ID [sampleID] corresponds to the individual sample for tracking purposes
        - For example, for "3DIPA_Au_Plate1_2.png" (Located [here](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-2024/blob/41b6656c9c1e08af4327d9a3535cf96991ccb853/Supplemental/Print%20Fidelity/BiocompatibilityResolutionCoupons/3DresynIPA/Autoclave/3DIPA_Au_Plate1_2.png)):
          - [resinType\]: "3DIPA" ([RESIN\] folder: "3DresynIPA")
          - [sterilizationType\]: "Au" ([STERILIZATION\] folder: "Autoclave")
          - [pID\]: 1
          - [sampleID\]: 2

*<a name="Table-2"></a>Table 2: Folder descriptions for the [RESIN] subdirectory*
| Folder Name  | Samples housed in the associated folder are composed of the following resin |
| ------------- | ------------- |
| 3DresynIPA | [3Dresyn Bioflex A10 MB Monomer Based](https://www.3dresyns.com/products/3dresyn-bioflex-a10-mb-monomer-based) postprocessed with 100% IPA  |
| 3DresynUNW2 | [3Dresyn Bioflex A10 MB Monomer Based](https://www.3dresyns.com/products/3dresyn-bioflex-a10-mb-monomer-based) postprocessed with the manufacturer's proprietary wash solution, [UNW2](https://www.3dresyns.com/products/cleaning-fluid-unw2-bio-ultra-non-whitening-biocompatible-cleaner-with-medium-viscosity-for-ultra-gloss-and-transparency-of-prints) |
| AsigaGuide | [Asiga DentaGUIDE](https://www.asiga.com/materials-dental/) |
| AsigaGum | [Asiga DentaGUM](https://www.asiga.com/materials-dental/) |
| FormlabsSiliconeIPA | [FormLabs Silicone 40A](https://formlabs.com/store/materials/silicone-40a-resin/) postprocessed with 100% IPA  |
| FormlabsSiliconMix | [FormLabs Silicone 40A](https://formlabs.com/store/materials/silicone-40a-resin/) postprocessed with an 80% IPA and 20% <i>n</i>-butyl acetate (BuOAc) solution  |
| LiqcreateBiomed | [Liqcreate Bio-Med Clear](https://www.liqcreate.com/product/bio-med-clear-biocompatible-resin/) |
| PhrozenAquaGrey | [Phrozen Aqua-Gray 8K](https://phrozen3d.com/products/aqua-8k-resin) |

*<a name="Table-3"></a>Table 3: Folder descriptions for the [STERILIZATION] subdirectory*
| Folder Name  | Images of samples housed in the associated folder underwent the following sterilization treatment prior to imaging|
| ------------- | ------------- |
| Autoclave | Autoclaving |
| Ethanol | Submerged in 70\% ethanol |
| NonSterile | No sterilization treatments prior to imaging |



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

