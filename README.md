# README for Dataset Described In:

**Resin Cytotoxic, Mechanical, and Print Fidelity Properties with a Low-Cost Printer for C2C12 Biohybrid Devicess**

AS Liao, K Dai, AB Irez, A Sun, MJ Bennington, S Schaffer, B Chopra, JM Seok, R Adams, YJ Zhang, VA Webster-Wood

2025

> [!NOTE]
> If using this dataset, please cite using the following DOI for this dataset:
> 
> [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14675757.svg)](https://doi.org/10.5281/zenodo.14675757)
>
> Liao, A. S., Dai, K., irez, A. B., Sun, A., Bennington, M. J., Schaffer, S., Chopra, B., Seok, J. M., Adams, R., Zhang, Y. J., & Webster-Wood, V. (2025). CMU-BORG/Commercial-Resins-for-Biohybrids-SciData [Data set]. Zenodo. https://doi.org/10.5281/zenodo.14675757

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

### Directory: [Cytotoxicity](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Cytotoxicity)
The Excel spreadsheets in these files are the raw data exported from the plate reader with some organization in the subsequent sheets. The Excel spreadsheets each have 2 sheets: 
- RawData
  - This sheet is directly exported from the Agilent Synergy H1 plate reader and the relevant metadata for the plate reader.
- TableFormat
  - This takes the data exported from the plate reader (in 'RawData') and reformats it into a table format with the metadata associated with the raw data

*<a name="Table-1"></a>Table 1: Columns and expected values for the 'TableFormat' sheet in cytotoxicity raw data spreadsheets ([CytotoxicityPlateReader_allExceptFL.xlsx](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Cytotoxicity/CytotoxicityPlateReader_allExceptFL.xlsx) and [CytotoxicityPlateReader_FollowUp_asl_20240402_resin_72h_followUpFL.xlsx](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Cytotoxicity/CytotoxicityPlateReader_FollowUp_asl_20240402_resin_72h_followUpFL.xlsx))*
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

- [CytotoxicityPlateReader_allExceptFL.xlsx](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Cytotoxicity/CytotoxicityPlateReader_allExceptFL.xlsx):
  - Raw plate reader data for the initial study, which included all resin samples except for the follow-up cytotoxicity analysis for Formlabs Silicone 40A (IPA/BuOAc)
  - See [Table 1](#Table-1) for descriptions on the 'TableFormat' sheet

- [CytotoxicityPlateReader_FollowUp_asl_20240402_resin_72h_followUpFL.xlsx](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Cytotoxicity/CytotoxicityPlateReader_FollowUp_asl_20240402_resin_72h_followUpFL.xlsx):
  - Raw plate reader data for the follow-up cytotoxicity analysis for Formlabs Silicone 40A (IPA/BuOAc)
  - See [Table 1](#Table-1) for descriptions on the 'TableFormat' sheet


### Directory: [CAD Files for 3D Printing Samples](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/CAD%20Files%20for%203D%20Printing%20Samples)
This directory includes the SolidWorks part files that house the designs for the samples used in the study.

#### File List

- [96wellspecimen.SLDPRT](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/CAD%20Files%20for%203D%20Printing%20Samples/96wellspecimen.SLDPRT)
  - Design for the cytotoxicity analysis samples - these printed parts were cultured with the cells in a 96 well plate

- [Compression_ASTM_D575_40percentscale.SLDPRT](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/CAD%20Files%20for%203D%20Printing%20Samples/Compression_ASTM_D575_40percentscale.SLDPRT)
  - Design for the compression test sample pucks for the elastomeric resins

- [Compression_ASTM_D695_20percentscale.SLDPRT](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/CAD%20Files%20for%203D%20Printing%20Samples/Compression_ASTM_D695_20percentscale.SLDPRT)
  - Design for the compression test sample pucks for the rigid resins

- [Dogbone_ASTM_D412C_halfsize.SLDPRT](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/CAD%20Files%20for%203D%20Printing%20Samples/Dogbone_ASTM_D412C_halfsize.SLDPRT)
  - Design for the tensile test dogbones for the elastomeric resins

- [Dogbone_ASTM_D638V_31percentscale.SLDPRT](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/CAD%20Files%20for%203D%20Printing%20Samples/Dogbone_ASTM_D412C_halfsize.SLDPRT)
  - Design for the tensile test dogbones for the rigid resins

- [resolutionCoupon_merged.SLDPRT](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/CAD%20Files%20for%203D%20Printing%20Samples/resolutionCoupon_merged.SLDPRT)
  - Design for the samples printed to assess print fidelity



### Directory: [Print Fidelity](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Print%20Fidelity)
This directory houses the raw data and the subsequent analyses for the print fidelity assessment. Images of the print fidelity resolution coupons can be found within the BiocompatibilityResolutionCoupons. These images were assessed manually via a Google Form (Print Fidelity Scoring for supplemental to paper - Google Forms.pdf). The raw data for the resulting assessment can be found in the Excel spreadsheet. Personally identifiable information has been recoded into a numerical value.

#### File List

- [printFidelity_rawData.xlsx](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Print%20Fidelity/printFidelity_rawData.xlsx)
  - Raw data from the manual assessment of print fidelity using the Google Form survey
  - The column names are associated with the questions in the questionnaire used for the manual scoring ([Print Fidelity Scoring for supplemental to paper - Google Forms.pdf](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Print%20Fidelity/Print%20Fidelity%20Scoring%20for%20supplemental%20to%20paper%20-%20Google%20Forms.pdf))
  - The 'Email' collected from the questionnaire was used as a scorer ID. In the raw data spreadsheet, it has been recoded as a number to protect scorer privacy (Column B)

- [Print Fidelity Scoring for supplemental to paper - Google Forms.pdf](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Print%20Fidelity/Print%20Fidelity%20Scoring%20for%20supplemental%20to%20paper%20-%20Google%20Forms.pdf)
  - PDF copy of the Google Form used as part of the manual assessment of print fidelity
  - For the resin type and sterilization method questions (Questions 1 and 2), these match the subdirectory names (see [Sub-Directory: BiocompatibilityResolutionCoupons](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData#sub-directory-biocompatibilityresolutioncoupons) for details)
  - The plate ID and the sample ID (Questions 4 and 5) are encoded in the file name of the microscope image
- [PrintFidelity_MEMonly.mpx](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/707753de7882992734af89946f87900972d393c4/Supplemental/Print%20Fidelity/PrintFidelity_MEMonly.mpx)
  - Mixed effects model for each scored feature from the manual assessment
  - The Worksheet matches the format in [printFidelity_rawData.xlsx](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Print%20Fidelity/printFidelity_rawData.xlsx), with the following adjustments:
    - Shortened column names instead of the full question
    - The responses for the rows containing the largest and smallest features (Columns C9 and C10) were recoded from the original text responses (A-G) to numeric (1-7) for the analysis. These recoded columns were added to the end of the Worksheet in C17 and C18

#### Sub-Directory: [BiocompatibilityResolutionCoupons](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Print%20Fidelity/BiocompatibilityResolutionCoupons)

- Original microscope images of each print fidelity resolution coupon used as part of the manual assessment of print fidelity

- Overview of structure: BiocompatibilityResolutionCoupons/[RESIN]/[STERILIZATION]/{images}

  - Sub-Directory: [RESIN]: This sub-directory houses more directories, each labeled for the specific resin sample images (See [Table 2](#Table-2) below)

    - Sub-Directory: [STERILIZATION]: Within each resin directory, there are sub-directories for each sterilization condition (See [Table 3](#Table-3) below)

      - {images} Within the sterilization folder are the images (*.png) of the samples. Each image is of a different, distinct sample and labeled as follows: [resinType\]_[sterilizationType\]_Plate[pID]\_[sampleID].png.
        - The resin type [resinType] is in a shortened format. Refer to the [RESIN] directory housing the image for the longer name
        - The sterilization type [sterilizationType] is in a shortened format. Refer to the [STERILIZATION] directory housing the image for the longer name
        - The plate ID [pID] corresponds to the printing round that sample is part of since multiple prints were needed to fabricate all of the samples in the study
        - The sample ID [sampleID] corresponds to the individual sample for tracking purposes
        - For example, for "3DIPA_Au_Plate1_2.png" (Located [here](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Print%20Fidelity/BiocompatibilityResolutionCoupons/3DresynIPA/Autoclave/3DIPA_Au_Plate1_2.png)):
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



### Directory: [Mechanical Testing Data and Modeling](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling)
The Excel document in this folder ("masterRawData_clean.xlsx") contains the mechanical data for all tensile and compressive tests. Excess data after specimens failed are removed in the Compession_Rem and Tensile_Rem tabs for the compression and tensile tests respectively (see [Table 3](#Table-3) for a full description). The MATLAB files contain all of the analysis code used to process the force-length data into stress-strain data, fit the mechanical models according to the methods laid out in the main manuscript, generate the mechanical analysis results figures, and calculate the numerical metrics reported in the manuscript. The .tgn file is used to generate a table of the Yeoh model parameters (table not shown in manuscript). The MATLAB Outputs folder contains the generated outputs of the MATLAB scripts. Within this folder, the two .mat files contain the results of the mechanical model fits, and the Excel file contains the compile scalar metrics and Yeoh model parameters. The Output Figures folder contains the rendered data plots generated by the mechanical analysis MATLAB scripts. Finally, the Beam Testing folder contains the data and MATLAB scripts to perform the mechanical analysis on cantilevered beam that was used to investigate the effects of PBS soaking on the material properties of the rigid resins. Within this folder is a set of MATLAB files that are used to solve the beam model equations referenced in the manuscript and fit the Young's modulus to the deflected beam geometries. The SummaryComparisons svg file shows a summary figure reporting the same data from Table 6. Additionally, there are two sub-folders: the Test Beam Geometry folder contains the files required to generate and print the beams that were tested, and the Test Beam Images folder contains the images taken during the experiment, the superimposed composite image used to digitize data, a .tar file containing the outputs from the WebPlotDigitizer session, and a csv file with the geometry data used to fit the beam model. These files are organized first by resin and then by sample ID.

To reproduce the model fitting (Yeoh model and Hookean model) for the tensile and compression tests performed here, you will need to run the "MaterialModelFitting_MJB.m" MATLAB script. Before doing so, create a new folder in the "MATLAB Outputs" folder and inside, create a folder called "Output Figures" as these will be used by the script. Then, within the script on line 7, change the variable "saveDir" to be equal to your new folder name. This is required because the script to set up to not overwrite the "Material Properties.xlsx" Excel file in an existing output folder, and the script will through an error if that file already exists in the folder to which you are trying to save data.

#### File List

- [masterRawData_clean.xlsx](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/masterRawData_clean.xlsx)
	- Contains all tensile and compression testing data for both rigid and elastomeric samples
 	- See [Table 3](#Table-3) for a full description of the contents

- [MaterialModelFitting_MJB.m](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/MaterialModelFitting_MJB.m):
	- MATLAB script that conducts all mechanical analysis and modeling fitting reported in the manuscript. 

- [Bootstrap_Model_Fit.m](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/Bootstrap_Model_Fit.m):
	- MATLAB function that conducts the bootstrapping parameter estimation described in the elastomeric resin modeling section

- [RoundToSigFigs.m](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/RoundToSigFigs.m):
	- MATLAB helper function that rounds numerical values to an appropriate number of significant figures

- Sub-Directory: [MATLAB Outputs](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/MATLAB%20Outputs)

  - Sub-Directory: [FinalAnalysis](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/MATLAB%20Outputs/FinalAnalysis) (regenerated data reporting the same data as [20240523](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/MATLAB%20Outputs/20240523))

    - Sub-Directory: Output Figures:
    	- All data plots generated from the mechanical analysis

    - Material Properties.xlsx:
    	- Calculated scalar metrics and model parameters for all tested samples.

    - RigidMaterialsModels.mat:
    	- MATLAB struct containing the mechanical models and scalar metrics for the rigid resin samples

    - SoftMaterialModels.mat:
    	- MATLAB struct containing the mechanical models and scalar metrics for the elastomeric resin samples

- Sub-Directory: [Beam Testing](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/Beam%20Testing])
	
	- Sub-Directory: [Test Beam Geometry](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/Beam%20Testing/Test%20Beam%20Geometry)
		
		- CenterLine_s_0_5.txt:
  			- Contains the geometry of the beam centerline used to generate the beams in SolidWorks
		
		- TestBeam1_s_0_5_10mm_3mm_with_weight_attachment.SLDPRT:
  			- SolidWorks file containing the geometry of the test beams
		
		- TestBeam1_s_0_5_10mm_3mm_Phrozen.ctb:
  			- print file with settings for the Phrozen AquaGray 8K test beam
		
		- TestBeam1_s_0_5_10mm_3mm_BiomedClear_v2.ctb:
  			- print file with settings for the Liqcreate Bio-Med Clear test beam
		
	- Sub-Directory: [Test Beam Images](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/Beam%20Testing/Test%20Beam%20Images)
		
		All following sub-directories contain the same file for each of the different beams with the names updated to reflect the current sample. The full list of files is laid out for Sample 1 of the Phrozen Beams, but can be applied to all other samples.
		
		- Sub-Directory: [Phrozen Beams](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/Beam%20Testing/Test%20Beam%20Images/Phrozen%20Beams)
			
			- Sub-Directory: [Sample 1](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/Beam%20Testing/Test%20Beam%20Images/Phrozen%20Beams/Sample%201)
				
				- ModelFit.png:
    					- Final model fit overlaid on top of the composite image of all three beam conditions (0g, 100g, 200g)
				
				- Overlay.png:
    					- Composite image with all three beam conditions used to digitize the beam positions
				
				- Overlay.svg:
    					- Inkscape file used to generate the overlay
				
				- Phrozen_Sample_1.svg:
    					- Output from the MATLAB optimization script
				
				- Phrozen_Sample1.tar:
    					- compressed folder containing all files from the WebPlotDigitizer session used to obtain the beam positions
				
				- PhrozenSample1_Data.csv:
    					- csv file containing all of the fiducial marker data. The first two columns give the location of two points on the ruler used to set the scale of the image. The remaining columns show the position of the beam for the different loading conditions. The last point for each loading condition corresponds to the origin of the beam.
				
				- Sample1_NoLoad.jpg:
    					- Image of the beam before any load is applied
				
				- Sample1_100g.jpg:
    					- Image of the beam with a 100g load applied
				
				- Sample1_200g.jpg:
    					- Image of the beam with a 200g load applied
				
			- Sub-Directory: [Sample 2](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/Beam%20Testing/Test%20Beam%20Images/Phrozen%20Beams/Sample%202)
			
			- Sub-Directory: [Sample 3](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/Beam%20Testing/Test%20Beam%20Images/Phrozen%20Beams/Sample%202)
			
			- Sub-Directory: [Sample 4](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/Beam%20Testing/Test%20Beam%20Images/Phrozen%20Beams/Sample%204)
			
		- Sub-Directory: [Liqcreate Beams](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/Beam%20Testing/Test%20Beam%20Images/Liqcreate%20Beams)
		
			- Sub-Directory: [Sample 2](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/Beam%20Testing/Test%20Beam%20Images/Liqcreate%20Beams/Sample%202)
			
			- Sub-Directory: [Sample 3](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/Beam%20Testing/Test%20Beam%20Images/Liqcreate%20Beams/Sample%203)
			
			- Sub-Directory: [Sample 4](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/Beam%20Testing/Test%20Beam%20Images/Liqcreate%20Beams/Sample%204)
			
			- Sub-Directory: [Sample 5](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/tree/main/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/Beam%20Testing/Test%20Beam%20Images/Liqcreate%20Beams/Sample%205)
			
	- BeamModel.m:
 		- Function that returns the deformed configuration of a geometrically exact Euler-Bernoulli beam under a dead weight.
	
	- D1_matrix:
 		- Function that generates the finite difference matrix for the first derivative
	
	- FittingBeamModel.m:
 		- Controller script that runs the full beam experiment analysis. To regenerate the results of the beam analysis, only this script needs to be run. This scipt utilized the other functions in this folder.
	
	- Jacobian.m:
 		- Function that calculates the Jacobian of a vector-valued function using a finite difference approach.
	
	- QuasistaticForceBalance.m:
 		- Function that returns the residuals of the governing equations of the beam model
	
	- ReturnSubset.m:
 		- Helper function that returns a subset of a vector-valued function
	
	- SummaryComparisons.svg:
 		- Figure generated by the FittingBeamModel.m script

*<a name="Table-3"></a>Table 3: Sheets, columns and expected values for the [masterRawData_clean.xlsx](https://github.com/CMU-BORG/Commercial-Resins-for-Biohybrids-SciData/blob/f351759a86b12f52eb7a917f0f6e4eb815e601a2/Supplemental/Mechanical%20Testing%20Data%20and%20Modeling/masterRawData_clean.xlsx). The 'Compression' and 'Tensile' sheets contain all of the raw data that was collected. The 'Compression_Rem' and the 'Tensile_Rem' sheets contain the raw data with specific data points removed based on the exclusion criteria.*
<!--header-->
<table>
  <thead>
    <th> Sheet Name </th>
    <th> Column Name </th>
    <th> Possible Values </th>
    <th> Description </th>
  </thead>
  <tbody>
    <tr>
      <td rowspan=19> SpecificGeometries </td>
      	<td> LookUp Helper Column </td>
      	  <td> String </td>
	    <td> Combined information from 'Test Type', 'Sterilization', 'Resin', 'Plate ID', and 'Specimen ID' in a single column </td>
    </tr>
    <tr>
      	<td> Test Type </td>
      	  <td> 'Compression' or 'Tensile' </td>
	    <td> Type of mechanical test associated with the given specimen and geometries </td>	    
    </tr>
    <tr>
      	<td> Sterilization </td>
      	  <td>'NS', 'EtOH', or 'Autoclave' </td>
	    <td> Sterilization treatment for the given specimen and geometries, where 'NS' indicates no sterilization treatment, 'EtOH' indicates the ethanol treatment, and 'Autoclave' indicates the autoclave treatment </td>	    
    </tr>
    <tr>
      	<td rowspan=8> Resin </td>
      	  <td> FormLabs Silicone IPA </td>
	    <td> <a href="https://formlabs.com/store/materials/silicone-40a-resin/">FormLabs Silicone 40A</a> postprocessed with an 100% IPA solution  </td>
    </tr>
    <tr>
      	  <td> FormLabs Silicone Mix </td>
	    <td> <a href="https://formlabs.com/store/materials/silicone-40a-resin/">FormLabs Silicone 40A</a> postprocessed with an 80% IPA and 20% <i>n</i>-butyl acetate (BuOAc) solution  </td>
    </tr>
    <tr>
      <td> 3DResyn Bioflex UNW2 </td>
      <td> <a href="https://www.3dresyns.com/products/3dresyn-bioflex-a10-mb-monomer-based">3Dresyn Bioflex A10 MB Monomer Based</a> postprocessed with the manufacturer's proprietary wash solution, <a href="https://www.3dresyns.com/products/cleaning-fluid-unw2-bio-ultra-non-whitening-biocompatible-cleaner-with-medium-viscosity-for-ultra-gloss-and-transparency-of-prints">UNW2</a> </td>
    </tr>
    <tr>
      <td> 3DResyn Bioflex IPA </td>
      <td> <a href="https://www.3dresyns.com/products/3dresyn-bioflex-a10-mb-monomer-based">3Dresyn Bioflex A10 MB Monomer Based</a> postprocessed with 100% IPA solution </td>
    </tr>
    <tr>
      <td> Asiga DentaGuide </td>
      <td> <a href="https://www.asiga.com/materials-dental/">Asiga DentaGUIDE</a> </td>
    </tr>
    <tr>
      <td> Asiga DentaGum </td>
      <td> <a href="https://www.asiga.com/materials-dental/">Asiga DentaGUM</a> </td>
    </tr>
    <tr>
      <td> Liqcreate Biomed </td>
      <td> <a href="https://www.liqcreate.com/product/bio-med-clear-biocompatible-resin/">Liqcreate Bio-Med Clear</a> </td>
    </tr>
    <tr>
      <td> Phrozen AquaGray </td>
      <td> <a href="https://phrozen3d.com/products/aqua-8k-resin">Phrozen Aqua-Gray 8K</a> </td>
    </tr>
    <tr>
      	<td> Plate ID </td>
      	  <td> Numeric Integer (1-3) </td>
	    <td> ID to indicate which parts of a given resin were printed in the same printing session </td>	    
    </tr>
    <tr>
      	<td> Specimen ID </td>
      	  <td> Numeric Integer (1-3) or String </td>
	    <td> ID to indicate the specific specimen. A few specimens have strings with the number, which was to indicate experimental loading issues </td>	    
    </tr>
    <tr>
      	<td> Thickness or Gage Length (mm) </td>
      	  <td> Numeric </td>
	    <td> Thickness or gage length of the sample in millimeters, as measured with a caliper prior to mechanical testing </td>	    
    </tr>
    <tr>
      	<td> Diameter or Tensile CSA Thickness (mm) </td>
      	  <td> Numeric </td>
	    <td> Diameter of compression sample puck or thickness cross-sectional area of tensile dogbone of the sample in millimeters, as measured with a caliper prior to mechanical testing </td>	    
    </tr>
    <tr>
      	<td> Tensile Width (mm) </td>
      	  <td> Numeric </td>
	    <td> If the sample was a dogbone specimen for tensile testing, this column has the width in millimeters, as measured with a caliper prior to mechanical testing </td>	    
    </tr>
    <tr>
      	<td> Area (mm2) </td>
      	  <td> Numeric </td>
	    <td> Cross-sectional area in square millimeters, as calculated based on measured values in the other columns </td>   
    </tr>
    <tr>
      	<td> Sample Broke (0 or 1) </td>
      	  <td> Numeric Integer (0-1) </td>
	    <td> Indicator on whether the sample broke (1) or remained intact (0) during the mechanical test </td>	    
    </tr>
    <tr>
      	<td> Use Sample (No for if Experimental Error) </td>
      	  <td> Numeric Integer (0-1) </td>
	    <td> Indicator on whether the sample was used in the computational mechanical model (1) or not (0) based on the exclusion criteria to account for experimental error </td>	    
    </tr>
    <tr>
      <td rowspan=8> 'Compression', 'Tensile', 'Compression_Rem', 'Tensile_Rem' </td>
      	<td> Resin </td>
      	  <td> String </td>
	    <td> Resin type (refer to description for the 'Resin' column in the 'SpecificGeometries' sheet above) </td>
    </tr>
    <tr>
      	<td> Type </td>
      	  <td> 'Soft' or 'Rigid' </td>
	    <td> Resin classification on whether it was considered 'Rigid' or 'Soft' (Elastomeric) </td>
    </tr>
    <tr>
      	<td> Sterilization </td>
      	  <td> 'NS', 'EtOH', 'Autoclave' </td>
	    <td> Sterilization treatment for the given specimen and geometries, where 'NS' indicates no sterilization treatment, 'EtOH' indicates the ethanol treatment, and 'Autoclave' indicates the autoclave treatment </td>
    </tr>
    <tr>
      	<td> Plate ID </td>
      	  <td> Numeric Integer (1-3) </td>
	    <td> ID to indicate which parts of a given resin were printed in the same printing session </td>	    
    </tr>
    <tr>
      	<td> Specimen ID </td>
      	  <td> Numeric Integer (1-3) or String </td>
	    <td> ID to indicate the specific specimen. A few specimens have strings with the number, which was to indicate experimental loading issues </td>	    
    </tr>
    <tr>
      	<td> Time (sec) </td>
      	  <td> Numeric </td>
	    <td> Time stamp, in seconds, during the experimental data collection </td>	    
    </tr>
    <tr>
      	<td> Crosshead (mm) </td>
      	  <td> Numeric </td>
	    <td> Crosshead position, in millimeters, during the experimental data collection associated with the corresponding time stamp in the same row </td>	    
    </tr>
    <tr>
      	<td> Load (N) </td>
      	  <td> Numeric </td>
	    <td> Load measured, in Newtons, during the experimental data collection associated with the corresponding time stamp in the same row </td>	    
    </tr>
  </tbody>
</table>

<!--/header-->
