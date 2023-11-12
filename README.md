# rewilding2021
The data we share here go with Downie et al. (2023), provisionally accepted at <i>Science Advances</i>.  They also are some of the data belonging with Oyesola et al. (2023); we have not yet provided the full data for this latter manuscript, as it is currently being revised.

The data are presented in comma-separated-value (CSV) files, except the microbiome data; these are in a tab-delimited file. The full description of the methods for producing these data and flow cytometry gating schemes are available in Downie et al. (2023), as well as Oyesola et al. (2023). Below we provide some description for each file. All NA entries in data are places where the quantity either was not measured or was not known; if you have particular questions about these, please contact us.

If you have particular questions about the software underlying the RFID data collection, please e-mail Jesse Saunders at 'jbs' at 'princeton.edu'

## Details on contents of data files
mice.metadata.2021.block1.GxE.csv: This file contains the metadata for each mouse used in the rewilding experiment. This includes mouse ID, genotype, the enclosure holding it when rewilded, the cage it was raised in, parentage (when known – on occasion this is not known due to issues with maintaining stable mouse IDs prior to release), whether the mouse was challenged with Trichuris muris, and information about its RFID and capture status during the experiment.

mice.metadata.2021.block2.GxE.csv: Same as above, but for the second block (repetition) of the experiment.

2021.block1.clean.GxE.csv: This file contains the RFID check-in information for each mouse across the duration of the rewilding period (five weeks). The columns are wedge (reader number), timestamp (to the second), RFID (the RFID recorded by the reader), and mouse.

2021.block2.clean.GxE.csv: Same as above, but for the second block (repetition) of the experiment.

2021.plasma.cyt.csv: This file contains the concentrations of different cytokines in the plasma at the time of sacrifice, as well as some relevant pieces of mouse metadata taken from the metadata files. The cytokines are IFNg, IL5, TNFa, IL6, IL17A, and IL22.

2021.MLN.cytokine.prod.csv: This file contains the concentrations of different cytokines after stimulation of mesenteric lymph node cells with various antigens in vitro to provoke cytokine responses. The columns are organized by antigen: each cytokine's production in response to PBS (the control), followed by CD3/CD28, LPS, Candida albicans, Clostridium perfringens, Bacteroides vulgatus, and Trichuris muris. See the manuscript methods for the precise details of the stimulation assay and the full list of the cytokines assayed.

2021.CBC.data.compliant.csv: This file contains the readouts of complete blood count analyses for each mouse at different time points: prior to release into enclosures, two weeks after release, and five week after release (at trapout). Not all timepoints are available for all mice due to mice not being captured at particular trapping sessions and trouble maintaining stable mouse IDs prior to release. Metadata from the metadata files is included, as well as age of mouse in weeks at time of sampling. Both total counts of different WBC types and percentages are included, along with data on RBC properties (unanalyzed in the manuscript). Details on method and machine for CBC analysis are in the manuscript methods.

2021.blood.flow.analysis.csv: This file contains flow cytometry results for analyses of lymphocyte (CD4 and CD8 T cell) populations in the peripheral blood, along with memory phenotypes (based on CD44 and CD62L expression). The data include some metadata drawn from the above sheets as well as both counts of cells and percentages (which is which is clearly marked in column names). Percentages for CD44 and CD62L expression phenotypes are given as percentages of the total cell subtype (CD4 T, CD8 T, or B220+ B cells). For details on flow cytometry procedures, gating, and cell type classification, please see Downie et al. (2023) and Oyesola et al. (2023). Naïve cells are CD62L+, CD44-; effector memory cells are CD62L+, CD44+; central memory cells are CD62L-, CD44+; double negative cells are CD62L-, CD44-.

2021.MLN.flow.analysis.csv: This file is much the same as the above, but for cells in the mesenteric lymph nodes instead of the peripheral blood. There are also a few additional cell types, including B220+ B cells (with CD44 and CD62L expression), NK cells, gd T cells, CD45+ cells, and CD45+ ILCs. All values are given strictly as percentages; percentages for CD44 and CD62L expression phenotypes are given as percentages of the total cell subtype (CD4 T, CD8 T, or B220+ B cells).


## Details on R scripts and their use
The code is in six R scripts; the titles provide rough descriptions of the material in each script. We ran it on R v4.1.2 with RStudio v2023.03.0+386. We also employed several further packages at different stages:

    dplyr v1.0.10
    tidyr v1.2.1
    ggplot2 v3.4.0
    ggridges v0.5.3
    brms v2.16.3
    MetBrewer v0.2.0
    vegan v2.6-4

The scripts are organized in a sequence laid out below; in most cases the use of one script will depend on the output of the prior one, although not always.

CI.data.preprocessing.R: This script takes the RFID check-in data and mouse metadata and provides additional information for each check-in, such as the location of the RFID reader, the Julian night of the experiment, etc.

individual.behavior.calculation.R: This script uses the RFID check-in data (after updating via the previous script) and mouse metadata to identify various aspects of individual behavior from the check-ins, as shown mostly in Fig. 1 of the manuscript. This includes number of check-ins per night, roaming entropy, and minimum distance traveled per night.

social.association.calculation.R: This script uses the RFID check-in data (after updating from the first script) and mouse metadata to determine the pairwise strength of spatiotemporal association between mice. We use a variety of different ways of defining an association, and as a result the script is very long (we never developed an efficient way to do all the calculation for individual definitions from only a few commands, but in theory it's probably not too difficult to streamline). This script produces the data frame referred to as full.2021.assoc.info, a very important data frame for subsequent analyses.

immune.similarity.calculation.R: This script uses the various immune data spreadsheets, along with the mouse metadata to calculate for each pair of mice the similarity of their immune phenotypes, for various aspects of immune phenotype. These pairwise metrics are needed for calculating the relationship between pairwise spatiotemporal association and immune similarity. The calculations are slightly different for different aspects of immune defense. This script also contains code for calculating microbiome similarity.

statistical.analyses.R: This script uses the outputs from the three preceding scripts to explore statistical relationships between metrics of spatiotemporal association and immune similarity, as well as predictors of variation in spatiotemporal association and individual behavior. The various models are a mix of standard linear regressions for some quantities and beta regressions for other quantities, particularly the pairwise similarities.

results.plotting.R: This script contains the code necessary to reproduce the figures in the manuscript. It relies on the model outputs in the preceding script, as well as the output from the other scripts when plotting, for example, distributions of check-ins across the experiment or variation in behavior metrics between individuals.
