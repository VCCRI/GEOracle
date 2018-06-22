# GEOracle: Mining the GEO database for perturbation data sets

#############################
 What is GEOracle?
#############################

GEOracle is a R Shiny app that greatly speeds up the identification and processing of large numbers of perturbation microarray gene expression data sets from GEO. It uses text mining of the GEO metadata along with machine learning techniques to automatically:

    Identify perturbation GSEs
    Cluster GSM samples
    Label clusters as control or perturbation
    Pair each perturbation cluster with it's control.
    Identify the molecule that is being perturbed and the direction of perturbation.


It provides this information via an interactive interface that allows the user to verify and change the details of the perturbation experiments. The work flow is described below.

#############################
Installing GEOracle
#############################

    Make sure you have an up to date version of R (3.3+) and the Shiny package installed.
    Use the following R commands to install georacle onto your computer:
    
    install.packages("devtools")
    library(devtools)
    devtools::install_github("VCCRI/Georacle")
    
    You can run georacle using the following commands:
    Georacle::Georacle()



GEOracle will then attempt to install all neccessary packages and download the full database of metadata from GEO (~5GB of space required). Once complete GEOracle will open as a Shiny app in a web browser or as a new R studio window. The user can then upload a list of GSE IDs and begin the analysis.



#############################
GEOracle analysis steps
#############################

    The user uploads a list of GSE IDs (simple .txt file, one GSE ID per line)
    The metadata for each GSE will then be analysed and clustering / labelling / matching will occur.
    The user decides whether to use the default filters for which GSEs to process, or manually set desired filters, then click "COMPUTE".
    GEOracle will now try to detect the names of the perturbed molecules and direction of perturbation for each GSE that passes filtering.
    Once the GSE IDs appear in the "Processed GSEs" table, the user clicks them one by one to modify / verify / remove the comparisons.
    Once the user is happy with the details of every GSE, they click "NEXT STEP" and enter a desired output directory.

GEOracle will then automatically perform differential expression using the limma method as implemented in NCBI's GEO2R. As well as differential expression results GEOracle will output a list of edges that can be used to instantly build a causal gene regulatory network based on the set of input perturbation experiments. This edge list can be loaded straight into Cytoscape or R's igraph package.


If you are having trouble installing or running GEOracle, please follow the below troubleshooting tips

    GEOracle works best on a relatively clean install of R (3.3+). You can perform a clean install of the latest version of R without removing your old versions.
    The most common problems occur from trying to use GEOracle with older versions of R (<3.3). Because GEOracle uses the latest features from Shiny, you will need to update your R and other core packages which Shiny depends on.
    Some users also encounter problems installing or updating packages due to permissions. Package install problems occur with R fairly frequently and are best managed by doing a clean install of the latest R with full administrator privilages.
    Some users report memory / lag issues after extended periods of time using the tool. If you have a long list of GSE IDs try to split them up into chunks of ~50 and process them individually. This will improve efficiency and make sure you don't lose any information in the case of a crash (R or web browsser).
    There are many cases where GSEs will not be represented in the final output. GEOracle will automatically remove GSEs which do not cluster confidently, or where a matched perturbation and control cluster can not be identified. Similarly many GSEs will not show significant results when using a standard differential expression approach. If you think some GSEs that should be present are missing, check the filters that were used and check for significant results using GEO2R on the GEO website for that particular GSE.
