################################################################################################################################################
# Set up
################################################################################################################################################

#this.dir <- dirname(parent.frame(2)$ofile)
#this.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
this.dir <- getwd()
#this.dir <- choose.dir(caption = "Please select the GEOracle folder")
#this.dir <- "/home/rstudio/ShinyApps/GEOracle_AWS"

#setwd(this.dir)

list.of.cran.packages <- c("cluster", "plyr", "e1071", "rockchalk", "RCurl", "shiny", "DT", "memisc", "igraph", "shinyBS","RSQLite","shinyjs","httr")
new.cran.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.cran.packages)) install.packages(new.cran.packages, dependencies = TRUE)

list.of.bioc.packages <- c("GEOmetadb", "limma", "Biobase", "GEOquery", "biomaRt")
new.bioc.packages <- list.of.bioc.packages[!(list.of.bioc.packages %in% installed.packages()[,"Package"])]
if(length(new.bioc.packages)){
  source("https://bioconductor.org/biocLite.R")
  biocLite(pkgs = new.bioc.packages, ask = FALSE)
  
} 

require(cluster)
require(plyr)
require(e1071)
require(rockchalk)
require(limma)
require(RCurl)
require(Biobase)
require(GEOquery)
require(biomaRt)
require(RSQLite)
require(RMySQL)
require(httr)

require(GEOmetadb)
require(shiny)
require(DT)
require(memisc) # for case where
require(igraph) # for graph
require(shinyBS)
require(shinyjs)

# Features.file <- "Features.txt"
extdata_path <- system.file("extdata",package = "Georacle")
Features.file <- paste(extdata_path,"Features.txt", sep="/")

outFolder <- "GEOracle_output"

#load("ClusterLabelModel.RData")
ClusterLabelModelFile <- paste(extdata_path,"ClusterLabelModel.RData", sep="/")

load(ClusterLabelModelFile)

##
load(paste(extdata_path,"GSEClassifierData/newFeatsToKeep.RData", sep="/"))
load(paste(extdata_path,"GSEClassifierData/PertModel.RData",sep="/"))


################################################################################################################################################
# GEOmetaDB functions
################################################################################################################################################

###
# This function over-ride the getSQLiteFile function in Geometadb 
###
getSQLiteFile<-function (destdir = getwd(), destfile = "GEOmetadb.sqlite.gz") 
{
    print("Retriving SQLite file");
    localfile <- file.path(destdir, destfile)
    url_geo_1 = "http://starbuck1.s3.amazonaws.com/sradb/GEOmetadb.sqlite.gz"
    url_geo_2 = "http://gbnci-abcc.ncifcrf.gov/geo/GEOmetadb.sqlite.gz"
    url_geo_3 = "http://watson.nci.nih.gov/~zhujack/GEOmetadb.sqlite.gz"
    if (!inherits(try(url(url_geo_1, open = "rb"), silent = TRUE), 
        "try-error")) {
        url_geo = url_geo_1
    }
    else if (!inherits(try(url(url_geo_2, open = "rb"), silent = TRUE), 
        "try-error")) {
        url_geo = url_geo_2
    }
    else {
        url_geo = url_geo_3
    }
    download.file(url_geo, destfile = localfile, mode = "wb")
    cat("Unzipping...\n")
    gunzip(localfile, overwrite = TRUE)
    unzippedlocalfile <- gsub("[.]gz$", "", localfile)
    con <- dbConnect(SQLite(), unzippedlocalfile)
    dat <- dbGetQuery(con, "select * from metaInfo")
    dbDisconnect(con)
    cat("Metadata associate with downloaded file:\n")
    print(dat)
    return(unzippedlocalfile)
}



# This function downloads and loads the metadata
loadMetadata <- function(default=TRUE){
  if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
  metadata <- dbConnect(SQLite(),'GEOmetadb.sqlite')
  #metadata <- dbConnect(MySQL(),'GEOmetadb.sql')
  #metadata <- dbConnect(SQLite(),host='http://georacle.victorchang.edu.au/', dbname='GEOmetadb.sqlite')
  #metadata@dbname <- "http://georacle.victorchang.edu.au/GEOmetadb.sqlite"
  return(metadata)
}

# load the metadata
metadata <- loadMetadata()

########################
########################

# This function grabs all the metadata related to a GSE (experiment) ID, includng the associated GSM (sample) IDs
grabAllMetadataGSE <- function(GSEid) {
  # Fetch the GSE metadata and make it a list
  GSEMeta <- as.list(dbGetQuery(metadata,paste("select * from gse where gse.gse =\"", GSEid,"\"", "\n", sep="")))
  # Get the GSM IDs contained within this GSE
  GSEtoGSM<- grabGSMids(GSEid)
  # Get the GSM metadata
  GSMMeta <- lapply(GSEtoGSM, grabAllMetadataGSM)
  # Add the GSM metadata to the GSE metadata
  if(length(unlist(GSMMeta)) == 0){return(NA)}
  
  GSEMeta$GSMMeta <- GSMMeta
  # Return the GSE metadata list
  return(GSEMeta)
}

# This function grabs all metadata related to a GSM (sample) ID
grabAllMetadataGSM <- function (gsmid){
  dbGetQuery(metadata, paste("select * from gsm where gsm.gsm =\"", gsmid,"\"\n", sep=""))
}

# This function finds the GSM (sample) IDs associated with a GSE (experiment) ID 
grabGSMids <- function(GSEid){
  gsm <- unlist(dbGetQuery(metadata,paste("select gse_gsm.gsm from gse_gsm where gse_gsm.gse =\"", GSEid,"\"", "\n", sep="")))
  return(gsm)
}

# This function finds the GPL (platform) ID associated with a GSM (sample) ID
grabGPL <- function (gsmid){
  gplid <- unlist(dbGetQuery(metadata,paste("select gsm.gpl from gsm where gsm.gsm =\"", gsmid,"\"", "\n", sep="")))
  return(gplid)
}

## select metadata from multiple GSE ids, returns a data.frame instead of a list
grabMultipleGSE <- function (gseids){ 
  dbGetQuery(metadata, paste("select * from gse where gse IN (\'" ,paste(gseids, collapse="\',\'"), "\')", sep="")) 
}

## select metadata from multiple GSM ids, returns a data.frame instead of a list
grabMultipleGSM <- function (gsmids){
  dbGetQuery(metadata, paste("select * from gsm where gsm.gsm IN (\'" ,paste(gsmids, collapse="\',\'"), "\')", sep="")) 
}


################################################################################################################################################
# GEORACLE FUNCTIONS
################################################################################################################################################



get.title.table <- function(GSM){
  title <- GSM$title
  o.title <- title
  if(!grepl("(wt)|(ko)(_)?\\d*?$",title, ignore.case = T) & !grepl("((\\s|\\_|\\-)[0-9]{1,2}(e|h|w|d|p|(pd)|(day)|(week)|(ed)|(hr)|(hour))$)|((\\s\\_\\-)(e|h|w|d|p|(pd)|(day)|(week)|(ed)|(hr)|(hour))[0-9]{1,2}$)", title, ignore.case = T)){ 
    title <- unlist(strsplit(gsub("[!a-z!0-9!A-Z][0-9a-zA-Z]{,2}$","",title), "\\,|\\_"))
  }
  if(identical(title, character(0))){ title <- o.title}
  title <- unlist(strsplit(gsub("\\s|\\-","_",title), "\\,|\\_"))
  title <- title[!title == ""]
  return(title)
}

get.characteristics.table <- function(GSM){
  characteristics <- unlist(strsplit(GSM$characteristics_ch1, ";\t|\\, "))
  characteristics.split <- do.call(rbind, lapply(characteristics, function(X){return(unlist(strsplit(X, ":\\s?")))}))
  
  if(as.numeric(ncol(characteristics.split)) > 1) {
    characteristics.table <- as.data.frame(t(characteristics.split[,2]))
    colnames(characteristics.table) <- characteristics.split[,1]
  } else {
    characteristics.table <- as.data.frame(characteristics.split)
    if(nrow(characteristics.table)>1) {
      characteristics.table <- as.data.frame(t(characteristics.table))
    }
  }
  return(characteristics.table)
}

#Calculates combined table
get.combined.table <- function(GSM) {
  characteristicsTable <- get.characteristics.table(GSM)
  titleTable <- get.title.table(GSM)
  combined.table <- cbind(t(titleTable),characteristicsTable)
  return(combined.table)
}

###########################################################
###########################################################


get.source <- function(GSM) {
  source <- GSM$source_name_ch1
  return(source)
}


## takes in GSM metadata
get.combined.table2 <- function(GSM) {
  characteristicsTable <- get.characteristics.table(GSM)
  titleTable <- get.title.table(GSM)
  souceInfo <- get.source(GSM)
  combined.table <- cbind(t(titleTable),characteristicsTable,souceInfo)
  return(combined.table)
}

get.min.diss <- function(Table){
  
  Table <- data.matrix(Table)
  Table[is.na(Table)] <- "a"
  Table <- data.frame(Table)
  diss.mat <- data.matrix(daisy(Table))
  
  min.dis <- lapply(1:nrow(diss.mat), function(X){
    return(unlist(which(diss.mat[X,] == min(diss.mat[X,]))))
  })
  return(min.dis)
}

TitleClustering <- function(gseMetadata) {
  
  print(paste0("Clustering ", gseMetadata$gse))
  TitleTable <- data.frame(do.call(rbind, lapply(gseMetadata$GSMMeta, get.title.table)))
  
  GSMtitles <- unlist(lapply(gseMetadata$GSMMeta, function(X){return(X$title)}))
  names(GSMtitles) <- do.call(rbind, lapply(gseMetadata$GSMMeta, function(x) {x$gsm}))
  if(sum(duplicated(GSMtitles)) == 0){
    rownames(TitleTable) <- GSMtitles
  } else {
    rownames(TitleTable) <- names(GSMtitles)
  }
  # finding least dissimilarity = best match
  title.min.dis <- get.min.diss(TitleTable)
  
  tNamesInClusters <- lapply(unique(title.min.dis), function(X){return(GSMtitles[X])})
  
  return(tNamesInClusters)  
}

#Group GSMs based on characteristics
CharacteristicsClustering <- function(gseMetadata) {
  
  CharTable <- data.frame(do.call(rbind.fill, lapply(gseMetadata$GSMMeta, get.characteristics.table)))
  
  GSMtitles <- unlist(lapply(gseMetadata$GSMMeta, function(X){return(X$title)}))
  names(GSMtitles) <- do.call(rbind, lapply(gseMetadata$GSMMeta, function(x) {x$gsm}))
  if(sum(duplicated(GSMtitles)) == 0){
    rownames(CharTable) <- GSMtitles
  } else {
    rownames(CharTable) <- names(GSMtitles)
  }# finding least dissimilarity = best match
  char.min.dis <- get.min.diss(CharTable)
  
  cNamesInClusters <- lapply(unique(char.min.dis), function(X){return(GSMtitles[X])})
  return(cNamesInClusters)  
}

#Group GSMs based on both title and characteristics
CombinedClustering <- function(gseMetadata) {
  combinedTable <- data.frame(do.call(rbind.fill, lapply(gseMetadata$GSMMeta, get.combined.table2)))
  GSMtitles <- unlist(lapply(gseMetadata$GSMMeta, function(X){return(X$title)}))
  names(GSMtitles) <- do.call(rbind, lapply(gseMetadata$GSMMeta, function(x) {x$gsm}))
  combined.min.dis <- get.min.diss(combinedTable)
  fNamesInClusters <- lapply(unique(combined.min.dis), function(X){return(GSMtitles[X])})
  return(fNamesInClusters)  
}

#Check if clustering is valid 
informative.clustering <- function(CL){
  if(length(CL) > 1 & length(CL) <  (length(unlist(CL)) - 1) ){
    return("TRUE")
  }
  return("FALSE")
}

#Check for invalid titles which are purely numeric
non.numeric.title <- function(GSM) {
  title <- GSM$title
  title <- gsub("[!a-z!0-9!A-Z][0-9a-zA-Z]{,2}$","",title)
  title <- gsub("\\s|\\-","",title)
  if(grepl("^[0-9]+$", title)) {
    return("FALSE")
  }
  return("TRUE")
}

#Decide whether to use title clustering or characteristics clustering 
GEOclustering <- function(gseMetadata){
  
  title.res <- TitleClustering(gseMetadata)
  char.res <- CharacteristicsClustering(gseMetadata)
  combined.res <- CombinedClustering(gseMetadata)
  
  names(title.res) <- unlist(lapply(title.res, function(X){ return( paste(names(X), collapse="_"))}))
  names(char.res) <- unlist(lapply(char.res, function(X){ return( paste(names(X), collapse="_"))}))
  names(combined.res) <- unlist(lapply(combined.res, function(X){ return( paste(names(X), collapse="_"))}))
  
  tinfo <- informative.clustering(title.res)
  non.numeric <- unlist(lapply(gseMetadata$GSMMeta, non.numeric.title))
  if("FALSE" %in% non.numeric) {
    tvalid <- "FALSE"
  } else {
    tvalid <- "TRUE"
  }
  
  cinfo <- informative.clustering(char.res)
  
  if (tinfo == "TRUE" && cinfo == "TRUE") {
    #Check whether title and characteristics clusterings match
    factor.char.res <- lapply(char.res, as.factor)
    level.char.res <- lapply(factor.char.res, function(X){
      combineLevels(X,c(X[1:length(X)]), newLabel = toString(X[1]))})
    #Unlist cluster 
    unlisted.char <- as.factor(unlist(level.char.res))
    tClustersInCharClusters <- lapply(title.res, function(X){
      return(unique(unlisted.char[match(X, as.factor(unlist(char.res)))]))
    })
    
    maxNumberOfCharPerT <- max(unlist(lapply(tClustersInCharClusters, length)))
    
    maxNumberOfTPerChar <- max(table(unlist(tClustersInCharClusters)))
    
    
    if(maxNumberOfCharPerT == 1 && maxNumberOfTPerChar == 1) {
      #Matching title and characteristics clustering
      return(title.res)
    } else if(tvalid == "FALSE") {
      #Numeric titles 
      return(char.res)
    } else {
      #Title and characteristics clustering not matching
      factor.combined.res <- lapply(combined.res, as.factor)
      level.combined.res <- lapply(factor.combined.res, function(X){
        combineLevels(X,c(X[1:length(X)]), newLabel = toString(X[1]))})
      
      unlisted.combined <- as.factor(unlist(level.combined.res))
      tClustersInCombinedClusters <- lapply(title.res, function(X){
        return(unique(unlisted.combined[match(X, as.factor(unlist(combined.res)))]))
      })
      
      charClustersInCombinedClusters <- lapply(char.res, function(X){
        return(unique(unlisted.combined[match(X, as.factor(unlist(combined.res)))]))
      })
      
      maxNumberOfCombinedPerT <- max(unlist(lapply(tClustersInCombinedClusters, length)))
      
      maxNumberOfTPerCombined <- max(table(unlist(tClustersInCombinedClusters)))
      
      maxNumberOfCombinedPerChar <- max(unlist(lapply(charClustersInCombinedClusters, length)))
      
      maxNumberOfCharPerCombined <- max(table(unlist(charClustersInCombinedClusters)))
      
      if((maxNumberOfCombinedPerT == 1 && maxNumberOfTPerCombined == 1)|
         (maxNumberOfCombinedPerChar == 1 && maxNumberOfCharPerCombined == 1)){
        return(combined.res)
      } else {
        return(combined.res)
      }
    }
  } else if(tinfo == "TRUE" && cinfo == "FALSE") {
    #Invalid characteristics clustering
    return(title.res)
  } else if(tinfo == "FALSE" && cinfo == "TRUE") {
    #Invalid title clustering 
    return(char.res)
  } else {
    return("Discard - Neither informative")
  }
}

###########################################################
###########################################################

#Extract title features for GSM
get.title.features <- function(cluster) {
  title <- paste(cluster, collapse = " ")
  clusterTitleFeatures <- levels(factor(unlist(strsplit(gsub("\\,|:|\\.|;|\\_|\\)|\\("," ",title), "\\s{1,3}"))))
  return(clusterTitleFeatures)
}

#Gets characteristics features for GSE
get.characteristics.features <- function(GSEclusters, CurGSEmetadata) {
  #Get clusters in GSM codes
  clusterCode <- lapply(GSEclusters,names)
  
  GSMChar <- lapply(lapply(clusterCode, function(x){
    lapply(x, function(y){
      lapply(CurGSEmetadata$GSMMeta, function(z) {
        z$characteristics_ch1[z$gsm == y]
      })
    })
  }), unlist)
  
  clusterChar <- lapply(GSMChar, function(x){ 
    paste(x, collapse = ";")
  })
  
  #Tokenise keywords of characteristics 
  clusterCharFeatures <- lapply(clusterChar, function(x) {
    x2 <- strsplit(x, ";")[[1]]
    x3 <- gsub(".*:(.*)","\\1",x2)
    levels(factor(unlist(strsplit(gsub("\\,|:|\\.|;|\\_|\\)|\\("," ",x3), "\\s+"))))
  })
  
  return(clusterCharFeatures)
}

#Return list with info on cluster GSMS, title keywords and characteristic keywords 
combine.cluster.features <- function(GSEclusters, CurGSEmetadata) {
  
  #Get features and pu into string
  clusterChar <- lapply(get.characteristics.features(GSEclusters, CurGSEmetadata), function(x){ 
    paste(x, collapse = " ")
  })
  
  clusterTitle <- lapply(lapply(GSEclusters, get.title.features), function(x){ 
    paste(x, collapse = " ")
  })
  
  clusterGSM <- lapply(lapply(GSEclusters, names), function(x){ 
    paste(x, collapse = " ")
  })
  
  combinedTable <- cbind(clusterGSM, clusterTitle, clusterChar)
  
  return(combinedTable)
}

############################################

#Load control and treatment features
get.features <- function(txtFile) {
  featureTxt <- read.table(txtFile)
  features <- as.vector(unlist(featureTxt))
}

############################################

#Check presence of features
checkPresence <- function(GSEFeatures) {
  
  combined.features <- get.features(Features.file)
  
  feat.matrix <- do.call(rbind, lapply(GSEFeatures, function(cluster){
    feat.vector <- unlist(lapply(combined.features, function(feat){
      grepl(feat, cluster, ignore.case = T)
    }))
  }))
  
  colnames(feat.matrix) <- combined.features
  return(feat.matrix)
}

#Produce logic matrix to inidicate presence of features in title and characteristics
get.logic <- function(GSEclusters, CurGSEmetadata) {
  
  GSEFeatures <- combine.cluster.features(GSEclusters, CurGSEmetadata)
  
  #Check presence of features in title
  titleFeatures <- checkPresence(GSEFeatures[,2])
  
  #Check presence of features in characteristics
  charFeatures <- checkPresence(GSEFeatures[,3])
  
  titleCharMatrix <- cbind(titleFeatures, charFeatures)
  return(titleCharMatrix)
}

###############


#Obtain matrix
get.svm.m <- function(GSEclusters, CurGSEmetadata) {
  
  m <- get.logic(GSEclusters, CurGSEmetadata)
  
  #Edit column names
  colnames(m) <- gsub("/","neg",colnames(m))
  colnames(m) <- gsub("-","\\.",colnames(m))
  colnames(m) <- gsub("\\?|\\|","\\.",colnames(m))
  
  #Distinguish between characteristics and title features
  colnames(m)[(ncol(m)/2+1):ncol(m)] <- paste(colnames(m)[(ncol(m)/2+1):ncol(m)],".char", sep='')
  
  return(m)
}

#Obtain predictions through SVM
get.preds <- function(data, model) {
  
  #Define dependent variable
  d <- data.frame(data)
  
  predictions <- predict(model, d, probability = T)
  return(predictions)
}




LabelClusters <- function(GSEclusters, CurGSEmetadata, model){
  data <- get.svm.m(GSEclusters, CurGSEmetadata)
  prediction <- get.preds(data, model)
  
  return(prediction)
}

#Reverse label
reverse.label <- function(oldLabel) {
  newLabel <- c("FALSE", "TRUE")[c("FALSE", "TRUE")!= as.vector(oldLabel)]
  return(newLabel)
}

FixLabels <- function(GSEclusters, predictions, analysis_type = "All experiments"){
  
  
  GSEClassInfo <- data.frame(predLabel = predictions, probs = attr(predictions, "probabilities"), row.names = lapply(lapply(GSEclusters, names), paste, collapse=" "))
  
  reverseRow <- vector()
  assessClass <- vector()
  threshold <- 0.83
  
  #Check if number of classes
  if(length(unique(GSEClassInfo$predLabel)) == 1) {
    #If one class:
    
    #Select probability values to work with 
    if(unique(GSEClassInfo$predLabel) == "FALSE") {
      GSEClassInfo$probability <- GSEClassInfo$probs.FALSE
    } else{
      GSEClassInfo$probability <- GSEClassInfo$probs.TRUE
    }
    
    #Check if all probability values equal 
    if (!all(GSEClassInfo$probability == GSEClassInfo$probability[1])) {
      
      #Check for probability values above threshold
      highConf <- rownames(GSEClassInfo)[GSEClassInfo$probability >= threshold]
      
      #Label high confidence clusters
      if(length(highConf) > 0 & length(highConf) < dim(GSEClassInfo)[1]) {
        for(i in 1:length(highConf)){
          assessClass[i] <- "High"
          names(assessClass)[i] <- highConf[i]
        }
      } else if (length(highConf) == nrow(GSEClassInfo)) {
        #Clear high confidence vector if all cluster probabilities above threshold
        length(highConf) <- 0
      }
      
      #Reverse labels with lowest confidence 
      reverseRow <- rownames(GSEClassInfo)[which(GSEClassInfo$probability == min(GSEClassInfo$probability))]
      GSEClassInfo[reverseRow,]$predLabel <- sapply(GSEClassInfo[reverseRow,]$predLabel,reverse.label)
      
      for(i in 1:length(reverseRow)){
        assessClass[i+length(highConf)] <- "Low"
        names(assessClass)[i+length(highConf)] <- reverseRow[i]
      }
      
      #Label everything neither high nor low as medium
      mediumClusters <- rownames(GSEClassInfo)[which(!rownames(GSEClassInfo) %in% names(assessClass))]
      j <- length(assessClass)
      if (length(mediumClusters) != 0) {
        for(i in 1:length(mediumClusters)){
          assessClass[i+j] <- "Medium"
          names(assessClass)[i+j] <- mediumClusters[i]
        }
      }
    } else {
# actually, lets always choose one to randomly be the "control" if we cant detect it...     
#      if(analysis_type == "All experiments"){
      
          ## if we want to include all GSE, we will randomly permute one cluster to be different
        assessClass <- rep("Low", nrow(GSEClassInfo))
        randomCluster <- sample(1:nrow(GSEClassInfo), 1)
        GSEClassInfo$predLabel[randomCluster] <- setdiff(c(T,F),GSEClassInfo$predLabel[randomCluster])
        
#      }else{
        #All probabilities equal - invalid clusters 
#        assessClass <- rep("Invalid", nrow(GSEClassInfo))
#      }
      names(assessClass) <- rownames(GSEClassInfo)
    }
  } else {
    #If both classes present:
    #Check for clusters that are high confidence
    highConfMut <- rownames(GSEClassInfo)[GSEClassInfo$probs.FALSE >= threshold]
    highConfWT <- rownames(GSEClassInfo)[GSEClassInfo$probs.TRUE >= threshold]
    highConf <- c(highConfMut,highConfWT)
    
    #Label high confidence clusters
    if(length(highConf) > 0) {
      for(i in 1:length(highConf)){
        assessClass[i] <- "High"
        names(assessClass)[i] <- highConf[i]
      }
    } 
    
    #Label medium confidence clusters  
    mediumClusters <- rownames(GSEClassInfo)[which(!rownames(GSEClassInfo) %in% names(assessClass))]
    j <- length(assessClass)
    if (length(mediumClusters) != 0) {
      for(i in 1:length(mediumClusters)){
        assessClass[i+j] <- "Medium"
        names(assessClass)[i+j] <- mediumClusters[i]
      }
    }
  }
  
  GSEClassInfo$confidence = assessClass[match(rownames(GSEClassInfo), names(assessClass))]
  
  return(GSEClassInfo)
}



###################################################################
#Get GSM clusters and labels
get.cluster.labels <- function(GSEId,class.res) {
  info.table <- class.res[grepl(GSEId,rownames(class.res)),][,1:2]
  colnames(info.table) <- c("GSM", "Label")
  
  return(info.table)
}

#check if multi-control analysis required
check.multi <- function(clusterGSMLabel) {
  if((nrow(clusterGSMLabel) == 2) | (sum(clusterGSMLabel[,2] == "FALSE") == 1)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}
######################################################
#Get label vector describing class or group
get.labels <- function(labelledClusters) {
  splitGSM <- lapply(labelledClusters,function(x){unlist(strsplit(as.character(x)," "))})
  GSMLabels <- sub("^([[:alpha:]]*).*", "\\1", names(unlist(splitGSM)))
  names(GSMLabels) <- unlist(splitGSM) 
  return(GSMLabels)
}

#Get table with title, characteristics and source keywords plus class and subgroup labels
get.label.table <- function(clusterGSMLabel,gseMetadata) {
  #Get combined table
  combinedTable <- data.frame(do.call(rbind.fill, lapply(gseMetadata$GSMMeta, get.combined.table2)))
  rownames(combinedTable) <- do.call(rbind, lapply(gseMetadata$GSMMeta, function(x) {x$gsm}))
  
  #Generating class vector for GSMs
  mC <- as.data.frame(clusterGSMLabel)
  classSplit <- split(mC$GSM,mC$Label)
  classLabel <- get.labels(classSplit)
  
  #Generating subgroup vector for GSMs
  ms <- clusterGSMLabel[,1]
  names(ms) <- letters[1:length(ms)]
  subgroupLabel <- get.labels(ms)
  
  labelT <- cbind(combinedTable,ClassLabels = classLabel[rownames(combinedTable)], SubgroupLabels = subgroupLabel[rownames(combinedTable)])
  
  return(labelT)
}

#Get dissimilarity matrix of all GSMs
get.diss.mat <- function(Table){ #Put in combinedTable
  Table <- data.matrix(Table)
  Table[is.na(Table)] <- "a"
  Table <- data.frame(Table)
  diss.m <- data.matrix(daisy(Table))
  return(diss.m)
}

#Input cluster of GSM titles as string and output list of matched pair
pair.clusters <- function(GSMCluster,labelTable, CurGSEmetadata) {
  clusterGSM <- strsplit(GSMCluster," ")[[1]]
  clusterLabel <- unique(labelTable[GSMCluster,]$predLabel)
  if(clusterLabel == "FALSE") {
    return(NA)
  }
  
  combinedTable <- data.frame(do.call(rbind.fill, lapply(CurGSEmetadata$GSMMeta, get.combined.table2)))
  rownames(combinedTable) <- do.call(rbind, lapply(CurGSEmetadata$GSMMeta, function(x) {x$gsm}))
  #Obtain dissimilarity matrix of cluster with potential GSM matches 
  diss.mat <- get.diss.mat(combinedTable)
  
  diss.mat.cluster <- diss.mat[rownames(diss.mat) %in% clusterGSM,, drop=FALSE]
  
  
  #Obtain GSM of potential matches
  pMatch <- rownames(labelTable)[labelTable$predLabel != clusterLabel]
  names(pMatch) <- pMatch
  pDisss <- lapply(pMatch, function(CL){
    GSMs <- unlist(strsplit(CL, " "))
    return(diss.mat.cluster[,GSMs, drop = FALSE])
  })
  
  pMatchMeanDiss <- lapply(pDisss, mean)
  
  matchedClusters <- pMatchMeanDiss[which(pMatchMeanDiss == min(as.numeric(pMatchMeanDiss)))]
  
  ##
  ## Currently in the case of a tie just return the first one
  matchedGSMs <- names(matchedClusters[1])
  
  #Assess confidence 
  if(length(matchedClusters) == 1) {
    confidence <- "High"
  } else {
    confidence <- "Low"
  }    
  
  result <- list(GSMCluster,matchedGSMs,confidence)
  
  names(result) <- c("Mut","WT","Confidence")
  
  #Return paired labelled clusters 
  return(result)
}

#Match clusters for simple GSEs
simple.pair <- function(GSMCluster,labelTable) {
  clusterLabel <- unique(labelTable[GSMCluster,]$predLabel)
  if(clusterLabel == "FALSE") {
    return(NA)
  }
  
  matchedGSMs <- rownames(labelTable)[labelTable$predLabel == "FALSE"]
  
  result <- list(GSMCluster,matchedGSMs,"High")
  names(result) <- c("Mut","WT","Confidence")
  #Return paired labelled clusters 
  return(result)
}

#save matched GSMs in text file 
save.match <- function(GSMres,GSEId) {
  fName <- paste(GSMres[["Mut"]],collapse = "-")
  res1 <- lapply(GSMres, function(x){
    paste(x,collapse = " ")
  })
  writeData <- as.data.frame(do.call(rbind,res1))
  write.table(writeData, file = paste(GSEId,"-",fName,".txt",sep =""),col.names = FALSE)
}

#Get match and confidence
get.match.outcome2 <- function(labelTable, CurGSEmetadata) { #Pass through matrix with labelled clusters 
  multi <- check.multi(labelTable)
  
  if(multi) {
    outcome <- lapply(rownames(labelTable[labelTable$predLabel == "TRUE",,drop = FALSE]),function(x) {
      pair.clusters(x,labelTable, CurGSEmetadata)
    })
  } else {
    outcome <- lapply(rownames(labelTable[labelTable$predLabel == "TRUE",,drop = FALSE]),function(x) {
      simple.pair(x,labelTable)
    })
  } 
  return(outcome)
}




MatchClusters <- function(clusterLabels, CurGSEmetadata){
  get.match.outcome2(clusterLabels, CurGSEmetadata)
}




#######################################################
#######################################################


#Retreive gpl 
get.gpl <- function(GSEId) {
  metadata <- grabAllMetadataGSE(GSEId)
  gplList <- lapply(metadata$GSMMeta, function(x){
    x$gpl
  })
  return(unique(gplList))
}

#Retreive gpl 
get.gpl2 <- function(metadata) {
  gplList <- lapply(metadata$GSMMeta, function(x){
    x$gpl
  })
  return(unlist(unique(gplList)))
}
#Retreive labels from targets 
## retreive labels from metadata / exprs(gset) colnames head(exprs(gset[[1]]))
#input matching results; sort output with exprs(gset) colnames order
get.pair.labels <- function(matchedPair,gSet) {
  #Pull out targets of use 
  
  
  mutCluster <- unlist(strsplit(matchedPair$Mut," "))
  wtCluster <- unlist(strsplit(matchedPair$WT," "))
  
  subgroupL <- c(rep("1",length(mutCluster)),rep("2",length(wtCluster)))
  names(subgroupL) <- c(mutCluster,wtCluster)
  
  #Order subgroup labels by 
  allGSM <- colnames(exprs(gSet))
  label <- as.character(subgroupL[allGSM])
  
  return(label)
}

getFtpList <- function(ftp){ 
  
  txt <- getURL(ftp) 
  
  dir <- read.table( textConnection(txt),as.is=TRUE) 
  out <- data.frame(Dir=ftp,Filename=dir[, ncol(dir)],Size=dir[ ,5], 
                    Month=dir[ ,6],Day=dir[ ,7],Time=dir[ 
                      ,8],stringsAsFactors=FALSE) 
  closeAllConnections() 
  return(out) 
} 

check.annotation <- function(GPLId){
  
  GEO <- toupper(GPLId)
  stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
  
  gplurl <- "ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/"
  myurl <- sprintf(gplurl, stub, GEO)
  
  contents <- getFtpList(myurl)
  
  return("annot" %in% contents$Filename)
  
}

#Get data and transform
get.transform.data <- function(GSEId) {
  
  GPLId <- unlist(get.gpl(GSEId))
  #Check annotation 
  if(!check.annotation(GPLId)) {
    print(paste0("No annotation file for this platform. Trying to use GEMMA to get annotation files: ", GSEId))
    
    GPL_anno <- get_GPL_annotations_GEMMA(GPLId)
    if(nrow(GPL_anno)>1){
      pdf <- data.frame(ID = GPL_anno$ProbeName, GT = GPL_anno$GeneNames, GS = GPL_anno$GeneSymbols, GI = GPL_anno$NCBIids)
      rownames(pdf) <- pdf$ID
      colnames(pdf) <- c("ID" ,"Gene title" ,"Gene symbol","Gene ID")
      
      apdf <- AnnotatedDataFrame()
      pData(apdf) <- pdf
      #return(NA)
      gset <- getGEO(GSEId, GSEMatrix =TRUE, getGPL = F) ## GSE41277 and GSE39553 gets stuck here (needs AnnotGPL= FALSE)
      if(is.na(gset)){return(NA)}
      
      if (length(gset) > 1) idx <- grep(GPLId, attr(gset, "names")) else idx <- 1 ##
      gset <- gset[[idx]]
      
      featureData(gset) <- apdf
      
      
    } else {
      print(paste0("No annotation for this platform on GEMMA. Trying user supplied annotations: ", GSEId))
      gset <- getGEO(GSEId, GSEMatrix =TRUE, getGPL = T)
      
      if(is.na(gset)){return(NA)}
      
      if (length(gset) > 1) idx <- grep(GPLId, attr(gset, "names")) else idx <- 1 ##
      gset <- gset[[idx]]
      
      pdf <- pData(featureData(gset))
      
      # This is to detect gene column names from non-standard annotations. This could be more sophisticated by layering the word gene with "ID" and "name".
      gene_columns <- grepl("^gene|ORF",colnames(pdf), ignore.case = T)
      
      if(sum(gene_columns)>0){
        
        pdf$GS <- gsub("(.{100}).*","\\1", pdf[,which(gene_columns)[1]])
        pdf$GT <- gsub("(.{100}).*","\\1", pdf[,which(gene_columns)[min(2,sum(gene_columns))]])
        pdf$GI <-gsub("(.{100}).*","\\1",  pdf[,which(gene_columns)[min(3,sum(gene_columns))]])
        
        
      } else {
        
        print("No Gene symbols could be found in any annotations.")
        
        if(SimpleGeneSymbolsOnly){
          print("Scrapping this GSE")
          return(NA)
        }
        
        print("Using probe ID only.")
        
        pdf$GT <- pdf$ID
        pdf$GS <- pdf$ID
        pdf$GI <- pdf$ID
        
      } 
      
      pdf <- pdf[,c("ID","GT","GS","GI")]
      colnames(pdf) <- c("ID" ,"Gene title" ,"Gene symbol","Gene ID")
      
      apdf <- AnnotatedDataFrame()
      pData(apdf) <- pdf
      
      featureData(gset) <- apdf
    }
    
  }else {
    
    gset <- getGEO(GSEId, GSEMatrix =TRUE, AnnotGPL=T) ## GSE41277 and GSE39553 gets stuck here (needs AnnotGPL= FALSE)
    if (length(gset) > 1) idx <- grep(GPLId, attr(gset, "names")) else idx <- 1 ##
    gset <- gset[[idx]]
  }  
  
  # make proper column names to match toptable 
  fvarLabels(gset) <- make.names(fvarLabels(gset))
  
  # log2 transform
  ex <- exprs(gset)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) ||
    max(ex) > 30
  if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }
  
  return(gset)
}

filter.toptable <- function(tt){
  
  na.filter <- !is.na(tt$AveExpr)
  
  tt <- tt[na.filter,]
  
  ## take the top 2/3 of expressed probes. 
  averages = tt$AveExpr
  
  a <- sort(averages)
  min.ave.exprs <- a[length(a)/3]
  exprs.filter <- averages>min.ave.exprs
  
  gSymbols <- tt$Gene.symbol
  no.symbol.filter <- !gSymbols == ""
  
  sags <- split(tt[,"AveExpr", drop=F], gSymbols)
  highest.probes.keep <- unlist(lapply(sags, function(X){return(row.names(X)[which.max(unlist(X))])}))
  
  highest.probe.filter <- row.names(tt) %in% highest.probes.keep
  
  combined.filter <- exprs.filter & no.symbol.filter & highest.probe.filter
  
  return(tt[combined.filter,])
}


#Get top table of matched clusters
# perform differential expression
get.match.top.table <- function(sml,gset) {
  
  # eliminate samples marked as "GNA"
  sel <- which(sml != "GNA")
  sml <- sml[sel]
  gset <- gset[ ,sel]
  
  fl <- as.factor(sml)
  gset$description <- fl
  ################
  
  design <- model.matrix(~ description + 0, gset)
  colnames(design) <- levels(fl)
  fit <- lmFit(gset, design)
  
  cont.matrix <- makeContrasts("G1-G2",levels=design) ##
  
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.01)
  
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
  
  #remove low expressed / redundant probes and genes without symbols
  tTF <- filter.toptable(tT)
  
  #readjust P value
  tTF$adj.P.Val <- p.adjust(tTF$P.Value, method = "BH")
  
  tT <- tTF[,c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title")]
  
  return(tT)
}

# These functions are not used any more
#Generate a list of valid genes from gset for comparison 
get.gset.genes <- function(gSet) {
  geneSymbols <- fData(gSet)$Gene.symbol
  
  gSStr <- paste(geneSymbols,collapse = " ") 
  
  gSList <- unlist(strsplit(gsub("\\/+"," ",gSStr)," "))
  
  #Find gene symbols with more than 3 characters 
  genesList <- unique(sort(gSList[nchar(gSList) >= 3])) 
  
  return(genesList)
}

#Figure out potential names by matching gene symbols in gset with GSM info
get.table.name <- function(mutGSM, gSet,combinedInfoTable){ #input GSMs obtained from top table names 
  
  mutGSM <- unlist(strsplit(mutGSM, " "))
  #Get title, characteristics and source data from combined table
  tCSInfo <- unique(unlist(combinedInfoTable[mutGSM,]))
  
  #Get list of gene symbols from gset 
  gS <- get.gset.genes(gSet)
  
  #Match gene symbols with metadata info for potential genes
  potGenes <- gS[which(sapply(gS,function(x){
    grepl(x,paste(tCSInfo,collapse = " "),ignore.case = T) #Flawed: GSE16740 "tnnt2a" exists in gS - how to match "tnnt2" in info?
  }))]
  
  return(potGenes)
}

#Process potential gene names into file name format 
get.file.name <- function(pNames) {
  maxlen <- max(sapply(pNames, nchar))
  fName <- pNames[nchar(pNames) == maxlen]
  
  if(length(fName) != 1) {
    fName <- paste(fName,collapse = "_")
  }
  return(fName)
}

#Assess confidence of file names(perturbed gene)
assess.confidence <- function(fNames) {
  #Identify number of genes in title 
  conf <- lapply(fNames, function(x) {
    #Account for empty strings
    if(length(x) == 0) {
      return("No result")
    }
    if(!grepl("\\_",x)){
      return("High")
    } else {
      return("Medium")
    }
  })
  
  #Identify duplicated titles 
  doubleUp <- unlist(unique(fNames[duplicated(fNames)]))
  
  conf <- unlist(conf)
  conf[fNames %in% doubleUp] <- "Low"
  
  return(conf)
}

#Pass in GSEId and all targets files 
get.GSE.top.tables <- function(GSEId,MatchedPairs,CurGSEmetadata) {
  gset <- get.transform.data(GSEId)
  
  
  
  #Check if gset exists
  if(is.na(gset)){
    return(NA)
  }
  names(MatchedPairs) <- lapply(MatchedPairs, '[[', 'Mut')
  
  labelList <- lapply(MatchedPairs, function(x) {
    paste("G", get.pair.labels(x,gset), sep="") # set group names
  })
  
  
  topTables <- lapply(labelList,function(x){
    get.match.top.table(x,gset)
  })
  
  ################
  ## BREAK HERE INTO SEPARATE NAMING FUNCTION
  ################
  
  combinedTable <- data.frame(do.call(rbind.fill, lapply(CurGSEmetadata$GSMMeta, get.combined.table2)))
  rownames(combinedTable) <- do.call(rbind, lapply(CurGSEmetadata$GSMMeta, function(x) {x$gsm}))
  
  #Get names of each topTable (==target filename)
  
  mutGSMClusters <- lapply(MatchedPairs, function(x) {
    return(x$Mut)
  })
  
  comparisons <- do.call(rbind, lapply(1:length(MatchedPairs), function(I) {
    return(data.frame(index = I, mut = paste(MatchedPairs[[I]]$Mut, collapse=" "), control = paste(MatchedPairs[[I]]$WT, collapse=" "), direction = MatchedPairs[[I]]$Perturbation))
  }))
  
  return(list(topTables = topTables, comparisons = comparisons))
}


writeTopTables <- function(GSEId, topTables, comparisons, outputFolder = "unnamed"){
  #Create folder
  
  if(!file.exists(paste(getwd(),"/",outFolder,sep = ""))){
    dir.create(paste(getwd(),"/",outFolder,sep = ""))
  }
  
  if(!file.exists(paste(getwd(),"/",outFolder,"/", outputFolder,sep = ""))){
    dir.create(paste(getwd(),"/",outFolder,"/", outputFolder,sep = ""))
  }
  
  folderDir <- paste(getwd(),"/",outFolder,"/", outputFolder,"/",GSEId,sep = "")
  dir.create(folderDir)
  #Write files 
  
  for(j in 1:length(topTables)){
    write.table(topTables[[j]], file = paste(folderDir,"/",j,"_",gsub("[[:punct:]]","-",names(topTables[j])),"_",comparisons$direction[j],".txt",sep =""), quote=FALSE, sep="\t", row.names = FALSE)
  }
  
}



writeComps <- function(GSEId, comparisons, topTables, outputFolder = "unnamed"){
  #Create folder
  
  
  if(!file.exists(paste(getwd(),"/",outFolder,sep = ""))){
    dir.create(paste(getwd(),"/",outFolder,sep = ""))
  }
  
  if(!file.exists(paste(getwd(),"/",outFolder,"/", outputFolder,sep = ""))){
    dir.create(paste(getwd(),"/",outFolder,"/", outputFolder,sep = ""))
  }
  
  folderDir <- paste(getwd(),"/",outFolder,"/", outputFolder,"/",GSEId,sep = "")
  dir.create(folderDir)
  #Write files 
  
  for(j in 1:length(topTables)){
    comparisons$verified_molecule[j] = names(topTables[j])
    comparisons$filename[j] = paste(j,"_",names(topTables[j]),"_",comparisons$direction[j],".txt",sep ="") 
  }
  
  write.table(comparisons, file = paste(folderDir,"/comparisons.txt",sep =""), quote=FALSE, sep="\t", row.names = FALSE)
  
}


rename.MatchedPairs <- function(GSE, MatchedPair, GSEmetadata){
  print(paste0("Renaming ", GSE))
  
  best.names <- best.name(GSE, MatchedPair, GSEmetadata)
  suppressWarnings(if(best.names == "Species not supported"){return(NA)})
  
  extract.best <- lapply(lapply(best.names, '[', 'Best'), function(X){paste(unlist(X), collapse=" ")})
  extract.conf <- lapply(best.names, '[', 'Confidence')
  
  pert.direction <- get.perturbation(GSEmetadata)
  
  names(MatchedPair) <- unlist(extract.best)
  
  MatchedPairwConf <- lapply(1:length(MatchedPair), function(I){
    m <- MatchedPair[[I]]
    m$NameConfidence <- unlist(extract.conf)[I]
    m$Perturbation <- pert.direction
    return(m)
  })
  
  names(MatchedPairwConf) <- names(MatchedPair)
  
  return(MatchedPairwConf)
}




#Get GSE title Genes
get.genes.in.text <- function(text, geneList){
  
  title.split <- strsplit(text, " |\\.|\\,|_|-")[[1]]
  geneMatches <- geneList[tolower(geneList)%in%tolower(title.split)]
  
  if(!length(geneMatches) > 0){
    
    geneMatches <- geneList[which(sapply(geneList, function(x){
      grepl(x, text, ignore.case = TRUE)
    }))]
  }
  
  if(!length(geneMatches) > 0){
    return("NOTHING")
  }
  
  return(geneMatches)
}

################

get.genes <- function(GSE, MatchedPair, GSEmetadata){ 
  species <- get.species(GSEmetadata)
  
  if(species %in% names(gene.List)){
    geneList <- gene.List[[species]]
  } else {
    print("Retrieving gene list from ensembl")
    geneList <- lapply(species, function(I){
      gene.List[[I]] <<- fix.list(get_ENSEMBL_symbol_map(I))
    })[[1]]
  }
  
  title.genes <- get.genes.in.text(GSEmetadata$title, geneList)
  summary.genes <- get.genes.in.text(GSEmetadata$summary, geneList)
  
  gsm.genes <- lapply(MatchedPair, function(x){
    mutGSM <- strsplit(x$Mut, " ")[[1]]
    mutTitles <- paste(sapply(GSEmetadata$GSMMeta[mutGSM],function(x){x$title}), collapse = " ")
    ghit <- get.genes.in.text(mutTitles, geneList)
  })
  
  
  
  gsm.gene.frequencies <- lapply(gsm.genes, table)
  
  return(list(title.summary.hits = table(c(title.genes, summary.genes)), gsm.hits = gsm.gene.frequencies))
  
}

check.names <- function(gene.counts){
  
  genes <- names(gene.counts)
  
  ## some common words should not be considered because they always match words like "embryo" "cell" "transcriptome" "myocardial"
  genes <- genes[(!tolower(genes) %in% c("cript","myoc","emb","cel","ell") )]
  gene.counts <- gene.counts[genes]
  
  if(length(genes) < 1){
    gene.counts <- 1
    names(gene.counts) <- "none"
    genes <- "none"
  }
  
  lengths <- unlist(lapply(genes, nchar))
  combined <- gene.counts * lengths
  best <- which(combined == max(combined))
  short <- lengths <= 3
  confidence <- rep("Low",length(gene.counts))
  names(confidence) <- genes
  confidence[!short] <- "Medium"
  
  if((!short[best]) & (gene.counts[best] > 1)){
    confidence[best] <- "High"
  } else if (short[best] & gene.counts[best] > 1){
    confidence[best] <- "Medium"
  }
  
  
  
  return(list(Confidence = confidence, Best = best))
}

best.name <- function(GSE, MatchedPair, GSEmetadata){
  GSEmetadata <- nameGSMs(GSEmetadata)
  
  supported <- get.species(GSEmetadata) %in% supported_species
  if(!supported){
    return("Species not supported")
  }
  
  genes <- get.genes(GSE, MatchedPair, GSEmetadata)
  
  checked <- suppressWarnings(check.names(genes$title.summary.hits))
  
  conf <- checked$Confidence
  best <- checked$Best
  best.g <- names(best)
  
  gsm.best <- suppressWarnings(lapply(genes$gsm.hits, function(X){
    
    checked_gsm <- check.names(X)
    
    b <- checked_gsm$Best
    c <- checked_gsm$Confidence
    g <- names(X[b])
    
    if(g == best.g){
      
      return.gene <- g
      return.conf <- "High"
      
    } else if (g%in%names(conf[conf%in%c("Medium","High")])){
      return.gene <- g
      return.conf <- "Medium"
      
    } else {
      return.gene <- best.g
      if (conf[best.g] == "High"){
        return.conf <- "Medium"
      } else {
        return.conf <- "Low"
      }
      
    }
    
    return(list(Best = return.gene, Confidence = return.conf))
  }))
  
  return(gsm.best)
}


get.perturbation <- function(GSEmetadata){
  GSEmetadata <- nameGSMs(GSEmetadata)
  
  up.feats <- c("overexp", "express", "transgen", "expos", "tg", "induc")
  down.feats <- c("knock", "null", "ko", "s[hi]rna", "delet", "reduc", "\\-\\/", "\\/\\-", "\\+\\/", "\\/\\+", "cre", "flox", "mut","defici")
  
  gse.text <- paste(GSEmetadata$title, GSEmetadata$summary, GSEmetadata$overall_design, collapse=" ")
  title.split <- strsplit(gse.text, " |\\.|\\,|_")[[1]]
  
  
  
  
  upMatches <- sapply(up.feats, function(x){
    length(grep(x, title.split, ignore.case = TRUE))
  })
  
  downMatches <- sapply(down.feats, function(x){
    length(grep(x, title.split, ignore.case = TRUE))
  })
  
  
  if(sum(upMatches) > sum(downMatches)){
    
    return("+")
  } else {
    
    return("-")
  }
}


#Rename all metadata with GSM Ids
nameGSMs <- function(GSEMetadata) {
  GSMId <- lapply(GSEMetadata$GSMMeta, function(x){
    x$gsm
  })
  
  names(GSEMetadata$GSMMeta) <- unlist(GSMId)
  
  return(GSEMetadata)
}

#Get species from Metadata and return emsembl format
get.species <- function(GSEmetadata){
  
  species <- GSEmetadata$GSMMeta[[1]]$organism_ch1
  
  if(is.null(species)){return(NULL)}
  
  eSpecies <- tolower(paste(substring(species, 1, 1),strsplit(species," ")[[1]][2],sep = ""))
  
  
  return(eSpecies)
  
}

####################################################################


parseGSEs <- function(GSEMeta){
  
  fields.to.use <- c("title","summary","overall_design","GSMMeta")
  
  gsmMeta <- paste(unique(unlist(GSEMeta$GSMMeta)), collapse=" ")
  GSEMeta$GSMMeta <- gsmMeta
  return(GSEMeta[fields.to.use])
}


checkPresenceGSEFieldsList <- function(data_and_fields) {
  
  GSEFeatures <- data_and_fields$data
  field.features <- data_and_fields$feats
  
  feat.list <- lapply(GSEFeatures, function(field){
    
    feat.vector <- unlist(lapply(field.features, function(feat){
      grepl(feat, field, ignore.case = T)
    }))
    
    names(feat.vector) <- field.features
    return(feat.vector)
  })
  
  
  return(feat.list)
}


get.GSE.feature.vector <- function(GSEid){
  GSEMeta <- grabAllMetadataGSE(GSEid)
  parsedGSEMeta <- parseGSEs(GSEMeta)
  
  
  feature.presence.list <- lapply(names(parsedGSEMeta), function(Fname){
    
    field.features <- gsub(paste(Fname, "\\.", sep=""), "", new.feats.to.keep[grepl(Fname, new.feats.to.keep)])
    
    data_and_fields <- list(data = parsedGSEMeta[[Fname]], feats = field.features)
    
    checkPresenceGSEFieldsList(data_and_fields)
  })
  
  names(feature.presence.list) <- names(parsedGSEMeta)
  
  final.vector <- unlist(feature.presence.list, use.names = T)
  return(final.vector)
}


############################################################################
biomart.ID <- "ENSEMBL_MART_ENSEMBL"
host.ID <-  "dec2015.archive.ensembl.org"

find_supported_datasets <- function(default=TRUE){
  print("Connecting to ENSEMBL...")
  ensembl <- useMart(biomart.ID, dataset='hsapiens_gene_ensembl', host=host.ID)
  sets <- listDatasets(ensembl)
  organisms <- gsub("([a-z]*)_.*","\\1",sets[,1],"_")
  return(organisms)
}

supported_species <- sort(find_supported_datasets())


##########################
# use ENSEMBL to generate mapping  between ENSEMBL ids and gene symbols
# input parameteres is the species string like 'hsapiens'
get_ENSEMBL_symbol_map <- function(species){
  require(biomaRt)
  dataset_name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomart.ID, dataset=dataset_name, host=host.ID)
  ENSEMBL_symbol_map <- getBM(attributes=c("external_gene_name"), mart = ensembl)
  return(ENSEMBL_symbol_map)
}

#Turn gene data frame into list and remove brackets etc
fix.list <- function(geneDF) {
  geneList <- geneDF[,1][nchar(geneDF[,1])>2] #Only gene symbols longer than 2 characters
  uniqueList <- unique(gsub("[[:space:]]\\([^\\]]*\\)", "",geneList)) #remove brackets and content eg. "TMEM151A (1 of 2)"
  return(uniqueList)
}


#Allow us to annotate outlier GPLs
get_GPL_annotations_GEMMA <- function(GPL){
  url <- paste0("https://gemma.msl.ubc.ca/rest/v2/platforms/",GPL,"/annotations/")
  raw.result <- GET(url = url)
  bin <- httr::content(raw.result, "raw")
  writeBin(bin, "myfile.gz")
  GPLannotations <- read.table("myfile.gz", header = T, sep="\t", comment.char = "#", quote="")
  return(GPLannotations)
}



################################################################################################################################################
# SHINY CODE
################################################################################################################################################



gene.List <<- list() 

checkedGSEs <- vector()

default.page.length <- 5

########################################################################
# Server functionality
########################################################################


shinyServer(function(input, output, session) {
  
  session$onSessionEnded(stopApp)
  
  hide(id = "loading-content", anim = TRUE, animType = "fade")
  
  CurTab <- reactiveValues(Tab = "Verify")
  
  #CurGSEs <- reactiveValues(GSEs = character(0))
  CurGSEs <- reactiveValues(GSEs = NULL)
  
  model <- reactive({
    if(exists("LabelModel")){
      return(LabelModel)
    } else {
      print("ERROR: No model for predicting labels of GSE")
      return(NULL)
    }
  })
  
  
  #########################################
  ## Add tool tips
  ## search tab
  addTooltip(session, id = "searchTerm", title = "You can use & or | for multiple keywords.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "uploadFile", title = "e.g. GSE14491 GSE28448 GSE42373...", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "analysis_select", title = "Which kind of experimental design are you interested in? 'Perturbations' streamline GRN construction. 'All experiments' will not restrict which GSE are processed.",  placement = "bottom", trigger = "hover")
  addTooltip(session, id = "species_select_search_div", title = "GEOracle can process data from one species at a time.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "microarray_checkbox_div", title = "Do you want to display only microarray results which can be rapidly processed in the next step? Recommended.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "keep_button_div", title = "Warning: This will remove all GSE that are not currently selected from the above list!", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "remove_button_div", title = "Warning: This will remove the currently selected GSE from the above list!", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "download_GSE_button_div", title = "For longer analyses and for reproducibility we suggest you download and save your list of GSE.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "analyse_GSE_button_div", title = "If you want to analyse the above list directly, you can. This will take you to the next tab and begin clustering samples. We suggest downloading this list of GSE first.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "GSEtext", title = "You can type or copy & paste your own list of GSE in here, separated by whitespace.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "upload_file_div", title = "You can upload a text file with a list of GSE separated by whitespace here.", placement = "bottom", trigger = "hover")
  
  #upload / filter tab
  addTooltip(session, id = "filter_select", title = "Select how to filter the list of GSE for your desired analysis.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "speciesUI", title = "By default the most common species in the list of GSE is selected.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "minClusterSize", title = "What is the minimum size of a cluster of GSM samples to consider? You need 3 for an inerpretable p-value.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "maxClusterSize", title = "What is the maximum size of a cluster of GSM samples to consider?", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "allbetween", title = "Do you need every cluster size to fit between the min and max? This effectively toggles cluster size filtering on or off.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "simpleOnly", title = "Only consider GSE with 1 comparison (likely to be paired correctly). This is useful for automated analyses.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "maxcomps", title = "The maximum number of comparisons within a GSE before it is filtered out. Some GSE have hundreds of irrelevent samples and comparisons which can slow down the analysis.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "singlePlatformOnly", title = "Only consider GSE with a single platform (technology). It is currently not trivial to process GSE with GSM from multiple platforms. Recommended.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "onechannel", title = "Filter out microarrays with more than one channel as they can not be automatically processed at this stage. Recommended.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "predpert", title = "Only consider GSE which are predicted to be perturbation experiments.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "gene_symbols_only", title = "Only include GSE where gene symbols can be easily automatically annotated. If not, the process is MUCH slower and some results may contain microarray probe IDs only, but you will process more GSE.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "see_all_GSE", title = "Just let all GSE through the filtering stage. Not recommended.", placement = "bottom", trigger = "hover")
  
# set DE 
  addTooltip(session, id = "folder", title = "The output folder name to save results in. This could represent your entire analysis. e.g. 'Human heart failure'.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "logfc", title = "The absolute Log2 fold change threshold to use for automatically detecting differentially expressed genes. The full results are also saved.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "pval", title =  "The Benajmini & Hochberg adjusted p-value threshold to use for automatically detecting differentially expressed genes. The full results are also saved.", placement = "bottom", trigger = "hover")
  addTooltip(session, id = "process", title =  "Calculate differential expression for all GSE. This can take a while.", placement = "bottom", trigger = "hover")

  #########################################
  ## search function inspired by ScanGEO implementation
  
  GSEsToInclude <- reactiveValues(include=1)
  
  output$analysis_select <- renderUI({ 
    selectInput('analysis_type', 'Experiment type', c("Perturbations","All experiments"), selected = "Perturbations", selectize = FALSE)
  })
  
  output$GEoracleLogo2 <- output$GEoracleLogo <- renderImage({
    
    list(src = "georacle.png",
         title = "GEOracle",
         contentType = 'image/png',
         alt = "GEOracle Logo",
         align = "left",
         width="90%")
  }, deleteFile = FALSE)
  
  datasetInput_1 <- eventReactive(input$find, {
    
    withProgress(message = "Searching GEO:",
                 detail = "this can be slow for generic search terms as we must retrieve GSM metadata to filter the species.", value = 0.2, {
                   
                   ### formatting search terms into SQL query with & and | functionality
                   
                   search_input <- toupper(input$searchTerm)
                   
                   if(grepl("\\&", search_input)&grepl("\\|", search_input)){
                     
                     showModal(modalDialog(
                       title = "Input Error",
                       "We can only handle one type of complex search symbol (& or |) at a time. Please try another search.",
                       easyClose = TRUE
                     ))
                     
                     Final_query <- "select title, gse, type, submission_date from gse where title like '%Input Error%'"
                     
                   } else if(grepl("\\&", search_input)){
                     
                     query_words <- trimws(strsplit(search_input,'\\&')[[1]])
                     
                     formatted_query_words <- paste("'%", query_words, "%'", sep="")
                     
                     combo <- lapply(formatted_query_words, function(query_word){
                       paste(c("title like","summary like","overall_design like"), query_word, collapse = " or ")
                     })
                     
                     AND_joint_query <- paste0("(", paste(combo, collapse=") and ("), ")")
                     
                     Final_query <- paste("select title, gse, type, submission_date from gse where ", AND_joint_query)
                     
                   } else if(grepl("\\|", search_input)){
                     
                     query_words <- trimws(strsplit(search_input,'\\|')[[1]])
                     
                     formatted_query_words <- paste("'%", query_words, "%'", sep="")
                     
                     combo <- lapply(formatted_query_words, function(query_word){
                       paste(c("title like","summary like","overall_design like"), query_word, collapse = " or ")
                     })
                     
                     OR_joint_query <- paste0("(", paste(combo, collapse=") | ("), ")")
                     
                     Final_query <- paste("select title, gse, type, submission_date from gse where ", OR_joint_query)
                   } else { 
                     query_words <- trimws(strsplit(search_input,'\\|')[[1]])
                     formatted_query_words <- paste("'%", query_words, "%'", sep="")
                     joint_query <- paste(c("title like","summary like","overall_design like"), formatted_query_words, collapse = " or ")
                     Final_query <- paste("select title, gse, type, submission_date from gse where ", joint_query, sep="")
                   }
                   
                   incProgress(amount = 0.2, detail = "Querying database: this can be slow for generic search terms.")
                   
                   ress <- dbGetQuery(metadata, Final_query)
                   
                   searchedGSEmetadata <- lapply(ress$gse, tryCatch({grabAllMetadataGSE}, error=function(e){
                     
                     cat("ERROR when getting Metadata: ",conditionMessage(e), "\n")
                     return(NA)
                   }))
                   
                   
                   
                   names(searchedGSEmetadata) <- ress$gse
                   searchedGSEmetadata <- searchedGSEmetadata[!is.na(searchedGSEmetadata)]
                   
                   incProgress(amount = 0.1, detail = "Figuring out species")
                   
                   species <- lapply(ress$gse, function(x){
                     
                     if(!x %in% names(searchedGSEmetadata)){
                       return("none") 
                     }
                     
                     X <- searchedGSEmetadata[[x]]
                     
                     if(nrow(X$GSMMeta[[1]]) < 1){return(NA)}
                     
                     gsespecies <- get.species(X)
                     
                     if(is.null(gsespecies)){return(NA)}
                     
                     return(tolower(gsub("(.).* (.*)","\\1\\2",gsespecies)))
                     
                   })
                   names(species) <- ress$gse
                   
                   ress$species <- species
                   
                 })
    
    
    if(input$GSE_platform_filter){
      filtered.datasetInput <- na.omit(ress[ress$species==input$selectSpecies,])
      filtered.datasetInput <- filtered.datasetInput[grepl("Expression profiling by array", filtered.datasetInput$type),]
    } else {
      filtered.datasetInput <- na.omit(ress[ress$species==input$selectSpecies,])
    }
    
    GSEsToInclude$include <- filtered.datasetInput$gse
    
    return(filtered.datasetInput)
  })
  
  
  output$Dim <- renderText(dim(datasetInput())[1])
  
  datasetInput <- reactive({
    datasetInput_1()[datasetInput_1()$gse%in%GSEsToInclude$include,]
  })
  
  output$Dim <- renderText(dim(datasetInput())[1])
  
  output$keep.button <-
    renderUI(expr = if (input$find & (dim(datasetInput())[1] > 0)) {
      actionButton("keepSelectedGSE", label = "Keep only selected GSE", icon("check"), style="color: #fff; background-color: green; border-color: #2e6da4")
    } else {
      NULL
    })
  
  output$remove.button <-
    renderUI(expr = if (input$find & (dim(datasetInput())[1] > 0)) {
      actionButton("removeSelectedGSE", label = "Remove selected GSE", icon("remove"), style="color: #fff; background-color: #ff6b30; border-color: #2e6da4")
    } else {
      NULL
    })
  
  observeEvent(
    {input[["keepSelectedGSE"]]},{
      
      GSEsToInclude$include <- datasetInput()$gse[input$GEOtable_rows_selected]
      
    }, ignoreNULL = TRUE, autoDestroy = TRUE, priority = 1)
  
  observeEvent(
    {input[["removeSelectedGSE"]]},{
      
      GSEsToInclude$include <- setdiff(GSEsToInclude$include, datasetInput()$gse[input$GEOtable_rows_selected])
      
    }, ignoreNULL = TRUE, autoDestroy = TRUE, priority = 1)
  
  output$download.button <-
    renderUI(expr = if (input$find & (dim(datasetInput())[1] > 0)) {
      tagList(
        renderText("For reproducibility we suggest you download the above list of GSE before proceeding.\n"),
        downloadButton("downloadGSE", label = "Download entire list of GSE")
      )
    } else {
      NULL
    })
  
  output$batch.download.button <-
    renderUI(expr = if (input$find & (dim(datasetInput())[1] > 20)) {
      tagList(
        hr(),
        renderText("The long list of GSE has been batched into groups of 20 for easier processing. Processing many GSE in batches will help ensure you don't lose too much if the app crashes. It also helps deal with some memory issues."),
        downloadButton("downloadGSEBatch", label = "Download pre-batched GSE")
      )
    } else {
      NULL
    })
  
  output$direct.analyse.button <-
    renderUI(expr = if (input$find & (dim(datasetInput())[1] > 0)) {
      tagList(
        hr(),
        #renderText("For longer analyses and for reproducibility we suggest you download and save your list of GSE. However if you want to analyse the above list now, you can."),
        actionButton("directAnalyseGSE", label = "Analyse this list of GSE", icon("check"), style="color: #fff; background-color: #ff6b30; border-color: #2e6da4")
      )
    } else {
      NULL
    })
  
  
  observeEvent(input$find, {
    
    if (input$searchTerm > 0){
      if (dim(datasetInput())[1] == 0) {
        output$datamessage <- renderText(paste(dim(datasetInput())[1],
                                               "studies found \n searching for \"", isolate(input$searchTerm), "\"",
                                               "\n in", isolate(input$selectSpecies),
                                               "\n Try modifying your search term"))        
      } else {
        output$datamessage <- renderText(paste(dim(datasetInput())[1],
                                               "studies found \n searching for \"", isolate(input$searchTerm), "\"",
                                               "\n in", isolate(input$selectSpecies)
        ))
      }
    } else {
      output$datamessage <- renderText("You didn't enter a search term")
    }
    
    output$searchStatus <- renderUI({
      verbatimTextOutput('datamessage')
    })
    
    output$GEOtable <- renderDataTable({
      
      Summary_tab <<- datasetInput()
      Summary_tab$gse <<- unlist(lapply(datasetInput()$gse, function(x){
        paste0("<a href='", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
               x, "' target='_blank'>", x,"</a>")}))
      return(Summary_tab)
    }, escape = FALSE)
    
    
  })
  
  output$downloadGSE <- downloadHandler(
    
    filename = function() { paste0("GEOracle_search_",input$selectSpecies, "_", input$searchTerm, '.txt') },
    content = function(file) {
      write.table(datasetInput()$gse, gsub("\\|","OR",file), row.names = F, col.names = F, quote= F, sep="\t")
    }
  )
  
  ###
  output$downloadGSEBatch <- downloadHandler(
    
    
    
    ###
    
    filename = paste0("GEOracle_search_",input$selectSpecies, "_", input$searchTerm, '.zip'),
    content = function(fname) {
      
      search_ID <- gsub("\\|","OR",paste0("GEOracle_search_",input$selectSpecies, "_", input$searchTerm))
      
      if(!file.exists(paste(getwd(),"/",outFolder,sep = ""))){
        dir.create(paste(getwd(),"/",outFolder,sep = ""))
      }
      
      if(!file.exists(paste(getwd(),"/",outFolder,"/", search_ID,sep = ""))){
        dir.create(paste(getwd(),"/",outFolder,"/", search_ID,sep = ""))
      }
      
      folderDir <- paste(getwd(),"/",outFolder,"/", search_ID,sep = "")
      
      gses_split <- split(datasetInput()$gse, ceiling(seq_along(datasetInput()$gse)/20))
      
      lapply(names(gses_split), function(S){
        write.table(gses_split[[S]], file=paste0(folderDir, "/", S, ".txt"), sep="\t", col.names = F, row.names = F, quote=F)
      })
      setwd(folderDir)
      
      fs <- list.files(path = folderDir)
      
      zipfile <- zip(zipfile=fname, files=fs)
      setwd(this.dir)
      
      return(zipfile)
    },
    contentType = "application/zip"
  )
  
  
  GDSmeta <- dbGetQuery(metadata, "select * from gds\n")
  GEOspecies <- sort(unique(tolower(gsub("(.).* ([a-z]*).*","\\1\\2",unique(GDSmeta$platform_organism)))))
  clean.species <- GEOspecies[nchar(GEOspecies) > 4]
  
  
  
  updateSelectInput(session, "selectSpecies", choices = clean.species,
                    selected = "hsapiens")
  
  #################
  
  obs_de <- observeEvent(
    input$submitGSETtext
    ,{
      CurGSEs$GSEs <- strsplit(input$GSEtext,"\\W")[[1]]
    })
  
  obs_uf <- observeEvent(
    input$uploadFile
    ,{
      file <- read.table(input$uploadFile$datapath, sep=" ", stringsAsFactors = F)
      if(ncol(file) > 1){
        reload.file <- read.table(input$uploadFile$datapath, sep=" ", header=TRUE, stringsAsFactors = F, row.names = 1)
        CurGSEs$GSEs = rownames(reload.file)[reload.file$Perturbation == 1]
      } else {
        CurGSEs$GSEs = unlist(file)
      }
    })
  
  obs_da <- observeEvent(
    input$directAnalyseGSE
    ,{
      CurGSEs$GSEs <- datasetInput()$gse
      updateTabsetPanel(session, "inTabset",
                        selected = "Process GSE list")
    })
  
  AllGSEmetadata <- eventReactive(CurGSEs$GSEs, {
    print("Obtaining Metadata")
    Metadata <- lapply(CurGSEs$GSEs, tryCatch({grabAllMetadataGSE}, error=function(e){
      cat("ERROR when getting Metadata: ",conditionMessage(e), "\n")
      return(NA)
    }))
    names(Metadata) <- CurGSEs$GSEs
    
    return(Metadata[!is.na(Metadata)])
    
  })
  
  
  
  
  GSEclusters <- eventReactive(AllGSEmetadata(),{
    
    print("Clustering Samples")
    
    
    output$filter_select <- renderUI({ 
      selectInput('filter', 'Strictness', c("Perturbations (Default)","All experiments (Loose)", "Custom"), selected = c("Perturbations (Default)","All experiments (Loose)")[(input$analysis_type == "All experiments") + 1], selectize = FALSE)
    })
    
    output$speciesUI <- renderUI({ 
      selectInput("species", "Species", names(table(unlist(GSEspecies()))), selected = names(table(unlist(GSEspecies())))[which.max(table(unlist(GSEspecies())))], width = 200, selectize = FALSE)
    })
    
    
    
    withProgress(message = "IN PROGRESS:",
                 detail = "Clustering: Initialising", value = 0.5, {
                   clusters <- lapply(AllGSEmetadata(), function(X){
                     incProgress(0.9/length(AllGSEmetadata()), detail = paste0("Clustering: ",X$gse))
                     return(GEOclustering(X))
                   })
                 })
    return(clusters)
  })
  
  
  MatchedPairs <- eventReactive(GSEclusters(), {
    
    print("Matching clusters")
    
    
    clusterLabels <- reactive({
      Predictions <- lapply(names(GSEclusters()), function(X){LabelClusters(GSEclusters()[[X]], AllGSEmetadata()[[X]], model())})
      names(Predictions) <- names(GSEclusters())
      return(Predictions)
    })
    
    FixedLabels <- reactive({
      
      FixedPredictions <- mapply(FixLabels, GSEclusters(), clusterLabels(), input$analysis_type, SIMPLIFY = F)
      
      invalid <- unlist(lapply(FixedPredictions, function(X){return(X$confidence[1] == "Invalid")}))
      return(FixedPredictions[!invalid])
    })
    
    withProgress(message = "IN PROGRESS:",
                 detail = "Matching Clusters: Initialising", value = 0.1, {
                   Pairs <- lapply(names(FixedLabels()), function(X){
                     incProgress(amount = 0.8/length(FixedLabels()), detail = paste0("Matching Clusters: ", X))
                     res <- MatchClusters(FixedLabels()[[X]], AllGSEmetadata()[[X]])
                     names(res) <- unlist(lapply(res, function(X){paste(X$Mut, X$WT, collapse="_")}))
                     return(res)
                   })
                   names(Pairs) <- names(FixedLabels())
                 })
    return(Pairs)
  })
  
  GSEplatforms <- reactive({
    res <- lapply(names(MatchedPairs()), function(x){
      X <- AllGSEmetadata()[[x]]
      gpls <- unique(unlist(lapply(X$GSMMeta, function(x){
        x$gpl
      })))
    })
    names(res) <- names(MatchedPairs())
    return(res)
  })
  
  GSEspecies <- reactive({
    res <- lapply(names(MatchedPairs()), function(x){
      X <- AllGSEmetadata()[[x]]
      species <- unique(unlist(lapply(X$GSMMeta, function(x){
        x$organism_ch1
      })))
    })
    names(res) <- names(MatchedPairs())
    return(res)
  })
  
  
  GSEchannels <- reactive({
    res <- lapply(names(MatchedPairs()), function(x){
      X <- AllGSEmetadata()[[x]]
      
      channels <- unique(unlist(lapply(X$GSMMeta, function(x){
        x$label_ch2
      })))
    })
    names(res) <- names(MatchedPairs())
    return(res)
  })
  
  
  output$numGSEs <- renderText(
    paste0('You have ', length(MatchedPairs()), ' GSE that match your selected experiment type.')
  )
  
  
  
  #########################################  
  ### when click the GO button
  obsr <- observeEvent(input$go, {
    MPCs <- reactiveValues(
      ## only those GSEs that pairing worked
      mpairs = isolate(MatchedPairs()[unlist(lapply(MatchedPairs(), length)) > 0])
    )
    
    # SETTING FILTERS
    include.pairs <- reactive({
      if(input$filter == "Perturbations (Default)"){
        ## set strict defaults    
        mincl <- 2 
        maxcl <- 6
        simple <- 0
        singlePlatform <- 1
        mainspecies <- names(table(unlist(GSEspecies())))[which.max(table(unlist(GSEspecies())))]
        allbetween <- 0
        maxcomps <- 10
        oncehannel <- 1
        predpert <- 1
        SimpleSimpleGeneSymbolsOnly <<- 1
      }else if(input$filter == "All experiments (Loose)"){
        ## set strict defaults    
        mincl <- 2 
        maxcl <- 10
        simple <- 0
        singlePlatform <- 1
        mainspecies <- names(table(unlist(GSEspecies())))[which.max(table(unlist(GSEspecies())))]
        allbetween <- 0
        maxcomps <- 20
        oncehannel <- 1
        predpert <- 0
        SimpleGeneSymbolsOnly <<- 0
      }else{
        mincl <- input$minClusterSize 
        maxcl <- input$maxClusterSize
        
        simple <- input$simpleOnly
        singlePlatform <- input$singlePlatformOnly
        mainspecies <- input$species
        allbetween <- input$allbetween
        maxcomps <- input$maxcomps
        onechannel <- input$onechannel
        predpert <- input$predpert
        SimpleGeneSymbolsOnly <<- input$gene_symbols_only
      }
      
      
      platforminclude <- unlist(lapply(names(MPCs$mpairs), function(X){
        numplatforms <- GSEplatforms()[[X]]
        channels <- GSEchannels()[[X]]
        return(((length(numplatforms) == 1) | (!singlePlatform)) & is.na(channels[1]))
      }))
      
      simpleinclude <- unlist(lapply(MPCs$mpairs, function(X){
        return(length(X) == 1)
      }))
      
      toobiginclude <- unlist(lapply(MPCs$mpairs, function(X){
        return(length(X) <= maxcomps)
      }))
      
      speciesinclude <- unlist(lapply(names(MPCs$mpairs), function(X){
        spec <- GSEspecies()[[X]]
        if(length(spec) > 1){ return(FALSE) } else {
          return(spec == mainspecies)
        }
      }))
      
      clustersizeinclude <- unlist(lapply(MPCs$mpairs, function(X){
        lengths <- lapply(X, function(Y){
          return(c(length(strsplit(Y$Mut," ")[[1]]), length(strsplit(Y$WT," ")[[1]])))
        })
        maxpass <- unlist(lengths) <= maxcl
        minpass <- unlist(lengths) >= mincl
        
        if(allbetween){
          return((sum(maxpass) == length(maxpass)) & (sum(minpass) == length(minpass)))
        } else {
          return((sum(maxpass) > 0) & (sum(minpass) > 0))
        }
      }))
      
      predpertinclude <- unlist(lapply(names(MPCs$mpairs), function(X){
        if(predpert == 0){
          return(TRUE)
        }
        return(as.logical(predict(PertModel, matrix(get.GSE.feature.vector(X), 1, length(new.feats.to.keep)))))
      }))
      
      includedf <- data.frame(si = (simpleinclude|!simple), ci = clustersizeinclude, pi = platforminclude, spi = speciesinclude, tbi= toobiginclude, ppi = predpertinclude)
      
      allinclude <- (simpleinclude|!simple) & clustersizeinclude & platforminclude & speciesinclude & toobiginclude & predpertinclude
      
      if(input$see_all_GSE){
        allinclude <- rep(T, length(allinclude))
      }
      
      return(allinclude)
    })
    
    ## perform filtering
    FMPCs <- reactive({
      m <- na.omit(MPCs$mpairs[include.pairs()])
      withProgress(message = "IN PROGRESS:",
                   detail = "Detecting names: Retrieving gene list from Ensembl", value = 0.1, {
                     
                     
                     m2 <- lapply(na.omit(names(m)), function(X){
                       renamed.ress <- rename.MatchedPairs(X, m[[X]], AllGSEmetadata()[[X]])
                       incProgress(0.9/length(m), detail = paste0("Detected names for: ", X))
                       return(renamed.ress)
                     })
                   })
      names(m2) <- na.omit(names(m))
      
      
      return(m2)
    })
    
    FMPC2 <<- reactiveValues(
      mpairs = isolate(FMPCs())
      #original_predicted_names = names(isolate(FMPCs())) - would be good to collect this information but cant yet
    )
    
    ## create list of filtered GSEs    
    processedGSEsDF <- reactive({
      valid.lengths <- unlist(lapply(FMPC2$mpairs, function(X){ 
        return(length(X) - sum(unlist(lapply(X, is.na))))
      }))
      
      
      
      df <- data.frame(Comps = valid.lengths, row.names = names(FMPC2$mpairs))
    })
    
    
    pages <- reactiveValues(
      page = 0
    )
    
    pagesync <- reactive({
      req(pages$page, input$processedGSEs_state)
      
      
      if(pages$page == floor(input$processedGSEs_state$start / input$processedGSEs_state$length)){
        return(TRUE)
      } else {
        return(FALSE)
      }
    })
    
    
    output$processedGSEs <- DT::renderDataTable({
      
      DT::datatable(processedGSEsDF(), selection = "single", options = list(lengthMenu = list(c(default.page.length, 3*default.page.length, -1), c(default.page.length, 3*default.page.length, 'All')),  pageLength = default.page.length, stateSave = TRUE
                                                                            #datatable(head(iris, 30), callback = JS('table.page("next").draw(false);'))
      ))
    })
    
    
    SelectedGSE <- reactive({
      req(input$processedGSEs_row_last_clicked)
      rownames(processedGSEsDF())[input$processedGSEs_row_last_clicked]
    })
    
    
    
    observe({
      invalidateLater(3000)
      
      
      if(!isolate(pagesync())){
        
        session$sendCustomMessage("pager", pages$page)
      }
      
    }, priority = -2)
    
    ##################################################################################
    
    
    output$matchedPairs <- renderUI({
      req(SelectedGSE())
      
      
      
      MPC.l <- FMPC2$mpairs
      req(MPC.l[[SelectedGSE()]])
      
      withProgress(message = "Loading...",
                   detail = "Just a sec.", value = 0.5, {
                     GSEmeta <- AllGSEmetadata()[[SelectedGSE()]]
                     GSMmeta <- GSEmeta$GSMMeta
                     matchoutputs <- lapply(1:length(MPC.l[[SelectedGSE()]]), function(Y){
                       X <-  MPC.l[[SelectedGSE()]][[Y]]
                       suppressWarnings(
                         if(is.na(X)) return(h3(paste("REMOVED COMPARISON ", Y), style = "color: #FF6200;", align = "center"))
                       )
                       
                       
                       
                       
                       muts <- unlist(strsplit(X$Mut, " "))
                       wts <- unlist(strsplit(X$WT, " "))
                       pert <- X$Perturbation
                       
                       GSMtitles <- data.frame(do.call(rbind, lapply(GSMmeta, function(XX){return(c(GSM = XX$gsm, Title = XX$title))})))
                       
                       muts.titles <- GSMtitles[GSMtitles$GSM%in%muts,]
                       wts.titles <- GSMtitles[GSMtitles$GSM%in%wts,]
                       
                       alloutputs <- list(
                         hr(style = "height: 10px; border: 0; box-shadow: 0 -10px 5px -10px #8c8b8b inset;")
                         ,
                         div(id="comp_name_box", style="display:inline-block",
                             if(nchar(names(MPC.l[[SelectedGSE()]])[Y])>20){
                               h4(renderText(paste0(substr(names(MPC.l[[SelectedGSE()]])[Y], 1, 20), "...")))
                             } else {
                               h4(renderText(names(MPC.l[[SelectedGSE()]])[Y]))
                             }
                         )
                         ,
                         bsTooltip(id = "comp_name_box", title = "The current name assigned to the comparison of samples below. Should identify the experimental condition. e.g. the gene that has been knocked down or disease name.",
                                   placement = "bottom", trigger = "hover")
                         ,
                         div(id ="comp_pert_box", style="display:inline-block",
                             h3(renderText(pert)) 
                         )
                         ,
                         bsTooltip(id = "comp_pert_box", title = "The current direction assigned to the comparison of samples below. Identifies whether the experimental condition represents the addition or subtraction of a perturbation agent or disease state.",
                                   placement = "bottom", trigger = "hover")
                         ,
                         
                         div(id="rename_comp_textbox", style="display:inline-block",
                             textInput(paste0("namebox_", SelectedGSE(), "_", Y), "Rename", width=100)
                         )
                         ,
                         bsTooltip(id = "rename_comp_textbox", title = "Rename the comparison here.",
                                   placement = "bottom", trigger = "hover")
                         ,
                         div(id="rename_comp_button", style="display:inline-block",
                             actionButton(paste0("rename_", SelectedGSE(), "_", Y),"", icon("pencil"), style="color: #000; background-color: #94F4FF; border-color: #000000")
                         )
                         ,
                         bsTooltip(id = "rename_comp_button", title = "Confirm renaming of this comparison.",
                                   placement = "bottom", trigger = "hover")
                         ,
                         div(id="comp_direction", style="display:inline-block",
                             selectInput(paste0("perturbation_", SelectedGSE(), "_", Y), "+ or -", choices = c("+","-"), selected = pert, width = 100, selectize = FALSE)
                         )
                         ,
                         bsTooltip(id = "comp_direction", title = "Select whether the perturbation agent is present (+) or absent (-). e.g. disease state (+) / gene knockdown (-).",
                                   placement = "bottom", trigger = "hover")
                         ,
                         h4(renderText("\nPerturbation samples"))
                         ,
                         Mutants <- DT::renderDataTable(as.data.frame(muts.titles), options = list(dom = 't'))
                         ,
                         br()
                         ,
                         h4(renderText("\nControl samples"))
                         ,
                         Wilds = DT::renderDataTable(as.data.frame(wts.titles), options = list(dom = 't'))
                         ,
                         div(id = "remove_comparison_button",
                             actionButton(paste0("remove_", SelectedGSE(), "_", Y),"REMOVE THIS COMPARISON", icon("remove"), style="color: #fff; background-color: #ff6b30; border-color: #2e6da4")
                         )
                         ,
                         bsTooltip(id = "remove_comparison_button", title = "Remove the above comparison from this GSE.",
                                   placement = "bottom", trigger = "hover")
                         ,
                         
                         hr(style = "height: 10px; border: 0; box-shadow: 0 10px 20px -10px #8c8b8b inset;")
                       )
                     })
                     
                     out <- list(
                       link = h4(a(href=paste0("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",SelectedGSE()), paste(SelectedGSE(), " - ", GSEmeta$title), target = "blank"))
                       ,
                       div(id="remove_GSE_button", style="display:inline-block",
                           actionButton(paste0("removeGSE_", SelectedGSE()),"REMOVE THIS ENTIRE GSE", icon("remove"), style="color: #fff; background-color: red; border-color: #2e6da4")
                       )
                       ,
                       bsTooltip(id = "remove_GSE_button", title = "Remove this GSE from the analysis.",
                                 placement = "bottom", trigger = "hover")
                       ,
                       div(id = "add_comparison_button", style="display:inline-block",
                           actionButton("addComp","ADD A COMPARISON", icon("plus"), style="color: #fff; background-color: green; border-color: #2e6da4")
                       )
                       ,
                       bsTooltip(id = "add_comparison_button", title = "Add a new / custom comparison of samples to this GSE.",
                                 placement = "bottom", trigger = "hover")
                       ,
                       summary = renderText(paste(substr(GSEmeta$summary,1,600)," ..."))
                       ,
                       matches = matchoutputs
                     )
                   })
      
      
      return(out)
      
    })
    
    
    ###########################################################
    
    
    output$MutantTable <- DT::renderDataTable({
      req(SelectedGSE())
      
      GSEmeta <- AllGSEmetadata()[[SelectedGSE()]]
      GSMmeta <- GSEmeta$GSMMeta
      GSMlist <- grabGSMids(SelectedGSE())
      
      MutantTable <- data.frame(Sample = GSMlist, Title = do.call(rbind, GSMmeta)[,2])
      
      DT::datatable(MutantTable)
      
    })
    
    output$ControlTable <- DT::renderDataTable({
      req(SelectedGSE())
      
      GSEmeta <- AllGSEmetadata()[[SelectedGSE()]]
      GSMmeta <- GSEmeta$GSMMeta
      GSMlist <- grabGSMids(SelectedGSE())
      
      ControlTable <- data.frame(Sample = GSMlist, Title = do.call(rbind, GSMmeta)[,2])
      
      DT::datatable(ControlTable)
      
    })
    
    output$addComp <- renderUI({
      req(SelectedGSE())
      
      GSEmeta <- AllGSEmetadata()[[SelectedGSE()]]
      GSMmeta <- GSEmeta$GSMMeta
      
      out <- list(
        link = h4(a(href=paste0("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",SelectedGSE()), paste(SelectedGSE(), " - ", GSEmeta$title), target = "blank"))
        ,
        message2 = h4(renderText("Create a new comparison here for the current GSE. Simply select those samples you want to designate as Mutants and Controls in the tables below. Make sure you have entered the name of the molecule / perturbation in the text box and selected the direction of perturbation. Then click 'ADD' at the bottom of the page."))
        ,
        hr()
        ,
        div(id="add_new_name_box", style="display:inline-block",
            textInput("addnamebox", "Comparison Name", value = "New Comparison", width=200)
        )
        ,
        bsTooltip(id = "add_new_name_box", title = "Name the comparison of samples below. Should identify the experimental condition. e.g. the gene that has been knocked down or disease name.",                  placement = "bottom", trigger = "hover")
        ,
        div(style="display:inline-block",
            selectInput("addperturbation", "+ or -", choices = c("+","-"), selected = "+", width = 100, selectize = FALSE)
        )
        ,
        bsTooltip(id = "addperturbation", title = "Select whether the experimental condition is present (+) or absent (-). e.g. disease state (+) / gene knockdown (-).",          placement = "bottom", trigger = "hover")
        ,
        
        message3 = h3(renderText("Select Perturbation samples:")),
        DT::dataTableOutput('MutantTable'),

        bsTooltip(id = "MutantTable", title = "Select the samples which represent the perturbation. e.g. gene knock down or disease state", placement = "bottom", trigger = "hover")
        ,
        hr(),
        message3 = h3(renderText("Select Control samples:")),
        DT::dataTableOutput('ControlTable'),
        bsTooltip(id = "ControlTable", title = "Select the samples which represent the control. e.g. normal samples or healthy patients", placement = "bottom", trigger = "hover")
        ,
        hr(),
        actionButton("finishAdd","ADD THIS COMPARISON", icon("addition"), style="color: #fff; background-color: green; border-color: #2e6da4")
        ,
        bsTooltip(id = "finishAdd", title = "Confirm and add the above comparison to the GSE", placement = "bottom", trigger = "hover")
        ,
        actionButton("cancelAdd","CANCEL", icon("step-backward"), style="color: #fff; background-color: red; border-color: #2e6da4")
        
      )
      
      
      
      
      
      return(out)
      
    })
    
    #addition observer
    observeEvent(
      {input$finishAdd},{
        GSEmeta <- AllGSEmetadata()[[SelectedGSE()]]
        GSMmeta <- GSEmeta$GSMMeta
        GSMTable <- do.call(rbind, GSMmeta)[,1:2]
        
        GSMids <- grabGSMids(SelectedGSE())
        
        tmpPair <- list(Mut = paste(GSMids[input$MutantTable_rows_selected], sep= " "), WT = paste(GSMids[input$ControlTable_rows_selected], sep= " "), Confidence = "High", NameConfidence = "High", Perturbation = input$addperturbation)
        
        FMPC2$mpairs[[SelectedGSE()]][[length(FMPC2$mpairs[[SelectedGSE()]])+1]] <- tmpPair
        names(FMPC2$mpairs[[SelectedGSE()]])[length(FMPC2$mpairs[[SelectedGSE()]])] <- input$addnamebox
        
        
        CurTab$Tab = "Verify"
        
      }, ignoreNULL = TRUE, autoDestroy = TRUE, priority = 1)
    
    
    observeEvent(
      {input$cancelAdd},{
        
        CurTab$Tab = "Verify"
        
      }, ignoreNULL = TRUE, autoDestroy = TRUE, priority = 1)
    
    
    ##############################################################################################################################
    ### start observerS
    
    
    removeMP.observers <- list()
    removeGSE.observers <- list()
    renameMP.observers <- list()
    pertMP.observers <- list()
    
    
    aco <- observeEvent(
      {input[["addComp"]]},{
        
        CurTab$Tab = "Add"
        
      }, ignoreNULL = TRUE, autoDestroy = TRUE, priority = 1)
    
    
    obsr2 <- observe({
      
      
      if(length(renameMP.observers) > 0){
        
        
        renameMP.observers <<- lapply(renameMP.observers, function(X){ X$destroy() })
        removeMP.observers <<- lapply(removeMP.observers, function(X){ X$destroy() })
        removeGSE.observers <<- lapply(removeGSE.observers, function(X){ X$destroy() })
        pertMP.observers <<- lapply(pertMP.observers, function(X){ X$destroy() })
        
      }
      
      
      removeGSE.observers <<- lapply(names(input)[grepl(paste0("removeGSE_",SelectedGSE()), names(input))], function(x) {
        
        
        observeEvent(
          {input[[x]]},{
            FMPC2$mpairs[[SelectedGSE()]] <<- NULL
            
            if(length(renameMP.observers) > 0){
              
              renameMP.observers <<- lapply(renameMP.observers, function(X){ X$destroy() })
              renameMP.observers <<- list()
              removeMP.observers <<- lapply(removeMP.observers, function(X){ X$destroy() })
              removeMP.observers <<- list()
              pertMP.observers <<- lapply(pertMP.observers, function(X){ X$destroy() })
              pertMP.observers <<- list()
            }
            
          }, ignoreNULL = TRUE, autoDestroy = TRUE, priority = 1)
        
      })
      
      removeMP.observers <<- lapply(names(input)[grepl(paste0("remove_",SelectedGSE()), names(input))], function(x) {
        
        
        observeEvent(
          {input[[x]]},{
            
            to.rem <- as.numeric(gsub(paste0("remove_",SelectedGSE(),"_"),"",x))
            suppressWarnings(if(!is.na(FMPC2$mpairs[[SelectedGSE()]][[to.rem]])) {
              FMPC2$mpairs[[SelectedGSE()]][[to.rem]] <<- NA
            })
          }, ignoreNULL = TRUE, autoDestroy = TRUE)
        
        
      })
      
      
      renameMP.observers <<- lapply(names(input)[grepl(paste0("rename_",SelectedGSE()), names(input))], function(x) {
        
        observeEvent(
          {input[[x]]},{
            
            to.ren <- as.numeric(gsub(paste0("rename_",SelectedGSE(),"_"),"",x))
            suppressWarnings(if(!is.na(FMPC2$mpairs[[SelectedGSE()]][[to.ren]])) {
              names(FMPC2$mpairs[[SelectedGSE()]])[to.ren] <<- input[[paste0("namebox_",SelectedGSE(),"_", to.ren)]]
            }) 
          }, ignoreNULL = TRUE, autoDestroy = TRUE)
        
      })
      
      
      pertMP.observers <<- lapply(names(input)[grepl(paste0("perturbation_",SelectedGSE()), names(input))], function(x) {
        observeEvent(
          {input[[x]]},{
            
            
            
            to.ren <- as.numeric(gsub(paste0("perturbation_",SelectedGSE(),"_"),"",x))
            suppressWarnings(if(!is.na(FMPC2$mpairs[[SelectedGSE()]][[to.ren]])) {
              if(!FMPC2$mpairs[[SelectedGSE()]][[to.ren]]["Perturbation"] == input[[paste0("perturbation_",SelectedGSE(),"_", to.ren)]]){
                
                FMPC2$mpairs[[SelectedGSE()]][[to.ren]]["Perturbation"] <<- input[[paste0("perturbation_",SelectedGSE(),"_", to.ren)]]
              }
            })
          }, ignoreNULL = TRUE, autoDestroy = TRUE)
        
      })
      
      
      checkedGSEs <<- unique(append(checkedGSEs, SelectedGSE()))
      
      if(SelectedGSE()%in%rownames(processedGSEsDF())){
        pages$page <<- floor(which(rownames(processedGSEsDF()) == SelectedGSE()) / default.page.length)
      }
      
      
    }, autoDestroy = T, priority = -1)
    
    #})       
    
    
    ##############################################################################################################################
    ### end sub observerS
    
    
    ###########################################################
    
    output$makeGRN <- renderUI({
      
      withProgress(message = "Working:",
                   detail = "Take a break...", value = 0, {
                     
                     ## remove NAs
                     MPC.g <- lapply(FMPC2$mpairs, function(X){return(X[which(!is.na(X))])})
                     
                     ## remove comparisons with <2 replicates in either class as these will fail differential expression
                     MPC.g.l2 <- lapply(MPC.g, function(X){
                       
                       Muts.AL.2 <- unlist(lapply(X, function(G){ return(length(unlist(strsplit(G$Mut, " "))))})) > 1
                       WT.AL.2 <- unlist(lapply(X, function(G){ return(length(unlist(strsplit(G$WT, " "))))})) > 1
                       
                       return(X[Muts.AL.2 & WT.AL.2])})
                     
                     #                     MPC.g.f <- MPC.g[unlist(lapply(MPC.g, length)) > 0]
                     
                     ## remove GSEs with nothing left after previous removal steps
                     MPC.g.f <- MPC.g.l2[unlist(lapply(MPC.g, length)) > 0]
                     req(MPC.g.f)
                     
                     ProcessedData <- suppressWarnings(lapply(names(MPC.g.f), function(X){ 
                       #progress feedback
                       incProgress(1/length(MPC.g.f), detail = paste0("Processing: ", X))
                       print(paste0("Processing: ", X))
                       
                       res <- tryCatch(get.GSE.top.tables(X, MPC.g.f[[X]], AllGSEmetadata()[[X]]), error=function(e) NA)
                       
                       if(!is.na(res)){
                         names(res$topTables) <- names(MPC.g.f[[X]])
                         #names(res$filenames) <- names(MPC.g.f[[X]])
                       }
                       return(res)
                       
                       
                     }))
                     names(ProcessedData) <- names(MPC.g.f)
                   })  
      
      errors <- is.na(ProcessedData)
      ProcessedData <- ProcessedData[!errors]
      
      
      withProgress(message = "Writing results to disk...",
                   detail = "This may take a while...", value = 0.5, {
                     lapply(names(ProcessedData), function(X){ writeTopTables(X, ProcessedData[[X]]$topTables, ProcessedData[[X]]$comparisons, input$folder) })
                     lapply(names(ProcessedData), function(X){ writeComps(X, ProcessedData[[X]]$comparisons, ProcessedData[[X]]$topTables, input$folder) })
                   })
      
      
      withProgress(message = "Calculating edges...",
                   detail = "This may take a while...", value = 0.5, {
                     
                     
                     ####################
                     ## THIS SHOULD BE A FUNCTION
                     ####################
                     
                     edges <- lapply(names(ProcessedData), function(X){
                       
                       
                       sub.edges <- lapply(names(ProcessedData[[X]]$topTables), function(Y){
                         tmptable <- ProcessedData[[X]]$topTables[[Y]]
                         
                         ### this needs to be revisited
                         pert = MPC.g.f[[X]][[Y]][["Perturbation"]]
                         
                         include <- tmptable$adj.P.Val < as.numeric(input$pval) & abs(tmptable$logFC) > as.numeric(input$logfc)
                         print(paste0(sum(na.omit(include)), " edges from ", X))
                         if(sum(na.omit(include)) == 0){
                           return(NA)
                         }
                         subtab <- na.omit(tmptable[include,])
                         
                         effectcol <- cases(
                           "+" = subtab$logFC > 0,
                           "-" = subtab$logFC < 0
                         )
                         
                         
                         
                         df <- data.frame(Regulator = Y, Target = subtab$Gene.symbol, Perturbation = pert, Effect = effectcol, Source = X)
                         
                         df$edgetype <- cases(
                           "Act" = (df$Perturbation == "+" & df$Effect == "+")|(df$Perturbation == "-" & df$Effect == "-"),
                           "Inh" = (df$Perturbation == "-" & df$Effect == "+")|(df$Perturbation == "+" & df$Effect == "-")
                         )
                         
                         return(df[!df$Target == "",])
                         
                       })
                       
                       names(sub.edges) <- names(ProcessedData[[X]]$topTables)
                       return(do.call(rbind, sub.edges))
                     })
                     
                     ####################
                     ## END FUNCTION
                     ####################
                     
                     
                     names(edges) <- names(ProcessedData)
                     
                     
                     all.edges <- data.frame(Regulator = "Empty",	Target = "Empty",	Perturbation = "Empty",	Effect = "Empty",	Source = "Empty",	edgetype = "Empty")
                     combined.edges <- unique(do.call(rbind, edges[ unlist(lapply(edges, function(X){return(sum(!is.na(X)) > 0 )})) ] ))     
                     
                     if(!is.null(combined.edges)){
                       all.edges <- combined.edges
                     }
                     
                     if(!file.exists(paste(getwd(),"/",outFolder,sep = ""))){
                       dir.create(paste(getwd(),"/",outFolder,sep = ""))
                     }
                     
                     if(!file.exists(paste(getwd(),"/",outFolder,"/", input$folder,sep = ""))){
                       dir.create(paste(getwd(),"/",outFolder,"/", input$folder,sep = ""))
                     }
                     
                     write.table(na.omit(all.edges), paste0(outFolder,"/",input$folder,"/","AllEdges.txt"), sep="\t", col.names = T, row.names = F, quote=FALSE)
                     #write.table("aaa", paste0(outFolder,"/",input$folder,"/","AllEdges.txt"), sep="\t", col.names = T, row.names = F, quote=FALSE)
                     
                     
                     
                   })
      
      ##### download data
      
      
      output$downloadData <- downloadHandler(
        filename = paste0(input$folder, ".zip"),
        content = function(fname) {
          #tmpdir <- tempdir()
          tmpdir <- paste0(outFolder,"/",input$folder,"/")
          setwd(tmpdir)
          print(tmpdir)
          
          fs <- list.files()
          
          zip(zipfile=fname, files=fs)
        },
        contentType = "application/zip"
      )
      
      
      
      withProgress(message = "Plotting Graph...",
                   detail = "This may take a while, especially for large networks...", value = 0.5, {
                     
                     if(length(edges) > 0){
                       gseids <- names(edges)
                     } else {
                       gseids <- NA
                     }
                     
                     sigedgesdf <- data.frame(GSE = gseids, UpReg = 0, DownReg = 0)
                     
                     if(!is.null(all.edges)){
                       
                       split.all.edges <- split(all.edges, all.edges$Source)
                       
                       lapply(names(split.all.edges), function(X){ 
                         
                         typecounts <- table(split.all.edges[[X]]$edgetype)
                         
                         if(!is.na(typecounts['Act'])) {
                           sigedgesdf[sigedgesdf$GSE == X, "UpReg"] <<- typecounts["Act"]
                         }
                         
                         if(!is.na(typecounts['Inh'])) {
                           sigedgesdf[sigedgesdf$GSE == X, "DownReg"] <<- typecounts["Inh"]
                         }
                         
                         return(TRUE)  
                       })
                       
                       
                       if(nrow(na.omit(all.edges)) > 0 ){
                         
                        
                         
                         
                         out <- list(
                           
                           message2 = h3(renderText("Processing is finished and the results are ready for download! Click the 'Download Results' button.\n\n Below you can see the number of significantly differentially expressed genes from each GSE using standard thresholds.")),
                           
                           message5 = downloadButton("downloadData", label = "Download Results"),
                           
                           message4 = renderTable(sigedgesdf),
                           
                           #                           message6 = p(HTML("<A HREF=\"javascript:history.go(0)\">Please restart GEOracle to perform another analysis.</A>"))
                           #message6 = p(HTML("<A HREF=\"javascript:close_window()\">Please restart GEOracle to perform another analysis.</A>"))
                           
                           
                           message6 = actionButton("close_app", label = "Close GEOracle")#, onclick= "setTimeout(function(){window.close();},500);")
                         )
                         
                         
                         
                       } else {
                         out <- list(
                           message2 = h3(renderText("Processing is finished and the results are ready for download! Click the 'Download Results' button.\n\n Unfortunately we couldn't detect any significantly differentially expressed genes from each GSE using standard thresholds.")),
                           
                           message5 = downloadButton("downloadData", label = "Download Results"),
                           
                           message6 = p(HTML("<A HREF=\"javascript:history.go(0)\">Please restart GEOracle to perform another analysis.</A>"))
                           
                           
                         )
                         
                       }
                       
                       
                     } else {
                       
                       out <- list(
                         
                         message2 = h3(renderText("Processing is finished and the results are ready for download! Click the 'Download Results' button.\n\n Unfortunately we couldn't detect any significantly differentially expressed genes from each GSE using standard thresholds.")),
                         
                         message5 = downloadButton("downloadData", label = "Download Results"),
                         
                         message6 = p(HTML("<A HREF=\"javascript:history.go(0)\">Please restart GEOracle to perform another analysis.</A>"))
                         
                       )
                     }
                     
                     
                   })
      return(out)
      
    })
    
    
    
    
  })
  
  ###########################################################
  
  
  
  output$curtab <- renderText({
    CurTab$Tab
  })
  
  output$curtabs <- renderUI({
    list(hr = hr(), message =  div(style="display:inline-block", h5(id="curtabmessage", paste0("Current stage is: "))), message2 =  div(style="display:inline-block", textOutput("curtab")))
  })
  
  observeEvent(input$finished, {
    CurTab$Tab = "GRN"  
  })
  
  observeEvent(input$goback, {
    CurTab$Tab = "Verify"  
  })
  
  observeEvent(input$close_app, {
    js$closeWindow()
    stopApp()
  })
  
})

