### R-script for the LCMS data processing and quality control for the lipidomics of blood samples for the ADNI1 cohort.
## Author Dinesh Kumar Barupal dinkumar@ucdavis.edu West Coast Metabolomics Center, UC Davis USA

### This script has five steps. Follow instructions on each step to reproduce the results. 

#Step 1. Extraction of peak heights from the extracted ion chromatogram files for targeted list of lipids.
#Step 2. LOESS normalization and batch effect removal.
#STEP 3. A rule based data Filtering to remove duplicate peaks
#STEP 4. Data visualization for the internal standards
#STEP 5: Combining POS and NEG mode data, additional filtering and data visualization

##################################################################################################################
### Step 1. Extraction of peak heights from the extracted ion chromatogram files for targeted list of lipids.  ###
##################################################################################################################

# You will need a windows computer. The script was not tested on a linux machine. 
# You need to have a basic understanding of R programming language to follow up the steps. 
# Read the Nature Scientific Data Descriptor paper to get familiar with the data generationa and pipeline. 
# Please download the "ADMC Lipidomics Extracted Ion Chromatograms" (adni1-lipidomics-eics.zip) file from ADNI repository. 
# Make a new folder "ADNI1_lipidomics" in D: or E: drive in your computer. Do not use C: drive to avoid permission errors. 
# Make two subfolders "ESIPOS_EIC" and "ESINEG_EIC" inside the ADNI1_lipidomics folder. 
# Extract the adni1-lipidomics-eics.zip to the ADNI1_lipidomics folder. It has two sub-directories one for the positive mode and another for the negeative mode EICs.
# Copy the csv files from ESI pos folder to the "ESIPOS_EIC" and the files from ESI neg folder to the "ESINEG_EIC" folders.
# Download the "wcmc-lipidomics-database.zip" from the Sage Bionetwork (https://www.synapse.org/#!Synapse:syn10208654). 
# Unzip the file. There are two txt files into it. Copy and paste the "esi_positive_mode_targets.txt" file into the "ESIPOS_EIC" folder and "esi_negative_mode_targets.txt" into the "ESINEG_EIC" folder. 
# Download the "positive_mode_metadata.txt" and "negative_mode_metadata.txt" files from the Sage Bionetwork (Exact page location). These files contain the subject id to raw LCMS data file mapping and the time-stamp of the data acquisition. Copy these files into the ADNI1_lipidomics folder and into the "ESIPOS_EIC" and "ESINEG_EIC" folders. 
# Download the duplicate_samples_ids.txt file from Synpase project page (https://www.synapse.org/#!Synapse:syn10208654).
# With this we are ready to proceeed. 

#  load Required R libraries for the step 1. 
library(pacman) # This is the R package manager for easy installation and loading of the R-packages. 
pacman::p_load(doSNOW, foreach,XLConnect)
options(java.parameters = "-Xmx4g") # Exporting Excel tables need more memory. 

##################### Processing of the Positive mode data #######################
# In this step we will extract the peak height for each target in the esi_positive_mode_targets.txt file. 

datadir <- "LOCATION of ESI Positive Mode CSV files" # provide the full location of the ESI positive mode EICs. 
datadirfiles <- dir(datadir)
#csvfiles <- datadirfiles[grep("[.]csv$",datadirfiles)]
#### Import of Global Databases
postargetdb <- read.delim("esi_positive_mode_targets.txt", header = T, stringsAsFactors = F, sep="\t")
isddb.pos <- postargetdb[which(postargetdb$Type=="ISD"),]
isddb.pos <- isddb.pos[-which(isddb.pos$Cpd%in%c("ISD_TG","ISD_PE")==TRUE),] 
targetdb.pos <- postargetdb[which(postargetdb$Type=="TargetCPD"),]
pos_metadata <- read.delim("positive_mode_metadata.txt", header = F, stringsAsFactors = F, sep="\t")

eiccsvfiles <- datadirfiles[grep("[.]CSV$",datadirfiles)]

cl<-makeCluster(6) ## number of cores goes here, for calculating faster.
#registerDoParallel(cl, cores = 6)
registerDoSNOW(cl)
foreach (i = 1:length(eiccsvfiles)) %dopar% {
  con1 <- file(paste0(datadir,gsub("[.]CSV","",eiccsvfiles[i]),"_ISD_eics.txt"),"w")
  con2 <- file(paste0(datadir,gsub("[.]CSV","",eiccsvfiles[i]),"_target_eics.txt"),"w")

  isintendf <- readLines(paste0(datadir,"/",eiccsvfiles[i]),n=-1L)
  isintendf <- c(isintendf,"+ESI EIC")

  eicanchors <- grep("+ESI EIC",isintendf)
  isdanchors.list <- lapply(1:(length(eicanchors)-1), function(x) {isintendf[eicanchors[x]:(eicanchors[x+1]-1)]  })
  targetnamevec <- sapply(isdanchors.list, function(x) { strsplit(x[1]," <|>:")[[1]][2]})
  targetnamevec <- gsub("ï¿½","",targetnamevec)

  ### First get the Internal Standards Data from a raw LCMS data file. 

  idf <- isddb.pos
  idf$samplert <- 1.0
  idf$sampleid <- gsub("[.]CSV","",eiccsvfiles[i])
  idf$deltart <- 0.0
  idf$sampleintensity <- 0.0
  for (k in 1:nrow(isddb.pos)) {
    isrt <- idf$RT[k]
    eicdf <- isdanchors.list[[which(targetnamevec==isddb.pos$Cpd[k])[1]]][-c(1:2)]
    eicdf <- as.data.frame(do.call(rbind,lapply(eicdf,function(x) {strsplit(x,",")[[1]]})),stringsAsFactors = F)
    eicdf$V2 <- as.numeric(eicdf$V2)
    eicdf$V3 <- as.numeric(eicdf$V3)
    rtrange <- which(eicdf$V2<(isrt+0.30) & eicdf$V2>(isrt-0.30)) #
    isrt.sample <- eicdf$V2[rtrange[which.max(eicdf$V3[rtrange])]]
    isinten.sample <- max(eicdf$V3[rtrange])
    idf$samplert[k] <- isrt.sample
    idf$deltart[k] <- isrt.sample-isrt
    idf$sampleintensity[k] <- isinten.sample
  }
  write.table(idf,con1,col.names = T,row.names = F,quote = F, sep="\t")
  ## Correct the retention time using a polynomial linear model. 

  rtmodel <- lm(samplert ~ poly(RT,2), idf)
  new.df <- data.frame(RT = targetdb.pos$RT) # predict the new RTs
  newRtVec <- predict(rtmodel,new.df)

  ## Extract the peak heights for all the targets in the "esi_positive_mode_targets.txt" file. 
  tdf <- targetdb.pos
  tdf$newRT <- newRtVec
  tdf$samplert <- 1.0
  tdf$sampleid <- gsub("[.]CSV","",eiccsvfiles[i])
  tdf$deltart <- 0.0
  tdf$sampleintensity <- 0.0
  tdf$integration <- 1.0
  for (k in 1:nrow(targetdb.pos)) {
    isrt <- tdf$newRT[k]
    eicdf <- isdanchors.list[[which(targetnamevec==targetdb.pos$Cpd[k])[1]]][-c(1:2)]
    eicdf <- as.data.frame(do.call(rbind,lapply(eicdf,function(x) {strsplit(x,",")[[1]]})),stringsAsFactors = F)
    eicdf$V2 <- as.numeric(eicdf$V2)
    eicdf$V3 <- as.numeric(eicdf$V3)
    rtrange <- which(eicdf$V2<(isrt+0.20) & eicdf$V2>(isrt-0.20))
    isrt.sample <- eicdf$V2[rtrange[which.max(eicdf$V3[rtrange])]]
    isinten.sample <- max(eicdf$V3[rtrange])
    if(isinten.sample==0) {
      rtrange <- which(eicdf$V2<(isrt+0.25) & eicdf$V2>(isrt-0.25)) # if the intensity of sample comes up as zero, we increases the threshold of the rt delta.
      isrt.sample <- eicdf$V2[rtrange[which.max(eicdf$V3[rtrange])]]
      isinten.sample <- max(eicdf$V3[rtrange])
      is.inten.integration <- max(eicdf$V3[rtrange])/median(eicdf$V3[rtrange])
      tdf$samplert[k] <- isrt.sample
      tdf$deltart[k] <- isrt.sample-isrt
      tdf$sampleintensity[k] <- isinten.sample
      tdf$integration[k] <-  is.inten.integration
    } else {
      is.inten.integration <- max(eicdf$V3[rtrange])/median(eicdf$V3[rtrange])
      tdf$samplert[k] <- isrt.sample
      tdf$deltart[k] <- isrt.sample-isrt
      tdf$sampleintensity[k] <- isinten.sample
      tdf$integration[k] <-  is.inten.integration
    }
  }
  write.table(tdf,con2,col.names = T,row.names = F,quote = F, sep="\t")
  print(i)
  close(con1)
}
stopCluster(cl)
closeAllConnections()

# Reading back the peak height values and combining them. 
datadirfiles <- dir(datadir)
target_eic_files <- datadirfiles[grep("_target_eics.txt$",datadirfiles)]

pos_inten_list_delta <-lapply(1:length(target_eic_files), function(x) {
  target_data <- read.delim(paste0(datadir,target_eic_files[x]), header = T, stringsAsFactors = F)
  res <- data.frame(intensity=target_data$sampleintensity, rtdelta =target_data$deltart,integration = target_data$integration)
  res
})

# Export the intensities values.
pos_inten_inten <- lapply(1:length(pos_inten_list_delta), function(x) {pos_inten_list_delta[[x]]$intensity   })
pos_inten_inten <- as.data.frame(do.call(rbind, pos_inten_inten), stringsAsFactors = F)
pos_inten_inten <- round(pos_inten_inten,digits = 0)

pos_inten_inten <- rbind(targetdb.pos$Cpd,pos_inten_inten)
pos_inten_inten <- cbind(data.frame(filename=c("",gsub("_target_eics.txt","",target_eic_files))),pos_inten_inten)
pos_inten_inten$filename <- gsub("[.]2$","",pos_inten_inten$filename)
pos_inten_inten <- pos_inten_inten[- (which(pos_inten_inten$filename%in%pos_metadata$V1==FALSE)[-1]),]

tdf <- rbind(t(pos_inten_inten)[1,],as.character(sapply(pos_inten_inten$filename, function(x) { pos_metadata$V2[which(pos_metadata$V1==x)]})),as.character(sapply(pos_inten_inten$filename, function(x) { pos_metadata$V3[which(pos_metadata$V1==x)]})),t(pos_inten_inten)[-1,])
tdf <- tdf[,-1]

etdf <- cbind(as.character(c("","",pos_inten_inten[1,])),c("","","InChiKeys",rep("",length(pos_inten_inten[1,-1]))) ,c("","SampleLabel","TimeStamp",rep("",length(pos_inten_inten[1,-1]))),tdf   )
write.table(etdf, file="etdf_positive_mode_target_intensities.txt",col.names = F, row.names = F, quote = F, sep="\t")

## export to excel directly.
exc <- XLConnect::loadWorkbook("adni-lipidomics-raw.xlsx", create = T) # this file will be used by the LOESS normalization script.
XLConnect::createSheet(exc,'sudo')
XLConnect::createSheet(exc,'POS')
XLConnect::writeWorksheet(exc, etdf, sheet = "POS", startRow = 1, startCol = 1,header = F)
#XLConnect::saveWorkbook(exc)

##################################################################################
##################### Processing of the Negetive mode data #######################
##################################################################################

# In this step we will extract the peak height for each target in the esi_negative_mode_targets.txt file. All the steps are same as for the ESI-positive mode data. 

datadir <- "LOCATION OF NEGATIVE MODE DATA" # switch to the ESI neg EIC directory. 
datadirfiles <- dir(datadir)
csvfiles <- datadirfiles[grep("[.]csv$",datadirfiles)]

#### Import of Global Databases
negtargetdb <- read.delim("esi_negative_mode_targets.txt", header = T, stringsAsFactors = F, sep="\t")
isddb.neg <- negtargetdb[which(negtargetdb$Type=="ISD"),]
targetdb.neg <- negtargetdb[which(negtargetdb$Type=="TargetCPD"),]
targetdb.neg$Cpd <- gsub("_a$|_b$","",targetdb.neg$Cpd)
neg_metadata <- read.delim("negative_mode_metadata.txt", header = F, stringsAsFactors = F, sep="\t")

eiccsvfiles <- datadirfiles[grep("[.]CSV$",datadirfiles)]

cl<-makeCluster(6) ## number of cores goes here, for calculating faster.
#registerDoParallel(cl, cores = 6)
registerDoSNOW(cl)
foreach (i = 1:length(eiccsvfiles)) %dopar% {
  con1 <- file(paste0(datadir,gsub("[.]CSV","",eiccsvfiles[i]),"_ISD_eics.txt"),"w")
  con2 <- file(paste0(datadir,gsub("[.]CSV","",eiccsvfiles[i]),"_target_eics.txt"),"w")

  isintendf <- readLines(paste0(datadir,"/",eiccsvfiles[i]),n=-1L)
  isintendf <- c(isintendf,"-ESI EIC")

  eicanchors <- grep("-ESI EIC",isintendf)
  isdanchors.list <- lapply(1:(length(eicanchors)-1), function(x) {isintendf[eicanchors[x]:(eicanchors[x+1]-1)]  })
  targetnamevec <- sapply(isdanchors.list, function(x) { strsplit(x[1]," <|>:")[[1]][2]})
  targetnamevec <- gsub("ï¿½","",targetnamevec) ## some names has these weird chars. Clean this char in the targeted DB file.

  ### Get the Internal Standards Data

  idf <- isddb.neg
  idf$samplert <- 1.0
  idf$sampleid <- gsub("[.]CSV","",eiccsvfiles[i])
  idf$deltart <- 0.0
  idf$sampleintensity <- 0.0
  for (k in 1:nrow(isddb.neg)) {
    isrt <- idf$RT[k]
    eicdf <- isdanchors.list[[which(targetnamevec==isddb.neg$Cpd[k])[1]]][-c(1:2)]
    eicdf <- as.data.frame(do.call(rbind,lapply(eicdf,function(x) {strsplit(x,",")[[1]]})),stringsAsFactors = F)
    eicdf$V2 <- as.numeric(eicdf$V2)
    eicdf$V3 <- as.numeric(eicdf$V3)
    rtrange <- which(eicdf$V2<(isrt+0.30) & eicdf$V2>(isrt-0.30))
    isrt.sample <- eicdf$V2[rtrange[which.max(eicdf$V3[rtrange])]]
    isinten.sample <- max(eicdf$V3[rtrange])
    idf$samplert[k] <- isrt.sample
    idf$deltart[k] <- isrt.sample-isrt
    idf$sampleintensity[k] <- isinten.sample
  }
  write.table(idf,con1,col.names = T,row.names = F,quote = F, sep="\t")
  ## Get the RT model

  rtmodel <- lm(samplert ~ poly(RT,2), idf)
  new.df <- data.frame(RT = targetdb.neg$RT) # predict the new RTs
  newRtVec <- predict(rtmodel,new.df)

  ## Get the Targeted Lipid Data

  tdf <- targetdb.neg
  tdf$newRT <- newRtVec
  tdf$samplert <- 1.0
  tdf$sampleid <- gsub("[.]CSV","",eiccsvfiles[i])
  tdf$deltart <- 0.0
  tdf$sampleintensity <- 0.0
  tdf$integration <- 1.0
  for (k in 1:nrow(targetdb.neg)) {
    #print(k)
    isrt <- tdf$newRT[k]
    eicdf <- isdanchors.list[[which(targetnamevec==targetdb.neg$Cpd[k])[1]]][-c(1:2)]
    eicdf <- as.data.frame(do.call(rbind,lapply(eicdf,function(x) {strsplit(x,",")[[1]]})),stringsAsFactors = F)
    eicdf$V2 <- as.numeric(eicdf$V2)
    eicdf$V3 <- as.numeric(eicdf$V3)
    rtrange <- which(eicdf$V2<(isrt+0.20) & eicdf$V2>(isrt-0.20))
    isrt.sample <- eicdf$V2[rtrange[which.max(eicdf$V3[rtrange])]]
    isinten.sample <- max(eicdf$V3[rtrange])
    if(isinten.sample==0) {
      rtrange <- which(eicdf$V2<(isrt+0.25) & eicdf$V2>(isrt-0.25)) #
      isrt.sample <- eicdf$V2[rtrange[which.max(eicdf$V3[rtrange])]]
      isinten.sample <- max(eicdf$V3[rtrange])
      is.inten.integration <- max(eicdf$V3[rtrange])/median(eicdf$V3[rtrange])
      tdf$samplert[k] <- isrt.sample
      tdf$deltart[k] <- isrt.sample-isrt
      tdf$sampleintensity[k] <- isinten.sample
      tdf$integration[k] <-  is.inten.integration
    } else {
      is.inten.integration <- max(eicdf$V3[rtrange])/median(eicdf$V3[rtrange])
      tdf$samplert[k] <- isrt.sample
      tdf$deltart[k] <- isrt.sample-isrt
      tdf$sampleintensity[k] <- isinten.sample
      tdf$integration[k] <-  is.inten.integration
    }
  }
  write.table(tdf,con2,col.names = T,row.names = F,quote = F, sep="\t")
  #print(i)
  close(con1)
}
closeAllConnections()

## Import the peak heights and combine them. 

datadirfiles <- dir(datadir)
target_eic_files <- datadirfiles[grep("_target_eics.txt$",datadirfiles)]

neg_inten_list_delta <-lapply(1:length(target_eic_files), function(x) {
  target_data <- read.delim(paste0(datadir,target_eic_files[x]), header = T, stringsAsFactors = F)
  res <- data.frame(intensity=target_data$sampleintensity, rtdelta =target_data$deltart, integration=target_data$integration)
  res
})

# Export the intensities values.
neg_inten_inten <- lapply(1:length(neg_inten_list_delta), function(x) {neg_inten_list_delta[[x]]$intensity   })
neg_inten_inten <- as.data.frame(do.call(rbind, neg_inten_inten),stringsAsFactors = F)
neg_inten_inten <- round(neg_inten_inten,digits = 0)

## These files need to be in same folder.
negtargetdb <- read.delim("esi_negative_mode_targets.txt", header = T, stringsAsFactors = F, sep="\t")
isddb.neg <- negtargetdb[which(negtargetdb$Type=="ISD"),]
targetdb.neg <- negtargetdb[which(negtargetdb$Type=="TargetCPD"),]

neg_inten_inten <- rbind(targetdb.neg$Cpd,neg_inten_inten)
neg_inten_inten <- cbind(data.frame(filename=c("",gsub("_target_eics.txt","",target_eic_files))),neg_inten_inten)
neg_inten_inten$filename <- gsub("[.]2$","",neg_inten_inten$filename)
neg_inten_inten <- neg_inten_inten[- (which(neg_inten_inten$filename%in%neg_metadata$V1==FALSE)[-1]),]

tdf <- rbind(t(neg_inten_inten)[1,],as.character(sapply(neg_inten_inten$filename, function(x) { neg_metadata$V2[which(neg_metadata$V1==x)]})),as.character(sapply(neg_inten_inten$filename, function(x) { neg_metadata$V3[which(neg_metadata$V1==x)]})),t(neg_inten_inten)[-1,])
tdf <- tdf[,-1]

etdf <- cbind(as.character(c("","",neg_inten_inten[1,])),c("","","InChiKeys",rep("",length(neg_inten_inten[1,-1]))) ,c("","SampleLabel","TimeStamp",rep("",length(neg_inten_inten[1,-1]))),tdf   )

write.table(etdf, file="etdf_negative_mode_target_intensities.txt",col.names = F, row.names = F, quote = F, sep="\t")

XLConnect::createSheet(exc,'NEG')
XLConnect::writeWorksheet(exc, etdf, sheet = "NEG", startRow = 1, startCol = 1,header = F)
XLConnect::saveWorkbook(exc)


#####################################################################
### Step 2. LOESS normalization and batch effect removal. ###
#####################################################################

# Copy the etdf_negative_mode_target_intensities.txt, etdf_positive_mode_target_intensities.txt and adni-lipidomics-raw.xlsx to the ADNI1_lipidomics folder. 
# Run the following block of script for the loess normalization. 

##########################################
### LOESS Normalization ##################
### Author Sili Fan slfan@ucdavis.edu ####
##########################################

# functions.
# ---- define batches ---- #
# Batches first defined as time intervals. Five batches defined. For each compound, if a t test shows a non-significant
# result between neighbor batches, then merge these two batches.
# Author Sili Fan, UCDavis May 2017

# BLOCK START - select from here to the line 521 and run. 

DefineBatches = function(e,f,p,
                         time.intval.batch = TRUE, time = 'TimeStamp',number_of_batches=4,
                         batch = NULL,
                         QC.index = p$QC.index,
                         auto.merge.batch = T, p.adjust.method = 'none'){
  # first define 5 batches according to time interval. We define 5 batches.
  if(time.intval.batch){
    time.diff = diff(p[[time]])
    batch.index = which(time.diff%in%sort(time.diff,decreasing=T)[1:(number_of_batches-1)])
    batch = rep(LETTERS[1:number_of_batches],times = diff(c(0,batch.index,nrow(p))))
  }
  # then if auto.merge.batch is true, for each compound, we use t test to merge two neighbor batches. Rule: if p > 0.05, we merge.
  # this makes differnt compounds have different batch.
  batches = matrix(NA,nrow=nrow(p),ncol=nrow(f))
  if(auto.merge.batch){
    for(i in 1:nrow(f)){
      test=pairwise.wilcox.test(as.numeric(e[i,QC.index]), batch[QC.index],p.adjust.method = p.adjust.method)
      test.p = diag(test$p.value)
      batch.opt = LETTERS[1:(length(which(test.p<0.05))+1)]
      batches[,i] = rep(batch.opt,times = diff(c(0,batch.index[which(test.p<0.05)],nrow(p))))
    }
  }else{
    batches = matrix(rep(batch,ncol(batches)),ncol = ncol(batches))
  }
  return(t(batches))
}

get_loess_para = function(x,y,loess.span.limit = 0.5){ # use leave-one-out to select the best span parameter.
  j = 0
  error = rep(0, length(c(seq(loess.span.limit,1.5,0.1),1.75,2,2.25,2.5)))
  if(class(x)=='numeric'|class(x)=='integer'){
    for(par in c(seq(loess.span.limit,1.5,0.1),1.75,2,2.25,2.5)){
      j = j+1
      for(i in 2:(length(y)-1)){ # if i from 1 or to length(y), the prediction would be NA
        o = loess(y[-i]~x[-i],span = par)
        if(sum(is.na(o))){
          error[j] = Inf
        }else{
          err = abs(predict(o,newdata = x[i]) - y[i])
          error[j] = sum(error[j],err,na.rm = T)
        }
      }
    }
  }else if(class(x)=='data.frame'){
    d = data.frame(y,x)
    for(par in c(seq(loess.span.limit,1.5,0.1),1.75,2,2.25,2.5)){
      j = j+1
      for(i in 2:(length(y)-1)){ # if i from 1 or to length(y), the prediction would be NA
        o = loess(y~.,d[-i,],span=par)
        # o = loess(y[-i]~x[-i,],span = par)
        if(sum(is.na(o))){
          error[j] = Inf
        }else{
          err = abs(predict(o,newdata = d[i,colnames(x)]) - y[i])
          error[j] = sum(error[j],err,na.rm = T)
        }
      }
    }
  }else{
    stop("Error: The class of x should be either numeric or matrix")
  }


  return(c(seq(loess.span.limit,1.5,0.1),1.75,2,2.25,2.5)[which.min(error)])
}
remove_outlier = function(v){ # sometimes, the outlier of QC would be disaster when fitting loess. So we have to remove them.
  out = boxplot(v,plot=F)$out
  return(list(value = v[!v%in%out],index = which(v%in%out)))
}
shiftData = function(ori,norm){
  # to make sure that after normalization the data still at the same scale as raw data.
  ori.min = apply(ori,1,min,na.rm=T)
  norm.min = apply(norm,1,min,na.rm=T)
  return(norm - c(norm.min - ori.min))
}
loessNormalization = function(e,f,p,
                              batch = DefineBatches(e,f,p),
                              QC.index = p$QC.index,
                              time = "TimeStamp",
                              span.para = 'auto',loess.span.limit=0.5,
                              divide = F){
  # batch is a matrix. columns are samples and rows are compounds.
  # batch = define_batch(e=e,p=p,f=f)
  # for each compound, calculate the loess line.
  library(parallel)
  cl = makeCluster(10)
  e = data.matrix(e)
  norms = parSapply(cl, X = 1:nrow(e), function(i,e,f,p,QC.index,batch,time,remove_outlier,span.para,get_loess_para,
                                                loess.span.limit){
    models = by(data.frame(v=e[i,QC.index],t=p[[time]][QC.index]),
                batch[i,QC.index],function(x){
                  # x = data.frame(v=e[i,QC.index],t=p[[time]][QC.index])[batch[i,QC.index]=="B",]
                  if(length(remove_outlier(x$v)[[2]])>0){# if outlier exists.
                    span = ifelse(span.para=='auto',
                                  get_loess_para(x=x$t[-remove_outlier(x$v)[[2]]],y=remove_outlier(x$v)[[1]],
                                                 loess.span.limit = loess.span.limit),span.para) # find a proper span.
                  }else{
                    span = ifelse(span.para=='auto',
                                  get_loess_para(x=x$t,y=x$v,
                                                 loess.span.limit = loess.span.limit),span.para) # find a proper span.

                  }
                  if(length(remove_outlier(x$v)[[2]])>0){
                    loess(v~t,data=x[-remove_outlier(x$v)[[2]],],span=span)
                  }else{
                    loess(v~t,data=x,span=span)
                  }
                })
    # predict using the models.
    norm = mapply(function(u,v){
      o = tryCatch({
        predict(u,newdata = v)
      },
      error = function(e){
        print(e)
        v
      })
    },models,by(p[[time]],batch[i,],function(x){x}))
    norm = unlist(norm)

    # replace NA with the closest value.
    if(length(which(is.na(norm)))>0){
      for(j in which(is.na(norm))){
        time_notNA = p[[time]][-which(is.na(norm))]
        closest_time = time_notNA[which.min(abs(time_notNA - p[[time]][j]))]
        norm[j] = norm[which(p[[time]]==closest_time)]
      }
    }
    return(norm)
  },e,f,p,QC.index,batch,time,remove_outlier,span.para,get_loess_para,loess.span.limit)
  norms = t(norms)
  e_norm = matrix(NA,nrow=nrow(f),ncol=nrow(p))
  if(divide){
    for(i in 1:nrow(f)){
      e_norm[i,] = e[i,] / (norms[i,] / median(e[i,]))
    }
  }else{
    for(i in 1:nrow(f)){
      e_norm[i,] = e[i,] - (norms[i,] - median(e[i,]))
    }
  }

  rownames(e_norm) = rownames(e)
  colnames(e_norm) = colnames(e)
  stopCluster(cl)
  return(list(e = shiftData(ori=e,norm = e_norm),f=f,p=p,normalize_line = norms))
}
# ---- read data ---- #
pacman::p_load(readxl, data.table)
# positive mode or negative mode
### Positive Mode Data

POSorNEG = "POS" ##
data = read_excel("adni-lipidomics-raw.xlsx",sheet = ifelse(POSorNEG=="POS",2,3))
f = data.table(data[3:nrow(data),c(1,2,3)]); colnames(f) = c("compound","RT","m.z")
p = data.table(t(data[c(1,2),])); colnames(p) = c("SampleLabel","TimeStamp"); p = p[-c(1,2,3)];
p$index = 1:nrow(p)
order = order(p$TimeStamp)
p$TimeStamp = as.numeric(p$TimeStamp)
p = p[order,]
p[,QC.index:=grepl("QC",SampleLabel)]
# define batches according to time.
timeInterval = which(diff(p[,TimeStamp])%in%sort(diff(p[,TimeStamp]),decreasing = T)[1:3])
batch = rep(LETTERS[1:4],diff(c(0,timeInterval,nrow(p))))
p[,batches:=batch]
e = as.matrix(data[3:nrow(data),4:ncol(data)]); e=apply(e,2,as.numeric)
e = e[,order]
batches=DefineBatches(e,f,p,auto.merge.batch = F)

p1 <- p[order(p$index),]

# loess normalization.
loess = loessNormalization(e,f,p) # model fitting
e_loess = loess$e # matrix of normalized data

p2 <- rbind(t(p1),e_loess[,order(p$index)])

f1 <- do.call(cbind,lapply(f,function(x) { c(rep("",(nrow(p2)-nrow(f))),x)})) ## adding some empty rows to the feature matrix.

p3 <- cbind(f1,p2)
# save normalized data.
fwrite(data.table(p3),paste0(POSorNEG,"_loessNormalized.csv"))

### Negative Modes Data
POSorNEG = "NEG" ##
data = read_excel("adni-lipidomics-raw.xlsx",sheet = ifelse(POSorNEG=="POS",2,3))
f = data.table(data[3:nrow(data),c(1,2,3)]); colnames(f) = c("compound","RT","m.z")
p = data.table(t(data[c(1,2),])); colnames(p) = c("SampleLabel","TimeStamp"); p = p[-c(1,2,3)];
p$index = 1:nrow(p)
order = order(p$TimeStamp)
p$TimeStamp = as.numeric(p$TimeStamp)
p = p[order,]
p[,QC.index:=grepl("QC",SampleLabel)]
# define batches according to time.
timeInterval = which(diff(p[,TimeStamp])%in%sort(diff(p[,TimeStamp]),decreasing = T)[1:3])
batch = rep(LETTERS[1:4],diff(c(0,timeInterval,nrow(p))))
p[,batches:=batch]
e = as.matrix(data[3:nrow(data),4:ncol(data)]); e=apply(e,2,as.numeric)
e = e[,order]
batches=DefineBatches(e,f,p,auto.merge.batch = F)

p1 <- p[order(p$index),]

# loess normalization.
loess = loessNormalization(e,f,p) # model fitting
e_loess = loess$e # matrix of normalized data

p2 <- rbind(t(p1),e_loess[,order(p$index)])

f1 <- do.call(cbind,lapply(f,function(x) { c(rep("",(nrow(p2)-nrow(f))),x)})) ## adding some empty rows to the feature matrix.

p3 <- cbind(f1,p2)
# save normalized data.
fwrite(data.table(p3),paste0(POSorNEG,"_loessNormalized.csv"))

# It has generated two files POS_loessNormalized.csv and NEG_loessNormalized.csv files. 
# Block END


#####################################################################
### STEP 3. A rule based data Filtering to remove duplicate peaks ###
#####################################################################

# We will apply a series of rules to remove the poor quality signals and merge the ESI Pos and Neg mode data into one final data matrix. 

###################
## Positive Mode ##
###################

postargetdb <- read.delim("esi_positive_mode_targets.txt", header = T, stringsAsFactors = F, sep="\t")
isddb.pos <- postargetdb[which(postargetdb$Type=="ISD"),]
isddb.pos <- isddb.pos[-which(isddb.pos$Cpd%in%c("ISD_TG","ISD_PE")==TRUE),] 
targetdb.pos <- postargetdb[which(postargetdb$Type=="TargetCPD"),]
pos_metadata <- read.delim("positive_mode_metadata.txt", header = F, stringsAsFactors = F, sep="\t")

postdf <- read.csv("POS_loessNormalized.csv",header = T,stringsAsFactors = F)
lipidids <- postdf[(6:nrow(postdf)),1]
sampleids <- postdf[1,4:ncol(postdf)]
ikvec <- as.character(sapply(lipidids, function(x) { strsplit(x,"_")[[1]][2]  }))

post.df.sb <- postdf[(6:nrow(postdf)),4:ncol(postdf)]
post.df.sb <- do.call(rbind,lapply(1:nrow(post.df.sb), function(x) { as.numeric(post.df.sb[x,])}))
post.df.sb <- round(as.data.frame(post.df.sb, stringsAsFactors = F),digits = 0)
post.df.qc <- post.df.sb[,grep("^QC",sampleids)]

## outliers length
outliervec <- lapply(1:nrow(post.df.qc), function(x) { boxplot.stats(as.numeric(post.df.qc[x,]),coef = 2)$out } )

## QC RSD calculation. 
qcrsdvec <- sapply(1:nrow(post.df.qc), function(x) {
  if(length(outliervec[[x]])>0 & length(outliervec[[x]]) < 9 ) {
    medvec <- post.df.qc[x,-which(post.df.qc[x,]%in%outliervec[[x]]==TRUE)];
    (sd(as.numeric(medvec))/mean(as.numeric(medvec)))*100
  } else {
    (sd(as.numeric(post.df.qc[x,]))/mean(as.numeric(post.df.qc[x,])))*100
  }
})

## QC Median
qcmedianvec <- sapply(1:nrow(post.df.qc), function(x) { median(as.numeric(post.df.qc[x,])) })

## Sample Median
samplemedian <- sapply(1:nrow(post.df.qc), function(x) { median(as.numeric(post.df.sb[x,-grep("^QC",sampleids)])) }  )
## 1000 Counts
sample1000counts <- sapply(1:nrow(post.df.qc), function(x) { length(which(as.numeric(post.df.sb[x,-grep("^QC",sampleids)])>=1000)) }  )

includevec <- sapply(1:length(ikvec), function(x) {
  len <- length(which(ikvec==ikvec[x]))
  qrsd <- qcrsdvec[x]
  qcmed <- qcmedianvec[x]
  samed <- samplemedian[x]
  c1000 <- sample1000counts[1]
  include <- TRUE
  q1 <- qcrsdvec[which(ikvec==ikvec[x])]
  if(is.na(qrsd)) {
    include <- FALSE
  } else {
    if ( qrsd > 25 ){ # QC RSD is over 25%, exclude it.
      include =FALSE
    } else {
      if(len>1){ 
        if(length(which(q1<25))==length(q1)){ # if all adducts have less than 20% CV then select the compounds with highest median value in QC
            med1 <- qcmedianvec[which(ikvec==ikvec[x])]
            if( which(ikvec==ikvec[x])[which.max(med1)]!=x){
              include = FALSE
            }
        } else {
            if(which(ikvec==ikvec[x])[which.min(q1)]!=x){
              include=FALSE
            }
          }

      }
    }
  }
  include
})

includevec[which(qcmedianvec<2000)] <- FALSE

pos.export.df <- post.df.sb[includevec,]
pos.export.df <- cbind(lipidids[includevec], pos.export.df)
pos.export.df[,1] <- as.character(pos.export.df[,1])
names(pos.export.df) <- c("FullName",as.character(sampleids))
prtvec <- as.numeric(sapply(pos.export.df[,1], function(x) { postargetdb$RT[which(postargetdb$Cpd==x)[1]] }))
pmzvec <- as.numeric(sapply(pos.export.df[,1], function(x) { postargetdb$Mass[which(postargetdb$Cpd==x)[1]] }))
inchikeyvec <- as.character(sapply(pos.export.df[,1], function(x) { strsplit(x,"_")[[1]][2] }))
cpdmdf <- do.call(rbind,lapply(pos.export.df[,1], function(x) { strsplit(x," [[]|[]]")[[1]][c(1:2)]}))
pos.export.df.final <- cbind(LipidName=cpdmdf,MZ=pmzvec,RT=prtvec,InchiKey=inchikeyvec,pos.export.df)
write.table(pos.export.df.final,"adni_pos_lipidomics_loess_with_isd.txt", col.names = T,row.names = F, quote = F, sep="\t")

###################
## Negative mode ##
###################

#### Import of Global Databases
negtargetdb <- read.delim("esi_negative_mode_targets.txt", header = T, stringsAsFactors = F, sep="\t")
isddb.neg <- negtargetdb[which(negtargetdb$Type=="ISD"),]
targetdb.neg <- negtargetdb[which(negtargetdb$Type=="TargetCPD"),]
targetdb.neg$Cpd <- gsub("_a$|_b$","",targetdb.neg$Cpd)
neg_metadata <- read.delim("negative_mode_metadata.txt", header = F, stringsAsFactors = F, sep="\t")


negtdf <- read.csv("NEG_loessNormalized.csv",header = T,stringsAsFactors = F)
lipidids <- negtdf[(6:nrow(negtdf)),1]
sampleids <- negtdf[1,4:ncol(negtdf)]

notinposindex <- which(sampleids%in%c("B-422","QC044","QC023")==TRUE) ## These samples are not present in the positive modes.
sampleids <- sampleids[-notinposindex]
negtdf <- negtdf[,-notinposindex]

ikvec <- as.character(sapply(lipidids, function(x) { strsplit(x,"_")[[1]][2]  }))
negt.df.sb <- negtdf[(6:nrow(negtdf)),4:ncol(negtdf)]
negt.df.sb <- do.call(rbind,lapply(1:nrow(negt.df.sb), function(x) { as.numeric(negt.df.sb[x,])}))
negt.df.sb <- round(as.data.frame(negt.df.sb, stringsAsFactors = F),digits = 0)
negt.df.qc <- negt.df.sb[,grep("^QC",sampleids)]

outliervec <- lapply(1:nrow(negt.df.qc), function(x) { boxplot.stats(as.numeric(negt.df.qc[x,]),coef = 2)$out } )
qcrsdvec <- sapply(1:nrow(negt.df.qc), function(x) {
  if(length(outliervec[[x]])>0 & length(outliervec[[x]]) < 9 ) {
    medvec <- negt.df.qc[x,-which(negt.df.qc[x,]%in%outliervec[[x]]==TRUE)];
    (sd(as.numeric(medvec))/mean(as.numeric(medvec)))*100
  } else {
    (sd(as.numeric(negt.df.qc[x,]))/mean(as.numeric(negt.df.qc[x,])))*100
  }
})

## QC Median
qcmedianvec <- sapply(1:nrow(negt.df.qc), function(x) { median(as.numeric(negt.df.qc[x,])) })

## Sample Median
samplemedian <- sapply(1:nrow(negt.df.qc), function(x) { median(as.numeric(negt.df.sb[x,-grep("^QC",sampleids)])) }  )
## 1000 Counts
sample1000counts <- sapply(1:nrow(negt.df.qc), function(x) { length(which(as.numeric(negt.df.sb[x,-grep("^QC",sampleids)])>=1000)) }  )

includevec <- sapply(1:length(ikvec), function(x) {
  len <- length(which(ikvec==ikvec[x]))
  qrsd <- qcrsdvec[x]
  qcmed <- qcmedianvec[x]
  samed <- samplemedian[x]
  c1000 <- sample1000counts[1]
  include <- TRUE
  q1 <- qcrsdvec[which(ikvec==ikvec[x])]
  if(is.na(qrsd)) {
    include <- FALSE
  } else {
    if ( qrsd > 25 ){ # QC RSD is over 25%, exclude it.
      include =FALSE
    } else {
      if(len>1){ 
        if(length(which(q1<25))==length(q1)){ # if all adducts have less than 20% CV then select the compounds with highest median value in QC
          med1 <- qcmedianvec[which(ikvec==ikvec[x])]
          if( which(ikvec==ikvec[x])[which.max(med1)]!=x){
            include = FALSE
          }
        } else {
          if(which(ikvec==ikvec[x])[which.min(q1)]!=x){
            include=FALSE
          }
        }

      }
    }
  }
  include
})

includevec[which(qcmedianvec<1000)] <- FALSE
neg.export.df <- negt.df.sb[includevec,]
neg.export.df <- cbind(lipidids[includevec], neg.export.df)
neg.export.df[,1] <- as.character(neg.export.df[,1])
names(neg.export.df) <- c("FullName",as.character(sampleids))
prtvec <- as.numeric(sapply(neg.export.df[,1], function(x) { negtargetdb$RT[which(negtargetdb$Cpd==x)[1]] }))
pmzvec <- as.numeric(sapply(neg.export.df[,1], function(x) { negtargetdb$Mass[which(negtargetdb$Cpd==x)[1]] }))
inchikeyvec <- as.character(sapply(neg.export.df[,1], function(x) { strsplit(x,"_")[[1]][2] }))
cpdmdf <- do.call(rbind,lapply(neg.export.df[,1], function(x) { strsplit(x," [[]|[]]")[[1]][c(1:2)]}))
neg.export.df.final <- cbind(LipidName=cpdmdf,MZ=pmzvec,RT=prtvec,InchiKey=inchikeyvec,neg.export.df)
write.table(neg.export.df.final,"adni_neg_lipidomics_loess_with_isd.txt", col.names = T,row.names = F, quote = F, sep="\t")


##############################################################
### STEP 4. Data visualization for the internal standards  ###
###############################################################

# Positive mode
pos_raw  <- read.delim("etdf_positive_mode_target_intensities.txt",header = T, stringsAsFactors = F)
pos_loess <- read.delim("adni_pos_lipidomics_loess_with_isd.txt", header = T, stringsAsFactors = F)
posisdvec <- grep("iSTD",pos_loess[,6])

## Only ISDs
pos_raw.isd <- pos_raw[sapply(posisdvec, function(x) { which(pos_raw$X==pos_loess[x,6])}),grep("QC",pos_raw[1,])]
pos_loess.isd <- pos_loess[posisdvec,grep("QC",names(pos_loess))]

pos_isd_rsd <- cbind(Name=pos_loess$LipidName.1[posisdvec], raw=round(sapply(1:nrow(pos_raw.isd), function(x) { sd(as.numeric(pos_raw.isd[x,]))/mean(as.numeric(pos_raw.isd[x,]))  })*100,digits = 1), loess= round(sapply(1:nrow(pos_loess.isd), function(x) { sd(as.numeric(pos_loess.isd[x,]))/mean(as.numeric(pos_loess.isd[x,]))  })*100,digits = 1) )
write.table(pos_isd_rsd,"adni_pos_isd_qc_rsd.txt",col.names = T, row.names = F,quote = F,sep="\t")

## All Compounds
pos_raw.isd <- pos_raw[sapply(posisdvec, function(x) { which(pos_raw$X==pos_loess[x,6])}),4:ncol(pos_raw)]
pos_loess.isd <- pos_loess[posisdvec,7:ncol(pos_loess)]

pos_isd_rsd <- cbind(Name=pos_loess$LipidName.1[posisdvec], raw=round(sapply(1:nrow(pos_raw.isd), function(x) { sd(as.numeric(pos_raw.isd[x,]))/mean(as.numeric(pos_raw.isd[x,]))  })*100,digits = 1), loess= round(sapply(1:nrow(pos_loess.isd), function(x) { sd(as.numeric(pos_loess.isd[x,]))/mean(as.numeric(pos_loess.isd[x,]))  })*100,digits = 1) )
write.table(pos_isd_rsd,"adni_pos_isd_all_rsd.txt",col.names = T, row.names = F,quote = F,sep="\t") ## RSD of ISD in all the samples.

# negative mode
neg_raw  <- read.delim("etdf_negative_mode_target_intensities.txt",header = T, stringsAsFactors = F)

neg_loess <- read.delim("adni_neg_lipidomics_loess_with_isd.txt", header = T, stringsAsFactors = F)
negisdvec <- grep("iSTD",neg_loess[,6])

neg_raw.isd <- neg_raw[sapply(negisdvec, function(x) { which(neg_raw$X==neg_loess[x,6])}),grep("QC",neg_raw[1,])]
neg_loess.isd <- neg_loess[negisdvec,grep("QC",names(neg_loess))]

neg_isd_rsd <- cbind(Name=gsub("iSTD iSTD ","",neg_loess$LipidName.1[negisdvec]), raw=round(sapply(1:nrow(neg_raw.isd), function(x) { sd(as.numeric(neg_raw.isd[x,]))/mean(as.numeric(neg_raw.isd[x,]))  })*100,digits = 1), loess= round(sapply(1:nrow(neg_loess.isd), function(x) { sd(as.numeric(neg_loess.isd[x,]))/mean(as.numeric(neg_loess.isd[x,]))  })*100,digits = 1) )
write.table(neg_isd_rsd,"adni_neg_isd_qc_rsd.txt",col.names = T, row.names = F,quote = F,sep="\t")

neg_raw.isd <- neg_raw[sapply(negisdvec, function(x) { which(neg_raw$X==neg_loess[x,6])}),4:ncol(neg_raw)]
neg_loess.isd <- neg_loess[negisdvec,7:ncol(neg_loess)]

neg_isd_rsd <- cbind(Name=gsub("iSTD iSTD ","",neg_loess$LipidName.1[negisdvec]), raw=round(sapply(1:nrow(neg_raw.isd), function(x) { sd(as.numeric(neg_raw.isd[x,]))/mean(as.numeric(neg_raw.isd[x,]))  })*100,digits = 1), loess= round(sapply(1:nrow(neg_loess.isd), function(x) { sd(as.numeric(neg_loess.isd[x,]))/mean(as.numeric(neg_loess.isd[x,]))  })*100,digits = 1) )
write.table(neg_isd_rsd,"adni_neg_isd_all_rsd.txt",col.names = T, row.names = F,quote = F,sep="\t")

#################################################################
##### STEP 5: Combining POS and NEG, additional filtering and data visualization####
#################################################################

if(length(table(names(neg_loess)==names(neg_loess)))==1) {print("Positive and Negative mode data have same order for lipids")}

pn_merged.df <- rbind(cbind(ESIMode="POS",pos_loess), cbind(ESIMode="NEG",neg_loess))
pn_merged.df <- pn_merged.df[-grep("iSTD",pn_merged.df$FullName),] ## remove the internal standards
qcindex <- grep("^QC",names(pn_merged.df))

cpdnamevec <- pn_merged.df$LipidName.1
cpd.qc.rsdvec <- sapply(1:nrow(pn_merged.df), function(x) { (sd(as.numeric(pn_merged.df[x,qcindex]))/mean(as.numeric(pn_merged.df[x,qcindex])))*100 })

dup.names <- sapply(1:length(cpdnamevec), function(x) {
  pn_merged.df$FullName[which(cpdnamevec==cpdnamevec[x])]
})

includevec <- sapply(1:length(cpdnamevec), function(x) {
  len <- length(which(cpdnamevec==cpdnamevec[x]))
  include <- TRUE
  if(len>1){
    q1 <- cpd.qc.rsdvec[which(cpdnamevec==cpdnamevec[x])]
    if(which(cpdnamevec==cpdnamevec[x])[which.min(q1)]!=x){
      include=FALSE
    }
  }
  include
})

pn_merged.df.filt <- pn_merged.df[-which(includevec==FALSE),]
write.table(pn_merged.df.filt,"adni_lipidomics_normalized_merged.txt", col.names = T, row.names = F, quote = F, sep = "\t")

######################
##### Final Export ###
######################

lipidnames <- paste0("UCD.Lipid.",c(1:nrow(pn_merged.df.filt)))
datadict.df <- cbind(lipidnames,pn_merged.df.filt[,1:7])
datadict.df$Type <- "Known"
datadict.df$Type[grep("^CSH",datadict.df$FullName)] <- "UnKnown"

datamat <- cbind(lipidnames,pn_merged.df.filt[,-c(1:7)])

#exclude sample
excludesamples <- c("B.280","B.305","B.538","B.289","B.670")
datamat <- datamat[,-which(names(datamat)%in%excludesamples==TRUE)]

# QC RSDs
datamat.rsdvec <- round(sapply(1:nrow(datamat), function(x) { sd(datamat[x,grep("QC",colnames(datamat))])/mean(as.numeric(datamat[x,grep("QC",colnames(datamat))]) )   })*100, digits = 1)
datadict.df$RSD <- datamat.rsdvec

names(datadict.df) <- c("LipidID","ESIMode","LipidName","Adduct","M/Z","RT","InChiKey","FullName","Type","RSD in QCs")

write.table(datadict.df,"adni-lipidomics-data-dict.csv",col.names = T, row.names = F, quote = F,sep = ",")
write.table(datamat,"adni1-lipidomics.csv",col.names = T, row.names = F, quote = F,sep = ",")
write.table(datamat,"adni_lipidomics_final_submit.txt",col.names = T, row.names = F, quote = F,sep = "\t")

####################################
### RSD calculation all compounds ##
####################################

datamat <- read.delim("adni_lipidomics_final_submit.txt", header = T, stringsAsFactors = F)
qcindex <- grep("QC",names(datamat))
post.df.qc <- datamat[,qcindex]
## outliers length
outliervec <- lapply(1:nrow(datamat), function(x) { boxplot.stats(as.numeric(datamat[x,qcindex]),coef = 2)$out } )

qcrsdvec <- sapply(1:nrow(datamat), function(x) {
  if(length(outliervec[[x]])>0 & length(outliervec[[x]]) < 9 ) {
    medvec <- post.df.qc[x,-which(post.df.qc[x,]%in%outliervec[[x]]==TRUE)];
    (sd(as.numeric(medvec))/mean(as.numeric(medvec)))*100
  } else {
    (sd(as.numeric(post.df.qc[x,]))/mean(as.numeric(post.df.qc[x,])))*100
  }
})

write.table(cbind(datamat[,1],qcrsdvec), file="qc_rsd_targets.txt", col.names = F, row.names = F, sep="\t", quote=F)


#########################################
### Principal Component Analysis ########
#########################################

stat_file <- readChar("adni_lipidomics_final_submit.txt", nchars = 100000000)
  library(pacman)
  pacman::p_load(mixOmics,ggplot2, ReporteRs,svglite)
  stat_file <- gsub("\r","",stat_file)
  cfile <- strsplit(stat_file,"\n")[[1]]
  df1 <- as.data.frame(do.call(rbind, lapply(cfile, function (x) { strsplit(x,"\t")[[1]]  } )),stringsAsFactors = F)
  samplenames <- df1[1,][-1]
  sampletype <- rep("Sample",length(samplenames))
  sampletype[grep("QC",samplenames)] <- "QC"
  cpdnames <- df1[,1][-1]
  df2 <- df1[2:nrow(df1),2:ncol(df1)]

  pcamat <- do.call("cbind",lapply(df2,as.numeric))

  pcamat[which(is.na(pcamat)==TRUE)] <- min(pcamat[which(is.na(pcamat)==FALSE)])
  pcamat[which(is.na(pcamat)==TRUE)] <- 100
  pcamat[which(pcamat==0)] <- 100 ## zero need to replaced with some number.
  pcamat <- round(log(pcamat),digits = 5)

  pca.res <- mixOmics::pca(log(pcamat), ncomp=2,max.iter=100,center = T, scale = F)
  pc_scores <- pca.res$loadings[[1]]

  data_bw <- data.frame(snames = as.character(samplenames), stype=sampletype, PC1 = pc_scores[,1], PC2=pc_scores[,2], stringsAsFactors = F)
  data_bw$PC1 <- as.numeric(pc_scores[,1])
  data_bw$PC2 <- as.numeric(pc_scores[,2])

  f2 <- ggplot(data_bw, aes(PC2,PC1, label=snames)) +
    scale_y_continuous(paste("PC2 - variance explained : ",signif(pca.res$explained_variance[2],digit=2)*100," % ",sep="")) +
    scale_x_continuous(paste("PC1 - variance explained : ",signif(pca.res$explained_variance[1],digit=2)*100," % ",sep="")) +
    theme_bw() +
    labs(title = "Principal component analysis (PCA)") +
    theme(
      plot.title = element_text(face="bold", size=20),
      axis.title.x = element_text(face="bold", size=20),
      axis.title.y = element_text(face="bold", size=20, angle=90),
      panel.grid.major = element_blank(), # switch off major gridlines
      panel.grid.minor = element_blank(), # switch off minor gridlines
      panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
      plot.background = element_rect(fill = "transparent",colour = NA),
      legend.position = c(0.9,0.9), # manually position the legend (numbers being from 0,0 at bottom left of whole plot to 1,1 at top right)
      legend.title = element_blank(), # switch off the legend title
      legend.text = element_text(size=20),
      legend.key.size = unit(1.5, "lines"),
      legend.key = element_blank(), # switch off the rectangle around symbols in the legend
      legend.spacing = unit(.05, "cm"),
      axis.text.x = element_text(size=10,angle = 90, hjust = 1),
      axis.text.y = element_text(size=10,angle = 0, hjust = 1)
    )
  f3 <- f2 + geom_point(aes(colour=stype), size=5) #+ geom_text(aes(label=snames),size=5,hjust=0, vjust=0)
  f3
  #f4 <- f3 +  stat_ellipse( aes_string(colour=colnames(data_bw)[2] ), linetype = 1 )
  #plot(f4)
  wbp <- pptx(template = "./data/aggieplot_template_v1.pptx")
  wbp <- addSlide( wbp, "lipidClust" )
  wbp <- addPlot( wbp, function( ) print( f3 ), offx =.1 , offy = .1, width = 10, height = 7 , vector.graphic = TRUE )
  wbp <- addParagraph( wbp, value = c("Principal component analysis shows that there is a clear difference in the two groups compared.") )
  writeDoc( wbp, file = paste0("adni_lipidomics_pca.pptx") )

### Paired ttest analysis of replicated samples
### Duplicate Samples Data Matrix.###

datamat <- read.delim("adni_lipidomics_final_submit.txt", header = T, stringsAsFactors = F)
dupids <- read.delim("duplicate_samples_ids.txt", header = T, stringsAsFactors = F)
dupids$SampleLabel <- gsub("-",".",dupids$SampleLabel)
dupdf <- as.data.frame(do.call(cbind,lapply(dupids$SampleLabel, function(x) { datamat[,which(colnames(datamat)==x)]})),stringsAsFactors =F )
dupdf <- rbind(dupids$Pairs,dupdf)
names(dupdf) <- dupids$SampleLabel
dupdf <- cbind(SampleLabel=c("Pairs",as.character(datamat[,1])), dupdf)
write.table(dupdf, "ADNI_duplicate_samples.txt",col.names = T, row.names = F, quote = F,sep="\t")

dupdf <- read.delim("ADNI_duplicate_samples.txt", header = F, stringsAsFactors = F)
dupdf.mdata <- dupdf[1:2,]
dupdf.mdata <- dupdf.mdata[,-1]
dupdf <- dupdf[-c(1:2),]
lipidnames <- dupdf[,1]
dupdf <- dupdf[,-1]
pairdf <- sapply(names(table(as.numeric(dupdf.mdata[2,]))), function(x) { which(dupdf.mdata[2,]==x)  })

ttest.pvals <- sapply( 1:nrow(dupdf), function(x) {
  res  <- t.test(as.numeric(dupdf[x,pairdf[1,]]),as.numeric(dupdf[x,pairdf[2,]]), paired = T)
  res$p.value
})
write.table(cbind(lipidnames,ttest.pvals), file="pairedpvalues.txt", col.names = F, row.names = F, sep="\t", quote=F)



