### ADNI lipidomics analysis - association modelling study
### Author Dinesh Kumar Barupal March 2019

library(XLConnect)
pacman::p_load(grid,rcdk, RJSONIO,RCurl, dynamicTreeCut,ape,ggplot2,ggrepel,XLConnect,phytools,plotrix,plotly, htmlwidgets,DT,extrafont,XLConnect)
pacman::p_load(officer,openxlsx,grid,RJSONIO,RCurl,dynamicTreeCut,ape,ggplot2,ggrepel,phytools,plotrix,plotly, htmlwidgets,DT,extrafont,devEMF,rvg,magrittr)# better idea will be to load the packages as needed.

adni_samples <- read.delim("sample_details_adni.txt", stringsAsFactors = F, header = T)
lipid_data <- read.delim("lipidomics_raw_data.txt",stringsAsFactors = F, header = T)
lipid_dict <- read.delim("lipid_data_dictionary.txt", stringsAsFactors = F, header = T)

## Added the fish oil intake data.
dr_df <- read.delim("adni_drug_use_data.txt", header = T, stringsAsFactors = F)
adni_samples$FISH <- as.character(sapply(adni_samples$RID, function(x) { dr_df$Fish.Oils[which(dr_df$RID==x)]  }))
adni_samples$FISH[which(adni_samples$FISH!="")] <- "Yes"
adni_samples$FISH[which(adni_samples$FISH=="")] <- "No"

adni_samples.sb <- adni_samples[!duplicated(adni_samples$RID),]
lipid_data <- lipid_data[!duplicated(adni_samples$RID),]

rem_ind <- grep("#",adni_samples.sb$DIAG)

adni_samples.sb <- adni_samples.sb[-rem_ind,]
lipid_data <- lipid_data[-rem_ind,]
adni_samples <- adni_samples.sb

## Added the MCI to AD conversion data
convertdata <- read.delim("conversion_subjects.txt", header = T, stringsAsFactors = F)
adni_samples$MCIConvert  <- sapply(1:nrow(adni_samples), function(x) { convertdata$DXCHANGE[which(convertdata$RID==adni_samples$RID[x])[1]] })



############################
## Calculation starts here #
############################
######### Regression modelling #######
#######################
## Raw models #########
#######################

xvarVec <- c("ABETA142","TAU","SPARE_AD","TOTAL13","APOE4","GENDER","BMI","AGE","DIAG","DIAG","FISH","MCIConvert","PTAU181P")
xvarVec.type <- c("linear","linear","linear","linear","logit","logit","linear","linear","linear","logit","logit","logit","linear")

con1 <- file("ADNI_lipidomics_all_models_raw_082018_v1.txt","w")
writeLines(paste(c("Variable","Starta","ModelType","LipidName","Raw (beta)","Raw (95% CI)","Raw (95% CI)","Raw (pvalue)"),collapse = "\t"), con1)

for(k in 1:length(xvarVec)) {
  for(i in (2:ncol(lipid_data))) {

    sdf <- data.frame(adni_samples[,c(xvarVec[k],"DIAG")],lipid_data[,i], stringsAsFactors = F)
    names(sdf) <- c(xvarVec[k],"DIAG","lipid")
    sdf <- sdf[grep("^[A-Z]",sdf$DIAG),] # If DIAG has no value, we delete it.

    #####################################
    ## Settings for the logistic model.##
    #####################################

    if(xvarVec.type[k]=="logit") {

      if((xvarVec[k]=="DIAG")) {

        sdf <- sdf[,c("DIAG","lipid")]
        sdf <- sdf[which(sdf$DIAG!="LMCI"),]

        sdf$DIAG <-as.factor(sdf$DIAG)
        sdf$DIAG <- relevel(sdf$DIAG, ref= "CN")
        levels(sdf$DIAG) <- c("CN","AD")

        #sdf <- sdf[grep("^[0-9]",sdf[,xvarVec[k]]),] # get the rows having numerical values in the y variable.

        # Get rid of outliers.
        bx.stats <- boxplot(sdf$lipid,range=2,plot=F)
        if( length(bx.stats$out) >0) {
          sdf <- sdf[which(sdf$lipid%in%bx.stats$out==FALSE),]
        }

        # Global model

        sdf.global.mod.raw <- glm(sdf$DIAG ~ scale(sdf$lipid), family=binomial(link='logit') )
        sdf.global.mod.raw.conf <- confint(sdf.global.mod.raw)
        sdf.global.mod.raw.summary <- summary(sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"Global\tLogistic",names(lipid_data)[i],sdf.global.mod.raw.summary$coefficients[2,1],sdf.global.mod.raw.conf[2,],sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)
      }

      if((xvarVec[k]=="MCIConvert")) {
        sdf <- sdf[which(sdf$DIAG=="LMCI"),]
        sdf$MCIConvert[which(sdf$MCIConvert=="Conversion: MCI to Dementia")] <- "Yes"
        sdf$MCIConvert[which(sdf$MCIConvert!="Yes")] <- "No"
        sdf$MCIConvert[is.na(sdf$MCIConvert)] <- "No"
        sdf$MCIConvert <-as.factor(sdf$MCIConvert)
        sdf$MCIConvert <- relevel(sdf$MCIConvert, ref= "No")
        levels(sdf$MCIConvert) <- c("No","Yes")

        bx.stats <- boxplot(sdf$lipid,range=2,plot=F)
        if( length(bx.stats$out) >0) {
          sdf <- sdf[which(sdf$lipid%in%bx.stats$out==FALSE),]
        }

        sdf.global.mod.raw <- glm(sdf$MCIConvert ~ scale(sdf$lipid), family=binomial(link='logit') )
        sdf.global.mod.raw.conf <- confint(sdf.global.mod.raw)
        sdf.global.mod.raw.summary <- summary(sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"Global\tLogistic",names(lipid_data)[i],sdf.global.mod.raw.summary$coefficients[2,1],sdf.global.mod.raw.conf[2,],sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)
      }

      if((xvarVec[k]=="APOE4")) {
        sdf$APOE4[which(sdf$APOE4!=0)] <- "Yes"
        sdf$APOE4[which(sdf$APOE4==0)] <- "No"
        sdf <- sdf[sdf$APOE4%in% c("Yes","No"),]
        sdf$APOE4 <-as.factor(sdf$APOE4)
        sdf$APOE4 <- relevel(sdf$APOE4, ref= "No")
        levels(sdf$APOE4) <- c("No","Yes")
      }

      if((xvarVec[k]=="GENDER")) {
        sdf <- sdf[sdf$GENDER%in% c("Male","Female"),]
        sdf$GENDER <-as.factor(sdf$GENDER)
        sdf$GENDER <- relevel(sdf$GENDER, ref= "Female")
        levels(sdf$GENDER) <- c("Female","Male")
      }

      if((xvarVec[k]=="FISH")) {
        sdf <- sdf[sdf$FISH%in% c("No","Yes"),]
        sdf$FISH <-as.factor(sdf$FISH)
        sdf$FISH <- relevel(sdf$FISH, ref= "No")
        levels(sdf$FISH) <- c("No","Yes")
      }



      if((xvarVec[k]!="DIAG" & xvarVec[k]!="MCIConvert")) {
        # Get rid of outliers.
        bx.stats <- boxplot(sdf$lipid,range=2,plot=F)
        if( length(bx.stats$out) >0) {
          sdf <- sdf[which(sdf$lipid%in%bx.stats$out==FALSE),]
        }
        # Global logistic model

        sdf.global.mod.raw <- glm(sdf[,xvarVec[k]] ~ scale(sdf$lipid ),family=binomial(link='logit'))
        sdf.global.mod.raw.conf <- confint(sdf.global.mod.raw)
        sdf.global.mod.raw.summary <- summary(sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"Global\tLogistic",names(lipid_data)[i],sdf.global.mod.raw.summary$coefficients[2,1],sdf.global.mod.raw.conf[2,],sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)

        # CN model
        cn.sdf  <- sdf[which(sdf$DIAG=="CN"),]
        cn.sdf.global.mod.raw <- glm(cn.sdf[,xvarVec[k]] ~ scale(cn.sdf$lipid ),family=binomial(link='logit'))
        cn.sdf.global.mod.raw.conf <- confint(cn.sdf.global.mod.raw)
        cn.sdf.global.mod.raw.summary <- summary(cn.sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"CN\tLogistic",names(lipid_data)[i],cn.sdf.global.mod.raw.summary$coefficients[2,1],cn.sdf.global.mod.raw.conf[2,],cn.sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)

        # LMCI model
        mci.sdf  <- sdf[which(sdf$DIAG=="LMCI"),]
        mci.sdf.global.mod.raw <- glm(mci.sdf[,xvarVec[k]] ~ scale(mci.sdf$lipid ),family=binomial(link='logit'))
        mci.sdf.global.mod.raw.conf <- confint(mci.sdf.global.mod.raw)
        mci.sdf.global.mod.raw.summary <- summary(mci.sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"LCMI\tLogistic",names(lipid_data)[i],mci.sdf.global.mod.raw.summary$coefficients[2,1],mci.sdf.global.mod.raw.conf[2,],mci.sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)

        # AD model
        ad.sdf  <- sdf[which(sdf$DIAG=="AD"),]
        ad.sdf.global.mod.raw <- glm(ad.sdf[,xvarVec[k]] ~ scale(ad.sdf$lipid ),family=binomial(link='logit'))
        ad.sdf.global.mod.raw.conf <- confint(ad.sdf.global.mod.raw)
        ad.sdf.global.mod.raw.summary <- summary(ad.sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"AD\tLogistic",names(lipid_data)[i],ad.sdf.global.mod.raw.summary$coefficients[2,1],ad.sdf.global.mod.raw.conf[2,],ad.sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)
      }
    }

    #####################################
    ## Settings for the linear model.####
    #####################################

    if(xvarVec.type[k]=="linear") {
      print(k)
      if((xvarVec[k]=="DIAG")) {
        sdf <- sdf[,c("DIAG","lipid")]

        sdf$DIAG[which(sdf$DIAG=="CN")] <- 0
        sdf$DIAG[which(sdf$DIAG=="LMCI")] <- 1
        sdf$DIAG[which(sdf$DIAG=="AD")] <- 2

        sdf <- sdf[sdf$DIAG%in% c(0,1,2),]

        # Get rid of outliers.
        bx.stats <- boxplot(sdf$lipid,range=2,plot=F)
        if( length(bx.stats$out) >0) {
          sdf <- sdf[which(sdf$lipid%in%bx.stats$out==FALSE),]
        }

        # Global model
        sdf.global.mod.raw <- lm(as.numeric(sdf[,xvarVec[k]]) ~ scale(sdf$lipid ))
        sdf.global.mod.raw.conf <- confint(sdf.global.mod.raw)
        sdf.global.mod.raw.summary <- summary(sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"Global\tLinear",names(lipid_data)[i],sdf.global.mod.raw.summary$coefficients[2,1],sdf.global.mod.raw.conf[2,],sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)
      } else {
        sdf <- sdf[grep("^[0-9]",sdf[,xvarVec[k]]),] # get the rows having numerical values in the y variable.

        # Get rid of outliers.
        bx.stats <- boxplot(sdf$lipid,range=2,plot=F)
        if( length(bx.stats$out) >0) {
          sdf <- sdf[which(sdf$lipid%in%bx.stats$out==FALSE),]
        }

        # Global model
        sdf.global.mod.raw <- lm(as.numeric(sdf[,xvarVec[k]]) ~ scale(sdf$lipid ))
        sdf.global.mod.raw.conf <- confint(sdf.global.mod.raw)
        sdf.global.mod.raw.summary <- summary(sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"Global\tLinear",names(lipid_data)[i],sdf.global.mod.raw.summary$coefficients[2,1],sdf.global.mod.raw.conf[2,],sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)

        # CN model
        cn.sdf  <- sdf[which(sdf$DIAG=="CN"),]
        cn.sdf.global.mod.raw <- lm(as.numeric(cn.sdf[,xvarVec[k]]) ~ scale(cn.sdf$lipid ))
        cn.sdf.global.mod.raw.conf <- confint(cn.sdf.global.mod.raw)
        cn.sdf.global.mod.raw.summary <- summary(cn.sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"CN\tLinear",names(lipid_data)[i],cn.sdf.global.mod.raw.summary$coefficients[2,1],cn.sdf.global.mod.raw.conf[2,],cn.sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)

        # LMCI model
        mci.sdf  <- sdf[which(sdf$DIAG=="LMCI"),]
        mci.sdf.global.mod.raw <- lm(as.numeric(mci.sdf[,xvarVec[k]]) ~ scale(mci.sdf$lipid ))
        mci.sdf.global.mod.raw.conf <- confint(mci.sdf.global.mod.raw)
        mci.sdf.global.mod.raw.summary <- summary(mci.sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"LCMI\tLinear",names(lipid_data)[i],mci.sdf.global.mod.raw.summary$coefficients[2,1],mci.sdf.global.mod.raw.conf[2,],mci.sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)

        # AD model
        ad.sdf  <- sdf[which(sdf$DIAG=="AD"),]
        ad.sdf.global.mod.raw <- lm(as.numeric(ad.sdf[,xvarVec[k]]) ~ scale(ad.sdf$lipid ))
        ad.sdf.global.mod.raw.conf <- confint(ad.sdf.global.mod.raw)
        ad.sdf.global.mod.raw.summary <- summary(ad.sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"AD\tLinear",names(lipid_data)[i],ad.sdf.global.mod.raw.summary$coefficients[2,1],ad.sdf.global.mod.raw.conf[2,],ad.sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)
      }
    }
  }
}
close(con1)

#########################
### Adjusted models #####
#########################

xvarVec <- c("PTAU181P","ABETA142","TAU","SPARE_AD","TOTAL13","APOE4","DIAG","DIAG","FISH")
xvarVec.type <- c("linear","linear","linear","linear","linear","logit","linear","logit","logit")

con1 <- file("ADNI_lipidomics_all_models_adj.txt","w")
writeLines(paste(c("Variable","Starta","ModelType","LipidName","Raw (beta)","Raw (95% CI)","Raw (95% CI)","Raw (pvalue)"),collapse = "\t"), con1)

for(k in 1:length(xvarVec)) {
  for(i in (2:ncol(lipid_data))) {

    sdf <- data.frame(adni_samples[,c(xvarVec[k],"DIAG","BMI","GENDER")],lipid_data[,i], stringsAsFactors = F)
    names(sdf) <- c(xvarVec[k],"DIAG","BMI","GENDER","lipid")
    sdf <- sdf[grep("^[A-Z]",sdf$DIAG),] # If DIAG has no value, we delete it.
    sdf <- sdf[grep("^[0-9]",sdf$BMI),] # If BMI has has no value, we delete it.
    sdf$BMI <- as.numeric(sdf$BMI)

    sdf <- sdf[sdf$GENDER%in% c("Male","Female"),]
    sdf$GENDER <-as.factor(sdf$GENDER)
    sdf$GENDER <- relevel(sdf$GENDER, ref= "Female")
    levels(sdf$GENDER) <- c("Female","Male")

    #####################################
    ## Settings for the logistic model.##
    #####################################

    if(xvarVec.type[k]=="logit") {

      if((xvarVec[k]=="DIAG")) {

        #sdf <- sdf[,c("DIAG","lipid")]
        sdf <- sdf[which(sdf$DIAG!="LMCI"),]

        sdf$DIAG <-as.factor(sdf$DIAG)
        sdf$DIAG <- relevel(sdf$DIAG, ref= "CN")
        levels(sdf$DIAG) <- c("CN","AD")

        #sdf <- sdf[grep("^[0-9]",sdf[,xvarVec[k]]),] # get the rows having numerical values in the y variable.

        # Get rid of outliers.
        bx.stats <- boxplot(sdf$lipid,range=2,plot=F)
        if( length(bx.stats$out) >0) {
          sdf <- sdf[which(sdf$lipid%in%bx.stats$out==FALSE),]
        }

        # Global model

        sdf.global.mod.raw <- glm(sdf$DIAG ~ scale(sdf$lipid)+sdf$GENDER+sdf$BMI, family=binomial(link='logit') )
        sdf.global.mod.raw.conf <- confint(sdf.global.mod.raw)
        sdf.global.mod.raw.summary <- summary(sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"Global\tLogistic",names(lipid_data)[i],sdf.global.mod.raw.summary$coefficients[2,1],sdf.global.mod.raw.conf[2,],sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)
      }

      if((xvarVec[k]=="APOE4")) {
        sdf$APOE4[which(sdf$APOE4!=0)] <- "Yes"
        sdf$APOE4[which(sdf$APOE4==0)] <- "No"
        sdf <- sdf[sdf$APOE4%in% c("Yes","No"),]
        sdf$APOE4 <-as.factor(sdf$APOE4)
        sdf$APOE4 <- relevel(sdf$APOE4, ref= "No")
        levels(sdf$APOE4) <- c("No","Yes")
      }

      if((xvarVec[k]=="FISH")) {
        sdf <- sdf[sdf$FISH%in% c("No","Yes"),]
        sdf$FISH <-as.factor(sdf$FISH)
        sdf$FISH <- relevel(sdf$FISH, ref= "No")
        levels(sdf$FISH) <- c("No","Yes")
      }

      if((xvarVec[k]!="DIAG")) {
        # Get rid of outliers.
        bx.stats <- boxplot(sdf$lipid,range=2,plot=F)
        if( length(bx.stats$out) >0) {
          sdf <- sdf[which(sdf$lipid%in%bx.stats$out==FALSE),]
        }
        # Global logistic model

        sdf.global.mod.raw <- glm(sdf[,xvarVec[k]] ~ scale(sdf$lipid ) +sdf$GENDER+sdf$BMI,family=binomial(link='logit'))
        sdf.global.mod.raw.conf <- confint(sdf.global.mod.raw)
        sdf.global.mod.raw.summary <- summary(sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"Global\tLogistic",names(lipid_data)[i],sdf.global.mod.raw.summary$coefficients[2,1],sdf.global.mod.raw.conf[2,],sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)

        # CN model
        cn.sdf  <- sdf[which(sdf$DIAG=="CN"),]
        cn.sdf.global.mod.raw <- glm(cn.sdf[,xvarVec[k]] ~ scale(cn.sdf$lipid )+cn.sdf$GENDER+cn.sdf$BMI,family=binomial(link='logit'))
        cn.sdf.global.mod.raw.conf <- confint(cn.sdf.global.mod.raw)
        cn.sdf.global.mod.raw.summary <- summary(cn.sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"CN\tLogistic",names(lipid_data)[i],cn.sdf.global.mod.raw.summary$coefficients[2,1],cn.sdf.global.mod.raw.conf[2,],cn.sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)

        # LMCI model
        mci.sdf  <- sdf[which(sdf$DIAG=="LMCI"),]
        mci.sdf.global.mod.raw <- glm(mci.sdf[,xvarVec[k]] ~ scale(mci.sdf$lipid )+mci.sdf$GENDER+mci.sdf$BMI,family=binomial(link='logit'))
        mci.sdf.global.mod.raw.conf <- confint(mci.sdf.global.mod.raw)
        mci.sdf.global.mod.raw.summary <- summary(mci.sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"LCMI\tLogistic",names(lipid_data)[i],mci.sdf.global.mod.raw.summary$coefficients[2,1],mci.sdf.global.mod.raw.conf[2,],mci.sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)

        # AD model
        ad.sdf  <- sdf[which(sdf$DIAG=="AD"),]
        ad.sdf.global.mod.raw <- glm(ad.sdf[,xvarVec[k]] ~ scale(ad.sdf$lipid )+ad.sdf$GENDER+ad.sdf$BMI,family=binomial(link='logit'))
        ad.sdf.global.mod.raw.conf <- confint(ad.sdf.global.mod.raw)
        ad.sdf.global.mod.raw.summary <- summary(ad.sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"AD\tLogistic",names(lipid_data)[i],ad.sdf.global.mod.raw.summary$coefficients[2,1],ad.sdf.global.mod.raw.conf[2,],ad.sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)
      }
    }

    #####################################
    ## Settings for the linear model.####
    #####################################

    if(xvarVec.type[k]=="linear") {
      print(k)
      if((xvarVec[k]=="DIAG")) {

        sdf$DIAG[which(sdf$DIAG=="CN")] <- 0
        sdf$DIAG[which(sdf$DIAG=="LMCI")] <- 1
        sdf$DIAG[which(sdf$DIAG=="AD")] <- 2

        sdf <- sdf[sdf$DIAG%in% c(0,1,2),]

        # Get rid of outliers.
        bx.stats <- boxplot(sdf$lipid,range=2,plot=F)
        if( length(bx.stats$out) >0) {
          sdf <- sdf[which(sdf$lipid%in%bx.stats$out==FALSE),]
        }

        # Global model
        sdf.global.mod.raw <- lm(as.numeric(sdf[,xvarVec[k]]) ~ scale(sdf$lipid ) +sdf$GENDER+sdf$BMI)
        sdf.global.mod.raw.conf <- confint(sdf.global.mod.raw)
        sdf.global.mod.raw.summary <- summary(sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"Global\tLinear",names(lipid_data)[i],sdf.global.mod.raw.summary$coefficients[2,1],sdf.global.mod.raw.conf[2,],sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)
      } else {
        sdf <- sdf[grep("^[0-9]",sdf[,xvarVec[k]]),] # get the rows having numerical values in the y variable.

        # Get rid of outliers.
        bx.stats <- boxplot(sdf$lipid,range=2,plot=F)
        if( length(bx.stats$out) >0) {
          sdf <- sdf[which(sdf$lipid%in%bx.stats$out==FALSE),]
        }

        # Global model
        sdf.global.mod.raw <- lm(as.numeric(sdf[,xvarVec[k]]) ~ scale(sdf$lipid )+sdf$GENDER+sdf$BMI)
        sdf.global.mod.raw.conf <- confint(sdf.global.mod.raw)
        sdf.global.mod.raw.summary <- summary(sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"Global\tLinear",names(lipid_data)[i],sdf.global.mod.raw.summary$coefficients[2,1],sdf.global.mod.raw.conf[2,],sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)

        # CN model
        cn.sdf  <- sdf[which(sdf$DIAG=="CN"),]
        cn.sdf.global.mod.raw <- lm(as.numeric(cn.sdf[,xvarVec[k]]) ~ scale(cn.sdf$lipid )+cn.sdf$GENDER+cn.sdf$BMI)
        cn.sdf.global.mod.raw.conf <- confint(cn.sdf.global.mod.raw)
        cn.sdf.global.mod.raw.summary <- summary(cn.sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"CN\tLinear",names(lipid_data)[i],cn.sdf.global.mod.raw.summary$coefficients[2,1],cn.sdf.global.mod.raw.conf[2,],cn.sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)

        # LMCI model
        mci.sdf  <- sdf[which(sdf$DIAG=="LMCI"),]
        mci.sdf.global.mod.raw <- lm(as.numeric(mci.sdf[,xvarVec[k]]) ~ scale(mci.sdf$lipid )+mci.sdf$GENDER+mci.sdf$BMI)
        mci.sdf.global.mod.raw.conf <- confint(mci.sdf.global.mod.raw)
        mci.sdf.global.mod.raw.summary <- summary(mci.sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"LCMI\tLinear",names(lipid_data)[i],mci.sdf.global.mod.raw.summary$coefficients[2,1],mci.sdf.global.mod.raw.conf[2,],mci.sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)

        # AD model
        ad.sdf  <- sdf[which(sdf$DIAG=="AD"),]
        ad.sdf.global.mod.raw <- lm(as.numeric(ad.sdf[,xvarVec[k]]) ~ scale(ad.sdf$lipid )+ad.sdf$GENDER+ad.sdf$BMI)
        ad.sdf.global.mod.raw.conf <- confint(ad.sdf.global.mod.raw)
        ad.sdf.global.mod.raw.summary <- summary(ad.sdf.global.mod.raw)
        writeLines(paste(c(xvarVec[k],"AD\tLinear",names(lipid_data)[i],ad.sdf.global.mod.raw.summary$coefficients[2,1],ad.sdf.global.mod.raw.conf[2,],ad.sdf.global.mod.raw.summary$coefficients[2,4]), collapse = "\t"), con1)
      }
    }
  }
}
close(con1)

adni_models_adj <- read.delim("ADNI_lipidomics_all_models_adj.txt", header = T, stringsAsFactors = F)

#######################
## Summary analysis####
#######################

adni_models_raw <- read.delim("ADNI_lipidomics_all_models_raw_082018_v1.txt", header = T, stringsAsFactors = F)
adni_models_adj <- read.delim("ADNI_lipidomics_all_models_adj.txt", header = T, stringsAsFactors = F)
df1 <- read.delim("spare_ad_results_for_chemrich.txt", header = T, stringsAsFactors = F)

#How many lipid were significant different in the logistic model between AD and CN model.
length(which(adni_models_raw$Raw..pvalue.[which(adni_models_raw$Variable=="DIAG" & adni_models_raw$ModelType=="Logistic")]<0.05)) # 73
length( which( adni_models_raw$Variable=="DIAG" & adni_models_raw$ModelType=="Logistic" & adni_models_raw$Raw..pvalue. <0.05 & adni_models_raw$Raw..beta. < 0 ) ) # 50

############################################
## Chemical Cluster Enrichment Analysis #
###########################################
varVec <- unique(adni_models_raw$Variable)
startavec <- unique(adni_models_raw$Starta)
modelVec <- unique(adni_models_raw$ModelType)

df1 <- read.delim("spare_ad_results_for_chemrich.txt", header = T, stringsAsFactors = F)

con1 <- file("chemical_cluster_enrichment_all_models.txt","w")
writeLines(paste(c("Variable","Starta","ModelType","ChemicalCluster","pvalue","AlteredMetabolites","Cluster size","AdjustedPvalue"), collapse = "\t"), con1)
for(i in 1:length(varVec)) {
 for(j in 1:length(startavec)) {
   for(k in 1:length(modelVec)) {
     sbdf <- adni_models_raw[which( adni_models_raw$Variable==varVec[i] & adni_models_raw$ModelType== modelVec[k] & adni_models_raw$Starta==startavec[j]),]
     if(nrow(sbdf) !=521) next
     df1$pvalue <- sapply(1:nrow(df1), function(x) {sbdf$Raw..pvalue.[which(sbdf$LipidName==df1$LipidName[x])]})
     df1$coefficient <- sapply(1:nrow(df1), function(x) {sbdf$Raw..beta.[which(sbdf$LipidName==df1$LipidName[x])]})
     df1$LOWCI <- sapply(1:nrow(df1), function(x) {sbdf$Raw..95..CI.[which(sbdf$LipidName==df1$LipidName[x])]})
     df1$UPCI <- sapply(1:nrow(df1), function(x) {sbdf$Raw..95..CI..1[which(sbdf$LipidName==df1$LipidName[x])]})

     clusterids <- unique(df1$chemclass)
     clusterids <- names(which(table(df1$chemclass)[clusterids]>2))
     clusterids <- clusterids[which(clusterids!="")] ## get rid of the empty cells.

     cluster.pvalues <- sapply(clusterids, function(x) { # pvalues were calculated if the set has at least 2 metabolites with less than 0.10 pvalue.
       cl.member <- which(df1$chemclass==x)
       if( length(which(df1$pvalue[cl.member]<.05)) >1 ){
         pval.cl.member <- df1$pvalue[cl.member]
         p.test.results <- ks.test(pval.cl.member,"punif",alternative="greater")
         p.test.results$p.value
       } else {
         1
       }
     })

     cluster.pvalues[which(cluster.pvalues==0)] <- 2.2e-20 ### All the zero are rounded to the double.eps pvalues.\
     clusterdf <- data.frame(name=clusterids,pvalues=cluster.pvalues, stringsAsFactors = F)
     clust_s_vec <- sapply(clusterdf$name, function (k) {
       length(which(df1$chemclass==k))
     })

     clusterdf$alteredMetabolites <- sapply(clusterdf$name, function (k) {length(which(df1$chemclass==k & df1$pvalue<0.05))})
     clusterdf$csize <- clust_s_vec
     clusterdf <- clusterdf[which(clusterdf$csize>2),]
     clusterdf$adjustedpvalue <- p.adjust(clusterdf$pvalues, method = "fdr")
     exportdf <- cbind(paste0(varVec[i],"\t", startavec[j],"\t",modelVec[k]), clusterdf)
     write.table(exportdf,con1, col.names = F, row.names = F, quote = F, sep="\t")
   }
 }
}
close(con1)

#############################################################
################ Correlation cluster enrichment analysis ####
#############################################################

df1 <- read.delim("spare_ad_results_for_chemrich.txt", header = T, stringsAsFactors = F)
con1 <- file("correlation_cluster_enrichment_all_models.txt","w")

writeLines(paste(c("Variable","Starta","ModelType","ChemicalCluster","pvalue","AlteredMetabolites","Up Metabolites","Cluster size","AdjustedPvalue"), collapse = "\t"), con1)
for(i in 1:length(varVec)) {
  for(j in 1:length(startavec)) {
    for(k in 1:length(modelVec)) {
      sbdf <- adni_models_raw[which( adni_models_raw$Variable==varVec[i] & adni_models_raw$ModelType== modelVec[k] & adni_models_raw$Starta==startavec[j]),]
      if(nrow(sbdf) !=521) next
      df1$pvalue <- sapply(1:nrow(df1), function(x) {sbdf$Raw..pvalue.[which(sbdf$LipidName==df1$LipidName[x])]})
      df1$coefficient <- sapply(1:nrow(df1), function(x) {sbdf$Raw..beta.[which(sbdf$LipidName==df1$LipidName[x])]})
      df1$LOWCI <- sapply(1:nrow(df1), function(x) {sbdf$Raw..95..CI.[which(sbdf$LipidName==df1$LipidName[x])]})
      df1$UPCI <- sapply(1:nrow(df1), function(x) {sbdf$Raw..95..CI..1[which(sbdf$LipidName==df1$LipidName[x])]})

      clusterids <- unique(df1$corrclass)
      clusterids <- names(which(table(df1$corrclass)[clusterids]>2))
      clusterids <- clusterids[which(clusterids!="")] ## get rid of the empty cells.

      cluster.pvalues <- sapply(clusterids, function(x) { # pvalues were calculated if the set has at least 2 metabolites with less than 0.10 pvalue.
        cl.member <- which(df1$corrclass==x)
        if( length(which(df1$pvalue[cl.member]<.05)) >1 ){
          pval.cl.member <- df1$pvalue[cl.member]
          p.test.results <- ks.test(pval.cl.member,"punif",alternative="greater")
          p.test.results$p.value
        } else {
          1
        }
      })

      cluster.pvalues[which(cluster.pvalues==0)] <- 2.2e-20 ### All the zero are rounded to the double.eps pvalues.\
      clusterdf <- data.frame(name=clusterids,pvalues=cluster.pvalues, stringsAsFactors = F)
      clust_s_vec <- sapply(clusterdf$name, function (x) {
        length(which(df1$corrclass==x))
      })

      clusterdf$alteredMetabolites <- sapply(clusterdf$name, function (x) {length(which(df1$corrclass==x & df1$pvalue<0.05))})
      clusterdf$upmetabolites <- sapply(clusterdf$name, function (x) {
        xdf <- df1[which(df1$corrclass==x),]
        length(which(xdf$coefficient > 0 & xdf$pvalue<0.05))
      })
      clusterdf$csize <- clust_s_vec
      clusterdf <- clusterdf[which(clusterdf$csize>2),]
      clusterdf$adjustedpvalue <- p.adjust(clusterdf$pvalues, method = "fdr")
      exportdf <- cbind(paste0(varVec[i],"\t", startavec[j],"\t",modelVec[k]), clusterdf)
      write.table(exportdf,con1, col.names = F, row.names = F, quote = F, sep="\t")
    }
  }
}
close(con1)
