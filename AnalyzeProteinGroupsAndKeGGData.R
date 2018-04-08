#options(warn=1)

library(gplots)
library(scales)
library(RColorBrewer)
library(VennDiagram)
library(ReactomePA)
library(corrplot)

source("/Users/dybasj/JoeRLib/ProteomicsAnalysisFunctions.R")
source("/Users/dybasj/JoeRLib/PlottingFunctions.R")


##enter data
##read ensembl to uniprot mapping
ProtIDdata <- read.table("./ProteinData_UniprotEnsemblGeneID.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors=FALSE)
#print(ProtIDdata)
#head(ProtIDdata)
#dim(ProtIDdata)

##read iBAQ quantification protein IDs
WCPiBAQProtIDdata <- read.table("./datafiles/ProteinGroupData_ProteinIdGeneIdProteinName_LFQ.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors=FALSE)
#print(WCPiBAQProtIDdata)
#head(WCPiBAQProtIDdata)
#dim(WCPiBAQProtIDdata)

##read whole cell proteome iBAQ Log2 transformed, normalized quantification data
#formatted with experiments as columns and samples as rows
WCPiBAQdata <- read.table("./datafiles/ProteinGroupData_iBAQ_LFQ_Log2Norm.txt", sep="\t", row.names=1, header=TRUE)
#print(WCPiBAQdata)
#head(WCPiBAQdata)
#dim(WCPiBAQdata)

##read whole cell proteome raw MSMS counts from LFQ experiments
#formatted with experiments as columns and samples as rows
WCPMSMScountdata <- read.table("./datafiles/ProteinGroupData_MSMSCounts_LFQ_Raw.txt", sep="\t", row.names=1, header=TRUE)
#change "NA" counts to 0
WCPMSMScountdata[WCPMSMScountdata=='NA'] <- NA
WCPMSMScountdata[is.na(WCPMSMScountdata)] <- 0
#change column names to specify MSMScount
WCPMSMScountdatacolnames <- colnames(WCPMSMScountdata) 
WCPMSMScountdatacolnamesadj <- paste("MSMScounts", WCPMSMScountdatacolnames, sep=".")
colnames(WCPMSMScountdata) <- WCPMSMScountdatacolnamesadj
#print(WCPMSMScountdata)
#head(WCPMSMScountdata)
#dim(WCPMSMScountdata)

##read SILAC quantification protein IDs
WCPSILAQProtIDdata <- read.table("./datafiles/ProteinGroupData_ProteinIdGeneIdProteinName_SILAC.txt", sep="\t", row.names=1, header=TRUE)
#print(WCPSILAQProtIDdata)
#head(WCPSILAQProtIDdata)
#dim(WCPSILAQProtIDdata)

##read whole cell proteome normHLratios Log2 transformed
#formatted with experiments as columns and samples as rows
WCPnormLHdata <- read.table("./datafiles/ProteinGroupData_RatioHLnorm_SILAC_Log2.txt", sep="\t", row.names=1, header=TRUE)
#since the MaxQuant ratios are H/L and L/H is informative for fold changes, multiply log2 scores by -1 (inverse)
WCPnormLHdata <- -(WCPnormLHdata)
#change column names to specify normLHratio
normLHratiodatacolnames <- colnames(WCPnormLHdata)
normLHratiodatacolnamesadj <- paste("normLHratio", normLHratiodatacolnames, sep=".")
colnames(WCPnormLHdata) <- normLHratiodatacolnamesadj
#print(WCPnormLHdata)
#head(WCPnormLHdata)
#dim(WCPnormLHdata)

##read whole cell proteome raw RatioHLcount from SILAC experiments
#formatted with experiments as columns and samples as rows
WCPLHcountdata <- read.table("./datafiles/ProteinGroupData_RatioHLcount_SILAC_Raw.txt", sep="\t", row.names=1, header=TRUE)
#change "NA" counts to 0
WCPLHcountdata[WCPLHcountdata=='NA'] <- NA
WCPLHcountdata[is.na(WCPLHcountdata)] <- 0
#change column names to specify HLcount
WCPLHcountdatacolnames <- colnames(WCPLHcountdata)
WCPLHcountdatacolnamesadj <- paste("LHcounts", WCPLHcountdatacolnames, sep=".")
colnames(WCPLHcountdata) <- WCPLHcountdatacolnamesadj
#print(WCPLHcountdata)
#head(WCPLHcountdata)
#dim(WCPLHcountdata)

##read whole cell proteome total intensity data
WCPnormIntensitydata <- read.table("./datafiles/ProteinGroupData_Intensity_SILAC_Log2Norm.txt", sep="\t", row.names=1, header=TRUE)
#change column names to specify Intensity
intensitydatacolnames <- colnames(WCPnormIntensitydata)
intensitydatacolnamesadj <- paste("Intensity", intensitydatacolnames, sep=".")
colnames(WCPnormIntensitydata) <- intensitydatacolnamesadj
#print(WCPnormIntensitydata)
#head(WCPnormIntensitydata)
#dim(WCPnormIntensitydata)

##read whole cell proteome light intensity data
WCPnormIntensityLdata <- read.table("./datafiles/ProteinGroupData_IntensityL_SILAC_Log2Norm.txt", sep="\t", row.names=1, header=TRUE)
#change column names to specify IntensityL
intensityLdatacolnames <- colnames(WCPnormIntensityLdata)
intensityLdatacolnamesadj <- paste("IntensityL", intensityLdatacolnames, sep=".")
colnames(WCPnormIntensityLdata) <- intensityLdatacolnamesadj
#print(WCPnormIntensityLdata)
#head(WCPnormIntensityLdata)
#dim(WCPnormIntensityLdata)

##read whole cell proteome heavy intensity data
WCPnormIntensityHdata <- read.table("./datafiles/ProteinGroupData_IntensityH_SILAC_Log2Norm.txt", sep="\t", row.names=1, header=TRUE)
#change column names to specify IntensityH
intensityHdatacolnames <- colnames(WCPnormIntensityHdata)
intensityHdatacolnamesadj <- paste("IntensityH", intensityHdatacolnames, sep=".")
colnames(WCPnormIntensityHdata) <- intensityHdatacolnamesadj
#print(WCPnormIntensityHdata)
#head(WCPnormIntensityHdata)
#dim(WCPnormIntensityHdata)

#combine dataframes
WCPstimQuant <- merge(WCPiBAQProtIDdata, WCPSILAQProtIDdata, by=c("row.names", "ProtName", "GeneId"), all=TRUE)
rownames(WCPstimQuant)=WCPstimQuant$Row.names
WCPstimQuant$Row.names <- NULL
dim(WCPstimQuant)
WCPstimQuant <- merge(WCPstimQuant, WCPnormIntensityLdata, by="row.names", all=TRUE)
rownames(WCPstimQuant)=WCPstimQuant$Row.names
WCPstimQuant$Row.names <- NULL
dim(WCPstimQuant)
WCPstimQuant <- merge(WCPstimQuant, WCPnormIntensityHdata, by="row.names", all=TRUE)
rownames(WCPstimQuant)=WCPstimQuant$Row.names
WCPstimQuant$Row.names <- NULL
dim(WCPstimQuant)
WCPstimQuant <- merge(WCPstimQuant, WCPnormIntensitydata, by="row.names", all=TRUE)
rownames(WCPstimQuant)=WCPstimQuant$Row.names
WCPstimQuant$Row.names <- NULL
dim(WCPstimQuant)
WCPstimQuant <- merge(WCPstimQuant, WCPnormLHdata, by="row.names", all=TRUE)
rownames(WCPstimQuant)=WCPstimQuant$Row.names
WCPstimQuant$Row.names <- NULL
dim(WCPstimQuant)
WCPstimQuant <- merge(WCPstimQuant, WCPLHcountdata, by="row.names", all=TRUE)
rownames(WCPstimQuant)=WCPstimQuant$Row.names
WCPstimQuant$Row.names <- NULL
dim(WCPstimQuant)
WCPstimQuant <- merge(WCPstimQuant, WCPiBAQdata, by="row.names", all=TRUE)
rownames(WCPstimQuant)=WCPstimQuant$Row.names
WCPstimQuant$Row.names <- NULL
dim(WCPstimQuant)
WCPstimQuant <- merge(WCPstimQuant, WCPMSMScountdata, by="row.names", all=TRUE)
rownames(WCPstimQuant)=WCPstimQuant$Row.names
WCPstimQuant$Row.names <- NULL
dim(WCPstimQuant)
WCPstimQuant <- merge(WCPstimQuant, ProtIDdata, by="row.names", all=TRUE)
rownames(WCPstimQuant)=WCPstimQuant$Row.names
WCPstimQuant$Row.names <- NULL
dim(WCPstimQuant)
#head(WCPstimQuant)
#str(WCPstimQuant)

###merging causes some added lines to have NA where there should be 0
##change NA to 0 where necessary
#changezero<-c("LHcounts.PO1063_HLM", "LHcounts.PO1063_HLP", "LHcounts.PO1075_K48P", "LHcounts.PO1075_PANM", "LHcounts.PO1075_WCPM", "LHcounts.PO1075_WCPP", "MSMScounts.PO972igd_0hr", "MSMScounts.PO972igd_1hr", "MSMScounts.PO972igd_4hr")
#WCPstimQuant[changezero][is.na(WCPstimQuant[changezero])] <- 0
### end change NA to 0 where necessary

#add uniprot ids from row names
uniprotids<-rownames(WCPstimQuant)
WCPstimQuant$UniprotID<-uniprotids
#dim(WCPstimQuant)
#head(WCPstimQuant)


###read KeGG protein-level Ub quantification
########################################KeGGnormLHdata <- read.table("./KeGGproteinData_PeptideBasedTotalProteinQuantification_PO1063-PO972.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors=FALSE)
####################KeGGnormLHdata <- read.table("./KeGGproteinData_PeptideBasedTotalProteinQuantification_PO1063-PO972_0810_NEW.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors=FALSE)
KeGGnormLHdata <- read.table("./KeGGproteinData_PeptideBasedTotalProteinQuantification_PO1063-PO972_0810_NewWeighting.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors=FALSE)
#head(KeGGnormLHdata)
#dim(KeGGnormLHdata)
KeGG_AllUbProtId_Unfiltered<-row.names(subset(KeGGnormLHdata, (KeGGnormLHdata$numsites.PO1063kgg_HLM>0 | KeGGnormLHdata$numsites.PO972kgg_0hr>0 | KeGGnormLHdata$numsites.PO972kgg_4hr>0)))
KeGG_AllUbProtId_Unfiltered_Num <-length(KeGG_AllUbProtId_Unfiltered)
print("all KeGG prots Unfiltered")
print(KeGG_AllUbProtId_Unfiltered_Num)

##remove some erroneous proteins...
print("remove proteins...")
RemoveProteins <- c("B1AT92")
#print(RemoveProteins)
WCPstimQuant <- WCPstimQuant[!row.names(WCPstimQuant)%in%RemoveProteins,]
#dim(WCPstimQuant)


##KeGG analysis average PO972 and PO1063 0-4hr ratios
KeGGnormLHdata$UbFcAvg<-rowMeans(KeGGnormLHdata[,c("intensityweightedmeanratio.PO1063kgg_HLM","intensityweightedmeanratio.PO972kgg_4diff0")], na.rm=TRUE)
#head(KeGGnormLHdata)
#dim(KeGGnormLHdata)


##MSMS counts
#min PO972 0hr-4hr MSMS count
WCPstimQuant$IGD0hr4hr_MinMSMScount <- apply(WCPstimQuant, MARGIN=1, FUN=function(x) as.numeric(min(x[["MSMScounts.PO972igd_0hr"]],x[["MSMScounts.PO972igd_4hr"]])))
#max PO972 0hr-4hr MSMS count
WCPstimQuant$IGD0hr4hr_MaxMSMScount <- apply(WCPstimQuant, MARGIN=1, FUN=function(x) as.numeric(max(x[["MSMScounts.PO972igd_0hr"]],x[["MSMScounts.PO972igd_4hr"]])))


##reactome interactions
#read reactome interaction data
WCPnetworkdata <- read.table("TCRUbDataset_ReactomeNetwork_NumInteractions.txt", sep="\t", row.names=1, header=TRUE)
#print(WCPnetworkdata)
#head(WCPnetworkdata)
#dim(WCPnetworkdata)


##RNAseq data input and merge
#
#
###SHOULD COMPUTE RPKM FOR THE RNAseq DATA
#
#
print("RNAseq analysis")
RnaData<-read.table("./RNAseq/wt1-3unstim_wt1-3stim_Results_I.txt", sep=",", row.names=1, header=TRUE, stringsAsFactors=FALSE)
ensemblids<-row.names(RnaData)
RnaData$EnsemblId<-ensemblids
#log10 transform of base mean count values
RnaData$log10baseMean <- apply(RnaData, MARGIN=1, FUN=function(x) log10(as.numeric(x[["baseMean"]])))
#change value of log2FoldChange to reflect intuitive way of reading stim vs unstim
RnaData$log2FoldChangeStimUnstim<-apply(RnaData, MARGIN=1, FUN=function(x) (as.numeric(x[["log2FoldChange"]])*-1))
#calc neg log10 pval for volcano plots
RnaData$neglog10padj<-apply(RnaData, MARGIN=1, FUN=function(x) (log10(as.numeric(x[["padj"]])))*-1)
RnaData$neglog10padj[RnaData$neglog10padj == "Inf"] <- NA
#head(RnaData)
#
##mearge RnaData with WCPstimQuant data frame, keeping only transcripts w/ Wcp protein expression
WCPstimQuant <- merge(WCPstimQuant, RnaData, by="EnsemblId", all.x=TRUE)
row.names(WCPstimQuant) <- WCPstimQuant$UniprotID
WCPstimQuant$UniprotID<-NULL
#head(WCPstimQuant)
#dim(WCPstimQuant)
##merge WCPstimQuant with RnaData dataframe, keeping all transcripts, not just the transcripts w/ protein expression
RnaData <- merge(RnaData, WCPstimQuant, by="EnsemblId", all.x=TRUE)
#row.names(RnaData)<-RnaData$EnsemblId
#RnaData$EnsemblId<-NULL
#head(RnaData)
#dim(RnaData)
##end RNAseq data input and merge

###
print("merge")
#mearge KeGG data with WCPstimQuant data frame
WCPstimQuant <- merge(WCPstimQuant, KeGGnormLHdata, by="row.names", all.x=TRUE)
rownames(WCPstimQuant)=WCPstimQuant$Row.names
WCPstimQuant$Row.names <- NULL
#head(WCPstimQuant)
#dim(WCPstimQuant)
###

print("intro NA")
##verify that all "NA" are NA
WCPstimQuant[WCPstimQuant=='NA'] <- NA
##


##merging causes some added lines to have NA where there should be 0
#change NA to 0 where necessary
changezero<-c("LHcounts.PO1063_HLM", "LHcounts.PO1063_HLP", "LHcounts.PO1075_K48P", "LHcounts.PO1075_PANM", "LHcounts.PO1075_WCPM", "LHcounts.PO1075_WCPP", "MSMScounts.PO972igd_0hr", "MSMScounts.PO972igd_1hr", "MSMScounts.PO972igd_4hr")
WCPstimQuant[changezero][is.na(WCPstimQuant[changezero])] <- 0
#head(WCPstimQuant)
## end change NA to 0 where necessary

print("zscore")
##z-scores for intensities
#experiment (column) means and standard deviations
WCPstimQuantSub<-subset(WCPstimQuant, select=c(PO972igd_0hr, PO972igd_1hr, PO972igd_4hr, Intensity.PO1063_HLM, IntensityH.PO1063_HLM, IntensityL.PO1063_HLM, Intensity.PO1075_WCPM, IntensityH.PO1075_WCPM, IntensityL.PO1075_WCPM))
WCPstimQuant_colmeans <- apply(WCPstimQuantSub, MARGIN=2, FUN=function(x) mean(x, na.rm=TRUE))
WCPstimQuant_colstdev <- apply(WCPstimQuantSub, 2, sd, na.rm=TRUE)
#kegg igd 0,1,4 hr
WCPstimQuant$IGD0hr_zscore <- sapply(WCPstimQuant$PO972igd_0hr, FUN=function(x) zscores(x,WCPstimQuant_colmeans[["PO972igd_0hr"]],WCPstimQuant_colstdev[["PO972igd_0hr"]]))
WCPstimQuant$IGD1hr_zscore <- sapply(WCPstimQuant$PO972igd_1hr, FUN=function(x) zscores(x,WCPstimQuant_colmeans[["PO972igd_1hr"]],WCPstimQuant_colstdev[["PO972igd_1hr"]]))
WCPstimQuant$IGD4hr_zscore <- sapply(WCPstimQuant$PO972igd_4hr, FUN=function(x) zscores(x,WCPstimQuant_colmeans[["PO972igd_4hr"]],WCPstimQuant_colstdev[["PO972igd_4hr"]]))
#nedd total, heavy and light (0-4hr) intensity
WCPstimQuant$HLMint_zscore <- sapply(WCPstimQuant$Intensity.PO1063_HLM, FUN=function(x) zscores(x,WCPstimQuant_colmeans[["Intensity.PO1063_HLM"]],WCPstimQuant_colstdev[["Intensity.PO1063_HLM"]]))
WCPstimQuant$HLMintH_zscore <- sapply(WCPstimQuant$IntensityH.PO1063_HLM, FUN=function(x) zscores(x,WCPstimQuant_colmeans[["IntensityH.PO1063_HLM"]],WCPstimQuant_colstdev[["IntensityH.PO1063_HLM"]]))
WCPstimQuant$HLMintL_zscore <- sapply(WCPstimQuant$IntensityL.PO1063_HLM, FUN=function(x) zscores(x,WCPstimQuant_colmeans[["IntensityL.PO1063_HLM"]],WCPstimQuant_colstdev[["IntensityL.PO1063_HLM"]]))
#tube total, heavy and light (0-4hr) intensity
WCPstimQuant$WCPMint_zscore <- sapply(WCPstimQuant$Intensity.PO1075_WCPM, FUN=function(x) zscores(x,WCPstimQuant_colmeans[["Intensity.PO1075_WCPM"]],WCPstimQuant_colstdev[["Intensity.PO1075_WCPM"]]))
WCPstimQuant$WCPMintH_zscore <- sapply(WCPstimQuant$IntensityH.PO1075_WCPM, FUN=function(x) zscores(x,WCPstimQuant_colmeans[["IntensityH.PO1075_WCPM"]],WCPstimQuant_colstdev[["IntensityH.PO1075_WCPM"]]))
WCPstimQuant$WCPMintL_zscore <- sapply(WCPstimQuant$IntensityL.PO1075_WCPM, FUN=function(x) zscores(x,WCPstimQuant_colmeans[["IntensityL.PO1075_WCPM"]],WCPstimQuant_colstdev[["IntensityL.PO1075_WCPM"]]))
#average intensity zscore 0hr
WCPstimQuant$AbundanceInt_zscore0hr <- apply(WCPstimQuant, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["IGD0hr_zscore"]],x[["HLMintH_zscore"]],x[["WCPMintH_zscore"]])), na.rm=TRUE))
#average intensity zscore 4hr
WCPstimQuant$AbundanceInt_zscore4hr <- apply(WCPstimQuant, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["IGD4hr_zscore"]],x[["HLMintL_zscore"]],x[["WCPMintL_zscore"]])), na.rm=TRUE))
#average intensities for igd 0,1,4hr and HLMint and WCPMing
WCPstimQuant$AbundanceInt_zscore <- apply(WCPstimQuant, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["IGD0hr_zscore"]],x[["IGD4hr_zscore"]],x[["HLMint_zscore"]],x[["WCPMint_zscore"]])), na.rm=TRUE))


##fold changes for whole cell proteome stim
print("fold changes")
WCPstimQuant$IGDStim1hr_fc<-mapply(foldchange, WCPstimQuant$PO972igd_0hr, WCPstimQuant$PO972igd_1hr)
WCPstimQuant$IGDStim4hr_fc<-mapply(foldchange, WCPstimQuant$PO972igd_0hr, WCPstimQuant$PO972igd_4hr)
WCPstimQuant$IGD1hr4hr_fc<-mapply(foldchange, WCPstimQuant$PO972igd_1hr, WCPstimQuant$PO972igd_4hr)
WCPstimQuant$HLMStim4hr_fc<-WCPstimQuant$normLHratio.PO1063_HLM
WCPstimQuant$HLPStim4hr_fc<-WCPstimQuant$normLHratio.PO1063_HLP
WCPstimQuant$WCPMStim4hr_fc<-WCPstimQuant$normLHratio.PO1075_WCPM
##average fold change
WCPstimQuant$Stim4hr_AvgFc<-rowMeans(WCPstimQuant[,c("IGDStim4hr_fc","HLMStim4hr_fc","WCPMStim4hr_fc")], na.rm=TRUE)

print("coeff var")
##fold change coefficient of variances
WCPstimQuant$Stim4hr_FcCovar <- apply(WCPstimQuant, MARGIN=1, FUN=function(x) coeffvar(as.numeric(c(x[["IGDStim4hr_fc"]],x[["HLMStim4hr_fc"]],x[["WCPMStim4hr_fc"]]))))
#head(WCPstimQuant)

print("nedd inhibitor diff reg")
##differential regulation of protein expression in nedd inhibitor dataset
WCPstimQuant$HLMHLPStim4hr_DiffReg <- apply(WCPstimQuant, MARGIN=1, FUN=function(x) (as.numeric(x[["HLPStim4hr_fc"]])-as.numeric(x[["HLMStim4hr_fc"]])))

print("ttests")
##t-tests for whole cell proteome stim
WCPstimQuant$Stim4hr_TtestPval<-apply(WCPstimQuant, MARGIN=1, FUN=function(x) onesetsigtest(as.numeric(c(x[["IGDStim4hr_fc"]],x[["HLMStim4hr_fc"]],x[["WCPMStim4hr_fc"]]))))
WCPstimQuant$Stim4hr_NegLog10Pval<-sapply(WCPstimQuant$Stim4hr_TtestPval, FUN=function(x) -log10(x))

#head(WCPstimQuant)


if ( FALSE )
{
##significance tests comparing distributions of fold changes based on counts
for ( i in 1:20 ){
	print(i)
	###for ( j in 0:10 ){
		###print(j)
		p<-i-1
		prevset<-WCPstimQuant$IGD1hr4hr_fc[WCPstimQuant$IGD0hr4hr_MaxMSMScount==p]
		prevset<-prevset[!is.na(prevset)]
		varprevset<-var(prevset)
		currentset<-WCPstimQuant$IGD1hr4hr_fc[WCPstimQuant$IGD0hr4hr_MaxMSMScount==i]
		currentset<-currentset[!is.na(currentset)]
		#iqrcurrenset<-IQR(currentset)
		#print(iqrcurrenset)
		varcurrentset<-var(currentset)
		#print(varcurrentset)
		vardiff<-varprevset-varcurrentset
		print(vardiff)

		#cutoffsetdist<-WCPstimQuant$IGD1hr4hr_fc[WCPstimQuant$IGD0hr4hr_MaxMSMScount<=i]
		#cutoffsetdist<-cutoffsetdist[!is.na(cutoffsetdist)]
		#cutoffsetsize<-length(cutoffsetdist)
		
		#####retainsetdist<-WCPstimQuant$IGD1hr4hr_fc[WCPstimQuant$IGD0hr4hr_MaxMSMScount>=i]
		#####retainsetdist<-WCPstimQuant$HLMStim4hr_fc[WCPstimQuant$LHcounts.PO1063_HLM>=i]
		
		#retainsetdist<-retainsetdist[!is.na(retainsetdist)]
		#retainsetsize<-length(retainsetdist)
		#retainsetdistsubset<-sample(retainsetdist, cutoffsetsize, replace = FALSE)

		#varcutoff<-var(cutoffsetdist)
		#varretained<-var(retainsetdistsubset)

		#print(varcutoff)
		#varretained<-var(retainsetdist, na.rm=TRUE)
		#print(varretained)

		#testres<-ks.test(zerocountdist,onecountdist)
		#print(testres)
		#xyTtest<-t.test(zerocountdist,onecountdist,paired=FALSE,var.equal=FALSE,na.action="na.omit")
		#xyTtestPval<-xyTtest$p.value
		#print(xyTtestPval)
		#xyFtest<-var.test(zerocountdist,onecountdist,na.action="na.omit")
		#xyFtestPval<-xyFtest$p.value
		#print(xyFtestPval)
	###}
	#cutoffvar<-var(zerocountdist, na.rm=TRUE)
	#compvar<-var(onecountdist, na.rm=TRUE)
	#print(cutoffvar)
	#print(compvar)

}

}

####
##mearge KeGG data with WCPstimQuant data frame
#WCPstimQuant <- merge(WCPstimQuant, KeGGnormLHdata, by="row.names", all.x=TRUE)
#rownames(WCPstimQuant)=WCPstimQuant$Row.names
#WCPstimQuant$Row.names <- NULL
##head(WCPstimQuant)
##dim(WCPstimQuant)
######
#
#
###verify that all "NA" are NA
#WCPstimQuant[WCPstimQuant=='NA'] <- NA




##set directory where analysis files are written
analysisdir<-"analysisoutput_current"


##write output file for complete data analysis data frame
write.table(WCPstimQuant, file=file.path(".", analysisdir, "WholeCellProteomeTCRstim_PO972-PO1075-PO1063_MassSpecRnaSeq_NotFiltered.txt"), sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
head(WCPstimQuant)


###set the significance threshold for the subsequent analysis
#ProtSigThresh<-0.1
ProtSigThresh<-0.05
RnaSigThresh<-0.01

print("set analysis")

###quantified set analysis for unfiltered data
#PO972 proteins identified at 0 hour only
PO972_0hr_Unique_ProtId_Unfiltered <- row.names(subset(WCPstimQuant, ((!is.na(WCPstimQuant$PO972igd_0hr) & WCPstimQuant$MSMScounts.PO972igd_0hr>0) & is.na(WCPstimQuant$IGDStim4hr_fc))))
PO972_0hr_Unique_ProtId_Unfiltered_Num <- length(PO972_0hr_Unique_ProtId_Unfiltered)
print("PO972 0hr unique (unfiltered):")
print(PO972_0hr_Unique_ProtId_Unfiltered_Num)
#
#PO972 proteins identified at 4 hour only
PO972_4hr_Unique_ProtId_Unfiltered <- row.names(subset(WCPstimQuant, ((!is.na(WCPstimQuant$PO972igd_4hr) & WCPstimQuant$MSMScounts.PO972igd_4hr>0) & is.na(WCPstimQuant$IGDStim4hr_fc))))
PO972_4hr_Unique_ProtId_Unfiltered_Num <- length(PO972_4hr_Unique_ProtId_Unfiltered)
print("PO972 4hr unique (unfiltered):")
print(PO972_4hr_Unique_ProtId_Unfiltered_Num)
#
#PO972 0-4hr stim quantified proteins. proteins identified in BOTH 0 AND 4 hours
PO972_0hr4hr_ProtId_Unfiltered <- row.names(subset(WCPstimQuant, !is.na(WCPstimQuant$IGDStim4hr_fc)))
PO972_0hr4hr_ProtId_Unfiltered_Num <- length(PO972_0hr4hr_ProtId_Unfiltered)
print("PO972 0-4hr quantified both timepoints (unfiltered):")
print(PO972_0hr4hr_ProtId_Unfiltered_Num)
#
#PO972 all 0hr proteins identified
PO972_0hr_All_ProtId_Unfiltered <- union(PO972_0hr_Unique_ProtId_Unfiltered, PO972_0hr4hr_ProtId_Unfiltered)
PO972_0hr_All_ProtId_Unfiltered_Num <- length(PO972_0hr_All_ProtId_Unfiltered)
print("PO972 0hr all (unfiltered):")
print(PO972_0hr_All_ProtId_Unfiltered_Num)
#
#PO972 all 4hr proteins identified
PO972_4hr_All_ProtId_Unfiltered <- union(PO972_4hr_Unique_ProtId_Unfiltered, PO972_0hr4hr_ProtId_Unfiltered)
PO972_4hr_All_ProtId_Unfiltered_Num <- length(PO972_4hr_All_ProtId_Unfiltered)
print("PO972 4hr all (unfiltered):")
print(PO972_4hr_All_ProtId_Unfiltered_Num)
#
#PO972 all proteins identified
PO972ProtId_Unfiltered <- Reduce(union, list(PO972_0hr_Unique_ProtId_Unfiltered, PO972_4hr_Unique_ProtId_Unfiltered, PO972_0hr4hr_ProtId_Unfiltered))
PO972ProtId_Unfiltered_Num <- length(PO972ProtId_Unfiltered)
print("PO972 all quantified; union of 0hr unique, 4hr unique, 0hr4hr identified (unfiltered):")
print(PO972ProtId_Unfiltered_Num)
write(PO972ProtId_Unfiltered, file = file.path(".", analysisdir, "PO972_AllIdentifiedProteins_Unfiltered.txt"))
#
#PO1063 0-4hr stim quantified proteins. proteins identified at 0 hour only
PO1063_0hr_Unique_ProtId_Unfiltered <- row.names(subset(WCPstimQuant, (!is.na(WCPstimQuant$IntensityH.PO1063_HLM) & is.na(WCPstimQuant$IntensityL.PO1063_HLM) & is.na(WCPstimQuant$normLHratio.PO1063_HLM))))
PO1063_0hr_Unique_ProtId_Unfiltered_Num <- length(PO1063_0hr_Unique_ProtId_Unfiltered)
print("PO1063 0hr unique (unfiltered):")
print(PO1063_0hr_Unique_ProtId_Unfiltered_Num)
#
#PO1063 0-4hr stim quantified proteins. proteins identified at 4 hour only
PO1063_4hr_Unique_ProtId_Unfiltered <- row.names(subset(WCPstimQuant, (!is.na(WCPstimQuant$IntensityL.PO1063_HLM) & is.na(WCPstimQuant$IntensityH.PO1063_HLM) & is.na(WCPstimQuant$normLHratio.PO1063_HLM))))
PO1063_4hr_Unique_ProtId_Unfiltered_Num <- length(PO1063_4hr_Unique_ProtId_Unfiltered)
print("PO1063 4hr unique (unfiltered):")
print(PO1063_4hr_Unique_ProtId_Unfiltered_Num)
#
#PO1063 0-4hr stim quantified proteins. proteins identified in BOTH 0 AND 4 hours
PO1063_0hr4hr_ProtId_Unfiltered <- row.names(subset(WCPstimQuant, !is.na(WCPstimQuant$HLMStim4hr_fc)))
PO1063_0hr4hr_ProtId_Unfiltered_Num <- length(PO1063_0hr4hr_ProtId_Unfiltered)
print("PO1063 0-4hr quantified both timepoints (unfiltered):")
print(PO1063_0hr4hr_ProtId_Unfiltered_Num)
#
#PO1063 all 0hr proteins identified
PO1063_0hr_All_ProtId_Unfiltered <- union(PO1063_0hr_Unique_ProtId_Unfiltered, PO1063_0hr4hr_ProtId_Unfiltered)
PO1063_0hr_All_ProtId_Unfiltered_Num <- length(PO1063_0hr_All_ProtId_Unfiltered)
print("PO1063 0hr all (unfiltered):")
print(PO1063_0hr_All_ProtId_Unfiltered_Num)
#
#PO1063 all 4hr proteins identified
PO1063_4hr_All_ProtId_Unfiltered <- union(PO1063_4hr_Unique_ProtId_Unfiltered, PO1063_0hr4hr_ProtId_Unfiltered)
PO1063_4hr_All_ProtId_Unfiltered_Num <- length(PO1063_4hr_All_ProtId_Unfiltered)
print("PO1063 4hr all (unfiltered):")
print(PO1063_4hr_All_ProtId_Unfiltered_Num)
#
#PO1063 all proteins identified
PO1063ProtId_Unfiltered <- Reduce(union, list(PO1063_0hr_Unique_ProtId_Unfiltered, PO1063_4hr_Unique_ProtId_Unfiltered, PO1063_0hr4hr_ProtId_Unfiltered))
PO1063ProtId_Unfiltered_Num <- length(PO1063ProtId_Unfiltered)
print("PO1063 all quantified; union of 0hr unique, 4hr unique, 0hr4hr identified (unfiltered):")
print(PO1063ProtId_Unfiltered_Num)
write(PO1063ProtId_Unfiltered, file = file.path(".", analysisdir, "PO1063_AllIdentifiedProteins_Unfiltered.txt"))
#
#PO1075 0-4hr stim quantified proteins. proteins identified at 0 hour only
PO1075_0hr_Unique_ProtId_Unfiltered <- row.names(subset(WCPstimQuant, (!is.na(WCPstimQuant$IntensityH.PO1075_WCPM) & is.na(WCPstimQuant$IntensityL.PO1075_WCPM) & is.na(WCPstimQuant$normLHratio.PO1075_WCPM))))
PO1075_0hr_Unique_ProtId_Unfiltered_Num <- length(PO1075_0hr_Unique_ProtId_Unfiltered)
print("PO1075 0hr unique (unfiltered):")
print(PO1075_0hr_Unique_ProtId_Unfiltered_Num)
#
#PO1075 0-4hr stim quantified proteins. proteins identified at 4 hour only
PO1075_4hr_Unique_ProtId_Unfiltered <- row.names(subset(WCPstimQuant, (!is.na(WCPstimQuant$IntensityL.PO1075_WCPM) & is.na(WCPstimQuant$IntensityH.PO1075_WCPM) & is.na(WCPstimQuant$normLHratio.PO1075_WCPM))))
PO1075_4hr_Unique_ProtId_Unfiltered_Num <- length(PO1075_4hr_Unique_ProtId_Unfiltered)
print("PO1075 4hr unique (unique):")
print(PO1075_4hr_Unique_ProtId_Unfiltered_Num)
#
#PO1075 0-4hr stim quantified proteins. proteins identified in BOTH 0 AND 4 hours
PO1075_0hr4hr_ProtId_Unfiltered <- row.names(subset(WCPstimQuant, !is.na(WCPstimQuant$WCPMStim4hr_fc)))
PO1075_0hr4hr_ProtId_Unfiltered_Num <- length(PO1075_0hr4hr_ProtId_Unfiltered)
print("PO1075 0-4hr quantified both timepoints (unique):")
print(PO1075_0hr4hr_ProtId_Unfiltered_Num)
#
#PO1075 all 0hr proteins identified
PO1075_0hr_All_ProtId_Unfiltered <- union(PO1075_0hr_Unique_ProtId_Unfiltered, PO1075_0hr4hr_ProtId_Unfiltered)
PO1075_0hr_All_ProtId_Unfiltered_Num <- length(PO1075_0hr_All_ProtId_Unfiltered)
print("PO1075 0hr all (unfiltered):")
print(PO1075_0hr_All_ProtId_Unfiltered_Num)
#
#PO1075 all 4hr proteins identified
PO1075_4hr_All_ProtId_Unfiltered <- union(PO1075_4hr_Unique_ProtId_Unfiltered, PO1075_0hr4hr_ProtId_Unfiltered)
PO1075_4hr_All_ProtId_Unfiltered_Num <- length(PO1075_4hr_All_ProtId_Unfiltered)
print("PO1075 4hr all (unfiltered):")
print(PO1075_4hr_All_ProtId_Unfiltered_Num)
#
#PO1075 all proteins identified
PO1075ProtId_Unfiltered <- Reduce(union, list(PO1075_0hr_Unique_ProtId_Unfiltered, PO1075_4hr_Unique_ProtId_Unfiltered, PO1075_0hr4hr_ProtId_Unfiltered))
PO1075ProtId_Unfiltered_Num <- length(PO1075ProtId_Unfiltered)
print("PO1075 all quantified; union of 0hr unique, 4hr unique, 0hr4hr identified (unfiltered):")
print(PO1075ProtId_Unfiltered_Num)
write(PO1075ProtId_Unfiltered, file = file.path(".", analysisdir, "PO1075_AllIdentifiedProteins_Unfiltered.txt"))
#
#PO972, PO1063, PO1075 union of all identified proteins
PO972PO1063PO1075_AllId_Union_ProtId_Unfiltered <- Reduce(union, list(PO972ProtId_Unfiltered, PO1063ProtId_Unfiltered, PO1075ProtId_Unfiltered))
PO972PO1063PO1075_AllId_Union_ProtId_Unfiltered_Num <- length(PO972PO1063PO1075_AllId_Union_ProtId_Unfiltered)
print("PO972-PO1063-PO1075 all identified (unfiltered):")
print(PO972PO1063PO1075_AllId_Union_ProtId_Unfiltered_Num)
#
#PO972, PO1063, PO1075 union of all 0hr identified proteins
PO972PO1063PO1075_0hr_Union_ProtId_Unfiltered <- Reduce(union, list(PO972_0hr_All_ProtId_Unfiltered, PO1063_0hr_All_ProtId_Unfiltered, PO1075_0hr_All_ProtId_Unfiltered))
PO972PO1063PO1075_0hr_Union_ProtId_Unfiltered_Num <- length(PO972PO1063PO1075_0hr_Union_ProtId_Unfiltered)
print("PO972-PO1063-PO1075 0hr identified (unfiltered):")
print(PO972PO1063PO1075_0hr_Union_ProtId_Unfiltered_Num)
#
#PO972, PO1063, PO1075 union of all 4hr identified proteins
PO972PO1063PO1075_4hr_Union_ProtId_Unfiltered <- Reduce(union, list(PO972_4hr_All_ProtId_Unfiltered, PO1063_4hr_All_ProtId_Unfiltered, PO1075_4hr_All_ProtId_Unfiltered))
PO972PO1063PO1075_4hr_Union_ProtId_Unfiltered_Num <- length(PO972PO1063PO1075_4hr_Union_ProtId_Unfiltered)
print("PO972-PO1063-PO1075 4hr identified (unfiltered):")
print(PO972PO1063PO1075_4hr_Union_ProtId_Unfiltered_Num)
#
#PO972, PO1063, PO1075 unique 0hr proteins
PO972PO1063PO1075_0hr_Unique_ProtId_Unfiltered <- setdiff(PO972PO1063PO1075_0hr_Union_ProtId_Unfiltered, PO972PO1063PO1075_4hr_Union_ProtId_Unfiltered)
PO972PO1063PO1075_0hr_Unique_ProtId_Unfiltered_Num <- length(PO972PO1063PO1075_0hr_Unique_ProtId_Unfiltered)
print("PO972-PO1063-PO1075 0hr unique identified (unfiltered):")
print(PO972PO1063PO1075_0hr_Unique_ProtId_Unfiltered_Num)
#
#PO972, PO1063, PO1075 unique 4hr proteins
PO972PO1063PO1075_4hr_Unique_ProtId_Unfiltered <- setdiff(PO972PO1063PO1075_4hr_Union_ProtId_Unfiltered, PO972PO1063PO1075_0hr_Union_ProtId_Unfiltered)
PO972PO1063PO1075_4hr_Unique_ProtId_Unfiltered_Num <- length(PO972PO1063PO1075_4hr_Unique_ProtId_Unfiltered)
print("PO972-PO1063-PO1075 4hr unique identified (unfiltered):")
print(PO972PO1063PO1075_4hr_Unique_ProtId_Unfiltered_Num)
#
#PO972, PO1063, PO1075 0hr4hr proteins
PO972PO1063PO1075_0hr4hr_Int_ProtId_Unfiltered <- intersect(PO972PO1063PO1075_0hr_Union_ProtId_Unfiltered, PO972PO1063PO1075_4hr_Union_ProtId_Unfiltered)
PO972PO1063PO1075_0hr4hr_Int_ProtId_Unfiltered_Num <- length(PO972PO1063PO1075_0hr4hr_Int_ProtId_Unfiltered)
print("PO972-PO1063-PO1075 0hr4hr intersection identified (unfiltered):")
print(PO972PO1063PO1075_0hr4hr_Int_ProtId_Unfiltered_Num)
#
###end quantified set analysis for unfiltered data
WCPstimQuantDataForSupplement <- WCPstimQuant[row.names(WCPstimQuant)%in%PO972PO1063PO1075_AllId_Union_ProtId_Unfiltered,]
write.table(WCPstimQuantDataForSupplement, file="WcpRnaKeGG_SupplementalTable_1.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)


### filter by number of identified peptides or LHcounts
#a copy of the unfiltered dataset
WCPstimQuantUnfiltered <- WCPstimQuant
#dim(WCPstimQuantUnfiltered)
#

##filter to remove proteins that have only one peptide identified in only one experiment
#remove if not identified in at least two experiments in rest or restim
#filterprots_neg <- row.names(subset(WCPstimQuant, (((is.na(WCPstimQuant$IntensityL.PO1063_HLM) & is.na(WCPstimQuant$IntensityL.PO1075_WCPM)) | 
#                                                   (is.na(WCPstimQuant$IntensityL.PO1063_HLM) & WCPstimQuant$MSMScounts.PO972igd_4hr==0) |
#												   (is.na(WCPstimQuant$IntensityL.PO1075_WCPM) & WCPstimQuant$MSMScounts.PO972igd_4hr==0)) &
#												   ((is.na(WCPstimQuant$IntensityH.PO1063_HLM) & is.na(WCPstimQuant$IntensityH.PO1075_WCPM)) |
#												   (is.na(WCPstimQuant$IntensityH.PO1063_HLM) & WCPstimQuant$MSMScounts.PO972igd_0hr==0) |
#												   (is.na(WCPstimQuant$IntensityH.PO1075_WCPM) & WCPstimQuant$MSMScounts.PO972igd_0hr==0))
#												  )
#												 ))
filterprots_neg <- row.names(subset(WCPstimQuant, (((is.na(WCPstimQuant$IntensityL.PO1063_HLM) & is.na(WCPstimQuant$IntensityL.PO1075_WCPM)) | 
                                                   (is.na(WCPstimQuant$IntensityL.PO1063_HLM) & is.na(WCPstimQuant$PO972igd_4hr)) |
												   (is.na(WCPstimQuant$IntensityL.PO1075_WCPM) & is.na(WCPstimQuant$PO972igd_4hr))) &
												   ((is.na(WCPstimQuant$IntensityH.PO1063_HLM) & is.na(WCPstimQuant$IntensityH.PO1075_WCPM)) |
												   (is.na(WCPstimQuant$IntensityH.PO1063_HLM) & is.na(WCPstimQuant$PO972igd_0hr)) |
												   (is.na(WCPstimQuant$IntensityH.PO1075_WCPM) & is.na(WCPstimQuant$PO972igd_0hr)))
												  )
												 ))
print ("filtered to remove proteins")
WCPstimQuant <- WCPstimQuant[!row.names(WCPstimQuant)%in%filterprots_neg,]
dim(WCPstimQuant)
### end filter by number of identified peptides or LHcounts


###quantified set analysis for filtered data
#PO972 0-4hr stim quantified proteins. proteins identified at 0 hour only
##########PO972_0hr_Unique_ProtId <- row.names(subset(WCPstimQuant, ((!is.na(WCPstimQuant$PO972igd_0hr) & WCPstimQuant$MSMScounts.PO972igd_0hr>0) & is.na(WCPstimQuant$IGDStim4hr_fc))))
PO972_0hr_Unique_ProtId <- row.names(subset(WCPstimQuant, (!is.na(WCPstimQuant$PO972igd_0hr) & is.na(WCPstimQuant$IGDStim4hr_fc))))
PO972_0hr_Unique_ProtId_Num <- length(PO972_0hr_Unique_ProtId)
print("PO972 0hr unique:")
print(PO972_0hr_Unique_ProtId_Num)
#
#PO972 0-4hr stim quantified proteins. proteins identified at 4 hour only
##########PO972_4hr_Unique_ProtId <- row.names(subset(WCPstimQuant, ((!is.na(WCPstimQuant$PO972igd_4hr) & WCPstimQuant$MSMScounts.PO972igd_4hr>0) & is.na(WCPstimQuant$IGDStim4hr_fc))))
PO972_4hr_Unique_ProtId <- row.names(subset(WCPstimQuant, (!is.na(WCPstimQuant$PO972igd_4hr) & is.na(WCPstimQuant$IGDStim4hr_fc))))
PO972_4hr_Unique_ProtId_Num <- length(PO972_4hr_Unique_ProtId)
print("PO972 4hr unique:")
print(PO972_4hr_Unique_ProtId_Num)
#
#PO972 0-4hr stim quantified proteins. proteins identified in BOTH 0 AND 4 hours
PO972_0hr4hr_ProtId <- row.names(subset(WCPstimQuant, !is.na(WCPstimQuant$IGDStim4hr_fc)))
PO972_0hr4hr_ProtId_Num <- length(PO972_0hr4hr_ProtId)
print("PO972 0-4hr quantified both timepoints:")
print(PO972_0hr4hr_ProtId_Num)
#
#PO972 all 0hr proteins identified
PO972_0hr_All_ProtId <- union(PO972_0hr_Unique_ProtId, PO972_0hr4hr_ProtId)
PO972_0hr_All_ProtId_Num <- length(PO972_0hr_All_ProtId)
print("PO972 0hr all:")
print(PO972_0hr_All_ProtId_Num)
#
#PO972 all 4hr proteins identified
PO972_4hr_All_ProtId <- union(PO972_4hr_Unique_ProtId, PO972_0hr4hr_ProtId)
PO972_4hr_All_ProtId_Num <- length(PO972_4hr_All_ProtId)
print("PO972 4hr all:")
print(PO972_4hr_All_ProtId_Num)
#
#PO972 all proteins identified
PO972ProtId <- Reduce(union, list(PO972_0hr_Unique_ProtId, PO972_4hr_Unique_ProtId, PO972_0hr4hr_ProtId))
PO972ProtId_Num <- length(PO972ProtId)
print("PO972 all quantified; union of 0hr unique, 4hr unique, 0hr4hr identified:")
print(PO972ProtId_Num)
write(PO972ProtId, file = file.path(".", analysisdir, "PO972_AllIdentifiedProteins.txt"))
#
#PO1063 0-4hr stim quantified proteins. proteins identified at 0 hour only
PO1063_0hr_Unique_ProtId <- row.names(subset(WCPstimQuant, (!is.na(WCPstimQuant$IntensityH.PO1063_HLM) & is.na(WCPstimQuant$IntensityL.PO1063_HLM) & is.na(WCPstimQuant$normLHratio.PO1063_HLM))))
PO1063_0hr_Unique_ProtId_Num <- length(PO1063_0hr_Unique_ProtId)
print("PO1063 0hr unique:")
print(PO1063_0hr_Unique_ProtId_Num)
#
#PO1063 0-4hr stim quantified proteins. proteins identified at 4 hour only
PO1063_4hr_Unique_ProtId <- row.names(subset(WCPstimQuant, (!is.na(WCPstimQuant$IntensityL.PO1063_HLM) & is.na(WCPstimQuant$IntensityH.PO1063_HLM) & is.na(WCPstimQuant$normLHratio.PO1063_HLM))))
PO1063_4hr_Unique_ProtId_Num <- length(PO1063_4hr_Unique_ProtId)
print("PO1063 4hr unique:")
print(PO1063_4hr_Unique_ProtId_Num)
#
#PO1063 0-4hr stim quantified proteins. proteins identified in BOTH 0 AND 4 hours
PO1063_0hr4hr_ProtId <- row.names(subset(WCPstimQuant, (!is.na(WCPstimQuant$HLMStim4hr_fc) | (!is.na(WCPstimQuant$IntensityL.PO1063_HLM) & !is.na(WCPstimQuant$IntensityH.PO1063_HLM)))))
PO1063_0hr4hr_ProtId_Num <- length(PO1063_0hr4hr_ProtId)
print("PO1063 0-4hr quantified both timepoints:")
print(PO1063_0hr4hr_ProtId_Num)
#
#PO1063 all 0hr proteins identified
PO1063_0hr_All_ProtId <- union(PO1063_0hr_Unique_ProtId, PO1063_0hr4hr_ProtId)
PO1063_0hr_All_ProtId_Num <- length(PO1063_0hr_All_ProtId)
print("PO1063 0hr all:")
print(PO1063_0hr_All_ProtId_Num)
#
#PO1063 all 4hr proteins identified
PO1063_4hr_All_ProtId <- union(PO1063_4hr_Unique_ProtId, PO1063_0hr4hr_ProtId)
PO1063_4hr_All_ProtId_Num <- length(PO1063_4hr_All_ProtId)
print("PO1063 4hr all:")
print(PO1063_4hr_All_ProtId_Num)
#
#PO1063 all proteins identified
PO1063ProtId <- Reduce(union, list(PO1063_0hr_Unique_ProtId, PO1063_4hr_Unique_ProtId, PO1063_0hr4hr_ProtId))
PO1063ProtId_Num <- length(PO1063ProtId)
print("PO1063 all quantified; union of 0hr unique, 4hr unique, 0hr4hr identified:")
print(PO1063ProtId_Num)
write(PO1063ProtId, file = file.path(".", analysisdir, "PO1063_AllIdentifiedProteins.txt"))
#
#PO1075 0-4hr stim quantified proteins. proteins identified at 0 hour only
PO1075_0hr_Unique_ProtId <- row.names(subset(WCPstimQuant, (!is.na(WCPstimQuant$IntensityH.PO1075_WCPM) & is.na(WCPstimQuant$IntensityL.PO1075_WCPM) & is.na(WCPstimQuant$normLHratio.PO1075_WCPM))))
PO1075_0hr_Unique_ProtId_Num <- length(PO1075_0hr_Unique_ProtId)
print("PO1075 0hr unique:")
print(PO1075_0hr_Unique_ProtId_Num)
#
#PO1075 0-4hr stim quantified proteins. proteins identified at 4 hour only
PO1075_4hr_Unique_ProtId <- row.names(subset(WCPstimQuant, (!is.na(WCPstimQuant$IntensityL.PO1075_WCPM) & is.na(WCPstimQuant$IntensityH.PO1075_WCPM) & is.na(WCPstimQuant$normLHratio.PO1075_WCPM))))
PO1075_4hr_Unique_ProtId_Num <- length(PO1075_4hr_Unique_ProtId)
print("PO1075 4hr unique:")
print(PO1075_4hr_Unique_ProtId_Num)
#
#PO1075 0-4hr stim quantified proteins. proteins identified in BOTH 0 AND 4 hours
PO1075_0hr4hr_ProtId <- row.names(subset(WCPstimQuant, (!is.na(WCPstimQuant$WCPMStim4hr_fc) | (!is.na(WCPstimQuant$IntensityL.PO1075_WCPM) & !is.na(WCPstimQuant$IntensityH.PO1075_WCPM)))))
PO1075_0hr4hr_ProtId_Num <- length(PO1075_0hr4hr_ProtId)
print("PO1075 0-4hr quantified both timepoints:")
print(PO1075_0hr4hr_ProtId_Num)
#
#PO1075 all 0hr proteins identified
PO1075_0hr_All_ProtId <- union(PO1075_0hr_Unique_ProtId, PO1075_0hr4hr_ProtId)
PO1075_0hr_All_ProtId_Num <- length(PO1075_0hr_All_ProtId)
print("PO1075 0hr all:")
print(PO1075_0hr_All_ProtId_Num)
#
#PO1075 all 4hr proteins identified
PO1075_4hr_All_ProtId <- union(PO1075_4hr_Unique_ProtId, PO1075_0hr4hr_ProtId)
PO1075_4hr_All_ProtId_Num <- length(PO1075_4hr_All_ProtId)
print("PO1075 4hr all:")
print(PO1075_4hr_All_ProtId_Num)
#
#PO1075 all proteins identified
PO1075ProtId <- Reduce(union, list(PO1075_0hr_Unique_ProtId, PO1075_4hr_Unique_ProtId, PO1075_0hr4hr_ProtId))
PO1075ProtId_Num <- length(PO1075ProtId)
print("PO1075 all quantified; union of 0hr unique, 4hr unique, 0hr4hr identified:")
print(PO1075ProtId_Num)
write(PO1075ProtId, file = file.path(".", analysisdir, "PO1075_AllIdentifiedProteins.txt"))
#
###end quantified set analysis for filtered data


##quantified set intersection analysis
#PO972, PO1063, PO1075 union of all identified proteins
PO972PO1063PO1075_AllId_Union_ProtId <- Reduce(union, list(PO972ProtId, PO1063ProtId, PO1075ProtId))
PO972PO1063PO1075_AllId_Union_ProtId_Num <- length(PO972PO1063PO1075_AllId_Union_ProtId)
print("PO972-PO1063-PO1075 all identified  union:")
print(PO972PO1063PO1075_AllId_Union_ProtId_Num)
write(PO972PO1063PO1075_AllId_Union_ProtId, file = file.path(".", analysisdir, "PO972PO1063PO1075_AllIdentifiedProteins_Union.txt"))
#
#PO972, PO1063, PO1075 intersection of all identified proteins
PO972PO1063PO1075_AllId_Intersect_ProtId <- Reduce(intersect, list(PO972ProtId, PO1063ProtId, PO1075ProtId))
PO972PO1063PO1075_AllId_Intersect_ProtId_Num <- length(PO972PO1063PO1075_AllId_Intersect_ProtId)
print("PO972-PO1063-PO1075 all identified intersect:")
print(PO972PO1063PO1075_AllId_Intersect_ProtId_Num)
write(PO972PO1063PO1075_AllId_Intersect_ProtId, file = file.path(".", analysisdir, "PO972PO1063PO1075_AllIdentifiedProteins_Intersection.txt"))
#
#PO972, PO1063, PO1075 union of all 0hr identified proteins
PO972PO1063PO1075_0hr_Union_ProtId <- Reduce(union, list(PO972_0hr_All_ProtId, PO1063_0hr_All_ProtId, PO1075_0hr_All_ProtId))
PO972PO1063PO1075_0hr_Union_ProtId_Num <- length(PO972PO1063PO1075_0hr_Union_ProtId)
print("PO972-PO1063-PO1075 0hr identified:")
print(PO972PO1063PO1075_0hr_Union_ProtId_Num)
#
#PO972, PO1063, PO1075 union of all 4hr identified proteins
PO972PO1063PO1075_4hr_Union_ProtId <- Reduce(union, list(PO972_4hr_All_ProtId, PO1063_4hr_All_ProtId, PO1075_4hr_All_ProtId))
PO972PO1063PO1075_4hr_Union_ProtId_Num <- length(PO972PO1063PO1075_4hr_Union_ProtId)
print("PO972-PO1063-PO1075 4hr identified:")
print(PO972PO1063PO1075_4hr_Union_ProtId_Num)
#
#PO972, PO1063, PO1075 unique 0hr proteins
PO972PO1063PO1075_0hr_Unique_ProtId <- setdiff(PO972PO1063PO1075_0hr_Union_ProtId, PO972PO1063PO1075_4hr_Union_ProtId)
PO972PO1063PO1075_0hr_Unique_ProtId_Num <- length(PO972PO1063PO1075_0hr_Unique_ProtId)
print("PO972-PO1063-PO1075 0hr unique identified:")
print(PO972PO1063PO1075_0hr_Unique_ProtId_Num)
#
#PO972, PO1063, PO1075 unique 4hr proteins
PO972PO1063PO1075_4hr_Unique_ProtId <- setdiff(PO972PO1063PO1075_4hr_Union_ProtId, PO972PO1063PO1075_0hr_Union_ProtId)
PO972PO1063PO1075_4hr_Unique_ProtId_Num <- length(PO972PO1063PO1075_4hr_Unique_ProtId)
print("PO972-PO1063-PO1075 4hr unique identified:")
print(PO972PO1063PO1075_4hr_Unique_ProtId_Num)
#
#PO972, PO1063, PO1075 0hr4hr proteins
PO972PO1063PO1075_0hr4hr_Int_ProtId <- intersect(PO972PO1063PO1075_0hr_Union_ProtId, PO972PO1063PO1075_4hr_Union_ProtId)
PO972PO1063PO1075_0hr4hr_Int_ProtId_Num <- length(PO972PO1063PO1075_0hr4hr_Int_ProtId)
print("PO972-PO1063-PO1075 0hr4hr intersection identified:")
print(PO972PO1063PO1075_0hr4hr_Int_ProtId_Num)
#
#PO972 0-4hr stim, PO1063 0-4hr stim all protein intersect
PO972PO1063_ProtId_Int <- intersect(PO972ProtId,PO1063ProtId)
PO972PO1063_ProtId_Int_Num <- length(PO972PO1063_ProtId_Int)
print("PO972-PO1063 all identified intersect:")
print(PO972PO1063_ProtId_Int_Num)
#PO972 0-4hr stim, PO1075 0-4hr stim quantified protein intersect
PO972PO1075_ProtId_Int <- intersect(PO972ProtId,PO1075ProtId)
PO972PO1075_ProtId_Int_Num <- length(PO972PO1075_ProtId_Int)
print("PO972-PO1075 all identified intersect:")
print(PO972PO1075_ProtId_Int_Num)
#PO1063 0-4hr stim, PO1075 0-4hr stim quantified protein intersect
PO1063PO1075_ProtId_Int <- intersect(PO1063ProtId,PO1075ProtId)
PO1063PO1075_ProtId_Int_Num <- length(PO1063PO1075_ProtId_Int)
print("PO1063-PO1075 all identified intersect:")
print(PO1063PO1075_ProtId_Int_Num)
#
#PO972, PO1063, PO1075 proteins with 0-4hr fold change union
PO972PO1063PO1075_All0hr4hrFc_Union_ProtId <- Reduce(union, list(PO972_0hr4hr_ProtId, PO1063_0hr4hr_ProtId, PO1075_0hr4hr_ProtId))
PO972PO1063PO1075_All0hr4hrFc_Union_ProtId_Num <- length(PO972PO1063PO1075_All0hr4hrFc_Union_ProtId)
print("PO972-PO1063-PO1075 0hr4hr fold change protein union:")
print(PO972PO1063PO1075_All0hr4hrFc_Union_ProtId_Num)
#
#PO972, PO1063, PO1075 proteins with 0-4hr fold change intersection
PO972PO1063PO1075_All0hr4hrFc_Intersect_ProtId <- Reduce(intersect, list(PO972_0hr4hr_ProtId, PO1063_0hr4hr_ProtId, PO1075_0hr4hr_ProtId))
PO972PO1063PO1075_All0hr4hrFc_Intersect_ProtId_Num <- length(PO972PO1063PO1075_All0hr4hrFc_Intersect_ProtId)
print("PO972-PO1063-PO1075 0hr4hr fold change protein intersect:")
print(PO972PO1063PO1075_All0hr4hrFc_Intersect_ProtId_Num)
#
##PO972 0hr unique, PO1063 0hr unique protein intersect
#PO972PO1063_0hr_Unique_Int <- intersect(PO972_0hr_Unique_ProtId, PO1063_0hr_Unique_ProtId)
#PO972PO1063_0hr_Unique_Int_Num <- length(PO972PO1063_0hr_Unique_Int)
#print("PO972-PO1063 0hr Unique intersection:")
#print(PO972PO1063_0hr_Unique_Int_Num)
##PO972 0hr unique, PO1075 0hr unique protein intersect
#PO972PO1075_0hr_Unique_Int <- intersect(PO972_0hr_Unique_ProtId, PO1075_0hr_Unique_ProtId)
#PO972PO1075_0hr_Unique_Int_Num <- length(PO972PO1075_0hr_Unique_Int)
#print("PO972-PO1075 0hr Unique intersection:")
#print(PO972PO1075_0hr_Unique_Int_Num)
##PO1063 0hr unique, PO1075 0hr unique protein intersect
#PO1063PO1075_0hr_Unique_Int <- intersect(PO1063_0hr_Unique_ProtId, PO1075_0hr_Unique_ProtId)
#PO1063PO1075_0hr_Unique_Int_Num <- length(PO1063PO1075_0hr_Unique_Int)
#print("PO1063-PO1075 0hr Unique intersection:")
#print(PO1063PO1075_0hr_Unique_Int_Num)
##PO972-PO1063-PO1075 0hr unique, appears unique in at least two experiments
#PO972PO1063PO1075_0hr_Unique_TwoExpId <- Reduce(union, list(PO972PO1063_0hr_Unique_Int, PO972PO1075_0hr_Unique_Int, PO1063PO1075_0hr_Unique_Int))
#PO972PO1063PO1075_0hr_Unique_TwoExpId_Num <- length(PO972PO1063PO1075_0hr_Unique_TwoExpId)
#print("PO972-PO1063-PO1075 0hr Unique in at least 2 experiments:")
#print(PO972PO1063PO1075_0hr_Unique_TwoExpId_Num)
##
##PO972 4hr unique, PO1063 4hr unique protein intersect
#PO972PO1063_4hr_Unique_Int <- intersect(PO972_4hr_Unique_ProtId, PO1063_4hr_Unique_ProtId)
#PO972PO1063_4hr_Unique_Int_Num <- length(PO972PO1063_4hr_Unique_Int)
#print("PO972-PO1063 4hr Unique intersection:")
#print(PO972PO1063_4hr_Unique_Int_Num)
##PO972 4hr unique, PO1075 4hr unique protein intersect
#PO972PO1075_4hr_Unique_Int <- intersect(PO972_4hr_Unique_ProtId, PO1075_4hr_Unique_ProtId)
#PO972PO1075_4hr_Unique_Int_Num <- length(PO972PO1075_4hr_Unique_Int)
#print("PO972-PO1075 4hr Unique intersection:")
#print(PO972PO1075_4hr_Unique_Int_Num)
##PO1063 4hr unique, PO1075 4hr unique protein intersect
#PO1063PO1075_4hr_Unique_Int <- intersect(PO1063_4hr_Unique_ProtId, PO1075_4hr_Unique_ProtId)
#PO1063PO1075_4hr_Unique_Int_Num <- length(PO1063PO1075_4hr_Unique_Int)
#print("PO1063-PO1075 4hr Unique intersection:")
#print(PO1063PO1075_4hr_Unique_Int_Num)
##PO972-PO1063-PO1075 4hr unique, appears unique in at least two experiments
#PO972PO1063PO1075_4hr_Unique_TwoExpId <- Reduce(union, list(PO972PO1063_4hr_Unique_Int, PO972PO1075_4hr_Unique_Int, PO1063PO1075_4hr_Unique_Int))
#PO972PO1063PO1075_4hr_Unique_TwoExpId_Num <- length(PO972PO1063PO1075_4hr_Unique_TwoExpId)
#print("PO972-PO1063-PO1075 4hr Unique in at least 2 experiments:")
#print(PO972PO1063PO1075_4hr_Unique_TwoExpId_Num)
##
##PO972 0hr unique, PO1063 0hr unique, PO1075 0hr unique quantified protein union
#PO972PO1063PO1075_0hr_Unique_Union <- Reduce(union, list(PO972_0hr_Unique_ProtId, PO1063_0hr_Unique_ProtId, PO1075_0hr_Unique_ProtId))
#PO972PO1063PO1075_0hr_Unique_Union_Num <- length(PO972PO1063PO1075_0hr_Unique_Union)
#print("PO972-PO1063-PO1075 0hr Unique union:")
#print(PO972PO1063PO1075_0hr_Unique_Union_Num)
##PO972 0hr unique, PO1063 0hr unique, PO1075 0hr unique quantified protein intersect
#PO972PO1063PO1075_0hr_Unique_Int <- Reduce(intersect, list(PO972_0hr_Unique_ProtId, PO1063_0hr_Unique_ProtId, PO1075_0hr_Unique_ProtId))
#PO972PO1063PO1075_0hr_Unique_Int_Num <- length(PO972PO1063PO1075_0hr_Unique_Int)
#print("PO972-PO1063-PO1075 0hr Unique intersection:")
#print(PO972PO1063PO1075_0hr_Unique_Int_Num)
##
##PO972 4hr unique, PO1063 4hr unique, PO1075 4hr unique quantified protein union
#PO972PO1063PO1075_4hr_Unique_Union <- Reduce(union, list(PO972_4hr_Unique_ProtId, PO1063_4hr_Unique_ProtId, PO1075_4hr_Unique_ProtId))
#PO972PO1063PO1075_4hr_Unique_Union_Num <- length(PO972PO1063PO1075_4hr_Unique_Union)
#print("PO972-PO1063-PO1075 4hr Unique union:")
#print(PO972PO1063PO1075_4hr_Unique_Union_Num)
##PO972 4hr unique, PO1063 4hr unique, PO1075 4hr unique quantified protein intersect
#PO972PO1063PO1075_4hr_Unique_Int <- Reduce(intersect, list(PO972_4hr_Unique_ProtId, PO1063_4hr_Unique_ProtId, PO1075_4hr_Unique_ProtId))
#PO972PO1063PO1075_4hr_Unique_Int_Num <- length(PO972PO1063PO1075_4hr_Unique_Int)
#print("PO972-PO1063-PO1075 4hr Unique intersection:")
#print(PO972PO1063PO1075_4hr_Unique_Int_Num)
##
###end quantified set analysis


##dataframe of only proteins with fold changes for all three experiments
WCPstimQuantIntersect <- WCPstimQuant[row.names(WCPstimQuant)%in%PO972PO1063PO1075_All0hr4hrFc_Intersect_ProtId,]
dim(WCPstimQuantIntersect)
WCPstimQuantIntersect_0hrIntensities <- WCPstimQuantIntersect[,c("IntensityH.PO1063_HLM", "IntensityH.PO1075_WCPM", "PO972igd_0hr")]
IntCorr_0hr<-cor(WCPstimQuantIntersect_0hrIntensities)
#col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
col <- colorRampPalette(c("#FFFFFF", "#6E9DCC", "#3F85CC"))
pdf(file.path(".", analysisdir, "PO972PO1063PO1075Wcp_0hr_Log2IntensityCorrelations_heatmap.pdf"))
corrplot(IntCorr_0hr, method="color", type="lower", addCoef.col = "grey55", number.cex=2.5, cl.lim = c(0.70, 1), is.corr = FALSE)
dev.off()

WCPstimQuantIntersect_4hrIntensities <- WCPstimQuantIntersect[,c("IntensityL.PO1063_HLM", "IntensityL.PO1075_WCPM", "PO972igd_4hr")]
IntCorr_4hr<-cor(WCPstimQuantIntersect_4hrIntensities)
pdf(file.path(".", analysisdir, "PO972PO1063PO1075Wcp_4hr_Log2IntensityCorrelations_heatmap.pdf"))
corrplot(IntCorr_4hr, method="color", type="lower", addCoef.col = "grey55", number.cex=2.5, cl.lim = c(0.70, 1), is.corr = FALSE)
dev.off()


if ( FALSE )
{
##calculate PCA plot
#make PCA data frame
PCAdata<-PCAdataTmp[,c("IntensityH.PO1063_HLM", "IntensityH.PO1075_WCPM", "PO972igd_0hr", "IntensityL.PO1063_HLM", "IntensityL.PO1075_WCPM", "PO972igd_4hr")]
#PCAdata<-PCAdataTmp[,c("HLMintH_zscore", "WCPMintH_zscore", "IGD0hr_zscore", "HLMintL_zscore", "WCPMintL_zscore", "IGD4hr_zscore")]
dim(PCAdata)
PCAdataFiltered<-na.omit(PCAdata)
#head(PCAdata)
#print(PCAdataFiltered)
dim(PCAdataFiltered)
#
##to subset the data
#proteindata<-proteindata[!(proteindata$proteinID=="mock_1" | proteindata$proteinID=="mock_2" | proteindata$proteinID=="mock_3"),]
#print(proteindata)
#dim(proteindata)
#
##separate protein ids and quantification data
##########PCAdata_proteins <- PCAdataFiltered[,1]
PCAdata_proteins <- row.names(PCAdataFiltered)
PCAdata_data <- PCAdataFiltered[,1:dim(PCAdataFiltered)[2]]
#
##perform pca using prcomp
PCAdata_pca <- prcomp(PCAdata_data, center=TRUE, scale=TRUE)
print(PCAdata_pca)
summary(PCAdata_pca)
#
test<-PCAdata_pca[["rotation"]]
head(test)

#xmin <- floor(min(test[,"PC1"])) 
#xmax <- ceiling(max(test[,"PC1"]))
#ymin <- floor(min(test[,"PC2"]))
#ymax <- ceiling(max(test[,"PC2"]))
xmin <- -0.42
xmax <- -0.4
ymin <- -0.6
ymax <- 0.6
pdf("PO972PO1063PO1075_Intensity_PCAplot.pdf")
par(bty="l")
plot(test[,"PC1"],test[,"PC2"],
	pch=c(16,16,16,18,18,18),
	cex=1,
	#col="gray50",
	#col=c("goldenrod1", "dodgerblue", "royalblue4", "goldenrod1", "dodgerblue", "royalblue4"),
	col=c("goldenrod1", "dodgerblue", "royalblue4", "red", "cyan", "purple"),
	xlab="",ylab="",
	xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	xaxt='n',yaxt='n'
)
#label and format axes
axis(1, c(seq(xmin,xmax,by=0.0025)), labels=F, col="black",cex.axis=1, tck=-0.01)
axis(1, c(seq(xmin,xmax,by=0.005)), labels=T, col="black",cex.axis=1)
mtext("PC1",1,line=2.75,cex=1.25)
axis(2, c(seq(ymin,ymax,by=0.1)), labels=F, col="black",cex.axis=1, tck=-0.01)
axis(2, c(seq(ymin,ymax,by=0.2)), labels=T, col="black",cex.axis=1)
mtext("PC2",2,line=2.75,cex=1.25)
dev.off()
#
}
##end calculate PCA plot


### overlap of RNAseq and proteomics
print("rna seq overlap...")
#all proteins identified in 0 and 4 hours, PO972PO1063PO1075_AllId_Union_ProtId from set analysis
#number of quantified RNA transcript genes in transcriptome (RNAseq dataset)
RnaIdentifiedGenes_InTranscriptome<-row.names(subset(RnaData, !is.na(RnaData$log2FoldChange.x)))
RnaIdentifiedGenes_InTranscriptome_Num<-length(RnaIdentifiedGenes_InTranscriptome)
print("rna IDs in transcriptome:")
print(RnaIdentifiedGenes_InTranscriptome_Num)
#write(RnaIdentifiedGenes_InTranscriptome, file = file.path(".", analysisdir, "RNAseq_IdentifiedProteinsInTranscriptome.txt"))
#number of quantified RNA transcript genes in proteome (proteome and RNAseq left join)
RnaIdentifiedGenes_InProteome<-row.names(subset(WCPstimQuant, !is.na(WCPstimQuant$log2FoldChangeStimUnstim)))
RnaIdentifiedGenes_InProteome_Num<-length(RnaIdentifiedGenes_InProteome)
print("rna IDs in proteome subet:")
print(RnaIdentifiedGenes_InProteome_Num)
#intersection of identified proteins and RNA genes
RnaIdentifiedGenes_InTranscriptome_EnsemblId <- RnaData[RnaIdentifiedGenes_InTranscriptome,"EnsemblId"]
PO972PO1063PO1075_AllId_Union_EnsemblId <- WCPstimQuant[PO972PO1063PO1075_AllId_Union_ProtId,"EnsemblId"]
#ProtRnaIdGenes_ProteomeProteome_Int<-intersect(PO972PO1063PO1075_AllId_Union_ProtId,RnaIdentifiedGenes_InTranscriptome)
ProtRnaIdGenes_ProteomeTranscriptome_Int<-intersect(PO972PO1063PO1075_AllId_Union_EnsemblId, RnaIdentifiedGenes_InTranscriptome_EnsemblId)
ProtRnaIdGenes_ProteomeTranscriptome_Int_Num<-length(ProtRnaIdGenes_ProteomeTranscriptome_Int)
print("protein AND rna IDs: proteome-proteome intersection")
print(ProtRnaIdGenes_ProteomeTranscriptome_Int_Num)
### end overlap of RNAseq and proteomics


### overlap between KeGG id and WCP id
KeGGProtWcpProtInt<-intersect(PO972PO1063PO1075_AllId_Union_ProtId,KeGG_AllUbProtId_Unfiltered)
KeGGProtWcpProtInt_Num<-length(KeGGProtWcpProtInt)
print("KeGGProtWcpProtInt_Num")
print(KeGGProtWcpProtInt_Num)
### end overlap between KeGG id and WCP id


###determine differentially regulated protein sets
print("differential expression...")
##determine up and down regulation based on average fold change > or < 0 and significant p-value
#proteins up regulated and significant based on PO972-PO1063-PO1075 average fold change and p-value
PO972PO1063PO1075_Wcp_Upreg_AvgFcPval <- row.names(subset(WCPstimQuant, (WCPstimQuant$Stim4hr_AvgFc > 0 & WCPstimQuant$Stim4hr_TtestPval<=ProtSigThresh)))
PO972PO1063PO1075_Wcp_Upreg_AvgFcPval_Num <- length(PO972PO1063PO1075_Wcp_Upreg_AvgFcPval)
PO972PO1063PO1075_Wcp_Upreg_AvgFcPval_Mean <- mean(WCPstimQuant[PO972PO1063PO1075_Wcp_Upreg_AvgFcPval,"Stim4hr_AvgFc"])
PO972PO1063PO1075_Wcp_Upreg_AvgFcPval_Sd <- sd(WCPstimQuant[PO972PO1063PO1075_Wcp_Upreg_AvgFcPval,"Stim4hr_AvgFc"])
print("PO972-PO1063-PO1075 significant upregulation by fc>0 and p-value<0.1 (num, mean, sd)")
print(PO972PO1063PO1075_Wcp_Upreg_AvgFcPval_Num)
print(PO972PO1063PO1075_Wcp_Upreg_AvgFcPval_Mean)
print(PO972PO1063PO1075_Wcp_Upreg_AvgFcPval_Sd)
write(PO972PO1063PO1075_Wcp_Upreg_AvgFcPval, file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Upreg_AvgFcPval_ProtId.txt"))
write(WCPstimQuant[PO972PO1063PO1075_Wcp_Upreg_AvgFcPval,"GeneId.x"], file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Upreg_AvgFcPval_GeneId.txt"))
#proteins down regulated and significant based on PO972-PO1063-PO1075 average fold change and p-value
PO972PO1063PO1075_Wcp_Downreg_AvgFcPval <- row.names(subset(WCPstimQuant, (WCPstimQuant$Stim4hr_AvgFc < 0 & WCPstimQuant$Stim4hr_TtestPval<=ProtSigThresh)))
PO972PO1063PO1075_Wcp_Downreg_AvgFcPval_Num <- length(PO972PO1063PO1075_Wcp_Downreg_AvgFcPval)
PO972PO1063PO1075_Wcp_Downreg_AvgFcPval_Mean <- mean(WCPstimQuant[PO972PO1063PO1075_Wcp_Downreg_AvgFcPval,"Stim4hr_AvgFc"])
PO972PO1063PO1075_Wcp_Downreg_AvgFcPval_Sd <- sd(WCPstimQuant[PO972PO1063PO1075_Wcp_Downreg_AvgFcPval,"Stim4hr_AvgFc"])
print("PO972-PO1063-PO1075 significant downregulation by fc<0 and p-value<0.1 (num, mean, sd)")
print(PO972PO1063PO1075_Wcp_Downreg_AvgFcPval_Num)
print(PO972PO1063PO1075_Wcp_Downreg_AvgFcPval_Mean)
print(PO972PO1063PO1075_Wcp_Downreg_AvgFcPval_Sd)
write(PO972PO1063PO1075_Wcp_Downreg_AvgFcPval, file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Downreg_AvgFcPval_ProtId.txt"))
write(WCPstimQuant[PO972PO1063PO1075_Wcp_Downreg_AvgFcPval,"GeneId.x"], file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Downreg_AvgFcPval_GeneId.txt"))
#proteins unchanged based on PO972-PO1063-PO1075 average fold change and p-value
PO972PO1063PO1075_Wcp_Unch_AvgFcPval <- row.names(subset(WCPstimQuant, (WCPstimQuant$Stim4hr_TtestPval>ProtSigThresh)))
PO972PO1063PO1075_Wcp_Unch_AvgFcPval_Num <- length(PO972PO1063PO1075_Wcp_Unch_AvgFcPval)
PO972PO1063PO1075_Wcp_Unch_AvgFcPval_Mean <- mean(WCPstimQuant[PO972PO1063PO1075_Wcp_Unch_AvgFcPval,"Stim4hr_AvgFc"])
PO972PO1063PO1075_Wcp_Unch_AvgFcPval_Sd <- sd(WCPstimQuant[PO972PO1063PO1075_Wcp_Unch_AvgFcPval,"Stim4hr_AvgFc"])
print("PO972-PO1063-PO1075 unch by p-value>0.1 (num, mean, sd)")
print(PO972PO1063PO1075_Wcp_Unch_AvgFcPval_Num)
print(PO972PO1063PO1075_Wcp_Unch_AvgFcPval_Mean)
print(PO972PO1063PO1075_Wcp_Unch_AvgFcPval_Sd)
write(PO972PO1063PO1075_Wcp_Unch_AvgFcPval, file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Unch_AvgFcPval_ProtId.txt"))
write(WCPstimQuant[PO972PO1063PO1075_Wcp_Unch_AvgFcPval,"GeneId.x"], file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Unch_AvgFcPval_GeneId.txt"))
#




##determine up and down regulation by > or < mean plus/minus standard deviation of PO972-PO1063-PO1075 average fold change distributions
#mean and standard deviation of average fold change
PO972PO1063PO1075_AvgFcMean<-mean(WCPstimQuant$Stim4hr_AvgFc, na.rm=TRUE)
PO972PO1063PO1075_AvgFcSd<-sd(WCPstimQuant$Stim4hr_AvgFc, na.rm=TRUE)
print("PO972-PO1063-PO1075 average fold change mean,sd")
print(PO972PO1063PO1075_AvgFcMean)
print(PO972PO1063PO1075_AvgFcSd)
#upregulation by PO972-PO1063-PO1075 avg fold change mean+sd
PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd<-row.names(subset(WCPstimQuant, WCPstimQuant$Stim4hr_AvgFc>=(PO972PO1063PO1075_AvgFcMean+PO972PO1063PO1075_AvgFcSd)))
PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd_Num<-length(PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd)
print("PO972-PO1063-PO1075 upregulation avg fold change > mean+sd")
print(PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd_Num)
write(PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd, file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd_ProtId.txt"))
write(WCPstimQuant[PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd,"GeneId.x"], file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd_GeneId.txt"))
#downregulation by PO972-PO1063-PO1075 avg fold change mean-sd
PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd<-row.names(subset(WCPstimQuant, WCPstimQuant$Stim4hr_AvgFc<=(PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd)))
PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd_Num<-length(PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd)
print("PO972-PO1063-PO1075 downregulation avg fold change < mean-sd")
print(PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd_Num)
write(PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd, file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd_ProtId.txt"))
write(WCPstimQuant[PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd,"GeneId.x"], file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd_GeneId.txt"))
#unchanged by PO972-PO1063-PO1075 avg fold change mean +/- sd
PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd<-row.names(subset(WCPstimQuant, (WCPstimQuant$Stim4hr_AvgFc<(PO972PO1063PO1075_AvgFcMean+PO972PO1063PO1075_AvgFcSd) & WCPstimQuant$Stim4hr_AvgFc>(PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd))))
PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd_Num<-length(PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd)
print("PO972-PO1063-PO1075 unchanged avg fold change mean +/- sd")
print(PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd_Num)
write(PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd, file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd_ProtId.txt"))
write(WCPstimQuant[PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd,"GeneId.x"], file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd_GeneId.txt"))
#

filterprots_1fc <- row.names(subset(WCPstimQuant, (((is.na(WCPstimQuant$IGDStim4hr_fc) & is.na(WCPstimQuant$HLMStim4hr_fc)) |
                                                    (is.na(WCPstimQuant$IGDStim4hr_fc) & is.na(WCPstimQuant$WCPMStim4hr_fc)) |
													(is.na(WCPstimQuant$HLMStim4hr_fc) & is.na(WCPstimQuant$WCPMStim4hr_fc))))))
print ("filtered to remove 1fc proteins")
WCPstimQuant2fc <- WCPstimQuant[!row.names(WCPstimQuant)%in%filterprots_1fc,]
dim(WCPstimQuant2fc)
### end filter by number of identified peptides or LHcounts
##determine up and down regulation by > or < mean plus/minus standard deviation of PO972-PO1063-PO1075 average fold change distributions
#mean and standard deviation of average fold change
PO972PO1063PO1075_AvgFcMean<-mean(WCPstimQuant2fc$Stim4hr_AvgFc, na.rm=TRUE)
PO972PO1063PO1075_AvgFcSd<-sd(WCPstimQuant2fc$Stim4hr_AvgFc, na.rm=TRUE)
print("PO972-PO1063-PO1075 average fold change mean,sd")
print(PO972PO1063PO1075_AvgFcMean)
print(PO972PO1063PO1075_AvgFcSd)
#upregulation by PO972-PO1063-PO1075 avg fold change mean+sd
PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd<-row.names(subset(WCPstimQuant2fc, WCPstimQuant2fc$Stim4hr_AvgFc>=(PO972PO1063PO1075_AvgFcMean+PO972PO1063PO1075_AvgFcSd)))
PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd_Num<-length(PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd)
print("PO972-PO1063-PO1075 upregulation avg fold change > mean+sd")
print(PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd_Num)
write(PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd, file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd_ProtId.txt"))
write(WCPstimQuant[PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd,"GeneId.x"], file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd_GeneId.txt"))
#downregulation by PO972-PO1063-PO1075 avg fold change mean-sd
PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd<-row.names(subset(WCPstimQuant2fc, WCPstimQuant2fc$Stim4hr_AvgFc<=(PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd)))
PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd_Num<-length(PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd)
print("PO972-PO1063-PO1075 downregulation avg fold change < mean-sd")
print(PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd_Num)
write(PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd, file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd_ProtId.txt"))
write(WCPstimQuant[PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd,"GeneId.x"], file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd_GeneId.txt"))
#unchanged by PO972-PO1063-PO1075 avg fold change mean +/- sd
PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd<-row.names(subset(WCPstimQuant2fc, (WCPstimQuant2fc$Stim4hr_AvgFc<(PO972PO1063PO1075_AvgFcMean+PO972PO1063PO1075_AvgFcSd) & WCPstimQuant2fc$Stim4hr_AvgFc>(PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd))))
PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd_Num<-length(PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd)
print("PO972-PO1063-PO1075 unchanged avg fold change mean +/- sd")
print(PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd_Num)
write(PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd, file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd_ProtId.txt"))
write(WCPstimQuant[PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd,"GeneId.x"], file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd_GeneId.txt"))
#


#intersection of upreg by pval and +/- standard deviation
PO972PO1063PO1075_Upreg_PvalStdev_Int <- intersect(PO972PO1063PO1075_Wcp_Upreg_AvgFcPval,PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd)
PO972PO1063PO1075_Upreg_PvalStdev_Int_Num <- length(PO972PO1063PO1075_Upreg_PvalStdev_Int)
print("PO972-PO1063-PO1075 upreg pval,stdev intersect")
print(PO972PO1063PO1075_Upreg_PvalStdev_Int_Num)
write(PO972PO1063PO1075_Upreg_PvalStdev_Int, file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Upreg_AvgFcPvalStdevInt_ProtId.txt"))
write(WCPstimQuant[PO972PO1063PO1075_Upreg_PvalStdev_Int,"GeneId.x"], file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Upreg_AvgFcPvalStdevInt_GeneId.txt"))
#intersection of downreg by pval and +/- standard deviation
PO972PO1063PO1075_Downreg_PvalStdev_Int <- intersect(PO972PO1063PO1075_Wcp_Downreg_AvgFcPval,PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd)
PO972PO1063PO1075_Downreg_PvalStdev_Int_Num <- length(PO972PO1063PO1075_Downreg_PvalStdev_Int)
print("PO972-PO1063-PO1075 downreg pval,stdev intersect")
print(PO972PO1063PO1075_Downreg_PvalStdev_Int_Num)
write(PO972PO1063PO1075_Downreg_PvalStdev_Int, file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Downreg_AvgFcPvalStdevInt_ProtId.txt"))
write(WCPstimQuant[PO972PO1063PO1075_Downreg_PvalStdev_Int,"GeneId.x"], file = file.path(".", analysisdir, "PO972PO1063PO1075_Wcp_Downreg_AvgFcPvalStdevInt_GeneId.txt"))
#
##determine up and down regulation by mean and standard deviation of individual fold change distributions
#PO1075 fold change distribution
PO1075FcMean<-mean(WCPstimQuant$WCPMStim4hr_fc, na.rm=TRUE)
PO1075FcSd<-sd(WCPstimQuant$WCPMStim4hr_fc, na.rm=TRUE)
print("PO1075 fold change mean,sd")
print(PO1075FcMean)
print(PO1075FcSd)
#PO1075 up regulation by mean + 1 standard deviation
PO1075_Upreg_MeanSd<-row.names(subset(WCPstimQuant, (WCPstimQuant$WCPMStim4hr_fc>=(PO1075FcMean+PO1075FcSd))))
PO1075_Upreg_MeanSd_Num<-length(PO1075_Upreg_MeanSd)
print("PO1075 upreg by fc mean+sd")
print(PO1075_Upreg_MeanSd_Num)
#PO1075 down regulation by mean - 1 standard deviation
PO1075_Downreg_MeanSd<-row.names(subset(WCPstimQuant, (WCPstimQuant$WCPMStim4hr_fc<=(PO1075FcMean-PO1075FcSd))))
PO1075_Downreg_MeanSd_Num<-length(PO1075_Downreg_MeanSd)
print("PO1075 downreg by fc mean+sd")
print(PO1075_Downreg_MeanSd_Num)
#PO1063 minus Nedd inhibitor fold change distribution
PO1063FcMean<-mean(WCPstimQuant$HLMStim4hr_fc, na.rm=TRUE)
PO1063FcSd<-sd(WCPstimQuant$HLMStim4hr_fc, na.rm=TRUE)
print("PO1063 fold change mean,sd")
print(PO1063FcMean)
print(PO1063FcSd)
#PO1063 minus Nedd inhibitor up regulation by mean + 1 standard deviation
PO1063_Upreg_MeanSd<-row.names(subset(WCPstimQuant, (WCPstimQuant$HLMStim4hr_fc>=(PO1063FcMean+PO1063FcSd))))
PO1063_Upreg_MeanSd_Num<-length(PO1063_Upreg_MeanSd)
print("PO1063 upreg by fc mean+sd")
print(PO1063_Upreg_MeanSd_Num)
#PO1063 minus Nedd inhibitor down regulation by mean - 1 standard deviation
PO1063_Downreg_MeanSd<-row.names(subset(WCPstimQuant, (WCPstimQuant$HLMStim4hr_fc<=(PO1063FcMean-PO1063FcSd))))
PO1063_Downreg_MeanSd_Num<-length(PO1063_Downreg_MeanSd)
print("PO1063 downreg by fc mean+sd")
print(PO1063_Downreg_MeanSd_Num)
#PO972 fold change distribution
PO972FcMean<-mean(WCPstimQuant$IGDStim4hr_fc, na.rm=TRUE)
PO972FcSd<-sd(WCPstimQuant$IGDStim4hr_fc, na.rm=TRUE)
print("PO972 fold change mean,sd")
print(PO972FcMean)
print(PO972FcSd)
#PO972 up regulation by mean + 1 standard deviation
PO972_Upreg_MeanSd<-row.names(subset(WCPstimQuant, (WCPstimQuant$IGDStim4hr_fc>=(PO972FcMean+PO972FcSd))))
PO972_Upreg_MeanSd_Num<-length(PO972_Upreg_MeanSd)
print("PO972 upreg by fc mean+sd")
print(PO972_Upreg_MeanSd_Num)
#PO972 down regulation by mean - 1 standard deviation
PO972_Downreg_MeanSd<-row.names(subset(WCPstimQuant, (WCPstimQuant$IGDStim4hr_fc<=(PO972FcMean-PO972FcSd))))
PO972_Downreg_MeanSd_Num<-length(PO972_Downreg_MeanSd)
print("PO972 downreg by fc mean+sd")
print(PO972_Downreg_MeanSd_Num)
#
#PO1075-PO1063-PO972 individual mean-sd upreg union
PO1075PO1063PO972_Union_MeanSt_Upreg <- Reduce(union, list(PO1075_Upreg_MeanSd, PO1063_Upreg_MeanSd, PO972_Upreg_MeanSd))
PO1075PO1063PO972_Union_MeanSt_Upreg_Num <- length(PO1075PO1063PO972_Union_MeanSt_Upreg)
print("PO972-PO1063-PO1075 union individual mean-sd up regulation")
print(PO1075PO1063PO972_Union_MeanSt_Upreg_Num)
##proteins up regulated and significant based on PO972, PO1075 and PO1063 mean + 1sd fold change intersection
#PO1075+PO1063 upreg
PO1075PO1063_MeanSd_Upreg<-intersect(PO1075_Upreg_MeanSd, PO1063_Upreg_MeanSd)
PO1075PO1063_MeanSd_Upreg_Num<-length(PO1075PO1063_MeanSd_Upreg)
#PO1075+PO972 upreg
PO1075PO972_MeanSd_Upreg<-intersect(PO1075_Upreg_MeanSd, PO972_Upreg_MeanSd)
PO1075PO972_MeanSd_Upreg_Num<-length(PO1075PO972_MeanSd_Upreg)
#PO0163+PO972 upreg
PO1063PO972_MeanSd_Upreg<-intersect(PO1063_Upreg_MeanSd, PO972_Upreg_MeanSd)
PO1063PO972_MeanSd_Upreg_Num<-length(PO1063PO972_MeanSd_Upreg)
#PO1075+PO1063 + PO1075+PO972 + PO1063+PO972 (union of intersection of two)
PO1075PO1063PO972_TwoInt_MeanSd_Upreg<-Reduce(union, list(PO1075PO1063_MeanSd_Upreg, PO1075PO972_MeanSd_Upreg, PO1063PO972_MeanSd_Upreg))
PO1075PO1063PO972_TwoInt_MeanSd_Upreg_Num<-length(PO1075PO1063PO972_TwoInt_MeanSd_Upreg)
print("PO972-PO1063-PO1075 at least two individual mean-sd up regulation")
print(PO1075PO1063PO972_TwoInt_MeanSd_Upreg_Num)
#PO1075+PO1063+PO972 upreg (intersection of three)
PO1075PO1063PO972_MeanSd_Upreg<-Reduce(intersect,  list(PO1075_Upreg_MeanSd, PO1063_Upreg_MeanSd, PO972_Upreg_MeanSd)) 
PO1075PO1063PO972_MeanSd_Upreg_Num<-length(PO1075PO1063PO972_MeanSd_Upreg)
print("PO972-PO1063-PO1075 intersection of individual mean-sd up regulation")
print(PO1075PO1063PO972_MeanSd_Upreg_Num)
#
#PO1075-PO1063-PO972 individual mean-sd upreg union
PO1075PO1063PO972_Union_MeanSt_Downreg <- Reduce(union, list(PO1075_Downreg_MeanSd, PO1063_Downreg_MeanSd, PO972_Downreg_MeanSd))
PO1075PO1063PO972_Union_MeanSt_Downreg_Num <- length(PO1075PO1063PO972_Union_MeanSt_Downreg)
print("PO972-PO1063-PO1075 union individual mean-sd down regulation")
print(PO1075PO1063PO972_Union_MeanSt_Downreg_Num)
##proteins down regulated and significant based on PO972, PO1075 and PO1063 mean - 1sd fold change intersection
#PO1075+PO1063 downreg
PO1075PO1063_MeanSd_Downreg<-intersect(PO1075_Downreg_MeanSd, PO1063_Downreg_MeanSd)
PO1075PO1063_MeanSd_Downreg_Num<-length(PO1075PO1063_MeanSd_Downreg)
#PO1075+PO972 downreg
PO1075PO972_MeanSd_Downreg<-intersect(PO1075_Downreg_MeanSd, PO972_Downreg_MeanSd)
PO1075PO972_MeanSd_Downreg_Num<-length(PO1075PO972_MeanSd_Downreg)
#PO1063+PO972 downreg
PO1063PO972_MeanSd_Downreg<-intersect(PO1063_Downreg_MeanSd, PO972_Downreg_MeanSd)
PO1063PO972_MeanSd_Downreg_Num<-length(PO1063PO972_MeanSd_Downreg)
#PO1075+PO1063 + PO1063+PO972 + PO1063+PO972 (union of intersection of two)
PO1075PO1063PO972_TwoInt_MeanSd_Downreg<-Reduce(union, list(PO1075PO1063_MeanSd_Downreg, PO1075PO972_MeanSd_Downreg, PO1063PO972_MeanSd_Downreg))
PO1075PO1063PO972_TwoInt_MeanSd_Downreg_Num<-length(PO1075PO1063PO972_TwoInt_MeanSd_Downreg)
print("PO972-PO1063-PO1075 at least two individual mean-sd down regulation")
print(PO1075PO1063PO972_TwoInt_MeanSd_Downreg_Num)
#PO1075+PO1063+PO972 downreg (intersection of three)
PO1075PO1063PO972_MeanSd_Downreg<-Reduce(intersect, list(PO1075_Downreg_MeanSd, PO1063_Downreg_MeanSd, PO972_Downreg_MeanSd))
PO1075PO1063PO972_MeanSd_Downreg_Num<-length(PO1075PO1063PO972_MeanSd_Downreg)
print("PO972-PO1063-PO1075 intersection of individual mean-sd down regulation")
print(PO1075PO1063PO972_MeanSd_Downreg_Num)
#
#PO972 1hr fold change distribution
PO972FcMean_1hr<-mean(WCPstimQuant$IGDStim1hr_fc, na.rm=TRUE)
PO972FcSd_1hr<-sd(WCPstimQuant$IGDStim1hr_fc, na.rm=TRUE)
print("PO972 1hr fold change mean,sd")
print(PO972FcMean_1hr)
print(PO972FcSd_1hr)
#PO972 1hr up regulation by mean + 1 standard deviation
PO972_Upreg_1hr_MeanSd<-row.names(subset(WCPstimQuant, (WCPstimQuant$IGDStim1hr_fc>=(PO972FcMean_1hr+PO972FcSd_1hr))))
PO972_Upreg_1hr_MeanSd_Num<-length(PO972_Upreg_1hr_MeanSd)
print("PO972 1hr upreg by fc mean+sd")
print(PO972_Upreg_1hr_MeanSd_Num)
#PO972 1hr down regulation by mean - 1 standard deviation
PO972_Downreg_1hr_MeanSd<-row.names(subset(WCPstimQuant, (WCPstimQuant$IGDStim1hr_fc<=(PO972FcMean_1hr-PO972FcSd_1hr))))
PO972_Downreg_1hr_MeanSd_Num<-length(PO972_Downreg_1hr_MeanSd)
print("PO972 1hr downreg by fc mean+sd")
print(PO972_Downreg_1hr_MeanSd_Num)
#
#proteins that are increased, moderately increased, decreased or moderatly decreased
wcpdeltathreshold <- log2(1.25)
PO972PO1063PO1075_IncDeProt_Num<-length(row.names(subset(WCPstimQuant, (WCPstimQuant$Stim4hr_AvgFc>=0.321928 & WCPstimQuant$Stim4hr_TtestPval<=ProtSigThresh))))
print("PO972-PO1063-PO1075 diff ex increased prots")
print(PO972PO1063PO1075_IncDeProt_Num)
PO972PO1063PO1075_ModIncDeProt_Num<-length(row.names(subset(WCPstimQuant, ((WCPstimQuant$Stim4hr_AvgFc<0.321928 & WCPstimQuant$Stim4hr_AvgFc>0) & WCPstimQuant$Stim4hr_TtestPval<=ProtSigThresh))))
print("PO972-PO1063-PO1075 diff ex moderately increased prots")
print(PO972PO1063PO1075_ModIncDeProt_Num)
PO972PO1063PO1075_DecDeProt_Num<-length(row.names(subset(WCPstimQuant, (WCPstimQuant$Stim4hr_AvgFc<= -0.321928 & WCPstimQuant$Stim4hr_TtestPval<=ProtSigThresh))))
print("PO972-PO1063-PO1075 diff ex decreased prots")
print(PO972PO1063PO1075_DecDeProt_Num)
PO972PO1063PO1075_ModDecDeProt_Num<-length(row.names(subset(WCPstimQuant, ((WCPstimQuant$Stim4hr_AvgFc>-0.321928 & WCPstimQuant$Stim4hr_AvgFc<0) & WCPstimQuant$Stim4hr_TtestPval<=ProtSigThresh))))
print("PO972-PO1063-PO1075 diff ex moderately decreased prots")
print(PO972PO1063PO1075_ModDecDeProt_Num)
###end determine differentially regulated protein sets


###determine RNAseq up and down regulation
print("RNAseq diff expression...")
#
#upregulation based on avg fold change and BH corrected p value
RNAseq_Upreg_AvgFcPval <- row.names(subset(WCPstimQuant, (WCPstimQuant$log2FoldChangeStimUnstim>0 & WCPstimQuant$padj<=RnaSigThresh)))
RNAseq_Upreg_AvgFcPval_Num <- length(RNAseq_Upreg_AvgFcPval)
print("RNAseq significant upregulation by fc>0 and p-value<0.01")
print(RNAseq_Upreg_AvgFcPval_Num)
#
#downregulation based on avg fold change and BH corrected p value
RNAseq_Downreg_AvgFcPval <- row.names(subset(WCPstimQuant, (WCPstimQuant$log2FoldChangeStimUnstim<0 & WCPstimQuant$padj<=RnaSigThresh)))
RNAseq_Downreg_AvgFcPval_Num <- length(RNAseq_Downreg_AvgFcPval)
print("RNAseq significant downregulation by fc<0 and p-value<0.01")
print(RNAseq_Downreg_AvgFcPval_Num)
#
#unchanged based on avg fold change and BH corrected p value
RNAseq_Unch_AvgFcPval <- row.names(subset(WCPstimQuant, (WCPstimQuant$padj>RnaSigThresh)))
RNAseq_Unch_AvgFcPval_Num <- length(RNAseq_Unch_AvgFcPval)
print("RNAseq unchanged by fc><0 and p-value<0.01")
print(RNAseq_Unch_AvgFcPval_Num)
#
#
#mean and standard deviation of average fold change
RNAseq_FcMean<-mean(WCPstimQuant$log2FoldChangeStimUnstim, na.rm=TRUE)
RNAseq_FcSd<-sd(WCPstimQuant$log2FoldChangeStimUnstim, na.rm=TRUE)
print("RNAseq average fold change mean,sd")
print(RNAseq_FcMean)
print(RNAseq_FcSd)
#
#upregulation by RNAseq avg fold change mean+sd
RNAseq_Upreg_AvgFcMeanSd<-row.names(subset(WCPstimQuant, WCPstimQuant$log2FoldChangeStimUnstim>=(RNAseq_FcMean+RNAseq_FcSd)))
RNAseq_Upreg_AvgFcMeanSd_Num<-length(RNAseq_Upreg_AvgFcMeanSd)
print("RNAseq upregulation avg fold change > mean+sd")
print(RNAseq_Upreg_AvgFcMeanSd_Num)
write(RNAseq_Upreg_AvgFcMeanSd, file = file.path(".", analysisdir, "RNAseq_Upreg_AvgFcMeanSd_ProtId.txt"))
write(WCPstimQuant[RNAseq_Upreg_AvgFcMeanSd,"EnsemblId"], file = file.path(".", analysisdir, "RNAseq_Upreg_AvgFcMeanSd_EnsemblId.txt"))
#
#downregulation by RNAseq avg fold change mean-sd
RNAseq_Downreg_AvgFcMeanSd<-row.names(subset(WCPstimQuant, WCPstimQuant$log2FoldChangeStimUnstim<=(RNAseq_FcMean-RNAseq_FcSd)))
RNAseq_Downreg_AvgFcMeanSd_Num<-length(RNAseq_Downreg_AvgFcMeanSd)
print("RNAseq downregulation avg fold change < mean-sd")
print(RNAseq_Downreg_AvgFcMeanSd_Num)
write(RNAseq_Downreg_AvgFcMeanSd, file = file.path(".", analysisdir, "RNAseq_Downreg_AvgFcMeanSd_ProtId.txt"))
write(WCPstimQuant[RNAseq_Downreg_AvgFcMeanSd,"EnsemblId"], file = file.path(".", analysisdir, "RNAseq_Downreg_AvgFcMeanSd_EnsemblId.txt"))
#
#unchanged by RNAseq avg fold change mean
RNAseq_Unch_AvgFcMeanSd<-row.names(subset(WCPstimQuant, (WCPstimQuant$log2FoldChangeStimUnstim<(RNAseq_FcMean+RNAseq_FcSd) & WCPstimQuant$log2FoldChangeStimUnstim>(RNAseq_FcMean-RNAseq_FcSd))))
RNAseq_Unch_AvgFcMeanSd_Num<-length(RNAseq_Unch_AvgFcMeanSd)
print("RNAseq unchanged avg fold change,sd")
print(RNAseq_Unch_AvgFcMeanSd_Num)
write(RNAseq_Unch_AvgFcMeanSd, file = file.path(".", analysisdir, "RNAseq_Unch_AvgFcMeanSd_ProtId.txt"))
write(WCPstimQuant[RNAseq_Unch_AvgFcMeanSd,"EnsemblId"], file = file.path(".", analysisdir, "RNAseq_Unch_AvgFcMeanSd_EnsemblId.txt"))
#
###end determine RNAseq up and down regulation


###identify KeGG protein ubiquitylation state and analyze increase or decrease during stim
##proteins identified by a KeGG in PO1063_HLM or PO972_0hr-4hr, requires KeGG identification in both time points)
print("KeGG analysis...")
#PO1063 HLM (minus Nedd inhibitor) KeGG identified proteins
PO1063HLMkgg_0hr4hr_UbProts<-row.names(subset(WCPstimQuant, !is.na(WCPstimQuant$intensityweightedmeanratio.PO1063kgg_HLM)))
PO1063HLMkgg_0hr4hr_UbProts_Num<-length(PO1063HLMkgg_0hr4hr_UbProts)
print("PO1063HLM 0hr4hr fold change Ubiquitylated proteins")
print(PO1063HLMkgg_0hr4hr_UbProts_Num)
#PO972 kgg (time course 0-1-4 hr) KeGG identified proteins (identified at 0 AND 4 hour)
PO972kgg_0hr4hr_UbProts<-row.names(subset(WCPstimQuant, !is.na(WCPstimQuant$intensityweightedmeanratio.PO972kgg_4diff0)))
PO972kgg_0hr4hr_UbProts_Num<-length(PO972kgg_0hr4hr_UbProts)
print("PO972 0hr4hr fold change Ubiquitylated proteins")
print(PO972kgg_0hr4hr_UbProts_Num)
#PO1063HLMkgg-PO972kgg 0hr-4hr fold change set union
PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts <- union(PO1063HLMkgg_0hr4hr_UbProts, PO972kgg_0hr4hr_UbProts)
PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts_Num <- length(PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts)
print("PO972kgg-PO1063HLM 0hr4hr fold change Ubiquitylated proteins (set union)")
print(PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts_Num)
write(PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts, file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbProteins_setUnion.txt"))
#
#PO1063HLMkgg-PO972kgg 0hr-4hr fold change set intersection
PO972kggPO1063HLMkgg_All0hr4hrFc_Int_UbProts <- intersect(PO1063HLMkgg_0hr4hr_UbProts, PO972kgg_0hr4hr_UbProts)
PO972kggPO1063HLMkgg_All0hr4hrFc_Int_UbProts_Num <- length(PO972kggPO1063HLMkgg_All0hr4hrFc_Int_UbProts)
print("PO972kgg-PO1063HLM 0hr4hr fold change Ubiquitylated proteins (set intersection)")
print(PO972kggPO1063HLMkgg_All0hr4hrFc_Int_UbProts_Num)
write(PO972kggPO1063HLMkgg_All0hr4hrFc_Int_UbProts, file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbProteins_setInt.txt"))
#
##ubiquitination at Light intensity (4hr) but without Heavy intensity (0hr) pair
##or ubiquitination at Heavy intensity (0hr) but without Light intensity (4hr)
#PO1063HLMkgg
#PO1063HLMkgg 4hr unique ubiquitin
PO1063HLMkgg_4hr_Unique_UbProts <- row.names(subset(WCPstimQuant, (is.na(intensityweightedmeanratio.PO1063kgg_HLM) & is.na(avgintenH.PO1063kgg_HLM) & !is.na(avgintenL.PO1063kgg_HLM))))
PO1063HLMkgg_4hr_Unique_UbProts_Num <- length(PO1063HLMkgg_4hr_Unique_UbProts)
print("PO1063HLMkgg 4hr unique")
print(PO1063HLMkgg_4hr_Unique_UbProts_Num)
#
#PO1063HLMkgg 0hr unique ubiquitin
PO1063HLMkgg_0hr_Unique_UbProts <- row.names(subset(WCPstimQuant, (is.na(intensityweightedmeanratio.PO1063kgg_HLM) & !is.na(avgintenH.PO1063kgg_HLM) & is.na(avgintenL.PO1063kgg_HLM))))
PO1063HLMkgg_0hr_Unique_UbProts_Num <- length(PO1063HLMkgg_0hr_Unique_UbProts)
print("PO1063HLMkgg 0hr unique")
print(PO1063HLMkgg_0hr_Unique_UbProts_Num)
#
#PO1063HLMkgg all identified
PO1063HLMkgg_AllUbId <- Reduce(union, list(PO1063HLMkgg_0hr_Unique_UbProts, PO1063HLMkgg_4hr_Unique_UbProts, PO1063HLMkgg_0hr4hr_UbProts))
PO1063HLMkgg_AllUbId_Num <- length(PO1063HLMkgg_AllUbId)
print("PO1063HLMkgg all Ub identified")
print(PO1063HLMkgg_AllUbId_Num)
#
#PO1063HLMkgg 0hr unique - 4hr unique ubiquitin intersection
PO1063HLMkgg_0hrUnique4hrUnique_Int <- intersect(PO1063HLMkgg_4hr_Unique_UbProts,PO1063HLMkgg_0hr_Unique_UbProts)
PO1063HLMkgg_0hrUnique4hrUnique_Int_Num <- length(PO1063HLMkgg_0hrUnique4hrUnique_Int)
print ("PO1063HLMkgg 0hr - 4hr unique intersection (should be zero)")
print(PO1063HLMkgg_0hrUnique4hrUnique_Int_Num)
if(PO1063HLMkgg_0hrUnique4hrUnique_Int_Num>0){
	stop("PO1063HLMkgg 0hr - 4hr unique intersection should be zero")
}
#
PO1063HLMkgg_0hrUnique_0hr4hr_Int <- intersect(PO1063HLMkgg_0hr_Unique_UbProts, PO1063HLMkgg_0hr4hr_UbProts)
PO1063HLMkgg_0hrUnique_0hr4hr_Int_Num <- length(PO1063HLMkgg_0hrUnique_0hr4hr_Int)
print ("PO1063HLMkgg 0hr unique, 0hr4hr intersection (should be zero)")
print(PO1063HLMkgg_0hrUnique_0hr4hr_Int_Num)
if(PO1063HLMkgg_0hrUnique_0hr4hr_Int_Num>0){
	stop("PO1063HLMkgg 0hr unique, 0hr4hr intersection should be zero")
}
#
PO1063HLMkgg_4hrUnique_0hr4hr_Int <- intersect(PO1063HLMkgg_4hr_Unique_UbProts, PO1063HLMkgg_0hr4hr_UbProts)
PO1063HLMkgg_4hrUnique_0hr4hr_Int_Num <- length(PO1063HLMkgg_4hrUnique_0hr4hr_Int)
print ("PO1063HLMkgg 4hr unique, 0hr4hr intersection (should be zero)")
print(PO1063HLMkgg_4hrUnique_0hr4hr_Int_Num)
if(PO1063HLMkgg_4hrUnique_0hr4hr_Int_Num>0){
	stop("PO1063HLMkgg 4hr unique, 0hr4hr intersection should be zero")
}
#
#PO972kgg
#PO972kgg 4hr unique ubiquitin
PO972kgg_4hr_Unique_UbProts <- row.names(subset(WCPstimQuant, (is.na(WCPstimQuant$avgintenL.PO972kgg_0hr) & (!is.na(WCPstimQuant$avgintenH.PO972kgg_4hr) & !is.na(WCPstimQuant$avgintenL.PO972kgg_4hr)))))
PO972kgg_4hr_Unique_UbProts_Num <- length(PO972kgg_4hr_Unique_UbProts)
print("PO972kgg 4hr unique")
print(PO972kgg_4hr_Unique_UbProts_Num)
#
#PO972kgg 0hr unique ubiquitin
PO972kgg_0hr_Unique_UbProts <- row.names(subset(WCPstimQuant, ((!is.na(WCPstimQuant$avgintenH.PO972kgg_0hr) & !is.na(WCPstimQuant$avgintenL.PO972kgg_0hr)) & is.na(WCPstimQuant$avgintenL.PO972kgg_4hr))))
PO972kgg_0hr_Unique_UbProts_Num <- length(PO972kgg_0hr_Unique_UbProts)
print("PO972kgg 0hr")
print(PO972kgg_0hr_Unique_UbProts_Num)
#
#PO972kgg all identified
PO972kgg_AllUbId <- Reduce(union, list(PO972kgg_0hr_Unique_UbProts, PO972kgg_4hr_Unique_UbProts, PO972kgg_0hr4hr_UbProts))
PO972kgg_AllUbId_Num <- length(PO972kgg_AllUbId)
print("PO972kgg all Ub identified")
print(PO972kgg_AllUbId_Num)
#
#PO972kgg 0hr unique - 4hr unique ubiquitin intersection
PO972kgg_0hrUnique4hrUnique_Int <- intersect(PO972kgg_4hr_Unique_UbProts,PO972kgg_0hr_Unique_UbProts)
PO1063HLMkgg_0hrUnique4hrUnique_Int_Num <- length(PO1063HLMkgg_0hrUnique4hrUnique_Int)
print ("PO972kgg 0hr - 4hr unique intersection (should be 0)")
print(PO1063HLMkgg_0hrUnique4hrUnique_Int_Num)
if(PO1063HLMkgg_0hrUnique4hrUnique_Int_Num>0){
	stop("PO972kgg 0hr - 4hr unique intersection should be 0")
}
#
PO972kgg_0hrUnique_0hr4hr_Int <- intersect(PO972kgg_0hr_Unique_UbProts, PO972kgg_0hr4hr_UbProts)
PO972kgg_0hrUnique_0hr4hr_Int_Num <- length(PO972kgg_0hrUnique_0hr4hr_Int)
print ("PO972kgg 0hr unique, 0hr4hr intersection (should be zero)")
print(PO972kgg_0hrUnique_0hr4hr_Int_Num)
if(PO972kgg_0hrUnique_0hr4hr_Int_Num>0){
	stop("PO972kgg 0hr unique, 0hr4hr intersection should be zero")
}
#
PO972kgg_4hrUnique_0hr4hr_Int <- intersect(PO972kgg_4hr_Unique_UbProts, PO972kgg_0hr4hr_UbProts)
PO972kgg_4hrUnique_0hr4hr_Int_Num <- length(PO972kgg_4hrUnique_0hr4hr_Int)
print ("PO972kgg 4hr unique, 0hr4hr intersection (should be zero)")
print(PO972kgg_4hrUnique_0hr4hr_Int_Num)
if(PO972kgg_4hrUnique_0hr4hr_Int_Num>0){
	stop("PO972kgg 4hr unique, 0hr4hr intersection should be zero")
}
#
#PO1063HLMkgg all 0hr
PO1063HLMkgg_0hr_All_UbProts <- union(PO1063HLMkgg_0hr_Unique_UbProts, PO1063HLMkgg_0hr4hr_UbProts)
PO1063HLMkgg_0hr_All_UbProts_Num <- length(PO1063HLMkgg_0hr_All_UbProts)
print("PO1063HLMkgg all 0hr Ub prots")
print(PO1063HLMkgg_0hr_All_UbProts_Num)
#
#PO972kgg all 0hr
PO972kgg_0hr_All_UbProts <- union(PO972kgg_0hr_Unique_UbProts, PO972kgg_0hr4hr_UbProts)
PO972kgg_0hr_All_UbProts_Num <- length(PO972kgg_0hr_All_UbProts)
print("PO972kgg all 0hr Ub prots")
print(PO972kgg_0hr_All_UbProts_Num)
#
#PO972kgg-PO1063HLMkgg 0hr union
PO972kggPO1063HLMkgg_0hr_All_UbProts <- union(PO972kgg_0hr_All_UbProts, PO1063HLMkgg_0hr_All_UbProts)
PO972kggPO1063HLMkgg_0hr_All_UbProts_Num <- length(PO972kggPO1063HLMkgg_0hr_All_UbProts)
print("PO972kgg-PO1063HLMkgg all 0hr Ub prots")
print(PO972kggPO1063HLMkgg_0hr_All_UbProts_Num)
#
#PO1063HLMkgg all 4hr
PO1063HLMkgg_4hr_All_UbProts <- union(PO1063HLMkgg_4hr_Unique_UbProts, PO1063HLMkgg_0hr4hr_UbProts)
PO1063HLMkgg_4hr_All_UbProts_Num <- length(PO1063HLMkgg_4hr_All_UbProts)
print("PO1063HLMkgg all 4hr Ub prots")
print(PO1063HLMkgg_4hr_All_UbProts_Num)
#
#PO972kgg all 4hr
PO972kgg_4hr_All_UbProts <- union(PO972kgg_4hr_Unique_UbProts, PO972kgg_0hr4hr_UbProts)
PO972kgg_4hr_All_UbProts_Num <- length(PO972kgg_4hr_All_UbProts)
print("PO972kgg all 4hr Ub prots")
print(PO972kgg_4hr_All_UbProts_Num)
#
#PO972kgg-PO1063HLMkgg 4hr union
PO972kggPO1063HLMkgg_4hr_All_UbProts <- union(PO972kgg_4hr_All_UbProts, PO1063HLMkgg_4hr_All_UbProts)
PO972kggPO1063HLMkgg_4hr_All_UbProts_Num <- length(PO972kggPO1063HLMkgg_4hr_All_UbProts)
print("PO972kgg-PO1063HLMkgg all 4hr Ub prots")
print(PO972kggPO1063HLMkgg_4hr_All_UbProts_Num)
#
#PO972kgg-PO1063HLMkgg 0hr unique Ub prots
PO972kggPO1063HLMkgg_0hr_Unique_UbProts <- setdiff(PO972kggPO1063HLMkgg_0hr_All_UbProts, PO972kggPO1063HLMkgg_4hr_All_UbProts)
PO972kggPO1063HLMkgg_0hr_Unique_UbProts_Num <- length(PO972kggPO1063HLMkgg_0hr_Unique_UbProts)
print("PO972kgg-PO1063HLMkgg 0hr unique Ub prots")
print(PO972kggPO1063HLMkgg_0hr_Unique_UbProts_Num);
#
#PO972kgg-PO1063HLMkgg 4hr unique Ub prots
PO972kggPO1063HLMkgg_4hr_Unique_UbProts <- setdiff(PO972kggPO1063HLMkgg_4hr_All_UbProts, PO972kggPO1063HLMkgg_0hr_All_UbProts)
PO972kggPO1063HLMkgg_4hr_Unique_UbProts_Num <- length(PO972kggPO1063HLMkgg_4hr_Unique_UbProts)
print("PO972kgg-PO1063HLMkgg 4hr unique Ub prots")
print(PO972kggPO1063HLMkgg_4hr_Unique_UbProts_Num);
#
#PO972kgg-PO1063HLMkgg 0hr4hr intersection Ub prots
PO972kggPO1063HLMkgg_0hr4hr_Int_UbProts <- intersect(PO972kggPO1063HLMkgg_0hr_All_UbProts, PO972kggPO1063HLMkgg_4hr_All_UbProts)
PO972kggPO1063HLMkgg_0hr4hr_Int_UbProts_Num <- length(PO972kggPO1063HLMkgg_0hr4hr_Int_UbProts)
print("PO972kgg-PO1063HLMkgg 0hr4hr intersect Ub prots")
print(PO972kggPO1063HLMkgg_0hr4hr_Int_UbProts_Num)
#
#PO972kgg-PO1063HLMkgg 0hr4hr union Ub prots
PO972kggPO1063HLMkgg_0hr4hr_Union_UbProts <- union(PO972kggPO1063HLMkgg_0hr_All_UbProts, PO972kggPO1063HLMkgg_4hr_All_UbProts)
PO972kggPO1063HLMkgg_0hr4hr_Union_UbProts_Num <- length(PO972kggPO1063HLMkgg_0hr4hr_Union_UbProts)
print("PO972kgg-PO1063HLMkgg 0hr4hr union Ub prots")
print(PO972kggPO1063HLMkgg_0hr4hr_Union_UbProts_Num)
#
#add column of 0hr unique and 4hr unique to data frame
WCPstimQuant$Ub_0hrUnique<-NA
WCPstimQuant$Ub_4hrUnique<-NA
WCPstimQuant[PO972kggPO1063HLMkgg_0hr_Unique_UbProts, "Ub_0hrUnique"]<-1
WCPstimQuant[PO972kggPO1063HLMkgg_4hr_Unique_UbProts, "Ub_4hrUnique"]<-1
#
PO972kggPO1063HLMkgg_AllUbId_Int<-intersect(PO972kgg_AllUbId, PO1063HLMkgg_AllUbId)
PO972kggPO1063HLMkgg_AllUbId_Int_Num<-length(PO972kggPO1063HLMkgg_AllUbId_Int)
print("PO972kgg-PO1063HLMkgg all identified intersection")
print(PO972kggPO1063HLMkgg_AllUbId_Int_Num)
#
#
##KeGG ubiquitination quantification analysis
#mean,sd of PO1063HLM-PO972kgg ubiquitination fold change
PO972kggPO1063HLMkgg_UbFc_Mean<-mean(WCPstimQuant$UbFcAvg, na.rm=TRUE)
PO972kggPO1063HLMkgg_UbFc_Sd<-sd(WCPstimQuant$UbFcAvg, na.rm=TRUE)
print("PO972kgg-PO1063HLMkgg protein ubiquitylation fold change (mean,sd)")
print(PO972kggPO1063HLMkgg_UbFc_Mean)
print(PO972kggPO1063HLMkgg_UbFc_Sd)
#protein Ub upregulation by KeGG avg fold change mean+sd
PO972kggPO1063HLMkgg_Ub_Upreg_AvgFcMeanSd<-row.names(subset(WCPstimQuant, WCPstimQuant$UbFcAvg>=(PO972kggPO1063HLMkgg_UbFc_Mean+PO972kggPO1063HLMkgg_UbFc_Sd)))
PO972kggPO1063HLMkgg_Ub_Upreg_AvgFcMeanSd_Num<-length(PO972kggPO1063HLMkgg_Ub_Upreg_AvgFcMeanSd)
print("protein Ub upregulation avg fold change > mean+sd")
print(PO972kggPO1063HLMkgg_Ub_Upreg_AvgFcMeanSd_Num)
#write(PO972kggPO1063HLMkgg_Ub_Upreg_AvgFcMeanSd, file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_Ub_Upreg_AvgFcMeanSd.txt"))
#write(WCPstimQuant[PO972kggPO1063HLMkgg_Ub_Upreg_AvgFcMeanSd,"GeneId"], file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_GeneID_Ub_Upreg_AvgFcMeanSd.txt"))
#protein Ub downregulation by KeGG avg fold change mean+sd
PO972kggPO1063HLMkgg_Ub_Downreg_AvgFcMeanSd<-row.names(subset(WCPstimQuant, WCPstimQuant$UbFcAvg<=(PO972kggPO1063HLMkgg_UbFc_Mean-PO972kggPO1063HLMkgg_UbFc_Sd)))
PO972kggPO1063HLMkgg_Ub_Downreg_AvgFcMeanSd_Num<-length(PO972kggPO1063HLMkgg_Ub_Downreg_AvgFcMeanSd)
print("protein Ub downregulation avg fold change < mean-sd")
print(PO972kggPO1063HLMkgg_Ub_Downreg_AvgFcMeanSd_Num)
#write(PO972kggPO1063HLMkgg_Ub_Downreg_AvgFcMeanSd, file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_Ub_Downreg_AvgFcMeanSd.txt"))
#write(WCPstimQuant[PO972kggPO1063HLMkgg_Ub_Downreg_AvgFcMeanSd,"GeneId"], file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_GeneID_Ub_Downreg_AvgFcMeanSd.txt"))
#protein Ub unchanged by KeGG avg fold change mean +/- sd
PO972kggPO1063HLMkgg_Ub_Unch_AvgFcMeanSd<-row.names(subset(WCPstimQuant, (WCPstimQuant$UbFcAvg>(PO972kggPO1063HLMkgg_UbFc_Mean-PO972kggPO1063HLMkgg_UbFc_Sd) & WCPstimQuant$UbFcAvg<(PO972kggPO1063HLMkgg_UbFc_Mean+PO972kggPO1063HLMkgg_UbFc_Sd))))
PO972kggPO1063HLMkgg_Ub_Unch_AvgFcMeanSd_Num<-length(PO972kggPO1063HLMkgg_Ub_Unch_AvgFcMeanSd)
print("protein Ub unch avg fold change >< mean+/-sd")
print(PO972kggPO1063HLMkgg_Ub_Unch_AvgFcMeanSd_Num)
#write(PO972kggPO1063HLMkgg_Ub_Unch_AvgFcMeanSd, file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_Ub_Unch_AvgFcMeanSd.txt"))
#write(WCPstimQuant[PO972kggPO1063HLMkgg_Ub_Unch_AvgFcMeanSd,"GeneId"], file = file.path(".", analysisdir, "PO972kggPO1063HLM_GeneID_Ub_Unch_AvgFcMeanSd.txt"))
#
#
##normalize Ub fold change by total protein abundance fold change
##########WCPstimQuant$UbFcAvg_4hr_NormWcpFc <- apply(WCPstimQuant, MARGIN=1, FUN=function(x) (as.numeric(x[["UbFcAvg"]])-as.numeric(x[["Stim4hr_AvgFc"]])))
WCPstimQuant$PO972kgg_NormWcpFc <- apply(WCPstimQuant, MARGIN=1, FUN=function(x) (as.numeric(x[["intensityweightedmeanratio.PO972kgg_4diff0"]])-as.numeric(x[["IGDStim4hr_fc"]])))
WCPstimQuant$PO1063HLMkgg_NormWcpFc <- apply(WCPstimQuant, MARGIN=1, FUN=function(x) (as.numeric(x[["intensityweightedmeanratio.PO1063kgg_HLM"]])-as.numeric(x[["HLMStim4hr_fc"]])))
WCPstimQuant$PO1063HLPkgg_NormWcpFc <- apply(WCPstimQuant, MARGIN=1, FUN=function(x) (as.numeric(x[["intensityweightedmeanratio.PO1063kgg_HLP"]])-as.numeric(x[["HLPStim4hr_fc"]])))
#
WCPstimQuant$UbFcAvg_4hr_NormWcpFc <- apply(WCPstimQuant, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["PO972kgg_NormWcpFc"]],x[["PO1063HLMkgg_NormWcpFc"]])),na.rm=TRUE))

#
#ubiquitination fold change norm to WCP upreg by fc >0 
PO972kggPO1063HLMkgg_UbNormFc_Upreg_FcGt0<-row.names(subset(WCPstimQuant, WCPstimQuant$UbFcAvg_4hr_NormWcpFc>0))
PO972kggPO1063HLMkgg_UbNormFc_Upreg_FcGt0_Num<-length(PO972kggPO1063HLMkgg_UbNormFc_Upreg_FcGt0)
print("PO972kgg-PO1063HLMkgg protein ubiquitylation fold change upreg norm to WCP (>0 norm)")
print(PO972kggPO1063HLMkgg_UbNormFc_Upreg_FcGt0_Num)
write(PO972kggPO1063HLMkgg_UbNormFc_Upreg_FcGt0, file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNorm_Upreg_FcGt0_ProtId.txt"))
write(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Upreg_FcGt0,"GeneId.x"], file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNorm_Upreg_FcGt0_GeneId.txt"))
#ubiquitination fold change norm to WCP downreg by fc < 0
PO972kggPO1063HLMkgg_UbNormFc_Downreg_FcGt0<-row.names(subset(WCPstimQuant, WCPstimQuant$UbFcAvg_4hr_NormWcpFc<0))
PO972kggPO1063HLMkgg_UbNormFc_Downreg_FcGt0_Num<-length(PO972kggPO1063HLMkgg_UbNormFc_Downreg_FcGt0)
print("PO972kgg-PO1063HLMkgg protein ubiquitylation fold change downreg norm to WCP (<0 norm)")
print(PO972kggPO1063HLMkgg_UbNormFc_Downreg_FcGt0_Num)
write(PO972kggPO1063HLMkgg_UbNormFc_Downreg_FcGt0, file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNorm_Downreg_FcGt0_ProtId.txt"))
write(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Downreg_FcGt0,"GeneId.x"], file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNorm_Downreg_FcGt0_GeneId.txt"))
#
#mean,sd of PO1063HLM-PO972kgg ubiquitination fold change norm to WCP
PO972kggPO1063HLMkgg_UbNormFc_Mean<-mean(WCPstimQuant$UbFcAvg_4hr_NormWcpFc, na.rm=TRUE)
PO972kggPO1063HLMkgg_UbNormFc_Sd<-sd(WCPstimQuant$UbFcAvg_4hr_NormWcpFc, na.rm=TRUE)
print("PO972kgg-PO1063HLMkgg protein ubiquitylation fold change norm to WCP (mean,sd)")
print(PO972kggPO1063HLMkgg_UbNormFc_Mean)
print(PO972kggPO1063HLMkgg_UbNormFc_Sd)
#protein Ub upregulation by KeGG avg fold change norm to WCP mean+sd
PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd<-row.names(subset(WCPstimQuant, WCPstimQuant$UbFcAvg_4hr_NormWcpFc>=(PO972kggPO1063HLMkgg_UbNormFc_Mean+PO972kggPO1063HLMkgg_UbNormFc_Sd)))
PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd_Num<-length(PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd)
print("protein Ub upregulation avg fold change norm to WCP > mean+sd")
print(PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd_Num)
write(PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd, file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd_ProtId.txt"))
write(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd,"GeneId.x"], file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd_GeneId.txt"))
#protein Ub downregulation by KeGG avg fold change norm to WCP mean+sd
PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd<-row.names(subset(WCPstimQuant, WCPstimQuant$UbFcAvg_4hr_NormWcpFc<=(PO972kggPO1063HLMkgg_UbNormFc_Mean-PO972kggPO1063HLMkgg_UbNormFc_Sd)))
PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd_Num<-length(PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd)
print("protein Ub downregulation avg fold change norm to WCP < mean-sd")
print(PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd_Num)
write(PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd, file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd_ProtId.txt"))
write(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd,"GeneId.x"], file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd_GeneId.txt"))
#protein Ub unchanged by KeGG avg fold change norm to WCP mean +/- sd
PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanSd<-row.names(subset(WCPstimQuant, (WCPstimQuant$UbFcAvg_4hr_NormWcpFc>(PO972kggPO1063HLMkgg_UbNormFc_Mean-PO972kggPO1063HLMkgg_UbNormFc_Sd) & WCPstimQuant$UbFcAvg_4hr_NormWcpFc<(PO972kggPO1063HLMkgg_UbNormFc_Mean+PO972kggPO1063HLMkgg_UbNormFc_Sd))))
PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanSd_Num<-length(PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanSd)
print("protein Ub unch avg fold change norm to WCP >< mean+/-sd")
print(PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanSd_Num)
write(PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanSd, file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanSd_ProtId.txt"))
write(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanSd,"GeneId.x"], file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanSd_GeneId.txt"))
#
ubdeltathreshold <- log2(1.25)
#protein Ub upregulation by KeGG avg fold change +/- threshold
PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanThreshold<-row.names(subset(WCPstimQuant, WCPstimQuant$UbFcAvg_4hr_NormWcpFc >= ubdeltathreshold))
PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanThreshold_Num<-length(PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanThreshold)
print("protein Ub upregulation avg fold change norm to WCP > threshold")
print(PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanThreshold_Num)
write(PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanThreshold, file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanThreshold_ProtID.txt"))
write(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanThreshold,"GeneId.x"], file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanThreshold_GeneId.txt"))
#protein Ub downregulation by KeGG avg fold change +/- threshold
PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanThreshold<-row.names(subset(WCPstimQuant, WCPstimQuant$UbFcAvg_4hr_NormWcpFc <= -1*ubdeltathreshold))
PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanThreshold_Num<-length(PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanThreshold)
print("protein Ub downregulation avg fold change norm to WCP < threshold")
print(PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanThreshold_Num)
write(PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanThreshold, file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanThreshold_ProtId.txt"))
write(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanThreshold,"GeneId.x"], file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanThreshold_GeneId.txt"))
#protein Ub unchanged by KeGG avg fold change norm to WCP mean +/- sd
PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanThreshold<-row.names(subset(WCPstimQuant, (WCPstimQuant$UbFcAvg_4hr_NormWcpFc > -1*ubdeltathreshold & WCPstimQuant$UbFcAvg_4hr_NormWcpFc < ubdeltathreshold)))
PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanThreshold_Num<-length(PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanThreshold)
print("protein Ub unch avg fold change norm to WCP >< threshold")
print(PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanThreshold_Num)
write(PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanThreshold, file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanThreshold_ProtId.txt"))
write(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanThreshold,"GeneId.x"], file = file.path(".", analysisdir, "PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanThreshold_GeneId.txt"))
#
#
####
####
# this section has not been checked
#ubiquitination at 1 hr
#PO972 kgg (time course 0-1-4 hr) KeGG identified proteins (identified at 0 AND 1 hour)
PO972kgg_UbProts_1hr<-row.names(subset(WCPstimQuant, (WCPstimQuant$numsites.PO972kgg_0hr>0 & WCPstimQuant$numsites.PO972kgg_1hr>0)))
PO972kgg_UbProts_1hr_Num<-length(PO972kgg_UbProts_1hr)
print("PO972 1hr Ubiquitylated proteins")
print(PO972kgg_UbProts_1hr_Num)
#normalize 1hr Ub fold change by total protein abundance fold change
WCPstimQuant$UbFcAvg_1hr_NormWcpFc <- apply(WCPstimQuant, MARGIN=1, FUN=function(x) (as.numeric(x[["intensityweightedmeanratio.PO972kgg_1diff0"]])-as.numeric(x[["IGDStim1hr_fc"]])))
PO972kgg_UbNormFc_1hr_Upreg_NormToWcp <- row.names(subset(WCPstimQuant, WCPstimQuant$UbFcAvg_1hr_NormWcpFc>0))
PO972kgg_UbNormFc_1hr_Upreg_NormToWcp_Num <- length(PO972kgg_UbNormFc_1hr_Upreg_NormToWcp)
print("PO972 1hr Ub upreg norm to Wcp")
print(PO972kgg_UbNormFc_1hr_Upreg_NormToWcp_Num)
#
#ubiquitination at 1-4hr
#PO972 kgg (time course 0-1-4 hr) KeGG identified proteins (identified at 1 AND 4 hour)
PO972kgg_UbProts_1hr4hr<-row.names(subset(WCPstimQuant, (WCPstimQuant$numsites.PO972kgg_1hr>0 & WCPstimQuant$numsites.PO972kgg_4hr>0)))
PO972kgg_UbProts_1hr4hr_Num<-length(PO972kgg_UbProts_1hr4hr)
print("PO972 1hr-4hr Ubiquitylated proteins")
print(PO972kgg_UbProts_1hr4hr_Num)
#normalize 1-4hr Ub fold change by total protein abundance fold change
WCPstimQuant$UbFcAvg_1hr4hr_NormWcpFc <- apply(WCPstimQuant, MARGIN=1, FUN=function(x) (as.numeric(x[["intensityweightedmeanratio.PO972kgg_4diff1"]])-as.numeric(x[["IGD1hr4hr_fc"]])))
PO972kgg_UbNormFc_1hr4hr_Upreg_NormToWcp <- row.names(subset(WCPstimQuant, WCPstimQuant$UbFcAvg_1hr4hr_NormWcpFc>0))
PO972kgg_UbNormFc_1hr4hr_Upreg_NormToWcp_Num <- length(PO972kgg_UbNormFc_1hr4hr_Upreg_NormToWcp)
print("PO972 1hr-4hr Ub upreg norm to Wcp")
print(PO972kgg_UbNormFc_1hr4hr_Upreg_NormToWcp_Num)
#
####
####
#
###end KeGG protein ubiquitylation state analysis


###potential Cullin substrate analysis
#
##proteins with decreased ubiquitination in + Nedd inhibitor vs - Nedd inhibitor experiments

predubsubs<-row.names(subset(WCPstimQuant, (WCPstimQuant$PO1063HLPkgg_NormWcpFc<WCPstimQuant$PO1063HLMkgg_NormWcpFc+1 & WCPstimQuant$PO1063HLPkgg_NormWcpFc>WCPstimQuant$PO1063HLMkgg_NormWcpFc-1)))
predubsubs_num<-length(predubsubs)
print("predubsubs_num")
print(predubsubs_num)
predneddsubs<-row.names(subset(WCPstimQuant, (WCPstimQuant$PO1063HLPkgg_NormWcpFc>WCPstimQuant$PO1063HLMkgg_NormWcpFc+1 | WCPstimQuant$PO1063HLPkgg_NormWcpFc<WCPstimQuant$PO1063HLMkgg_NormWcpFc-1)))
predneddsubs_num<-length(predneddsubs)
print("predneddsubs_num")
print(predneddsubs_num)

WCPstimQuant$PO1063HLMHLPkgg_NormUbFc_DiffReg <- apply(WCPstimQuant, MARGIN=1, FUN=function(x) (as.numeric(x[["PO1063HLPkgg_NormWcpFc"]])-as.numeric(x[["PO1063HLMkgg_NormWcpFc"]])))
PO1063HLMHLPkgg_NormUbFc_DiffReg_Mean <- mean(WCPstimQuant$PO1063HLMHLPkgg_NormUbFc_DiffReg, na.rm=TRUE)
PO1063HLMHLPkgg_NormUbFc_DiffReg_Sd <- sd(WCPstimQuant$PO1063HLMHLPkgg_NormUbFc_DiffReg, na.rm=TRUE)
print("PO1063HLP-PO1063HLM Norm Ub expression fold change (mean,sd)")
print(PO1063HLMHLPkgg_NormUbFc_DiffReg_Mean)
print(PO1063HLMHLPkgg_NormUbFc_DiffReg_Sd)
#differential regulation decreased > mean - sd
PO1063HLMHLPkgg_UbDiffReg<-row.names(subset(WCPstimQuant, (WCPstimQuant$PO1063HLMHLPkgg_NormUbFc_DiffReg<=PO1063HLMHLPkgg_NormUbFc_DiffReg_Mean-PO1063HLMHLPkgg_NormUbFc_DiffReg_Sd)))
PO1063HLMHLPkgg_UbDiffReg_Num<-length(PO1063HLMHLPkgg_UbDiffReg)
print("protein Ub downregulated w/ Nedd inhibitor treatment")
print(PO1063HLMHLPkgg_UbDiffReg_Num)
write.table(WCPstimQuant[PO1063HLMHLPkgg_UbDiffReg,], file=file.path(".", analysisdir, "NeddylationInhibitorDownregProteins.txt"), sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#
###########proteins with increased ubiquitination in + Nedd inhibitor vs - Nedd inhibitor experiments
##########PO1063HLMHLPkgg_UbDiffUpReg<-row.names(subset(WCPstimQuant, (intensityweightedmeanratio.PO1063kgg_HLM-intensityweightedmeanratio.PO1063kgg_HLP) < -1))
##########PO1063HLMHLPkgg_UbDiffUpReg_Num<-length(PO1063HLMHLPkgg_UbDiffUpReg)
##########print("protein Ub upregulated w/ Nedd inhibitor treatment")
##########print(PO1063HLMHLPkgg_UbDiffUpReg_Num)
##########NeddUpreg<-WCPstimQuant[row.names(WCPstimQuant)%in%PO1063HLMHLPkgg_UbDiffUpReg,]
##########write.table(NeddUpreg, file="NeddylationInhibitorUpregProteins.txt", sep="\t", quote=FALSE, na="NA", dec=".")
#
#proteins with increased wcp expression with nedd inhibitor based on mean +/- 1 sd of log2 fold change difference
HLMHLPStim4hr_DiffReg_Mean<-mean(WCPstimQuant$HLMHLPStim4hr_DiffReg, na.rm=TRUE)
HLMHLPStim4hr_DiffReg_Sd<-sd(WCPstimQuant$HLMHLPStim4hr_DiffReg, na.rm=TRUE)
print("PO1063HLP-PO1063HLM protein expression fold change (mean,sd)")
print(HLMHLPStim4hr_DiffReg_Mean)
print(HLMHLPStim4hr_DiffReg_Sd)
#HLP-HLM protein expression upregulation by avg fold change mean+sd
PO1063HLPPO1063HLM_Upreg_AvgFcMeanSd<-row.names(subset(WCPstimQuant, WCPstimQuant$HLMHLPStim4hr_DiffReg>=(HLMHLPStim4hr_DiffReg_Mean+(1*HLMHLPStim4hr_DiffReg_Sd))))
PO1063HLPPO1063HLM_Upreg_AvgFcMeanSd_Num<-length(PO1063HLPPO1063HLM_Upreg_AvgFcMeanSd)
print("PO1063HLP-PO1063HLM protein expression upregulation avg fold change > mean+sd")
print(PO1063HLPPO1063HLM_Upreg_AvgFcMeanSd_Num)
#HLP-HLM protein expression downregulation by avg fold change mean+sd
PO1063HLPPO1063HLM_Downreg_AvgFcMeanSd<-row.names(subset(WCPstimQuant, WCPstimQuant$HLMHLPStim4hr_DiffReg<=(HLMHLPStim4hr_DiffReg_Mean-(1*HLMHLPStim4hr_DiffReg_Sd))))
PO1063HLPPO1063HLM_Downreg_AvgFcMeanSd_Num<-length(PO1063HLPPO1063HLM_Downreg_AvgFcMeanSd)
print("PO1063HLP-PO1063HLM protein expression downregulation avg fold change < mean+sd")
print(PO1063HLPPO1063HLM_Downreg_AvgFcMeanSd_Num)
#HLP-HLM protein expression unchanged by avg fold change mean+sd
PO1063HLPPO1063HLM_Unch_AvgFcMeanSd<-row.names(subset(WCPstimQuant, (WCPstimQuant$HLMHLPStim4hr_DiffReg>(HLMHLPStim4hr_DiffReg_Mean-(1*HLMHLPStim4hr_DiffReg_Sd)) & WCPstimQuant$HLMHLPStim4hr_DiffReg<(HLMHLPStim4hr_DiffReg_Mean+(1*HLMHLPStim4hr_DiffReg_Sd)))))
PO1063HLPPO1063HLM_Unch_AvgFcMeanSd_Num<-length(PO1063HLPPO1063HLM_Unch_AvgFcMeanSd)
print("PO1063HLP-PO1063HLM protein expression unchanged avg fold change < mean+sd")
print(PO1063HLPPO1063HLM_Unch_AvgFcMeanSd_Num)
#
#proteins with upregulated wcp expression and downregulated ubiquitination
HLMHLP_NormUbDownReg_WcpUpReg <- intersect(PO1063HLPPO1063HLM_Upreg_AvgFcMeanSd,PO1063HLMHLPkgg_UbDiffReg)
print(HLMHLP_NormUbDownReg_WcpUpReg)
#
###end potential Cullin substrate analysis


##data analysis to predict degradative vs. non-degradative ubiquitination
#increased ub proteins
IncreaseUb <- PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanThreshold
##########IncreaseUb <- PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd
IncreaseUb_Num <- length(IncreaseUb)
print("substrates w/ inc Ub identified at least once")
print(IncreaseUb_Num)
IncreaseUb_PO972kggPO1063HLMkggId <- intersect(PO972kggPO1063HLMkgg_All0hr4hrFc_Int_UbProts, IncreaseUb)
IncreaseUb_PO972kggPO1063HLMkggId_Num <- length(IncreaseUb_PO972kggPO1063HLMkggId)
print("substrates w/ inc Ub and identified in both PO972kgg and PO1063kgg")
print(IncreaseUb_PO972kggPO1063HLMkggId_Num)
#
#decreased ub proteins
DecreaseUb <- PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanThreshold
##########DecreaseUb <- PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd
DecreaseUb_Num <- length(DecreaseUb)
print("substrates w/ dec Ub identified at least once")
print(DecreaseUb_Num)
DecreaseUb_PO972kggPO1063HLMkggId <- intersect(PO972kggPO1063HLMkgg_All0hr4hrFc_Int_UbProts, DecreaseUb)
DecreaseUb_PO972kggPO1063HLMkggId_Num <- length(DecreaseUb_PO972kggPO1063HLMkggId)
print("substrates w/ dec Ub and identified in both PO972kgg and PO1063kgg")
print(DecreaseUb_PO972kggPO1063HLMkggId_Num)
#
#unchanged ub proteins
UnchangeUb <- PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanThreshold
##########UnchangeUb <- PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanSd
UnchangeUb_Num <- length(UnchangeUb)
print("substrates w/ unch Ub identified at least once")
print(UnchangeUb_Num)
UnchangeUb_PO972kggPO1063HLMkggId <- intersect(PO972kggPO1063HLMkgg_All0hr4hrFc_Int_UbProts, UnchangeUb)
UnchangeUb_PO972kggPO1063HLMkggId_Num <- length(UnchangeUb_PO972kggPO1063HLMkggId)
print("substrates w/ unch Ub and identified in both PO972kgg and PO1063kgg")
print(UnchangeUb_PO972kggPO1063HLMkggId_Num)
#
#increased/decreased RNA and WCP
UbDegRegion_RnaInc<-RNAseq_Upreg_AvgFcPval
UbDegRegion_RnaDec<-RNAseq_Downreg_AvgFcPval
UbDegRegion_RnaUch<-RNAseq_Unch_AvgFcPval
#UbDegRegion_WcpInc<-PO972PO1063PO1075_Wcp_Upreg_AvgFcPval
#UbDegRegion_WcpDec<-PO972PO1063PO1075_Wcp_Downreg_AvgFcPval
#UbDegRegion_WcpUch<-PO972PO1063PO1075_Wcp_Unch_AvgFcPval
#
#UbDegRegion_RnaInc<-RNAseq_Upreg_AvgFcMeanSd
#UbDegRegion_RnaDec<-RNAseq_Downreg_AvgFcMeanSd
#UbDegRegion_RnaUch<-RNAseq_Unch_AvgFcMeanSd
UbDegRegion_WcpInc<-PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd
UbDegRegion_WcpDec<-PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd
UbDegRegion_WcpUch<-PO972PO1063PO1075_Wcp_Unch_AvgFcMeanSd
#
##divide protein and RNA differential expression into regions, overlay KeGG data
#dprotein downreg, RNAseq upreg, ubiquitylated. "1"
ProtDownregRnaUpreg<-intersect(UbDegRegion_WcpDec, UbDegRegion_RnaInc)
ProtDownregRnaUpreg_Num<-length(ProtDownregRnaUpreg)
ProtDownregRnaUpregUbUpreg<-intersect(ProtDownregRnaUpreg, IncreaseUb)
ProtDownregRnaUpregUbUpreg_Num<-length(ProtDownregRnaUpregUbUpreg)
ProtDownregRnaUpregUbDownreg<-intersect(ProtDownregRnaUpreg, DecreaseUb)
ProtDownregRnaUpregUbDownreg_Num<-length(ProtDownregRnaUpregUbDownreg)
ProtDownregRnaUpregUbUnch<-intersect(ProtDownregRnaUpreg, UnchangeUb)
ProtDownregRnaUpregUbUnch_Num<-length(ProtDownregRnaUpregUbUnch)
print("protein downreg, RNAseq upreg, ubiquitylated upreg/downreg :1")
print(ProtDownregRnaUpreg_Num)
print(ProtDownregRnaUpregUbUpreg_Num)
print(ProtDownregRnaUpregUbDownreg_Num)
print(ProtDownregRnaUpregUbUnch_Num)
#####write(ProtDownregRnaUpreg, file = file.path(".", analysisdir, "ProtDownRnaUp_1_ProtId.txt"))
write(ProtDownregRnaUpregUbUpreg, file = file.path(".", analysisdir, "ProtDownRnaUpUbUp_1_ProtId.txt"))
#####write(WCPstimQuant[ProtDownregRnaUpreg,"GeneId.x"], file = file.path(".", analysisdir, "ProtDownRnaUp_1_GeneId.txt"))
write(WCPstimQuant[ProtDownregRnaUpregUbUpreg,"GeneId.x"], file = file.path(".", analysisdir, "ProtDownRnaUpUbUp_1_GeneId.txt"))
#
#protein downreg, RNAseq unch, ubiquitylated. "2"
ProtDownregRnaUnch<-intersect(UbDegRegion_WcpDec, UbDegRegion_RnaUch)
ProtDownregRnaUnch_Num<-length(ProtDownregRnaUnch)
ProtDownregRnaUnchUbUpreg<-intersect(ProtDownregRnaUnch, IncreaseUb)
ProtDownregRnaUnchUbUpreg_Num<-length(ProtDownregRnaUnchUbUpreg)
ProtDownregRnaUnchUbDownreg<-intersect(ProtDownregRnaUnch, DecreaseUb)
ProtDownregRnaUnchUbDownreg_Num<-length(ProtDownregRnaUnchUbDownreg)
ProtDownregRnaUnchUbUnch<-intersect(ProtDownregRnaUnch, UnchangeUb)
ProtDownregRnaUnchUbUnch_Num<-length(ProtDownregRnaUnchUbUnch)
print("protein downreg, RNAseq unch, ubiquitylated upreg/downreg :2")
print(ProtDownregRnaUnch_Num)
print(ProtDownregRnaUnchUbUpreg_Num)
print(ProtDownregRnaUnchUbDownreg_Num)
print(ProtDownregRnaUnchUbUnch_Num)
#####write(ProtDownregRnaUnch, file = file.path(".", analysisdir, "ProtDownRnaUnch_2_ProtId.txt"))
write(ProtDownregRnaUnchUbUpreg, file = file.path(".", analysisdir, "ProtDownRnaUnchUbUp_2_ProtId.txt"))
#####write(WCPstimQuant[ProtDownregRnaUnch,"GeneId.x"], file = file.path(".", analysisdir, "ProtDownRnaUnch_2_GeneId.txt"))
write(WCPstimQuant[ProtDownregRnaUnchUbUpreg,"GeneId.x"], file = file.path(".", analysisdir, "ProtDownRnaUnchUbUp_2_GeneId.txt"))
#
#protein downreg, RNAseq downreg, ubiquitylated. "3"
ProtDownregRnaDownreg<-intersect(UbDegRegion_WcpDec, UbDegRegion_RnaDec)
ProtDownregRnaDownreg_Num<-length(ProtDownregRnaDownreg)
ProtDownregRnaDownregUbUpreg<-intersect(ProtDownregRnaDownreg, IncreaseUb)
ProtDownregRnaDownregUbUpreg_Num<-length(ProtDownregRnaDownregUbUpreg)
ProtDownregRnaDownregUbDownreg<-intersect(ProtDownregRnaDownreg, DecreaseUb)
ProtDownregRnaDownregUbDownreg_Num<-length(ProtDownregRnaDownregUbDownreg)
ProtDownregRnaDownregUbUnch<-intersect(ProtDownregRnaDownreg, UnchangeUb)
ProtDownregRnaDownregUbUnch_Num<-length(ProtDownregRnaDownregUbUnch)
print("protein downreg, RNAseq downreg, ubiquitylated upreg/downreg :3")
print(ProtDownregRnaDownreg_Num)
print(ProtDownregRnaDownregUbUpreg_Num)
print(ProtDownregRnaDownregUbDownreg_Num)
print(ProtDownregRnaDownregUbUnch_Num)
#####write(ProtDownregRnaDownreg, file = file.path(".", analysisdir, "ProtDownRnaDown_3_ProtId.txt"))
write(ProtDownregRnaDownregUbUpreg, file = file.path(".", analysisdir, "ProtDownRnaDownUbUp_3_ProtId.txt"))
#####write(WCPstimQuant[ProtDownregRnaDownreg,"GeneId.x"], file = file.path(".", analysisdir, "ProtDownRnaDown_3_GeneId.txt"))
write(WCPstimQuant[ProtDownregRnaDownregUbUpreg,"GeneId.x"], file = file.path(".", analysisdir, "ProtDownRnaDownUbUp_3_GeneId.txt"))
#
##protein upreg, RNAseq upreg, ubiquitylated. "4"
ProtUpregRnaUpreg<-intersect(UbDegRegion_WcpInc, UbDegRegion_RnaInc)
ProtUpregRnaUpreg_Num<-length(ProtUpregRnaUpreg)
ProtUpregRnaUpregUbUpreg<-intersect(ProtUpregRnaUpreg, IncreaseUb)
ProtUpregRnaUpregUbUpreg_Num<-length(ProtUpregRnaUpregUbUpreg)
ProtUpregRnaUpregUbDownreg<-intersect(ProtUpregRnaUpreg, DecreaseUb)
ProtUpregRnaUpregUbDownreg_Num<-length(ProtUpregRnaUpregUbDownreg)
ProtUpregRnaUpregUbUnch<-intersect(ProtUpregRnaUpreg, UnchangeUb)
ProtUpregRnaUpregUbUnch_Num<-length(ProtUpregRnaUpregUbUnch)
print("protein upreg, RNAseq upreg, ubiquitylated upreg/downreg :4")
print(ProtUpregRnaUpreg_Num)
print(ProtUpregRnaUpregUbUpreg_Num)
print(ProtUpregRnaUpregUbDownreg_Num)
print(ProtUpregRnaUpregUbUnch_Num)
#####write(ProtUpregRnaUpreg, file = file.path(".", analysisdir, "ProtUpRnaUp_4_ProtId.txt"))
write(ProtUpregRnaUpregUbUpreg, file = file.path(".", analysisdir, "ProtUpRnaUpUbUp_4_ProtId.txt"))
#####write(WCPstimQuant[ProtUpregRnaUpreg,"GeneId.x"], file = file.path(".", analysisdir, "ProtUpRnaUp_4_GeneId.txt"))
write(WCPstimQuant[ProtUpregRnaUpregUbUpreg,"GeneId.x"], file = file.path(".", analysisdir, "ProtUpRnaUpUbUp_4_GeneId.txt"))
#
#protein upreg, RNAseq unch, ubiquitylated. "5"
ProtUpregRnaUnch<-intersect(UbDegRegion_WcpInc, UbDegRegion_RnaUch)
ProtUpregRnaUnch_Num<-length(ProtUpregRnaUnch)
ProtUpregRnaUnchUbUpreg<-intersect(ProtUpregRnaUnch, IncreaseUb)
ProtUpregRnaUnchUbUpreg_Num<-length(ProtUpregRnaUnchUbUpreg)
ProtUpregRnaUnchUbDownreg<-intersect(ProtUpregRnaUnch, DecreaseUb)
ProtUpregRnaUnchUbDownreg_Num<-length(ProtUpregRnaUnchUbDownreg)
ProtUpregRnaUnchUbUnch<-intersect(ProtUpregRnaUnch, UnchangeUb)
ProtUpregRnaUnchUbUnch_Num<-length(ProtUpregRnaUnchUbUnch)
print("protein upreg, RNAseq unch, ubiquitylated upreg :5")
print(ProtUpregRnaUnch_Num)
print(ProtUpregRnaUnchUbUpreg_Num)
print(ProtUpregRnaUnchUbDownreg_Num)
print(ProtUpregRnaUnchUbUnch_Num)
#####write(ProtUpregRnaUnch, file = file.path(".", analysisdir, "ProtUpRnaUnch_5_ProtId.txt"))
write(ProtUpregRnaUnchUbUpreg, file = file.path(".", analysisdir, "ProtUpRnaUnchUbUp_5_ProtId.txt"))
#####write(WCPstimQuant[ProtUpregRnaUnch,"GeneId.x"], file = file.path(".", analysisdir, "ProtUpRnaUnch_5_GeneId.txt"))
write(WCPstimQuant[ProtUpregRnaUnchUbUpreg,"GeneId.x"], file = file.path(".", analysisdir, "ProtUpRnaUnchUbUp_5_GeneId.txt"))
#
#protein upreg, RNAseq downreg, ubiquitylated. "6"
ProtUpregRnaDownreg<-intersect(UbDegRegion_WcpInc, UbDegRegion_RnaDec)
ProtUpregRnaDownreg_Num<-length(ProtUpregRnaDownreg)
ProtUpregRnaDownregUbUpreg<-intersect(ProtUpregRnaDownreg, IncreaseUb)
ProtUpregRnaDownregUbUpreg_Num<-length(ProtUpregRnaDownregUbUpreg)
ProtUpregRnaDownregUbDownreg<-intersect(ProtUpregRnaDownreg, DecreaseUb)
ProtUpregRnaDownregUbDownreg_Num<-length(ProtUpregRnaDownregUbDownreg)
ProtUpregRnaDownregUbUnch<-intersect(ProtUpregRnaDownreg, UnchangeUb)
ProtUpregRnaDownregUbUnch_Num<-length(ProtUpregRnaDownregUbUnch)
print("protein upreg, RNAseq downreg, ubiquitylated upreg :6")
print(ProtUpregRnaDownreg_Num)
print(ProtUpregRnaDownregUbUpreg_Num)
print(ProtUpregRnaDownregUbDownreg_Num)
print(ProtUpregRnaDownregUbUnch_Num)
#####write(ProtUpregRnaDownreg, file = file.path(".", analysisdir, "ProtUpRnaDown_6_ProtId.txt"))
write(ProtUpregRnaDownregUbUpreg, file = file.path(".", analysisdir, "ProtUpRnaDownUbUp_6_ProtId.txt"))
#####write(WCPstimQuant[ProtUpregRnaDownreg,"GeneId.x"], file = file.path(".", analysisdir, "ProtUpRnaDown_6_GeneId.txt"))
write(WCPstimQuant[ProtUpregRnaDownregUbUpreg,"GeneId.x"], file = file.path(".", analysisdir, "ProtUpRnaDownUbUp_6_GeneId.txt"))
#
#protein unch, RNAseq upreg, ubiquitylated. "7"
ProtUnchRnaUpreg<-intersect(UbDegRegion_WcpUch, UbDegRegion_RnaInc)
ProtUnchRnaUpreg_Num<-length(ProtUnchRnaUpreg)
ProtUnchRnaUpregUbUpreg<-intersect(ProtUnchRnaUpreg, IncreaseUb)
ProtUnchRnaUpregUbUpreg_Num<-length(ProtUnchRnaUpregUbUpreg)
ProtUnchRnaUpregUbDownreg<-intersect(ProtUnchRnaUpreg, DecreaseUb)
ProtUnchRnaUpregUbDownreg_Num<-length(ProtUnchRnaUpregUbDownreg)
ProtUnchRnaUpregUbUnch<-intersect(ProtUnchRnaUpreg, UnchangeUb)
ProtUnchRnaUpregUbUnch_Num<-length(ProtUnchRnaUpregUbUnch)
print("protein unch, RNAseq upreg, ubiquitylated upreg :7")
print(ProtUnchRnaUpreg_Num)
print(ProtUnchRnaUpregUbUpreg_Num)
print(ProtUnchRnaUpregUbDownreg_Num)
print(ProtUnchRnaUpregUbUnch_Num)
#####write(ProtUnchRnaUpreg, file = file.path(".", analysisdir, "ProtUnchRnaUp_7_ProtId.txt"))
write(ProtUnchRnaUpregUbUpreg, file = file.path(".", analysisdir, "ProtUnchRnaUpUbUp_7_ProtId.txt"))
#####write(WCPstimQuant[ProtUnchRnaUpreg,"GeneId.x"], file = file.path(".", analysisdir, "ProtUnchRnaUp_7_GeneId.txt"))
write(WCPstimQuant[ProtUnchRnaUpregUbUpreg,"GeneId.x"], file = file.path(".", analysisdir, "ProtUnchRnaUpUbUp_7_GeneId.txt"))
#
#protein unch, RNAseq unch, ubiquitylated. "8"
ProtUnchRnaUnch<-intersect(UbDegRegion_WcpUch, UbDegRegion_RnaUch)
ProtUnchRnaUnch_Num<-length(ProtUnchRnaUnch)
ProtUnchRnaUnchUbUpreg<-intersect(ProtUnchRnaUnch, IncreaseUb)
ProtUnchRnaUnchUbUpreg_Num<-length(ProtUnchRnaUnchUbUpreg)
ProtUnchRnaUnchUbDownreg<-intersect(ProtUnchRnaUnch, DecreaseUb)
ProtUnchRnaUnchUbDownreg_Num<-length(ProtUnchRnaUnchUbDownreg)
ProtUnchRnaUnchUbUnch<-intersect(ProtUnchRnaUnch, UnchangeUb)
ProtUnchRnaUnchUbUnch_Num<-length(ProtUnchRnaUnchUbUnch)
print("protein unch, RNAseq unch, ubiquitylated upreg :8")
print(ProtUnchRnaUnch_Num)
print(ProtUnchRnaUnchUbUpreg_Num)
print(ProtUnchRnaUnchUbDownreg_Num)
print(ProtUnchRnaUnchUbUnch_Num)
#####write(ProtUnchRnaUnch, file = file.path(".", analysisdir, "ProtUnchRnaUnch_8_ProtId.txt"))
write(ProtUnchRnaUnchUbUpreg, file = file.path(".", analysisdir, "ProtUnchRnaUnchUbUp_8_ProtId.txt"))
#####write(WCPstimQuant[ProtUnchRnaUnch,"GeneId.x"], file = file.path(".", analysisdir, "ProtUnchRnaUnch_8_GeneId.txt"))
write(WCPstimQuant[ProtUnchRnaUnchUbUpreg,"GeneId.x"], file = file.path(".", analysisdir, "ProtUnchRnaUnchUbUp_8_GeneId.txt"))
#
#protein unch, RNAseq downreg, ubiquitylated. "9"
ProtUnchRnaDownreg<-intersect(UbDegRegion_WcpUch, UbDegRegion_RnaDec)
ProtUnchRnaDownreg_Num<-length(ProtUnchRnaDownreg)
ProtUnchRnaDownregUbUpreg<-intersect(ProtUnchRnaDownreg, IncreaseUb)
ProtUnchRnaDownregUbUpreg_Num<-length(ProtUnchRnaDownregUbUpreg)
ProtUnchRnaDownregUbDownreg<-intersect(ProtUnchRnaDownreg, DecreaseUb)
ProtUnchRnaDownregUbDownreg_Num<-length(ProtUnchRnaDownregUbDownreg)
ProtUnchRnaDownregUbUnch<-intersect(ProtUnchRnaDownreg, UnchangeUb)
ProtUnchRnaDownregUbUnch_Num<-length(ProtUnchRnaDownregUbUnch)
print("protein unch, RNAseq downreg, ubiquitylated upreg :9")
print(ProtUnchRnaDownreg_Num)
print(ProtUnchRnaDownregUbUpreg_Num)
print(ProtUnchRnaDownregUbDownreg_Num)
print(ProtUnchRnaDownregUbUnch_Num)
#####write(ProtUnchRnaDownreg, file = file.path(".", analysisdir, "ProtUnchRnaDown_9_ProtId.txt"))
write(ProtUnchRnaDownregUbUpreg, file = file.path(".", analysisdir, "ProtUnchRnaDownUbUp_9_ProtId.txt"))
#####write(WCPstimQuant[ProtUnchRnaDownreg,"GeneId.x"], file = file.path(".", analysisdir, "ProtUnchRnaDown_9_GeneId.txt"))
write(WCPstimQuant[ProtUnchRnaDownregUbUpreg,"GeneId.x"], file = file.path(".", analysisdir, "ProtUnchRnaDownUbUp_9_GeneId.txt"))
#

##combine data to get all modeled degraded or non degraded proteins
#increased ubiquitination
ModelDeg_UbInc<-c(ProtDownregRnaUpregUbUpreg,ProtDownregRnaUnchUbUpreg,ProtUnchRnaUpregUbUpreg)
ModelDeg_UbInc_Num<-length(ModelDeg_UbInc)
print("ModelDeg_UbInc_Num")
print(ModelDeg_UbInc_Num)
ModelNonDeg_UbInc<-c(ProtDownregRnaDownregUbUpreg,ProtUpregRnaUpregUbUpreg,ProtUpregRnaUnchUbUpreg,ProtUpregRnaDownregUbUpreg,ProtUnchRnaUnchUbUpreg,ProtUnchRnaDownregUbUpreg)
ModelNonDeg_UbInc_Num<-length(ModelNonDeg_UbInc)
print("ModelNonDeg_UbInc_Num")
print(ModelNonDeg_UbInc_Num)
#decreased ubiquitination
ModelDeg_UbDec<-c(ProtUpregRnaDownregUbDownreg,ProtUpregRnaUnchUbDownreg,ProtUnchRnaDownregUbDownreg)
ModelDeg_UbDec_Num<-length(ModelDeg_UbDec)
print("ModelDeg_UbDec_Num")
print(ModelDeg_UbDec_Num)
ModelNonDeg_UbDec<-c(ProtUnchRnaUpregUbDownreg,ProtDownregRnaUpregUbDownreg,ProtDownregRnaUnchUbDownreg,ProtDownregRnaDownregUbDownreg,ProtUpregRnaUpregUbDownreg,ProtUnchRnaUnchUbDownreg)
ModelNonDeg_UbDec_Num<-length(ModelNonDeg_UbDec)
print("ModelNonDeg_UbDec_Num")
print(ModelNonDeg_UbDec_Num)
#unchanged ubiquitination
ModelDeg_UbUnch<-c(ProtDownregRnaUpregUbUnch,ProtDownregRnaUnchUbUnch,ProtUnchRnaUpregUbUnch)
ModelDeg_UbUnch_Num<-length(ModelDeg_UbUnch)
print("ModelDeg_UbUnch_Num")
print(ModelDeg_UbUnch_Num)
ModelNonDeg_UbUnch<-c(ProtDownregRnaDownregUbUnch,ProtUpregRnaUpregUbUnch,ProtUpregRnaUnchUbUnch,ProtUpregRnaDownregUbUnch,ProtUnchRnaUnchUbUnch,ProtUnchRnaDownregUbUnch)
ModelNonDeg_UbUnch_Num<-length(ModelNonDeg_UbUnch)
print("ModelNonDeg_UbUnch_Num")
print(ModelNonDeg_UbUnch_Num)


#write files
write(ModelDeg_UbInc, file = file.path(".", analysisdir, "PredictedDegradativeIncUbProteins_ProtId.txt"))
write(WCPstimQuant[ModelDeg_UbInc,"GeneId.x"], file = file.path(".", analysisdir, "PredictedDegradativeIncUbProteins_GeneId.txt"))
write(ModelNonDeg_UbInc, file = file.path(".", analysisdir, "PredictedNonDegradativeIncUbProteins_ProtId.txt"))
write(WCPstimQuant[ModelNonDeg_UbInc,"GeneId.x"], file = file.path(".", analysisdir, "PredictedNonDegradativeIncUbProteins_GeneId.txt"))

write(ModelDeg_UbUnch, file = file.path(".", analysisdir, "PredictedDegradativeUnchUbProteins_ProtId.txt"))
write(WCPstimQuant[ModelDeg_UbUnch,"GeneId.x"], file = file.path(".", analysisdir, "PredictedDegradativeUnchUbProteins_GeneId.txt"))
write(ModelNonDeg_UbUnch, file = file.path(".", analysisdir, "PredictedNonDegradativeUnchUbProteins_ProtId.txt"))
write(WCPstimQuant[ModelNonDeg_UbUnch,"GeneId.x"], file = file.path(".", analysisdir, "PredictedNonDegradativeUnchUbProteins_GeneId.txt"))

#
#pie chart of deg/non-deg for increased ub
pdf(file.path(".", analysisdir, "DegradationModel_IncreasedUb_piechart.pdf"))
pie(c(ModelDeg_UbInc_Num,ModelNonDeg_UbInc_Num), labels = c(NA,NA), col=c("grey90", "red"))
dev.off()
#pie chart of deg/non-deg for decreased ub
pdf(file.path(".", analysisdir, "DegradationModel_DecreasedUb_piechart.pdf"))
pie(c(ModelDeg_UbDec_Num,ModelNonDeg_UbDec_Num), labels = c(NA,NA), col=c("grey90", "red"))
dev.off()
#pie chart of deg/non-deg for unchanged ub
pdf(file.path(".", analysisdir, "DegradationModel_UnchangedUb_piechart.pdf"))
pie(c(ModelDeg_UbUnch_Num,ModelNonDeg_UbUnch_Num), labels = c(NA,NA), col=c("grey90", "red"))
dev.off()
###end divide protein and RNA differential expression into regions, overlay KeGG data


##correlations between wcp experiments
#
print("PO972igd_0hr IntensityH.PO1063_HLM intensity correlations")
WCPstimQuant_PO972igd0hrPO1063HLMIntH_pearson<-cor.test(WCPstimQuant$PO972igd_0hr, WCPstimQuant$IntensityH.PO1063_HLM, alternative="two.sided", method="pearson", conf.level = 0.95)
print(WCPstimQuant_PO972igd0hrPO1063HLMIntH_pearson)
WCPstimQuant_PO972igd0hrPO1063HLMIntH_spearman<-cor.test(WCPstimQuant$PO972igd_0hr, WCPstimQuant$IntensityH.PO1063_HLM, alternative="two.sided", method="spearman", conf.level = 0.95)
print(WCPstimQuant_PO972igd0hrPO1063HLMIntH_spearman)
#
print("PO972igd_0hr IntensityH.PO1075_WCPM intensity correlations")
WCPstimQuant_PO972igd0hrPO1075WCPMIntH_pearson<-cor.test(WCPstimQuant$PO972igd_0hr, WCPstimQuant$IntensityH.PO1075_WCPM, alternative="two.sided", method="pearson", conf.level = 0.95)
print(WCPstimQuant_PO972igd0hrPO1075WCPMIntH_pearson)
WCPstimQuant_PO972igd0hrPO1075WCPMIntH_spearman<-cor.test(WCPstimQuant$PO972igd_0hr, WCPstimQuant$IntensityH.PO1075_WCPM, alternative="two.sided", method="spearman", conf.level = 0.95)
print(WCPstimQuant_PO972igd0hrPO1075WCPMIntH_spearman)
#
print("IntensityH.PO1063_HLM IntensityH.PO1075_WCPM intensity correlations")
WCPstimQuant_PO1063HLMIntHPO1075WCPMIntH_pearson<-cor.test(WCPstimQuant$IntensityH.PO1063_HLM, WCPstimQuant$IntensityH.PO1075_WCPM, alternative="two.sided", method="pearson", conf.level = 0.95)
print(WCPstimQuant_PO1063HLMIntHPO1075WCPMIntH_pearson)
WCPstimQuant_PO1063HLMIntHPO1075WCPMIntH_spearman<-cor.test(WCPstimQuant$IntensityH.PO1063_HLM, WCPstimQuant$IntensityH.PO1075_WCPM, alternative="two.sided", method="spearman", conf.level = 0.95)
print(WCPstimQuant_PO1063HLMIntHPO1075WCPMIntH_spearman)
#
print("PO972igd_4hr IntensityL.PO1063_HLM intensity correlations")
WCPstimQuant_PO972igd4hrPO1063HLMIntL_pearson<-cor.test(WCPstimQuant$PO972igd_4hr, WCPstimQuant$IntensityL.PO1063_HLM, alternative="two.sided", method="pearson", conf.level = 0.95)
print(WCPstimQuant_PO972igd0hrPO1063HLMIntH_pearson)
WCPstimQuant_PO972igd4hrPO1063HLMIntL_spearman<-cor.test(WCPstimQuant$PO972igd_4hr, WCPstimQuant$IntensityL.PO1063_HLM, alternative="two.sided", method="spearman", conf.level = 0.95)
print(WCPstimQuant_PO972igd4hrPO1063HLMIntL_spearman)
#
print("PO972igd_4hr IntensityL.PO1075_WCPM intensity correlations")
WCPstimQuant_PO972igd4hrPO1075WCPMIntL_pearson<-cor.test(WCPstimQuant$PO972igd_4hr, WCPstimQuant$IntensityL.PO1075_WCPM, alternative="two.sided", method="pearson", conf.level = 0.95)
print(WCPstimQuant_PO972igd4hrPO1075WCPMIntL_pearson)
WCPstimQuant_PO972igd4hrPO1075WCPMIntL_spearman<-cor.test(WCPstimQuant$PO972igd_4hr, WCPstimQuant$IntensityL.PO1075_WCPM, alternative="two.sided", method="spearman", conf.level = 0.95)
print(WCPstimQuant_PO972igd4hrPO1075WCPMIntL_spearman)
#
print("IntensityL.PO1063_HLM IntensityL.PO1075_WCPM intensity correlations")
WCPstimQuant_PO1063HLMIntLPO1075WCPMIntL_pearson<-cor.test(WCPstimQuant$IntensityL.PO1063_HLM, WCPstimQuant$IntensityL.PO1075_WCPM, alternative="two.sided", method="pearson", conf.level = 0.95)
print(WCPstimQuant_PO1063HLMIntLPO1075WCPMIntL_pearson)
WCPstimQuant_PO1063HLMIntLPO1075WCPMIntL_spearman<-cor.test(WCPstimQuant$IntensityL.PO1063_HLM, WCPstimQuant$IntensityL.PO1075_WCPM, alternative="two.sided", method="spearman", conf.level = 0.95)
print(WCPstimQuant_PO1063HLMIntLPO1075WCPMIntL_spearman)
##

##correlations for nedd inhibitor experiment
#
print("PO1063HLMkgg_NormWcpFc vs PO1063HLPkgg_NormWcpFc fold change correlations")
WCPstimQuant_PO1063HLMHLPNormWcpFc_pearson<-cor.test(WCPstimQuant$PO1063HLMkgg_NormWcpFc, WCPstimQuant$PO1063HLPkgg_NormWcpFc, alternative="two.sided", method="pearson", conf.level = 0.95)
print(WCPstimQuant_PO1063HLMHLPNormWcpFc_pearson)
WCPstimQuant_PO1063HLMHLPNormWcpFc_spearman<-cor.test(WCPstimQuant$PO1063HLMkgg_NormWcpFc, WCPstimQuant$PO1063HLPkgg_NormWcpFc, alternative="two.sided", method="spearman", conf.level = 0.95)
print(WCPstimQuant_PO1063HLMHLPNormWcpFc_spearman)

##correlations between RNA, protein, ubiquitination abundance and fold change
#
print("RNA log2 fold change vs WCP log2 fold change correlations")
RNAseqWcpProts <- subset(WCPstimQuant, (!is.na(WCPstimQuant$Stim4hr_AvgFc) & !is.na(WCPstimQuant$log2FoldChangeStimUnstim)))
rnaprotfoldchangecortest_allexp_pearson<-cor.test(RNAseqWcpProts$Stim4hr_AvgFc, RNAseqWcpProts$log2FoldChangeStimUnstim, alternative="two.sided", method="pearson", conf.level = 0.95)
print(rnaprotfoldchangecortest_allexp_pearson)
rnaprotfoldchangecortest_allexp_spearman<-cor.test(RNAseqWcpProts$Stim4hr_AvgFc, RNAseqWcpProts$log2FoldChangeStimUnstim, alternative="two.sided", method="spearman", conf.level = 0.95)
print(rnaprotfoldchangecortest_allexp_spearman)
#
print("RNA log2 fold change vs WCP log2 fold change for differentially expressed WCP proteins correlations")
RNAseqWcpDiffExWcpProts <- subset(WCPstimQuant, (!is.na(WCPstimQuant$Stim4hr_AvgFc) & !is.na(WCPstimQuant$log2FoldChangeStimUnstim) & WCPstimQuant$Stim4hr_TtestPval<ProtSigThresh))
dim(RNAseqWcpDiffExWcpProts)
rnaprotfoldchangecortest_diffexp_pearson<-cor.test(RNAseqWcpDiffExWcpProts$Stim4hr_AvgFc, RNAseqWcpDiffExWcpProts$log2FoldChangeStimUnstim, alternative="two.sided", method="pearson", conf.level = 0.95)
print(rnaprotfoldchangecortest_diffexp_pearson)
rnaprotfoldchangecortest_diffexp_spearman<-cor.test(RNAseqWcpDiffExWcpProts$Stim4hr_AvgFc, RNAseqWcpDiffExWcpProts$log2FoldChangeStimUnstim, alternative="two.sided", method="spearman", conf.level = 0.95)
print(rnaprotfoldchangecortest_diffexp_spearman)
#
print("RNA abundance vs WCP abundance correlations")
rnaprotcortest_allexp_pearson<-cor.test(WCPstimQuant$log10baseMean, WCPstimQuant$AbundanceInt_zscore, alternative="two.sided", method="pearson", conf.level = 0.95)
print(rnaprotcortest_allexp_pearson)
rnaprotcortest_allexp_spearman<-cor.test(WCPstimQuant$log10baseMean, WCPstimQuant$AbundanceInt_zscore, alternative="two.sided", method="spearman", conf.level = 0.95)
print(rnaprotcortest_allexp_spearman)
#
print("RNA abundance vs WCP abundance correlations for Ubiquitinated proteins")
RNAseqUbProts <- WCPstimQuant[row.names(WCPstimQuant)%in%PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts,]
rnaprotcortest_Ubexp_pearson<-cor.test(RNAseqUbProts$log10baseMean, RNAseqUbProts$AbundanceInt_zscore, alternative="two.sided", method="pearson", conf.level = 0.95)
print(rnaprotcortest_Ubexp_pearson)
rnaprotcortest_Ubexp_spearman<-cor.test(RNAseqUbProts$log10baseMean, RNAseqUbProts$AbundanceInt_zscore, alternative="two.sided", method="spearman", conf.level = 0.95)
print(rnaprotcortest_Ubexp_spearman)
#
print("Ubiquitination Norm FC PO972-PO1063 correlations")
UbProtsInt <- WCPstimQuant[row.names(WCPstimQuant)%in%PO972kggPO1063HLMkgg_All0hr4hrFc_Int_UbProts,]
ubprotcortest_UbFc_pearson<-cor.test(UbProtsInt$PO972kgg_NormWcpFc, UbProtsInt$PO1063HLMkgg_NormWcpFc, alternative="two.sided", method="pearson", conf.level = 0.95)
print(ubprotcortest_UbFc_pearson)
ubprotcortest_UbFc_spearman<-cor.test(UbProtsInt$PO972kgg_NormWcpFc, UbProtsInt$PO1063HLMkgg_NormWcpFc, alternative="two.sided", method="spearman", conf.level = 0.95)
print(ubprotcortest_UbFc_spearman)
#
##end correlations between RNA, protein, ubiquitination abundance and fold change


##percentage of protein expression for RNA count abundance
#Merge filtered WCP data with RnaData dataframe, keeping all transcripts, not just the transcripts w/ protein expression
#
#
#
#
#
#
#
#
#
#
# check this merge section
RnaData <- merge(RnaData, WCPstimQuant, by="EnsemblId", all.x=TRUE)
int<-0.05
bins<-seq(0,ceiling(max(RnaData$log10baseMean.x,na.rm=TRUE)),by=int)
abdpcts<-vector()
for ( i in 1:length(bins) ){
	min<-bins[i]
	max<-min+int
	#print(min)
	TmpRnaData<-subset(RnaData, RnaData$log10baseMean.x>=min & RnaData$log10baseMean.x<max)
	#dim(TmpRnaData)
	wcpabd<-TmpRnaData$AbundanceInt_zscore
	wcpabd_num<-length(wcpabd)
	#print("wcpabd_num")
	#print(wcpabd_num)
	wcpabd_exp_num<-sum(!is.na(wcpabd))
	#print("wcpabd_exp_num")
	#print(wcpabd_exp_num)
	wcpabd_exp_pct<-wcpabd_exp_num/wcpabd_num
	#print("wcpabd_exp_pct")
	#print(wcpabd_exp_pct)
	abdpcts<-c(abdpcts,wcpabd_exp_pct)
}
#
#
#
#
#
#
#
# check this merge section

###write output file for complete data analysis data frame
write.table(WCPstimQuant, file=file.path(".", analysisdir, "WholeCellProteomeTCRstim_PO972-PO1075-PO1063_iBAQ-SILAC_FilteredByCount_analyzeddata.txt"), sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)
#head(WCPstimQuant)


##establish a data frame of only ubiquitinated proteins
UbProts<-c(PO972kggPO1063HLMkgg_0hr_Unique_UbProts, PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd, PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanSd, PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd, PO972kggPO1063HLMkgg_4hr_Unique_UbProts)
UbProtsData<-subset(WCPstimQuant, row.names(WCPstimQuant)%in%UbProts)
#dim(UbProtsData)
#head(UbProtsData)
##set normalized fold change of unique ubiquitinated proteins to min/max fold change
for( i in 1:length(row.names(UbProtsData)) ){
	protID<-row.names(UbProtsData)[i]
	if( !is.na(UbProtsData[protID,"Ub_0hrUnique"]) ){
		UbProtsData[protID,"UbFcAvg_4hr_NormWcpFc"]=min(UbProtsData$UbFcAvg_4hr_NormWcpFc,na.rm=TRUE)
	}
	if( !is.na(UbProtsData[protID,"Ub_4hrUnique"]) ){
		UbProtsData[protID,"UbFcAvg_4hr_NormWcpFc"]=max(UbProtsData$UbFcAvg_4hr_NormWcpFc,na.rm=TRUE)
	}
}
#dim(UbProtsData)

##analyze TCR pathway protein ubiquitination
#read TCR pathway protein list
TCRprots <- scan("AllTcrProts_PsmProtsRemovedMHCProtsRemoved.txt", what = "character")
print("num TCR pathway proteins compiled")
length(TCRprots)
#subset TCRprotData data frame of TCR pathway proteins
TCRprotsData<-subset(WCPstimQuant, WCPstimQuant$GeneId.x%in%TCRprots)
dim(TCRprotsData)
#subset TCRprotsUbData data frame of ubiquitinated TCR pathway proteins
#####TCRprotsUbData<-subset(TCRprotsData, !is.na(TCRprotsData$UbFcAvg_4hr_NormWcpFc))
TCRprotsUbData<-subset(UbProtsData, UbProtsData$GeneId.x%in%TCRprots)
#order TCRprotsUbData
TCRprotsUbData <- TCRprotsUbData[order(TCRprotsUbData$UbFcAvg_4hr_NormWcpFc, decreasing=TRUE),]
print("dim TCRprotsUbData")
dim(TCRprotsUbData)
#increased Ub TCR prots
TCRprots_IncUb <- row.names(subset(TCRprotsUbData, TCRprotsUbData$UbFcAvg_4hr_NormWcpFc>0))
TCRprots_IncUb_Num <- length(TCRprots_IncUb)
print("TCRprots_IncUb_Num")
print(TCRprots_IncUb_Num)
#decreased Ub TCR prots
TCRprots_DecUb <- row.names(subset(TCRprotsUbData, TCRprotsUbData$UbFcAvg_4hr_NormWcpFc<0))
TCRprots_DecUb_Num <- length(TCRprots_DecUb)
print("TCRprots_DecUb_Num")
print(TCRprots_DecUb_Num)
#increased Ub 25%
TCRprots_SigIncUb <- row.names(subset(TCRprotsUbData, TCRprotsUbData$UbFcAvg_4hr_NormWcpFc > 0.584963))
TCRprots_SigIncUb_Num <- length(TCRprots_SigIncUb)
print("TCRprots_SigIncUb_Num")
print(TCRprots_SigIncUb_Num)
#decreased Ub 24%
TCRprots_SigDecUb <- row.names(subset(TCRprotsUbData, TCRprotsUbData$UbFcAvg_4hr_NormWcpFc < -0.584963))
TCRprots_SigDecUb_Num <- length(TCRprots_SigDecUb)
print("TCRprots_SigDecUb_Num")
print(TCRprots_SigDecUb_Num)
#Ub increased, prot change
TcrUbInc_ProtInc <- intersect(TCRprots_IncUb,PO972PO1063PO1075_Wcp_Upreg_AvgFcPval)
TcrUbInc_ProtInc_Num <- length(TcrUbInc_ProtInc)
print("TcrUbInc_ProtInc_Num")
print(TcrUbInc_ProtInc_Num)
print(TcrUbInc_ProtInc)
TcrUbInc_ProtDec <- intersect(TCRprots_IncUb,PO972PO1063PO1075_Wcp_Downreg_AvgFcPval)
TcrUbInc_ProtDec_Num <- length(TcrUbInc_ProtDec)
print("TcrUbInc_ProtDec_Num")
print(TcrUbInc_ProtDec_Num)
TcrUbInc_ProtUnch <- intersect(TCRprots_IncUb,PO972PO1063PO1075_Wcp_Unch_AvgFcPval)
TcrUbInc_ProtUnch_Num <- length(TcrUbInc_ProtUnch)
print("TcrUbInc_ProtUnch_Num")
print(TcrUbInc_ProtUnch_Num)
#Ub decreased, prot change
TcrUbDec_ProtInc <- intersect(TCRprots_DecUb,PO972PO1063PO1075_Wcp_Upreg_AvgFcPval)
TcrUbDec_ProtInc_Num <- length(TcrUbDec_ProtInc)
print("TcrUbDec_ProtInc_Num")
print(TcrUbDec_ProtInc_Num)
TcrUbDec_ProtDec <- intersect(TCRprots_DecUb,PO972PO1063PO1075_Wcp_Downreg_AvgFcPval)
TcrUbDec_ProtDec_Num <- length(TcrUbDec_ProtDec)
print("TcrUbDec_ProtDec_Num")
print(TcrUbDec_ProtDec_Num)
TcrUbDec_ProtUnch <- intersect(TCRprots_DecUb,PO972PO1063PO1075_Wcp_Unch_AvgFcPval)
TcrUbDec_ProtUnch_Num <- length(TcrUbDec_ProtUnch)
print("TcrUbDec_ProtUnch_Num")
print(TcrUbDec_ProtUnch_Num)

TcrUbInc_RnaInc <- intersect(TCRprots_IncUb,RNAseq_Upreg_AvgFcPval)
TcrUbInc_RnaInc_Num <- length(TcrUbInc_RnaInc)
print("TcrUbInc_RnaInc_Num")
print(TcrUbInc_RnaInc_Num)
print(TcrUbInc_RnaInc)
TcrUbInc_RnaDec <- intersect(TCRprots_IncUb,RNAseq_Downreg_AvgFcPval)
TcrUbInc_RnaDec_Num <- length(TcrUbInc_RnaDec)
print("TcrUbInc_RnaDec_Num")
print(TcrUbInc_RnaDec_Num)
print(TcrUbInc_RnaDec)
TcrUbInc_RnaUnch <- intersect(TCRprots_IncUb,RNAseq_Unch_AvgFcPval)
TcrUbInc_RnaUnch_Num <- length(TcrUbInc_RnaUnch)
print("TcrUbInc_RnaUnch_Num")
print(TcrUbInc_RnaUnch_Num)
print(TcrUbInc_RnaUnch)

TcrUbDec_RnaInc <- intersect(TCRprots_DecUb,RNAseq_Upreg_AvgFcPval)
TcrUbDec_RnaInc_Num <- length(TcrUbDec_RnaInc)
print("TcrUbDec_RnaInc_Num")
print(TcrUbDec_RnaInc_Num)
TcrUbDec_RnaDec <- intersect(TCRprots_DecUb,RNAseq_Downreg_AvgFcPval)
TcrUbDec_RnaDec_Num <- length(TcrUbDec_RnaDec)
print("TcrUbDec_RnaDec_Num")
print(TcrUbDec_RnaDec_Num)
TcrUbDec_RnaUnch <- intersect(TCRprots_DecUb,RNAseq_Unch_AvgFcPval)
TcrUbDec_RnaUnch_Num <- length(TcrUbDec_RnaUnch)
print("TcrUbDec_RnaUnch_Num")
print(TcrUbDec_RnaUnch_Num)


###establish a data frame of only ubiquitinated proteins
#UbProts<-c(PO972kggPO1063HLMkgg_0hr_Unique_UbProts, PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd, PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanSd, PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd, PO972kggPO1063HLMkgg_4hr_Unique_UbProts)
#UbProtsData<-subset(WCPstimQuant, row.names(WCPstimQuant)%in%UbProts)
##dim(UbProtsData)
##head(UbProtsData)
###set normalized fold change of unique ubiquitinated proteins to min/max fold change
#for( i in 1:length(row.names(UbProtsData)) ){
#	protID<-row.names(UbProtsData)[i]
#	if( !is.na(UbProtsData[protID,"Ub_0hrUnique"]) ){
#		UbProtsData[protID,"UbFcAvg_4hr_NormWcpFc"]=min(UbProtsData$UbFcAvg_4hr_NormWcpFc,na.rm=TRUE)
#	}
#	if( !is.na(UbProtsData[protID,"Ub_4hrUnique"]) ){
#		UbProtsData[protID,"UbFcAvg_4hr_NormWcpFc"]=max(UbProtsData$UbFcAvg_4hr_NormWcpFc,na.rm=TRUE)
#	}
#}
##dim(UbProtsData)


###write output file for TCR proteins data frame
write.table(TCRprotsUbData, file=file.path(".", analysisdir, "WholeCellProteomeTCRstim_PO972-PO1075-PO1063_iBAQ-SILAC_FilteredByCount_analyzeddata_TCRProtsUb.txt"), sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)


#### plotting ####

#plot UpSet R plot of unfiltered identified proteins for each dataset
library(UpSetR)
listInput <- list("exp 1" = PO972ProtId_Unfiltered, 
                  "exp 2" = PO1063ProtId_Unfiltered, 
				  "exp 3" = PO1075ProtId_Unfiltered)
sapply(listInput, length)
upsetfile<-file.path(".", analysisdir, "PO972PO1063PO1075_UpsetR_UnfilteredProteinOverlap.pdf")
pdf(upsetfile, width=7, height=5)
upset(fromList(listInput), sets=c("exp 3","exp 2","exp 1"), keep.order=TRUE, 
      mainbar.y.label = "Num Common Proteins", sets.x.label = "Num Id Proteins",
	  mb.ratio = c(0.75,0.25),
	  text.scale = c(2, 2, 2, 1.75, 2, 2),
	  point.size = 2.5,
	  line.size = 1.25, main.bar.color = c("black", "black", "black", "black", "grey70", "grey70", "grey70"), sets.bar.color = c("royalblue","dodgerblue","skyblue"),
	  order.by = c("degree"), empty.intersections = "on")
dev.off()
#textscale c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
#end plot UpSet R plot of filtered identified proteins for each dataset


##venn diagram of number of unfiltered protein quantifications for 0hr, 4hr, 0hr4hr for PO972
plotpairwisevenn(PO972_0hr_All_ProtId_Unfiltered_Num,
                 PO972_4hr_All_ProtId_Unfiltered_Num,
				 PO972_0hr4hr_ProtId_Unfiltered_Num,
				 "", "", "skyblue", "royalblue", "skyblue", "royalblue", TRUE,
				 file.path(".", analysisdir, "PO972_NumQuantifiedProteins0hr4hr_Unfiltered.pdf"))
##end venn diagram of number of unfiltered protein quantifications for 0hr, 4hr, 0hr4hr for PO972


##venn diagram of number of unfiltered protein quantifications for 0hr, 4hr, 0hr4hr for PO1063
plotpairwisevenn(PO1063_0hr_All_ProtId_Unfiltered_Num,
                 PO1063_4hr_All_ProtId_Unfiltered_Num,
				 PO1063_0hr4hr_ProtId_Unfiltered_Num,
				 "", "", "skyblue", "royalblue", "skyblue", "royalblue", TRUE,
				 file.path(".", analysisdir, "PO1063_NumQuantifiedProteins0hr4hr_Unfiltered.pdf"))
##end venn diagram of number of unfiltered protein quantifications for 0hr, 4hr, 0hr4hr for PO972


##venn diagram of number of unfiltered protein quantifications for 0hr, 4hr, 0hr4hr for PO1063
plotpairwisevenn(PO1075_0hr_All_ProtId_Unfiltered_Num,
                 PO1075_4hr_All_ProtId_Unfiltered_Num,
				 PO1075_0hr4hr_ProtId_Unfiltered_Num,
				 "", "", "skyblue", "royalblue", "skyblue", "royalblue", TRUE,
				 file.path(".", analysisdir, "PO1075_NumQuantifiedProteins0hr4hr_Unfiltered.pdf"))
##end venn diagram of number of unfiltered protein quantifications for 0hr, 4hr, 0hr4hr for PO972


##venn diagram of number of unfiltered protein quantifications for 0hr, 4hr, 0hr4hr common
plotpairwisevenn(PO972PO1063PO1075_0hr_Union_ProtId_Unfiltered_Num,
                 PO972PO1063PO1075_4hr_Union_ProtId_Unfiltered_Num,
				 PO972PO1063PO1075_0hr4hr_Int_ProtId_Unfiltered_Num,
				 "", "", "blue", "red", "dodgerblue", "salmon", TRUE,
				 file.path(".", analysisdir, "PO972-PO1063-PO1075_NumQuantifiedProteins0hr4hr_Unfiltered.pdf"))
##end venn diagram of number of unfiltered protein quantifications for 0hr, 4hr, 0hr4hr common


##venn diagram of number of quantified proteins in RNA and proteomics and intersection
plotpairwisevenn(PO972PO1063PO1075_AllId_Union_ProtId_Num,
				 RnaIdentifiedGenes_InTranscriptome_Num,
				 ProtRnaIdGenes_ProteomeTranscriptome_Int_Num,
				 "", "", "steelblue4","gray40", "steelblue2","gray90", TRUE,
				 file.path(".", analysisdir, "ProteomicsRNAseq_ProteomeTranscriptome_VennDiagram.pdf"))
##end venn diagram of number of quantified proteins in RNA and proteomics and intersection


##venn diagram of number of identified ubiquitinated proteins at 0hr, 4hr and 0hr4hr common
plotpairwisevenn(PO972kggPO1063HLMkgg_0hr_All_UbProts_Num,
                 PO972kggPO1063HLMkgg_4hr_All_UbProts_Num,
				 PO972kggPO1063HLMkgg_0hr4hr_Int_UbProts_Num,
				 "", "", "blue", "red", "dodgerblue", "salmon", TRUE,
				 file.path(".", analysisdir, "PO972kgg-PO1063HLMkgg_NumQuantifiedUbProteins0hr4hr.pdf"))
##end venn diagram of number of identified ubiquitinated proteins at 0hr, 4hr and 0hr4hr common


##venn diagram of intersection between KeGG proteins in PO1063 and PO972
plotpairwisevenn(PO1063HLMkgg_AllUbId_Num,
                 PO972kgg_AllUbId_Num,
				 PO972kggPO1063HLMkgg_AllUbId_Int_Num,
				 "", "", "steelblue", "steelblue1", "steelblue", "steelblue1", TRUE,
				 file.path(".", analysisdir, "PO972kgg-PO1063HLMkgg_NumKeGGUbProteinsId_VennDiagram.pdf"))
##end venn diagram of intersection between KeGG proteins in PO1063 and PO972


##venn diagram of intersection of all identified KeGG and Wcp proteins
plotpairwisevenn(PO972PO1063PO1075_AllId_Union_ProtId_Num,
                 KeGG_AllUbProtId_Unfiltered_Num,
				 KeGGProtWcpProtInt_Num,
				 "", "", "grey60", "royalblue1", "grey60", "royalblue1", TRUE,
				 file.path(".", analysisdir, "UbProtWholeCellProt_NumProteinsId_VennDiagram.pdf"))
##end venn diagram of intersection of all identified KeGG and Wcp proteins


if ( FALSE )
{
#plot distributions of intensity values for 0hr or 4hr unique proteins and for all proteins to compare
plotnormdistribution(WCPstimQuant[PO972_4hr_Unique_ProtId,"PO972igd_4hr"], "NA", "NA", 1, "NA", "NA", "NA", "NA", "log2 norm intensity", "norm freq", file.path(".", analysisdir, "PO972igd_4hrUnique_Distribution.pdf"))
plotnormdistribution(WCPstimQuant[PO1063_4hr_Unique_ProtId,"IntensityL.PO1063_HLM"], "NA", "NA", 1, "NA", "NA", "NA", "NA", "log2 norm intensity", "norm freq", file.path(".", analysisdir, "PO1063HLM_4hrUnique_Distribution.pdf"))
plotnormdistribution(WCPstimQuant[PO1075_4hr_Unique_ProtId,"IntensityL.PO1075_WCPM"], "NA", "NA", 1, "NA", "NA", "NA", "NA", "log2 norm intensity", "norm freq", file.path(".", analysisdir, "PO1075WCPM_4hrUnique_Distribution.pdf"))
plotnormdistribution(WCPstimQuant[,"PO972igd_4hr"], "NA", "NA", 1, "NA", "NA", "NA", "NA", "log2 norm intensity", "norm freq", file.path(".", analysisdir, "PO972igd_4hr_AllProts_Distribution.pdf"))
plotnormdistribution(WCPstimQuant[,"IntensityL.PO1063_HLM"], "NA", "NA", 1, "NA", "NA", "NA", "NA", "log2 norm intensity", "norm freq", file.path(".", analysisdir, "PO1063HLM_4hr_AllProts_Distribution.pdf"))
plotnormdistribution(WCPstimQuant[,"IntensityL.PO1075_WCPM"], "NA", "NA", 1, "NA", "NA", "NA", "NA", "log2 norm intensity", "norm freq", file.path(".", analysisdir, "PO1075WCPM_4hr_AllProts_Distribution.pdf"))
plotnormdistribution(WCPstimQuant[PO972_0hr_Unique_ProtId,"PO972igd_0hr"], "NA", "NA", 1, "NA", "NA", "NA", "NA", "log2 norm intensity", "norm freq", file.path(".", analysisdir, "PO972igd_0hrUnique_Distribution.pdf"))
plotnormdistribution(WCPstimQuant[PO1063_0hr_Unique_ProtId,"IntensityH.PO1063_HLM"], "NA", "NA", 1, "NA", "NA", "NA", "NA", "log2 norm intensity", "norm freq", file.path(".", analysisdir, "PO1063HLM_0hrUnique_Distribution.pdf"))
plotnormdistribution(WCPstimQuant[PO1075_0hr_Unique_ProtId,"IntensityH.PO1075_WCPM"], "NA", "NA", 1, "NA", "NA", "NA", "NA", "log2 norm intensity", "norm freq", file.path(".", analysisdir, "PO1075WCPM_0hrUnique_Distribution.pdf"))
plotnormdistribution(WCPstimQuant[,"PO972igd_0hr"], "NA", "NA", 1, "NA", "NA", "NA", "NA", "log2 norm intensity", "norm freq", file.path(".", analysisdir, "PO972igd_0hr_AllProts_Distribution.pdf"))
plotnormdistribution(WCPstimQuant[,"IntensityH.PO1063_HLM"], "NA", "NA", 1, "NA", "NA", "NA", "NA", "log2 norm intensity", "norm freq", file.path(".", analysisdir, "PO1063HLM_0hr_AllProts_Distribution.pdf"))
plotnormdistribution(WCPstimQuant[,"IntensityH.PO1075_WCPM"], "NA", "NA", 1, "NA", "NA", "NA", "NA", "log2 norm intensity", "norm freq", file.path(".", analysisdir, "PO1075WCPM_ohr_AllProts_Distribution.pdf"))
}

##plot MSMS counts at 0 and 4 hours for TCR activation markers
#must use unfiltered data frame as these proteins are generally removed after filtering
TCRStimMarkers<-c("Cd40lg", "Egr1", "Il4", "Lta", "Icos", "Myc", "Relb", "Ifng", "Jun", "Tnfrsf9", "Egr2", "Cd69", "Tnfsf8", "Tnfrsf4", "Junb", "Irf8")
TCRStimMarkersNames<-c("CD40L", "EGR1", "IL4", "TNFbeta", "ICOS", "c-MYC", "RELB", "IFNgamma", "AP1", "CD137", "EGR2", "CD69", "CD30L", "OX40", "JUNB", "IRF8")
MSMS0hr<-vector()
MSMS4hr<-vector()
for( i in 1:length(TCRStimMarkers) ){
	currtcrstimgene<-TCRStimMarkers[i]
	MScount0hr<-WCPstimQuantUnfiltered[which(WCPstimQuantUnfiltered$GeneId.x==currtcrstimgene),"MSMScounts.PO972igd_0hr"]
	MScount4hr<-WCPstimQuantUnfiltered[which(WCPstimQuantUnfiltered$GeneId.x==currtcrstimgene),"MSMScounts.PO972igd_4hr"]
	MSMS0hr<-c(MSMS0hr,MScount0hr)
	MSMS4hr<-c(MSMS4hr,MScount4hr)
}
ymax<-30
TCRStimMarkersData <- rbind(t(MSMS0hr), t(MSMS4hr))
#print(TCRStimMarkersData)
pdf(file.path(".", analysisdir, "TCRstimMarkers_WcpMSMScounts_barplot.pdf"))
par(mar=c(8,5,4,2)+0.1)
barplot(TCRStimMarkersData, beside=T, space=c(0,0.5), ylab="", names=TCRStimMarkersNames, cex.names=1.5, las=2, cex.axis=1.5, width=c(5,5), ylim=c(0,ymax), col=c("grey70","red"))
box(bty="l", lwd=3)
mtext("MS/MS counts", side=2, line=2.75, cex=1.75)
dev.off()
##end plot MSMS counts at 0 and 4 hours for TCR activation markers


##Ub fold change number bar graph
UbNums<-c(PO972kggPO1063HLMkgg_0hr_Unique_UbProts_Num, PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanThreshold_Num, PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanThreshold_Num, PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanThreshold_Num, PO972kggPO1063HLMkgg_4hr_Unique_UbProts_Num)
pdf(file.path(".", analysisdir, "UbProts_NumUpDownReg_barchart.pdf"), width=4, height=4)
barplot(UbNums, beside=TRUE, names.arg=c("unique", "increased", "unchanged", "increased", "unique"), xlab="", ylab="", ylim=c(0,500), width=0.5, space=0.1, cex.axis=1, cex.names=1, las=1, col=c("steelblue4", "steelblue1", "grey90", "salmon", "red3"))
box(bty="l", lwd=3)
#width=c(2)
#space=c(1,10)
mtext("proteins identified", side=2, line=2.75, cex=1)
dev.off()
##end Ub fold change number bar graph


##RNAseq volcano plot
#CD4markers <- c("Cd4","Cd69","Cd40lg","Cd28","Cd81","Cd96","Cd27","Cd44","Cd5")
#alltcr <- c("Cd247","Ptprc","Cd3g","Lat","Lck","Itk","Grb2","Sos1","Zap70","Prkcq","Ikbkg","Nfkb1","Cblb","Cd3e","Rhoa","Nras","Kras","Ppp3r1")
#labelgenes <- c("Cul1", "Cul2", "Cul3", "Cul4a", "Cul4b", "Cul5")
#labelgenes <- c("Cul4b", "Ddb1", "Rbx1", "Crbn", "Ddb2", "Dcaf5", "Dcaf7", "Dcaf8", "Dcaf11", "Dcaf13")
labelgenes <- TCRStimMarkers
pdf(file.path(".", analysisdir, "RNAseq_0hr4hrStimFoldChange_VolcanoPlot.pdf"))
par(bty='l')
#xmin <- floor(min(WCPstimQuantUnfiltered$log2FoldChangeStimUnstim, na.rm=TRUE))
#xmax <- ceiling(max(WCPstimQuantUnfiltered$log2FoldChangeStimUnstim, na.rm=TRUE))
ymin <- 0
ymax <- ceiling(max(WCPstimQuantUnfiltered$neglog10padj, na.rm=TRUE))
xmin <- -5
xmax <- 8
#ymin <- 0
#ymax <- 120
#make an axes
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
#add horizontal line
par(new=T)
##abline(h=2, col="gray50", lwd=2, lty=2)
#add points
par(new=T)
plot(WCPstimQuantUnfiltered$log2FoldChangeStimUnstim,
     WCPstimQuantUnfiltered$neglog10padj,
     pch=19, cex=0.75,
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 col="grey70",
	 xlab="",ylab="", xaxt='n', yaxt='n'
)
#add points for genes in vector
for ( i in 1:length(labelgenes) )
{
	par(new=T)
	labelgene<-labelgenes[i]
	x<-WCPstimQuantUnfiltered[which(WCPstimQuantUnfiltered$GeneId.x==labelgene),"log2FoldChangeStimUnstim"]
	y<-WCPstimQuantUnfiltered[which(WCPstimQuantUnfiltered$GeneId.x==labelgene),"neglog10padj"]
	points(WCPstimQuantUnfiltered[which(WCPstimQuantUnfiltered$GeneId.x==labelgene),"log2FoldChangeStimUnstim"],
	       WCPstimQuantUnfiltered[which(WCPstimQuantUnfiltered$GeneId.x==labelgene),"neglog10padj"],
		   pch=19,
		   cex=2,
		   col="red"
          )
	#text(WCPstimQuantUnfiltered[which(WCPstimQuantUnfiltered$GeneId.x==labelgene),"log2FoldChangeStimUnstim"],
	 #    WCPstimQuantUnfiltered[which(WCPstimQuantUnfiltered$GeneId.x==labelgene),"neglog10padj"],
	#	 labels=labelgene,
	#	 cex=0.75,
	#	 pos=2,
	#	 col="red",
	#	 font=2
	#	)
}
#label axes
axis(1,at=c(seq(xmin,xmax,by=0.5)),labels=F,col="black",cex.axis=1,tck=-0.01)
axis(1,at=c(seq(xmin,xmax,by=1)),col="black",lwd=2,cex.axis=1.2)
mtext("log2 fold change (restim/rest)",1,line=2.75,cex=1.5)
axis(2,at=c(seq(ymin,ymax,by=10)),labels=F,col="black",cex.axis=1,tck=-0.01)
axis(2,at=c(seq(ymin,ymax,by=20)),col="black",lwd=2,cex.axis=1.2)
mtext("-log10 p-value",2,line=2.75,cex=1.5)
dev.off()
##end RNAseq volcano plot


##Wcp volcano plot
pdf(file.path(".", analysisdir, "Wcp_0hr4hrStimFoldChange_VolcanoPlot.pdf"))
par(bty='l')
#xmin <- floor(min(WCPstimQuant$Stim4hr_AvgFc, na.rm=TRUE))
#xmax <- ceiling(max(WCPstimQuant$Stim4hr_AvgFc, na.rm=TRUE))
#ymin <- 0
#ymax <- ceiling(max(WCPstimQuant$Stim4hr_NegLog10Pval, na.rm=TRUE))
xmin <- -3
xmax <- 5
ymin <- 0
ymax <- 3.5
#make an axes
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
#add horizontal line
par(new=T)
#abline(h=1, col="gray50", lwd=2, lty=2)
#add points
par(new=T)
plot(WCPstimQuant$Stim4hr_AvgFc,
     WCPstimQuant$Stim4hr_NegLog10Pval,
     pch=19, cex=1,
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 #col="grey50",
	 #col=ifelse((WCPstimQuant$Stim4hr_AvgFc>=0.321928 & WCPstimQuant$Stim4hr_NegLog10Pval>1.3), "red3",
	 #    ifelse((WCPstimQuant$Stim4hr_AvgFc>0 & WCPstimQuant$Stim4hr_AvgFc<0.321928 & WCPstimQuant$Stim4hr_NegLog10Pval>1.3), "lightsalmon",
	 #    ifelse((WCPstimQuant$Stim4hr_AvgFc<=-0.321928 & WCPstimQuant$Stim4hr_NegLog10Pval>1.3), "blue", 
	 #	  ifelse((WCPstimQuant$Stim4hr_AvgFc<0 & WCPstimQuant$Stim4hr_AvgFc>-0.321928 & WCPstimQuant$Stim4hr_NegLog10Pval>1.3), "lightblue", "grey50")))),
	 col=ifelse((WCPstimQuant$Stim4hr_AvgFc>=0 & WCPstimQuant$Stim4hr_NegLog10Pval>1.3), "red", ifelse((WCPstimQuant$Stim4hr_AvgFc<0 & WCPstimQuant$Stim4hr_NegLog10Pval>1.3), "blue", "grey50")),
	 xlab="", ylab="", xaxt='n', yaxt='n'
)
box(bty="l", lwd=3)
#add points for genes in vector
#for ( i in 1:length(labelgenes) )
#{
#	par(new=T)
#	labelgene<-labelgenes[i]
#	print(labelgene)
#	x<-WCPstimQuant[which(WCPstimQuant$GeneId.x==labelgene),"Stim4hr_AvgFc"]
#	y<-WCPstimQuant[which(WCPstimQuant$GeneId.x==labelgene),"Stim4hr_NegLog10Pval"]
#	print(x)
#	print(y)
#	points(WCPstimQuant[which(WCPstimQuant$GeneId.x==labelgene),"Stim4hr_AvgFc"],
#	       WCPstimQuant[which(WCPstimQuant$GeneId.x==labelgene),"Stim4hr_NegLog10Pval"],
#		   pch=19,
#		   cex=2,
#		   col="red"
#         )
#	#text(WCPstimQuant[which(WCPstimQuant$GeneId.x==labelgene),"Stim4hr_AvgFc"],
#	 #    WCPstimQuant[which(WCPstimQuant$GeneId.x==labelgene),"Stim4hr_NegLog10Pval"],
#	#	 labels=labelgene,
#	#	 cex=0.75,
#	#	 pos=2,
#	#	 col="red",
#	#	 font=2
#	#	)
#}
#label axes
axis(1,at=c(seq(xmin,xmax,by=0.5)),labels=F,col="black",cex.axis=1.25,tck=-0.01)
axis(1,at=c(seq(xmin,xmax,by=1)),col="black", lwd=2, cex.axis=1.5)
mtext("log2 fold change (restim/rest)",1,line=2.75,cex=1.75)
axis(2,at=c(seq(ymin,ymax,by=0.25)),labels=F,col="black",cex.axis=1.25,tck=-0.01)
axis(2,at=c(seq(ymin,ymax,by=0.5)),col="black", lwd=2, cex.axis=1.5)
mtext("-log10 p-value",2,line=2.75,cex=1.75)
dev.off()
##end Wcp volcano plot


##plot proteomics results vs RNAseq results
pdf(file.path(".", analysisdir, "ProteomicsRNAseq_Log2FoldChangeComp.pdf"))
par(bty="l")
xmin <- -3
xmax <- 5
ymin <- -5
ymax <- 7
#xmin <- floor(min(WCPstimQuant$Stim4hr_AvgFc, na.rm=TRUE))
#xmax <- ceiling(max(WCPstimQuant$Stim4hr_AvgFc, na.rm=TRUE))
#ymin <- floor(min(WCPstimQuant$log2FoldChangeStimUnstim, na.rm=TRUE))
#ymax <- ceiling(max(WCPstimQuant$log2FoldChangeStimUnstim, na.rm=TRUE))
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
par(new=T)
abline(v=0, col="gray70", lwd=2, lty=1)
par(new=T)
abline(h=0, col="gray70", lwd=2, lty=1)
par(new=T)
plot(WCPstimQuant$Stim4hr_AvgFc, WCPstimQuant$log2FoldChangeStimUnstim,
     pch=19,
	 cex=0.75,
	 #cex=4*(abs(WCPstimQuant$UbFcAvg_4hr_NormWcpFc)),
	 col="grey70",
	 #col=ifelse(WCPstimQuant$Stim4hr_TtestPval<ProtSigThresh, "red", "grey70"),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
)
diffexprots<-row.names(subset(WCPstimQuant, WCPstimQuant$Stim4hr_TtestPval<ProtSigThresh))
for ( i in 1:length(diffexprots) )
{
	par(new=T)
	diffprot<-diffexprots[i]
	points(WCPstimQuant[diffprot,"Stim4hr_AvgFc"],
	       WCPstimQuant[diffprot,"log2FoldChangeStimUnstim"],
		   pch=19,
		   cex=0.75,
		   col="red"
		  )
}
axis(1,c(seq(xmin,xmax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=1)),col="black",lwd=2,cex.axis=1.5)
mtext("WCP log2 fold change",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)),col="black",lwd=2,cex.axis=1.5)
mtext("RNAseq log2 fold change",2,line=2.75,cex=1.5)
dev.off()


if ( FALSE )
{
##plot heatmap of Cullin 0-4hr fold changes
cullins <- c("Cul1", "Cul2", "Cul3", "Cul4a", "Cul4b", "Cul5")
#cullins <- c("Cul4b", "Ddb1", "Rbx1", "Crbn", "Ddb2", "Dcaf5", "Dcaf7", "Dcaf8", "Dcaf11", "Dcaf13")
cullinsId <- vector()
for ( i in 1:length(cullins) ){
	appendgene <- cullins[i]
	cullinsId <- c(cullinsId, which(WCPstimQuant$GeneId.x == appendgene))
}
CullinData<-WCPstimQuant[cullinsId,]
CullinFoldChangeMat <- as.matrix(CullinData[,c("IGD0hr_zscore", "HLMintH_zscore", "WCPMintH_zscore", "IGD4hr_zscore", "HLMintL_zscore", "WCPMintL_zscore")])
#CullinFoldChangeMat <- as.matrix(CullinData[,c("IGDStim4hr_fc", "HLMStim4hr_fc","WCPMStim4hr_fc")])
#CullinFoldChangeMat[CullinFoldChangeMat=="Inf"] <- NA
#CullinFoldChangeMat <- as.matrix(CullinData[,c("IGDStim1hr_fc","Stim4hr_AvgFc")])
#CullinFoldChangeMat <- as.matrix(CullinData[,c("UbFcAvg_1hr_NormWcpFc","UbFcAvg_4hr_NormWcpFc","UbFcAvg_1hr4hr_NormWcpFc")])
#CullinFoldChangeMat <- as.matrix(CullinData[,c("PO972kgg_NormWcpFc","PO1063HLMkgg_NormWcpFc")])
CullinFoldChangeMat<-apply(CullinFoldChangeMat,2,as.numeric)
pdf(file.path(".", analysisdir, "TCRproteins_WcpZscore_heatmap.pdf"))
col_palette <- colorRampPalette(c("blue", "grey90", "red"))(n = 98)
#col_breaks = c(seq(-0.75,-0.01,length=33), # for red
#               seq(-0.009,0.009,length=33),  # for yellow
#               seq(0.01,0.75,length=33)) # for green
heatmap.2(CullinFoldChangeMat,
          notecol="black",      # change font color of cell labels to black
		  density.info="none",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  dendrogram="row",    # only draw a row dendrogram
		  #Rowv="NA",          # turn off row clustering
		  Colv="NA",            # turn off column clustering
		  #col=bluered(100),     #col_palette or col=bluered(100),
		  col=col_palette,
		  #breaks=col_breaks,
		  na.color = "white",
		  scale="none",
		  labRow=CullinData$GeneId.x,
		  labCol=c("exp1 0hr (LFQ)", "exp2 0hr (SILAC)", "exp3 0hr (SILAC)", "exp1 4hr (LFQ)", "exp2 4hr (SILAC)", "exp3 4hr (SILAC)"),
		  margins=c(22,24),
		  cexRow=0.8,
		  cexCol=0.8,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1.1,
		  sepwidth=c(0.005, 0.005),  # width of the borders
		  sepcolor='grey80',
		  colsep=0:ncol(CullinFoldChangeMat),
		  rowsep=0:nrow(CullinFoldChangeMat)
		 )
dev.off()
##end Cullin heatmap
}


##plot clustered heatmap of differentially expressed wcp fold change replicates
#DiffExWcp <- c(PO972PO1063PO1075_Wcp_Upreg_AvgFcPval,PO972PO1063PO1075_Wcp_Downreg_AvgFcPval)
DiffExWcp <- c(PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd,PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd)
WCPDiffExQuant <- WCPstimQuant[DiffExWcp,]
WCPDiffExQuantMat <- as.matrix(WCPDiffExQuant[,c("IGDStim4hr_fc","HLMStim4hr_fc","WCPMStim4hr_fc","Stim4hr_AvgFc")])
WCPDiffExQuantMat <- apply(WCPDiffExQuantMat,2,as.numeric)
WCPDiffExQuantMatOrdered <- WCPDiffExQuantMat[order(WCPDiffExQuantMat[,4], decreasing=TRUE), ]
WCPDiffExQuantMat <- as.matrix(WCPDiffExQuantMatOrdered[,c("IGDStim4hr_fc","HLMStim4hr_fc","WCPMStim4hr_fc")])
WCPDiffExQuantMat <- apply(WCPDiffExQuantMat,2,as.numeric)
pdf(file.path(".", analysisdir, "WcpDiffExProteins_PO972PO1063PO1075FoldChangesClustered_heatmap.pdf"))
col_palette <- colorRampPalette(c("blue", "grey90", "red"))(n = 59)
#col_breaks = c(seq(-5.00, -0.32, length=20), # for red
#               seq(-0.31, 0.31, length=20),  # for yellow
#			   seq(0.32, 5.00, length=20)) # for green
col_breaks = c(seq(-5.00, -0.3, length=20), # for red
               seq(-0.29, 0.37, length=20),
			   seq(0.38, 5.00, length=20)) #for green
heatmap.2(WCPDiffExQuantMat,
          notecol="black",      # change font color of cell labels to black
		  density.info="none",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  ###dendrogram="row",    # only draw a row dendrogram
		  dendrogram="none",
		  ###Rowv=TRUE,          # turn off row clustering
		  Rowv="NA",
		  Colv="NA",            # turn off column clustering
		  #col=bluered(100),     #col_palette or col=bluered(100),
		  col=col_palette,
		  breaks=col_breaks,
		  na.color = "white",
		  scale="none",
		  #labRow=UbProtsData$GeneId.x,
		  #labCol=c("WCP 4hr log2FC", "RNAseq 4hr log2FC"),
		  #margins=c(2,28),
		  cexRow=1,
		  cexCol=0.4,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1,
		  sepwidth=c(0.005, 0.00),  # width of the borders
		  #sepcolor='grey80',
		  colsep=0:ncol(WCPDiffExQuantMat),
		  #rowsep=0:nrow(WcpUbProtMat),
		  symkey=F
		 )
dev.off()


##plot clustered heatmap of Ub and Wcp regulation
print("ub-wcp heatmap")
if ( FALSE )
{
WcpUbProtMat <- as.matrix(UbProtsData[,c("UbFcAvg_4hr_NormWcpFc","Stim4hr_AvgFc")])
WcpUbProtMat <- apply(WcpUbProtMat,2,as.numeric)
WcpUbProtMatOrdered <- WcpUbProtMat[order(WcpUbProtMat[,1], decreasing=TRUE), ]
pdf(file.path(".", analysisdir, "UbProteins_UbFcWcpFcNonClustered_heatmap.pdf"))
col_palette <- colorRampPalette(c("blue", "grey90", "red"))(n = 59)
col_breaks = c(seq(-5.00, -0.32, length=20), # for red
               seq(-0.31, 0.31, length=20),  # for yellow
			   seq(0.32, 5.00, length=20)) # for green
heatmap.2(WcpUbProtMatOrdered,
          notecol="black",      # change font color of cell labels to black
		  density.info="none",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  ###dendrogram="row",    # only draw a row dendrogram
		  dendrogram="none",
		  ###Rowv=TRUE,          # turn off row clustering
		  Rowv="NA",
		  Colv="NA",            # turn off column clustering
		  #col=bluered(100),     #col_palette or col=bluered(100),
		  col=col_palette,
		  breaks=col_breaks,
		  na.color = "white",
		  scale="none",
		  #labRow=UbProtsData$GeneId.x,
		  #labCol=c("WCP 4hr log2FC", "RNAseq 4hr log2FC"),
		  #margins=c(2,28),
		  cexRow=1,
		  cexCol=0.4,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1,
		  sepwidth=c(0.005, 0.00),  # width of the borders
		  #sepcolor='grey80',
		  colsep=0:ncol(WcpUbProtMat),
		  #rowsep=0:nrow(WcpUbProtMat),
		  symkey=F
		 )
dev.off()
}


##plot heatmap of TCR proteins 0-4hr fold changes
print("tcr prot heatmap")
#####TCRprotsUbData <- UbProtsData
TCRprotsUbData <- subset(UbProtsData, UbProtsData$GeneId.x%in%TCRprots)
#####TCRprotsUbData <- TCRprotsUbData[-c("Cul1", "H2−T23", "H2−K1", "H2−D1"), ]
#####TCRprotsUbData <- TCRprotsUbData[!rownames(TCRprotsUbData) %in% c("Cul1", "H2−T23", "H2−K1", "H2−D1"), ]
TCRprotsUbData <- subset(TCRprotsUbData, !(TCRprotsUbData$GeneId.x%in%c("Cul1", "H2−T23", "H2−K1", "H2−D1")))
TCRprotsUbData <- subset(TCRprotsUbData, !(row.names(TCRprotsUbData)%in%c("Q3V014", "Q792Z7", "Q7JJ15")))
print(TCRprotsUbData)
TCRprotsUbData <- TCRprotsUbData[rev(order(TCRprotsUbData$UbFcAvg_4hr_NormWcpFc)),]
TCRprotUbFoldChangeMat <- as.matrix(TCRprotsUbData[,c("UbFcAvg_4hr_NormWcpFc", "Stim4hr_AvgFc",  "log2FoldChangeStimUnstim")])
TCRprotUbFoldChangeMat <- apply(TCRprotUbFoldChangeMat,2,as.numeric)
#TCRprotUbFoldChangeMatOrdered <- TCRprotUbFoldChangeMat[order(TCRprotUbFoldChangeMat[,"UbFcAvg_4hr_NormWcpFc"], decreasing=TRUE), ]
TCRprotUbFoldChangeMatOrdered <- TCRprotUbFoldChangeMat
#RNAseq coloring
pdf(file.path(".", analysisdir, "TCRUbproteins_RnaIntensity_heatmap.pdf"))
col_palette <- colorRampPalette(c("blue", "grey90", "red"))(n = 98)
col_breaks_rna = c(seq(-6.00,-1.15,length=33), # for red
                   seq(-1.1499,1.4699,length=33),  # for yellow
				   seq(1.47,6,length=33)) # for green
#col_breaks_rna = c(seq(-6.00,-1.15,length=33), seq(-1.1499,1.4699,length=33), seq(1.47,6,length=33))
#col_breaks_rna = c(seq(-6.00,-1.5,length=33), seq(-1.49,1.99,length=33), seq(2.0,6,length=33))
###TCRprotUbFoldChangeMatRna <- TCRprotUbFoldChangeMatOrdered
###TCRprotUbFoldChangeMatRna[,2] <- NA
###TCRprotUbFoldChangeMatRna[,3] <- NA
heatmap.2(TCRprotUbFoldChangeMatOrdered,
          notecol="black",      # change font color of cell labels to black
		  density.info="none",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  dendrogram="none",    # only draw a row dendrogram
		  Rowv="NA",          # turn off row clustering
		  Colv="NA",            # turn off column clustering
		  #col=bluered(100),     #col_palette or col=bluered(100),
		  col=col_palette,
		  breaks=col_breaks_rna,
		  na.color = "white",
		  scale="none",
		  labRow=TCRprotsUbData$GeneId.x,
		  #labCol=c("WCP 4hr log2FC", "RNAseq 4hr log2FC"),
		  margins=c(2,28),
		  cexRow=1,
		  cexCol=0.4,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1,
		  sepwidth=c(0.10, 0.00),  # width of the borders
		  #sepcolor='grey80',
		  colsep=0:ncol(TCRprotUbFoldChangeMatOrdered),
		  #rowsep=0:nrow(TCRprotUbFoldChangeMatOrdered),
		  symkey=F
		 )
dev.off()
#WCP coloring
pdf(file.path(".", analysisdir, "TCRUbproteins_WcpIntensity_heatmap.pdf"))
col_palette <- colorRampPalette(c("blue", "grey90", "red"))(n = 98)
col_breaks_wcp = c(seq(-5.00,-0.305,length=33), # for red
                   seq(-0.3049,0.3749,length=33),  # for yellow
				   seq(0.375,5,length=33)) # for green
heatmap.2(TCRprotUbFoldChangeMatOrdered,
          notecol="black",      # change font color of cell labels to black
		  density.info="none",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  dendrogram="none",    # only draw a row dendrogram
		  Rowv="NA",          # turn off row clustering
		  Colv="NA",            # turn off column clustering
		  #col=bluered(100),     #col_palette or col=bluered(100),
		  col=col_palette,
		  breaks=col_breaks_wcp,
		  na.color = "white",
		  scale="none",
		  labRow=TCRprotsUbData$GeneId.x,
		  #labCol=c("WCP 4hr log2FC", "RNAseq 4hr log2FC"),
		  margins=c(2,28),
		  cexRow=1,
		  cexCol=0.4,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1,
		  sepwidth=c(0.10, 0.00),  # width of the borders
		  #sepcolor='grey80',
		  colsep=0:ncol(TCRprotUbFoldChangeMatOrdered),
		  #rowsep=0:nrow(TCRprotUbFoldChangeMatOrdered),
		  symkey=F
		 )
dev.off()
#Ub coloring
pdf(file.path(".", analysisdir, "TCRUbproteins_UbIntensity_heatmap.pdf"))
col_palette <- colorRampPalette(c("blue", "grey90", "red"))(n = 98)
col_breaks_ub = c(seq(-5.00, -0.32, length=33), # for red
               seq(-0.31, 0.31, length=33),  # for yellow
               seq(0.32, 5.00, length=33)) # for green
heatmap.2(TCRprotUbFoldChangeMatOrdered,
          notecol="black",      # change font color of cell labels to black
		  density.info="none",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  dendrogram="none",    # only draw a row dendrogram
		  Rowv="NA",          # turn off row clustering
		  Colv="NA",            # turn off column clustering
		  #col=bluered(100),     #col_palette or col=bluered(100),
		  col=col_palette,
		  breaks=col_breaks_ub,
		  na.color = "white",
		  scale="none",
		  labRow=TCRprotsUbData$GeneId.x,
		  #labCol=c("WCP 4hr log2FC", "RNAseq 4hr log2FC"),
		  margins=c(2,28),
		  cexRow=1,
		  cexCol=0.4,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1,
		  sepwidth=c(0.10, 0.00),  # width of the borders
		  #sepcolor='grey80',
		  colsep=0:ncol(TCRprotUbFoldChangeMatOrdered),
		  #rowsep=0:nrow(TCRprotUbFoldChangeMatOrdered),
		  #symkey=F
		 )
dev.off()
##end TCP prot heatmap


##plot heatmap of Ub proteins 0-4hr fold changes
if ( FALSE ){
print("Ub non-norm Fc and Ub norm Fc heatmap")
UbData<-WCPstimQuant[PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts,]
UbFoldChangeMat <- as.matrix(UbData[,c("UbFcAvg","UbFcAvg_4hr_NormWcpFc")])
UbFoldChangeMat <- apply(UbFoldChangeMat,2,as.numeric)
UbFoldChangeMatOrdered <- UbFoldChangeMat[order(UbFoldChangeMat[, dim(UbFoldChangeMat)[2]]), ]
pdf(file.path(".", analysisdir, "Ubproteins_UbFoldChange_heatmap.pdf"))
col_palette <- colorRampPalette(c("blue", "grey90", "red"))(n = 98)
col_breaks = c(seq(-8.00,-0.7,length=33), # for red
               seq(-0.699,0.699,length=33),  # for yellow
               seq(0.7,8,length=33)) # for green
heatmap.2(UbFoldChangeMatOrdered,
          notecol="black",      # change font color of cell labels to black
		  density.info="none",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  dendrogram="none",    # only draw a row dendrogram
		  Rowv="NA",          # turn off row clustering
		  Colv="NA",            # turn off column clustering
		  #col=bluered(100),     #col_palette or col=bluered(100),
		  col=col_palette,
		  breaks=col_breaks,
		  na.color = "white",
		  scale="none",
		  #labRow=UbData$GeneId.x,
		  labCol=c("4hr log2FC", "4hr log2FC (norm)"),
		  margins=c(2,30),
		  cexRow=0.4,
		  cexCol=0.4,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1,
		  #sepwidth=c(0.005, 0.005),  # width of the borders
		  #sepcolor='grey80',
		  #colsep=0:ncol(CullinFoldChangeMat),
		  #rowsep=0:nrow(CullinFoldChangeMat)
		 )
dev.off()
}
##end Ub proteins heatmap

if ( TRUE ){
##plot heatmap of differentially expressed wcp proteins 0-4hr fold changes
print("RNAseq and WCP diff ex heatmap")
#WcpDiffEqIncData <- WCPstimQuant[c(ProtUpregRnaUpreg, ProtUpregRnaUnch, ProtUpregRnaDownreg),]
#WcpDiffEqIncMat <- as.matrix(WcpDiffEqIncData[,c("log2FoldChangeStimUnstim","Stim4hr_AvgFc")])
#WcpDiffEqIncMat <- apply(WcpDiffEqIncMat,2,as.numeric)
#WcpDiffEqIncMatOrdered <- WcpDiffEqIncMat[order(WcpDiffEqIncMat[,1], decreasing=TRUE), ]
#WcpDiffEqDecData <- WCPstimQuant[c(ProtDownregRnaDownreg, ProtDownregRnaUnch, ProtDownregRnaUpreg),]
#WcpDiffEqDecMat <- as.matrix(WcpDiffEqDecData[,c("log2FoldChangeStimUnstim","Stim4hr_AvgFc")])
#WcpDiffEqDecMat <- apply(WcpDiffEqDecMat,2,as.numeric)
#WcpDiffEqDecMatOrdered <- WcpDiffEqDecMat[order(WcpDiffEqDecMat[,1]), ]
#WcpDiffEqMat<-rbind(WcpDiffEqIncMatOrdered,WcpDiffEqDecMatOrdered)
#
WcpDiffEqData <- WCPstimQuant[c(PO972PO1063PO1075_Wcp_Upreg_AvgFcMeanSd, PO972PO1063PO1075_Wcp_Downreg_AvgFcMeanSd),]
WcpDiffEqMat <- as.matrix(WcpDiffEqData[,c("log2FoldChangeStimUnstim","Stim4hr_AvgFc")])
WcpDiffEqMat <- apply(WcpDiffEqMat,2,as.numeric)
WcpDiffEqMatOrdered <- WcpDiffEqMat[order(WcpDiffEqMat[,1], decreasing=TRUE), ]
WcpDiffEqMat <- WcpDiffEqMatOrdered
#
#colored for wcp change
pdf(file.path(".", analysisdir, "RNAseqWCP_WcpDiffExproteins_WcpIntensity_heatmap.pdf"))
col_palette <- colorRampPalette(c("blue", "grey90", "red"))(n = 98)
col_breaks_wcp = c(seq(-5.00,-0.305,length=33), # for red
               seq(-0.3049,0.3749,length=33),  # for yellow
               seq(0.375,5,length=33)) # for green
heatmap.2(WcpDiffEqMat,
          notecol="black",      # change font color of cell labels to black
		  density.info="none",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  dendrogram="none",    # only draw a row dendrogram
		  Rowv="NA",          # turn off row clustering
		  Colv="NA",            # turn off column clustering
		  col=col_palette,
		  #col=bluered(100),
		  #breaks=col_breaks_wcp,
		  na.color = "white",
		  scale="none",
		  #labRow=UbData$GeneId.x,
		  #labCol=c("RNA", "WCP"),
		  margins=c(2,20),
		  cexRow=0.075,
		  cexCol=0.4,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=2,
		  sepwidth=c(0.005, 0.00),  # width of the borders
		  #sepcolor='grey80',
		  colsep=0:ncol(WcpDiffEqMat),
		  #rowsep=0:nrow(WcpDiffEqMat),
		  symkey=F
		 )
dev.off()
#colored for rna change
pdf(file.path(".", analysisdir, "RNAseqWCP_WcpDiffExproteins_RnaIntensity_heatmap.pdf"))
col_palette <- colorRampPalette(c("blue", "grey90", "red"))(n = 98)
col_breaks_rna = c(seq(-6.00,-1.5,length=33), 
               seq(-1.49,1.99,length=33),
			   seq(2.0,6,length=33))
heatmap.2(WcpDiffEqMat,
          notecol="black",      # change font color of cell labels to black
		  density.info="none",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  dendrogram="none",    # only draw a row dendrogram
		  Rowv="NA",          # turn off row clustering
		  Colv="NA",            # turn off column clustering
		  col=col_palette,
		  breaks=col_breaks_rna,
		  na.color = "white",
		  scale="none",
		  #labRow=UbData$GeneId.x,
		  #labCol=c("RNA", "WCP"),
		  margins=c(2,20),
		  cexRow=0.075,
		  cexCol=0.4,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1,
		  sepwidth=c(0.005, 0.00),  # width of the borders
		  #sepcolor='grey80',
		  colsep=0:ncol(WcpDiffEqMat),
		  #rowsep=0:nrow(WcpDiffEqMat),
		  symkey=F
		 )
dev.off()
#clustered by RNAseq and WCP fold change
pdf(file.path(".", analysisdir, "RNAseqWCP_WcpDiffExproteins_ClusteredRnaWcpFc_heatmap.pdf"))
col_palette <- colorRampPalette(c("blue", "grey90", "red"))(n = 98)
col_breaks_wcp = c(seq(-5.00,-0.34,length=33), # for red
                   seq(-0.33,0.33,length=33),  # for yellow
				   seq(0.34,5,length=33)) # for green
heatmap.2(WcpDiffEqMat,
          notecol="black",      # change font color of cell labels to black
		  density.info="none",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  dendrogram="row",    # only draw a row dendrogram
		  #Rowv="TRUE",
		  Colv="NA",            # turn off column clustering
		  #col=bluered(100),
		  col=col_palette,
		  breaks=col_breaks_wcp,
		  na.color = "white",
		  scale="none",
		  #labRow=UbProtsData$GeneId.x,
		  #labCol=c("RNA", "WCP"),
		  margins=c(2,28),
		  cexRow=0.075,
		  cexCol=0.4,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1,
		  sepwidth=c(0.005, 0.00),  # width of the borders
		  #sepcolor='grey80',
		  colsep=0:ncol(WcpDiffEqMat),
		  #rowsep=0:nrow(WcpUbProtMat),
		  symkey=F
		 )
dev.off()
##end WCPdiffEq proteins heatmap
}


if ( FALSE )
{
##plot number of protein quantifications for each exp
pdf(file.path(".", analysisdir, "CullinAdaptors_AbundanceZscores.pdf"), width=4, height=7)
quants <- vector()
for ( i in 1:length(cullins) ){
	appendgene <- cullins[i]
	quants <- c(quants, WCPstimQuant[which(WCPstimQuant$GeneId.x == appendgene),"AbundanceInt_zscore"])
}
barplot(quants, beside=TRUE, names.arg=cullins, xlab="", ylab="", ylim=c(-1,2), col=c("grey70"), width=c(2), space=c(1,10), cex.axis=1, cex.names=1)
mtext("abundance Z-score", side=2, line=2.75, cex=1)
dev.off()
}


##plot distribution of counts for each experiment
plotnormdistribution(WCPstimQuant$MSMScounts.PO972igd_0hr, "NA", "NA", 1, "NA", "NA", "NA", "NA", "MSMS counts", "norm freq", file.path(".", analysisdir, "PO972igd_0hrMSMScounts_Distribution.pdf"))
plotnormdistribution(WCPstimQuant$MSMScounts.PO972igd_1hr, "NA", "NA", 1, "NA", "NA", "NA", "NA", "MSMS counts", "norm freq", file.path(".", analysisdir, "PO972igd_1hrMSMScounts_Distribution.pdf"))
plotnormdistribution(WCPstimQuant$MSMScounts.PO972igd_4hr, "NA", "NA", 1, "NA", "NA", "NA", "NA", "MSMS counts", "norm freq", file.path(".", analysisdir, "PO972igd_4hrMSMScounts_Distribution.pdf"))
plotnormdistribution(WCPstimQuant$IGD0hr4hr_MinMSMScount, "NA", "NA", 1, "NA", "NA", "NA", "NA", "min MSMS counts", "norm freq", file.path(".", analysisdir, "PO972igd_0hr4hrMinMSMScounts_Distribution.pdf"))
plotnormdistribution(WCPstimQuant$IGD0hr4hr_MaxMSMScount, "NA", "NA", 1, "NA", "NA", "NA", "NA", "max MSMS counts", "norm freq", file.path(".", analysisdir, "PO972igd_0hr4hrMaxMSMScounts_Distribution.pdf"))
plotnormdistribution(WCPstimQuant$LHcounts.PO1063_HLM, "NA", "NA", 1, "NA", "NA", "NA", "NA", "MSMS counts", "norm freq", file.path(".", analysisdir, "PO1063HLM_MSMScounts_Distribution.pdf"))
plotnormdistribution(WCPstimQuant$LHcounts.PO1075_WCPM, "NA", "NA", 1, "NA", "NA", "NA", "NA", "MSMS counts", "norm freq", file.path(".", analysisdir, "PO1075WCPM_MSMScounts_Distribution.pdf"))


##plot relationship of KeGG 0hr and KeGG 4hr MSMS counts
pdf(file.path(".", analysisdir, "Kegg_0hr4hr_MSMScounts_Boxplot.pdf"))
boxplot(WCPstimQuant$MSMScounts.PO972igd_4hr~WCPstimQuant$MSMScounts.PO972igd_0hr, xlim=c(0,20), ylim=c(0,50))
dev.off()


##comparison of counts across experiments
makescatterplot(WCPstimQuant$IGD0hr4hr_MaxMSMScount, WCPstimQuant$LHcounts.PO1063_HLM, "NA", "NA", "NA", "NA", "PO972 0-4hr max counts", "PO1063 LHratio counts", file.path(".", analysisdir, "PO972PO1063_CountComp.pdf"))
makescatterplot(WCPstimQuant$IGD0hr4hr_MaxMSMScount, WCPstimQuant$LHcounts.PO1075_WCPM, "NA", "NA", "NA", "NA", "PO972 0-4hr max counts", "PO1075 LHratio counts", file.path(".", analysisdir, "PO972PO1075_CountComp.pdf"))
makescatterplot(WCPstimQuant$LHcounts.PO1063_HLM, WCPstimQuant$LHcounts.PO1075_WCPM, "NA", "NA", "NA", "NA", "PO1063 LHratio counts", "PO1075 LHratio counts", file.path(".", analysisdir, "PO1063PO1075_CountComp.pdf"))


##plot relationships between counts and diff eq
pdf(file.path(".", analysisdir, "PO972_MinCountDiffExComp.pdf"))
boxplot(WCPstimQuant$IGDStim4hr_fc~WCPstimQuant$IGD0hr4hr_MinMSMScount, xlim=c(0,20), xlab="PO972 0hr-4hr min MSMS count", ylab="PO972 0hr-4hr Log2 fc", cex.axis=0.75)
dev.off()
pdf(file.path(".", analysisdir, "PO972_MaxCountDiffExComp.pdf"))
boxplot(WCPstimQuant$IGDStim4hr_fc~WCPstimQuant$IGD0hr4hr_MaxMSMScount, xlim=c(0,20), xlab="PO972 0hr-4hr max MSMS count", ylab="PO972 0hr-4hr Log2 fc", cex.axis=0.75)
dev.off()
pdf(file.path(".", analysisdir, "PO1063_CountDiffExComp.pdf"))
boxplot(WCPstimQuant$normLHratio.PO1063_HLM~WCPstimQuant$LHcounts.PO1063_HLM, xlim=c(0,20), xlab="PO1063 0hr-4hr LHcount", ylab="PO1063 0hr-4hr Log2 normLHratio", cex.axis=0.75)
dev.off()
pdf(file.path(".", analysisdir, "PO1075_CountDiffExComp.pdf"))
boxplot(WCPstimQuant$normLHratio.PO1075_WCPM~WCPstimQuant$LHcounts.PO1075_WCPM, xlim=c(0,20), xlab="PO1075 0hr-4hr LHcount", ylab="PO1075 0hr-4hr Log2 normLHratio", cex.axis=0.75)
dev.off()


##plot intensity comparisons
makescatterplot(WCPstimQuant$PO972igd_0hr, WCPstimQuant$IntensityH.PO1063_HLM, "NA", "NA", "NA", "NA", "PO972 iBAQ 0hr intensity", "PO1063 SILAC Heavy (0hr) intensity", file.path(".", analysisdir, "PO972PO1063_0hr_IntensityComp.pdf"))
makescatterplot(WCPstimQuant$PO972igd_0hr, WCPstimQuant$IntensityH.PO1075_WCPM, "NA", "NA", "NA", "NA", "PO972 iBAQ 0hr intensity", "PO1075 SILAC Heavy (0hr) intensity", file.path(".", analysisdir, "PO972PO1075_0hr_IntensityComp.pdf"))
makescatterplot(WCPstimQuant$IntensityH.PO1063_HLM, WCPstimQuant$IntensityH.PO1075_WCPM, "NA", "NA", "NA", "NA", "PO1063 SILAC Heavy (0hr) intensity", "PO1075 SILAC Heavy (0hr) intensity", file.path(".", analysisdir, "PO1063PO1075_0hr_IntensityComp.pdf"))
makescatterplot(WCPstimQuant$PO972igd_4hr, WCPstimQuant$IntensityL.PO1063_HLM, "NA", "NA", "NA", "NA", "PO972 iBAQ 4hr intensity", "PO1063 SILAC Light (4hr) intensity", file.path(".", analysisdir, "PO972PO1063_4hr_IntensityComp.pdf"))
makescatterplot(WCPstimQuant$PO972igd_4hr, WCPstimQuant$IntensityL.PO1075_WCPM, "NA", "NA", "NA", "NA", "PO972 iBAQ 4hr intensity", "PO1075 SILAC Light (4hr) intensity", file.path(".", analysisdir, "PO972PO1075_4hr_IntensityComp.pdf"))
makescatterplot(WCPstimQuant$IntensityL.PO1063_HLM, WCPstimQuant$IntensityL.PO1075_WCPM, "NA", "NA", "NA", "NA", "PO1063 SILAC Light (4hr) intensity", "PO1075 SILAC Light (4hr) intensity", file.path(".", analysisdir, "PO1063PO1075_4hr_IntensityComp.pdf"))

makedensitycoloredscatterplot(WCPstimQuant$PO972igd_0hr, WCPstimQuant$IntensityH.PO1063_HLM, -12, 12, -12, 12, "PO972 iBAQ 0hr intensity", "PO1063 SILAC Heavy (0hr) intensity", file.path(".", analysisdir, "PO972PO1063_0hr_IntensityComp_DensityPlot.pdf"))
makedensitycoloredscatterplot(WCPstimQuant$PO972igd_0hr, WCPstimQuant$IntensityH.PO1075_WCPM, -12, 12, -12, 12, "PO972 iBAQ 0hr intensity", "PO1075 SILAC Heavy (0hr) intensity", file.path(".", analysisdir, "PO972PO1075_0hr_IntensityComp_DensityPlot.pdf"))
makedensitycoloredscatterplot(WCPstimQuant$IntensityH.PO1063_HLM, WCPstimQuant$IntensityH.PO1075_WCPM, -12, 12, -12, 12, "PO1063 SILAC Heavy (0hr) intensity", "PO1075 SILAC Heavy (0hr) intensity", file.path(".", analysisdir, "PO1063PO1075_0hr_IntensityComp_DensityPlot.pdf"))
makedensitycoloredscatterplot(WCPstimQuant$PO972igd_4hr, WCPstimQuant$IntensityL.PO1063_HLM, -12, 12, -12, 12, "PO972 iBAQ 4hr intensity", "PO1063 SILAC Light (4hr) intensity", file.path(".", analysisdir, "PO972PO1063_4hr_IntensityComp_DensityPlot.pdf"))
makedensitycoloredscatterplot(WCPstimQuant$PO972igd_4hr, WCPstimQuant$IntensityL.PO1075_WCPM, -12, 12, -12, 12, "PO972 iBAQ 4hr intensity", "PO1075 SILAC Light (4hr) intensity", file.path(".", analysisdir, "PO972PO1075_4hr_IntensityComp_DensityPlot.pdf"))
makedensitycoloredscatterplot(WCPstimQuant$IntensityL.PO1063_HLM, WCPstimQuant$IntensityL.PO1075_WCPM, -12, 12, -12, 12, "PO1063 SILAC Light (4hr) intensity", "PO1075 SILAC Light (4hr) intensity", file.path(".", analysisdir, "PO1063PO1075_4hr_IntensityComp_DensityPlot.pdf"))

#tmpdat <- WCPstimQuant[,c("PO972igd_0hr","MSMScounts.PO972igd_0hr","IntensityH.PO1063_HLM","LHcounts.PO1063_HLM")]
#makescatterplotcolbypeps(tmpdat, "NA", "NA", "NA", "NA", "PO972 iBAQ 0hr intensity", "PO1063 SILAC Heavy (0hr) intensity", file.path(".", analysisdir, "PO972PO1063_0hr_IntensityComp_ColByPep.pdf"))


##plot z-score comparisons
makescatterplot(WCPstimQuant$IGD0hr_zscore, WCPstimQuant$HLMintH_zscore, "NA", "NA", "NA", "NA", "PO972 iBAQ 0hr z-score", "PO1063 SILAC IntH (0hr) z-score", file.path(".", analysisdir, "PO972PO1063_0hr_ZscoreComp.pdf"))
makescatterplot(WCPstimQuant$IGD0hr_zscore, WCPstimQuant$WCPMintH_zscore, "NA", "NA", "NA", "NA", "PO972 iBAQ 0hr z-score", "PO1075 SILAC IntH (0hr) z-score", file.path(".", analysisdir, "PO972PO1075_0hr_ZscoreComp.pdf"))
makescatterplot(WCPstimQuant$HLMintH_zscore, WCPstimQuant$WCPMintH_zscore, "NA", "NA", "NA", "NA", "PO1063 SILAC IntH (0hr) z-score", "PO1075 SILAC IntH (0hr) z-score", file.path(".", analysisdir, "PO1063PO1075_0hr_ZscoreComp.pdf"))
makescatterplot(WCPstimQuant$IGD4hr_zscore, WCPstimQuant$HLMintL_zscore, "NA", "NA", "NA", "NA", "PO972 iBAQ 4hr z-score", "PO1063 SILAC IntL (4hr) z-score", file.path(".", analysisdir, "PO972PO1063_4hr_ZscoreComp.pdf"))
makescatterplot(WCPstimQuant$IGD4hr_zscore, WCPstimQuant$WCPMintL_zscore, "NA", "NA", "NA", "NA", "PO972 iBAQ 4hr z-score", "PO1075 SILAC IntL (4hr) z-score", file.path(".", analysisdir, "PO972PO1075_4hr_ZscoreComp.pdf"))
makescatterplot(WCPstimQuant$HLMintL_zscore, WCPstimQuant$WCPMintL_zscore, "NA", "NA", "NA", "NA", "PO1063 SILAC IntL (4hr) z-score", "PO1075 SILAC IntL (4hr) z-score", file.path(".", analysisdir, "PO1063PO1075_4hr_ZscoreComp.pdf"))
makedensitycoloredscatterplot(WCPstimQuant$IGD0hr_zscore, WCPstimQuant$HLMintH_zscore, "NA", "NA", "NA", "NA", "PO972 iBAQ 0hr z-score", "PO1063 SILAC IntH (0hr) z-score", file.path(".", analysisdir, "PO972PO1063_0hr_ZscoreComp_DensityPlot.pdf"))
makedensitycoloredscatterplot(WCPstimQuant$IGD0hr_zscore, WCPstimQuant$WCPMintH_zscore, "NA", "NA", "NA", "NA", "PO972 iBAQ 0hr z-score", "PO1075 SILAC IntH (0hr) z-score", file.path(".", analysisdir, "PO972PO1075_0hr_ZscoreComp_DensityPlot.pdf"))
makedensitycoloredscatterplot(WCPstimQuant$HLMintH_zscore, WCPstimQuant$WCPMintH_zscore, "NA", "NA", "NA", "NA", "PO1063 SILAC IntH (0hr) z-score", "PO1075 SILAC IntH (0hr) z-score", file.path(".", analysisdir, "PO1063PO1075_0hr_ZscoreComp_DensityPlot.pdf"))
makedensitycoloredscatterplot(WCPstimQuant$IGD4hr_zscore, WCPstimQuant$HLMintL_zscore, "NA", "NA", "NA", "NA", "PO972 iBAQ 4hr z-score", "PO1063 SILAC IntL (4hr) z-score", file.path(".", analysisdir, "PO972PO1063_4hr_ZscoreComp_DensityPlot.pdf"))
makedensitycoloredscatterplot(WCPstimQuant$IGD4hr_zscore, WCPstimQuant$WCPMintL_zscore, "NA", "NA", "NA", "NA", "PO972 iBAQ 4hr z-score", "PO1075 SILAC IntL (4hr) z-score", file.path(".", analysisdir, "PO972PO1075_4hr_ZscoreComp_DensityPlot.pdf"))
makedensitycoloredscatterplot(WCPstimQuant$HLMintL_zscore, WCPstimQuant$WCPMintL_zscore, "NA", "NA", "NA", "NA", "PO1063 SILAC IntL (4hr) z-score", "PO1075 SILAC IntL (4hr) z-score", file.path(".", analysisdir, "PO1063PO1075_4hr_ZscoreComp_DensityPlot.pdf"))
##end plot z-score comparisons


##plot RNAseq counts vs abundance expression percentage
filename<-file.path(".", analysisdir, "RNAseqCountsProteinExpProb.pdf")
pdf(filename)
par(bty="l")
xmin<-floor(min(bins,na.rm=TRUE))
xmax<-ceiling(max(bins,na.rm=TRUE))
ymin<-floor(min(abdpcts,na.rm=TRUE))
ymax<-ceiling(max(abdpcts,na.rm=TRUE))
#make an axes
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
par(new=T)
plot(bins, abdpcts,
     type="b",
	 pch=19,
	 cex=1,
	 col="gray50",
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n',yaxt='n'
	)
#label and format axes
axis(1, c(seq(xmin,xmax,by=0.5)), labels=F, col="black",cex.axis=1, tck=-0.01)
axis(1, c(seq(xmin,xmax,by=1.0)), labels=T, col="black",cex.axis=2)
mtext("RNAseq counts (log10)",1,line=2.75,cex=3)
axis(2, c(seq(ymin,ymax,by=0.05)), labels=F, col="black",cex.axis=1, tck=-0.01)
axis(2, c(seq(ymin,ymax,by=0.1)), labels=T, col="black",cex.axis=2)
mtext("protein wcp expression prob",2,line=2.75,cex=3)
dev.off()


##plot fold change comparisons
#scatter plot
makescatterplot(WCPstimQuant$IGDStim4hr_fc, WCPstimQuant$HLMStim4hr_fc, "NA", "NA", "NA", "NA", "PO972 4hr TCR stim iBAQ log2 fc", "PO1063 4hr TCR stim SILAC intL-H log2 fc", file.path(".", analysisdir, "PO972PO1063_0hr4hrStim_FoldChangeComp.pdf"))
makescatterplot(WCPstimQuant$IGDStim4hr_fc, WCPstimQuant$WCPMStim4hr_fc, "NA", "NA", "NA", "NA", "PO972 4hr TCR stim iBAQ log2 fc", "PO1075 4hr TCR stim SILAC intL-H log2 fc", file.path(".", analysisdir, "PO972PO1075_0hr4hrStim_FoldChangeComp.pdf"))
makescatterplot(WCPstimQuant$HLMStim4hr_fc, WCPstimQuant$WCPMStim4hr_fc, "NA", "NA", "NA", "NA", "PO1063 4hr TCR stim SILAC intL-H log2 fc", "PO1075 4hr TCR stim SILAC intL-H log2 fc", file.path(".", analysisdir, "PO1063PO1075_0hr4hrStim_FoldChangeComp.pdf"))
makescatterplot(WCPstimQuant$IGDStim1hr_fc, WCPstimQuant$IGDStim4hr_fc, "NA", "NA", "NA", "NA", "PO972 1hr TCR stim iBAQ log2 fc", "PO972 4hr TCR stim iBAQ log2 fc", file.path(".", analysisdir, "PO9721hrPO9724hr_1hr4hrStim_FoldChangeComp.pdf"))
#density colored scatter plot
makedensitycoloredscatterplot(WCPstimQuant$IGDStim4hr_fc, WCPstimQuant$HLMStim4hr_fc, "NA", "NA", "NA", "NA", "PO972 4hr TCR stim iBAQ log2 fc", "PO1063 4hr TCR stim SILAC intL-H log2 fc", file.path(".", analysisdir, "PO972PO1063_0hr4hrStim_FoldChangeComp_DensityPlot.pdf"))
makedensitycoloredscatterplot(WCPstimQuant$IGDStim4hr_fc, WCPstimQuant$WCPMStim4hr_fc, "NA", "NA", "NA", "NA", "PO972 4hr TCR stim iBAQ log2 fc", "PO1075 4hr TCR stim SILAC intL-H log2 fc", file.path(".", analysisdir, "PO972PO1075_0hr4hrStim_FoldChangeComp_DensityPlot.pdf"))
makedensitycoloredscatterplot(WCPstimQuant$HLMStim4hr_fc, WCPstimQuant$WCPMStim4hr_fc, "NA", "NA", "NA", "NA", "PO1063 4hr TCR stim SILAC intL-H log2 fc", "PO1075 4hr TCR stim SILAC intL-H log2 fc", file.path(".", analysisdir, "PO1063PO1075_0hr4hrStim_FoldChangeComp_DensityPlot.pdf"))
#msms count colored scatter plot
tmpdat <- WCPstimQuant[,c("IGDStim4hr_fc","IGD0hr4hr_MinMSMScount","HLMStim4hr_fc","LHcounts.PO1063_HLM")]
makescatterplotcolbypeps(tmpdat, "NA", "NA", "NA", "NA", "PO972 4hr TCR stim iBAQ log2 fc", "PO1063 4hr TCR stim SILAC intL-H log2 fc", file.path(".", analysisdir, "PO972PO1063_0hr4hrStim_FoldChangeComp_ColByPep.pdf"))
#comparisons for only proteins that were found on all three experiments
tmpdat <- WCPstimQuant[row.names(WCPstimQuant)%in%PO972PO1063PO1075_AllId_Intersect_ProtId,]
makedensitycoloredscatterplot(tmpdat$IGDStim4hr_fc, tmpdat$HLMStim4hr_fc, "NA", "NA", "NA", "NA", "PO972 4hr TCR stim iBAQ log2 fc", "PO1063 4hr TCR stim SILAC intL-H log2 fc", file.path(".", analysisdir, "PO972PO1063_0hr4hrStim_ConsProts_FoldChangeComp_DensityPlot.pdf"))
makedensitycoloredscatterplot(tmpdat$IGDStim4hr_fc, tmpdat$WCPMStim4hr_fc, "NA", "NA", "NA", "NA", "PO972 4hr TCR stim iBAQ log2 fc", "PO1075 4hr TCR stim SILAC intL-H log2 fc", file.path(".", analysisdir, "PO972PO1075_0hr4hrStim_ConsProts_FoldChangeComp_DensityPlot.pdf"))
makedensitycoloredscatterplot(tmpdat$HLMStim4hr_fc, tmpdat$WCPMStim4hr_fc, "NA", "NA", "NA", "NA", "PO1063 4hr TCR stim SILAC intL-H log2 fc", "PO1075 4hr TCR stim SILAC intL-H log2 fc", file.path(".", analysisdir, "PO1063PO1075_0hr4hrStim_ConsProts_FoldChangeComp_DensityPlot.pdf"))


#plot fold change distributions
plotnormdistribution(WCPstimQuant$IGDStim4hr_fc, "NA", "NA", 0.1, "NA", "NA", "NA", "NA", "PO972 0-4hr TCR stim log2 fold change", "norm freq", file.path(".", analysisdir, "PO972_FoldChangeDistribution.pdf"))
plotnormdistribution(WCPstimQuant$HLMStim4hr_fc, "NA", "NA", 0.1, "NA", "NA", "NA", "NA", "PO1063 0-4hr TCR stim log2 fold change", "norm freq", file.path(".", analysisdir, "PO1063_FoldChangeDistribution.pdf"))
plotnormdistribution(WCPstimQuant$WCPMStim4hr_fc, "NA", "NA", 0.1, "NA", "NA", "NA", "NA", "PO1075 0-4hr TCR stim log2 fold change", "norm freq", file.path(".", analysisdir, "PO1075_FoldChangeDistribution.pdf"))
plotnormdistribution(WCPstimQuant$Stim4hr_AvgFc, "NA", "NA", 0.05, -2, 2, "NA", 0.12, "log2 fold change", "normalized frequency", file.path(".", analysisdir, "PO972-PO1063-PO1075_AvgFoldChangeDistribution.pdf"))
plotnormdistribution(WCPstimQuant$log2FoldChangeStimUnstim, "NA", "NA", 0.25, -5, 5, "NA", "NA", "log2 fold change", "normalized frequency", file.path(".", analysisdir, "RNAseq_AvgFoldChangeDistribution.pdf"))
plotnormdistribution(RnaData$log10baseMean.x, "NA", "NA", 0.1, "NA", "NA", "NA", 0.06, "log10 RNAseq count", "normalized frequency", file.path(".", analysisdir, "RNAseq_MeanCountDistribution.pdf"))


#plot ubiquitylation fold change distributions
plotnormdistribution(WCPstimQuant$intensityweightedmeanratio.PO972kgg_4diff0, "NA", "NA", 0.2, "NA", "NA", "NA", "NA", "log2 fold change", "norm freq", file.path(".", analysisdir, "PO972kgg_KeGGfoldchange_Distribution.pdf"))
plotnormdistribution(WCPstimQuant$intensityweightedmeanratio.PO1063kgg_HLM, "NA", "NA", 0.2, "NA", "NA", "NA", "NA", "log2 fold change", "norm freq", file.path(".", analysisdir, "PO1063HLM_KeGGfoldchange_Distribution.pdf"))
plotnormdistribution(WCPstimQuant$UbFcAvg, "NA", "NA", 0.1, "NA", "NA", "NA", 0.16, "log2 fold change", "norm freq", file.path(".", analysisdir, "PO972kggPO1063HLM_AvgKeGGfoldchange_Distribution.pdf"))
plotnormdistribution(WCPstimQuant$UbFcAvg_4hr_NormWcpFc, "NA", "NA", 0.1, "NA", "NA", "NA", "NA", "log2 fold change norm to wcp", "norm freq", file.path(".", analysisdir, "PO972kggPO1063HLM_AvgNormKeGGfoldchange_Distribution.pdf"))


##plot KeGG protein quantification comparisons
#scatterplot
makescatterplot(KeGGnormLHdata$intensityweightedmeanratio.PO972kgg_4diff0, KeGGnormLHdata$intensityweightedmeanratio.PO1063kgg_HLM, -3, 3, -3, 3, "Ub (exp 1) log2 fold change", "Ub (exp 2) log2 fold change", file.path(".", analysisdir, "PO972PO1063_0hr4hrStim_KeGGQuantFoldChangeComp.pdf"))
makescatterplot(WCPstimQuant$PO972kgg_NormWcpFc, WCPstimQuant$PO1063HLMkgg_NormWcpFc, -6, 4, -5, 3, "Ub (exp 1) log2 fold change", "Ub (exp 2) log2 fold change", file.path(".", analysisdir, "PO972PO1063_0hr4hrStim_NormKeGGQuantFoldChangeComp.pdf"))
#density colored scatterplot
makedensitycoloredscatterplot(KeGGnormLHdata$intensityweightedmeanratio.PO972kgg_4diff0, KeGGnormLHdata$intensityweightedmeanratio.PO1063kgg_HLM, "NA", "NA", "NA", "NA", "protein Ub (exp 1) log2 fold change", "protein Ub (exp 2) log2 fold change", file.path(".", analysisdir, "PO972PO1063_0hr4hrStim_KeGGQuantFoldChangeComp_DensityPlot.pdf"))
makedensitycoloredscatterplot(WCPstimQuant$PO972kgg_NormWcpFc, WCPstimQuant$PO1063HLMkgg_NormWcpFc, -6, 4, -5, 3, "protein Ub (exp 1) log2 fold change", "protein Ub (exp 2) log2 fold change", file.path(".", analysisdir, "PO972PO1063_0hr4hrStim_NormKeGGQuantFoldChangeComp_DensityPlot.pdf"))


##plot number of protein quantifications for each exp
pdf(file.path(".", analysisdir, "PO972PO1063_NumKeGGUbProteins0hr4hr.pdf"), width=4, height=7)
quants<-c(PO1063HLMkgg_0hr4hr_UbProts_Num, PO972kgg_0hr4hr_UbProts_Num)
barplot(quants, beside=TRUE, names.arg=c("exp 1","exp 2"), xlab="", ylab="", ylim=c(0,1000), col=c("steelblue","steelblue1"), cex.axis=1.25, cex.names=1.25)
mtext("num KeGG identified proteins", side=2, line=2.75, cex=1.75)
dev.off()



if ( FALSE )
{
##### where are these variables
##plot number of protein quantifications for each exp
pdf(file.path(".", analysisdir, "PO972-PO1063-PO1075_NumQuantifiedProteins0hr4hr.pdf"), width=5, height=7)
quants<-c(PO972ProtId_Num,PO1063ProtIdNum,PO1075ProtIdNum)
barplot(quants, beside=TRUE, names.arg=c("LFQ", "SILAC (I)", "SILAC (II)"), xlab="", ylab="", ylim=c(0,7000), col=c("steelblue4","dodgerblue1","dodgerblue2"), cex.axis=1.25, cex.names=1.25)
mtext("number of identified proteins", side=2, line=2.75, cex=1.75)
dev.off()

##venn diagram of intersection between quantified proteins in each exp
pdf(file.path(".", analysisdir, "PO972-PO1063-PO1075_NumQuantifiedProteins0hr4hr_VennDiagram.pdf"))
draw.triple.venn(PO972ProtId_Num,
                 PO1063ProtIdNum,
				 PO1075ProtIdNum,
				 PO972PO1063_ProtId_Int_Num,
				 PO1063PO1075_ProtId_Int_Num,
				 PO972PO1075_ProtId_Int_Num,
				 PO972PO1063PO1075_AllId_Intersect_ProtId_Num,
				 #category=c("KeGG","Nedd", "TUBE"),
				 overrideTriple=1, euler.d=FALSE, scaled=FALSE,
				 lwd=rep(6,3), col=c("steelblue4","dodgerblue1","dodgerblue2"), fill=c("steelblue2","skyblue1","lightskyblue1"), alpha=rep(0.25,3),
				 label.col=rep("gray40",7), fontface=rep("plain",7), cex=c(1.5,1.5,1.5,1.5,2,1.5,1.5)
				)
dev.off()

##venn diagram of intersection between quantified proteins in each exp
pdf(file.path(".", analysisdir, "Wcp-Ub-RNAseq_NumQuantifiedProteins0hr4hr_VennDiagram.pdf"))
draw.triple.venn(length(PO972PO1063PO1075_AllId_Union_ProtId),
                 length(PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts),
				 length(RnaIdentifiedGenes_InTranscriptome),
				 length(intersect(PO972PO1063PO1075_AllId_Union_ProtId,PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts)),
				 length(intersect(PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts,RnaIdentifiedGenes_InProteome)),
				 length(intersect(PO972PO1063PO1075_AllId_Union_ProtId,RnaIdentifiedGenes_InProteome)),
				 length(Reduce(intersect, list(PO972PO1063PO1075_AllId_Union_ProtId, RnaIdentifiedGenes_InProteome, PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts))),
				 category=c("Wcp","Ub", "RNAseq"),
				 overrideTriple=1, euler.d=FALSE, scaled=FALSE,
				 lwd=rep(6,3), col=c("goldenrod","dodgerblue1","olivedrab4"), fill=c("gold","skyblue1","olivedrab3"), alpha=rep(0.25,3),
				 label.col=rep("gray40",7), fontface=c("plain","plain","plain","plain","bold","plain","plain"), cex=c(1.5,1.5,1.5,1.5,2,1.5,1.5)
				)
dev.off()
} #end venn diagram FALSE section

##venn diagram of intersection upregulated proteins in each experiment
pdf(file.path(".", analysisdir, "PO972PO1063PO1075_UpregMeanSd0hr4hr_VennDiagram.pdf"))
draw.triple.venn(PO972_Upreg_MeanSd_Num,
                 PO1063_Upreg_MeanSd_Num,
                 PO1075_Upreg_MeanSd_Num,
                 PO1063PO972_MeanSd_Upreg_Num,
                 PO1075PO1063_MeanSd_Upreg_Num,
                 PO1075PO972_MeanSd_Upreg_Num,
                 PO1075PO1063PO972_MeanSd_Upreg_Num,
                 category=c("PO972","PO1063", "PO1075"),
                 overrideTriple=1, euler.d=TRUE, scaled=TRUE,
                 lwd=rep(4,3), col=c("dodgerblue4","skyblue4","cadetblue4"), fill=c("dodgerblue1","skyblue1","cadetblue1"),
                 label.col=rep("gray40",7), fontface=rep("plain",7), cex=rep(2,7)
                )
dev.off()


##venn diagram of intersection of downregulated proteins in each experiment
pdf(file.path(".", analysisdir, "PO972PO1063PO1075_DownregMeanSd0hr4hr_VennDiagram.pdf"))
draw.triple.venn(PO972_Downreg_MeanSd_Num,
                 PO1063_Downreg_MeanSd_Num,
                 PO1075_Downreg_MeanSd_Num,
                 PO1063PO972_MeanSd_Downreg_Num,
                 PO1075PO1063_MeanSd_Downreg_Num,
                 PO1075PO972_MeanSd_Downreg_Num,
                 PO1075PO1063PO972_MeanSd_Downreg_Num,
                 category=c("PO972","PO1063", "PO1075"),
                 overrideTriple=1, euler.d=TRUE, scaled=TRUE,
                 lwd=rep(4,3), col=c("gold3","darkgoldenrod4","goldenrod3"), fill=c("gold1","darkgoldenrod1","goldenrod1"),
                 label.col=rep("gray40",7), fontface=rep("plain",7), cex=rep(2,7)
                )
dev.off()


##scatterplot of WCP vs RNAseq for proteins increased in ubiqutination
UbUpregData<-WCPstimQuant[row.names(WCPstimQuant)%in%IncreaseUb,]
ModelDegData<-WCPstimQuant[row.names(WCPstimQuant)%in%ModelDeg_UbInc,]
ModelNonDegData<-WCPstimQuant[row.names(WCPstimQuant)%in%ModelNonDeg_UbInc,]
ModelNonDegData<-ModelNonDegData[rownames(ModelNonDegData) != "Q5F2A7", ]
#print(ModelNonDegData)
pdf(file.path(".", analysisdir, "ProteomicsRNAseq_Log2FoldChangeComp_UbUpregProts.pdf"))
par(bty="l")
xmin <- floor(min(UbUpregData$Stim4hr_AvgFc, na.rm=TRUE))
#####xmax <- ceiling(max(UbUpregData$Stim4hr_AvgFc, na.rm=TRUE))
xmax <- 1.5
ymin <- floor(min(UbUpregData$log2FoldChangeStimUnstim, na.rm=TRUE))
#####ymax <- ceiling(max(UbUpregData$log2FoldChangeStimUnstim, na.rm=TRUE))
ymax <- 3.5
#make axes
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
#add horizontal and vertical lines at 0
par(new=T)
abline(v=0, col="gray70", lwd=2, lty=1)
par(new=T)
abline(h=0, col="gray70", lwd=2, lty=1)
par(new=T)
#add vertical lines at x-axis proteomics mean +/- standard deviation
par(new=T)
#####abline(v=PO972PO1063PO1075_AvgFcMean+PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=1, lty=2)
par(new=T)
#####abline(v=PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=1, lty=2)
#add horizontal lines at y-axis RNAseq mean +/- standard deviation
par(new=T)
#####abline(h=RNAseq_FcMean+RNAseq_FcSd, col="gray70", lwd=1, lty=2)
par(new=T)
#####abline(h=RNAseq_FcMean-RNAseq_FcSd, col="gray70", lwd=1, lty=2)
#plot points
par(new=T)
plot(ModelNonDegData$Stim4hr_AvgFc,
     ModelNonDegData$log2FoldChangeStimUnstim,
	 pch=1,
	 cex=3*(abs(ModelNonDegData$UbFcAvg_4hr_NormWcpFc)),
	 col="red",
	 lwd=2.5,
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
	)
par(new=T)
plot(ModelDegData$Stim4hr_AvgFc,
	 ModelDegData$log2FoldChangeStimUnstim,
	 pch=1,
	 cex=3*(abs(ModelDegData$UbFcAvg_4hr_NormWcpFc)),
	 col="grey50",
	 lwd=2.5,
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
	)
#text(ModelDegData$Stim4hr_AvgFc,
#     ModelDegData$log2FoldChangeStimUnstim,
#	 labels=ModelDegData$GeneId.x,
#	 cex=0.5, pos=2, offset=0.2, col="red", font=2
#   )
#testgenes<-c("Zap70", "Lat", "Snx18", "Sin3b", "Prkcq", "Grap", "Mycbp2")
testgenes<-c("P43404", "A0A0U1RP03", "Q8C788", "Q62141-1", "A6H667", "Q9CX99", "F6SMY7")
#testgenes<-c("P25799", "A2RTT4", "P43404", "A0A0U1RP03", "Q9CX99", "Q4VAE6", "P60766", "K7Q7T7", "A6H667", "P24161", "Q3U4Y3", "A6H6M1")
TestGenesData<-WCPstimQuant[row.names(WCPstimQuant)%in%testgenes,]
TestGenesData["F6SMY7","UbFcAvg_4hr_NormWcpFc"]<-max(TestGenesData$UbFcAvg_4hr_NormWcpFc, na.rm=TRUE)
par(new=T)
plot(TestGenesData$Stim4hr_AvgFc,
	 TestGenesData$log2FoldChangeStimUnstim,
	 pch=21,
	 cex=3*(abs(TestGenesData$UbFcAvg_4hr_NormWcpFc)),
	 #col=ifelse((TestGenesData$GeneId.x=="Nfkb1" | TestGenesData$GeneId.x=="Ube2n" | TestGenesData$GeneId.x=="Zap70" | TestGenesData$GeneId.x=="Lat"), "grey50", "red"),
	 col=ifelse((TestGenesData$GeneId.x=="Prkcq" | TestGenesData$GeneId.x=="Grap"), "red", "grey50"),
	 bg=alpha("skyblue", 0.5),
	 lwd=2.5,
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
	)
#text(TestGenesData$Stim4hr_AvgFc,
#     TestGenesData$log2FoldChangeStimUnstim,
#	 labels=ModelDegData$GeneId.x,
#	 cex=0.5, pos=2, offset=0.2, col="red", font=2
#   )
axis(1,c(seq(xmin,xmax,by=0.25)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=0.5)),col="black",lwd=2,cex.axis=1.50)
mtext("wcp log2 fold change",1,line=2.75,cex=1.75)
axis(2,c(seq(ymin,ymax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)),col="black",lwd=2,cex.axis=1.50)
mtext("RNAseq log2 fold change",2,line=2.75,cex=1.75)
dev.off()
##end scatterplot of WCP vs RNAseq for proteins increased in ubiqutination


if ( FALSE )
{
##plot proteomics results vs RNAseq results
#tcrub<-c("Lat","Zap70","Prkcq","Cd3e","Rhoa","Grap","Pi4k2a","Pip5k1a")
#tcrub<-c("Slc11a2", "Slc12a6", "Slc12a7", "Slc15a4", "Slc16a1", "Slc16a10", "Slc16a3", "Slc1a5", "Slc20a1", "Slc25a1", "Slc25a10", "Slc25a11", "Slc25a12", "Slc25a13", "Slc25a17", "Slc25a19", "Slc25a20", "Slc25a22", "Slc25a24", "Slc25a3", "Slc25a32", "Slc25a4", "Slc25a40", "Slc25a42", "Slc25a46", "Slc25a5", "Slc25a51", "Slc27a4", "Slc29a1", "Slc29a1", "Slc2a1", "Slc2a3", "Slc30a5", "Slc30a6", "Slc30a7", "Slc33a1", "Slc35a1", "Slc35a2", "Slc35a3", "Slc35a4", "Slc35b2", "Slc35e1", "Slc35f2", "Slc38a1", "Slc38a10", "Slc38a2", "Slc39a7", "Slc3a2", "Slc43a3", "Slc44a2", "Slc4a1ap", "Slc4a1ap", "Slc4a2", "Slc4a7", "Slc50a1", "Slc6a11", "Slc6a6", "Slc7a1", "Slc7a5", "Slc7a6", "Slc7a6os", "Slc9a3r1", "Slco1a1")
#SlcTcellProts<-c("Slc38a1", "Slc38a2", "Slc6a6", "Slc1a5", "Slc7a5", "Slc7a1")
#SlcUbProts<-c("Slc11a2", "Slc1a5", "Slc2a1", "Slc2a3", "Slc38a2", "Slc3a2", "Slc43a3", "Slc6a6", "Slc7a1")
sub <- c(ProtDownregRnaUpreg,ProtDownregRnaUnch,ProtDownregRnaDownreg,
         ProtUpregRnaUpreg, ProtUpregRnaUnch, ProtUpregRnaDownreg,
		 ProtUnchRnaUpreg, ProtUnchRnaUnch, ProtUnchRnaDownreg)
#SubData<-WCPstimQuant[WCPstimQuant$GeneId.x%in%sub,]
SubData<-WCPstimQuant[row.names(WCPstimQuant)%in%sub,]
#SubDataUbUp<-WCPstimQuant[row.names(WCPstimQuant)%in%PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd,]
#SubDataUbDown<-WCPstimQuant[row.names(WCPstimQuant)%in%PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd,]
#SubData<-WCPstimQuant
#pdf(file.path(".", analysisdir, "ProteomicsRNAseq_Log2FoldChangeComp_UbNonNorm_SlcProts_TEST.pdf"))
pdf(file.path(".", analysisdir, "ProteomicsRNAseq_Log2FoldChangeComp_UbUpregProts_wLabels.pdf"))
par(bty="l")
#xmin <- floor(min(WCPstimQuant$Stim4hr_AvgFc, na.rm=TRUE))
#xmax <- ceiling(max(WCPstimQuant$Stim4hr_AvgFc, na.rm=TRUE))
#ymin <- floor(min(WCPstimQuant$log2FoldChangeStimUnstim, na.rm=TRUE))
#ymax <- ceiling(max(WCPstimQuant$log2FoldChangeStimUnstim, na.rm=TRUE))
xmin <- -2
xmax <- 1.5
ymin <- -4
ymax <- 4
#make axes
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
#add horizontal and vertical lines at 0
par(new=T)
abline(v=0, col="gray70", lwd=2, lty=1)
par(new=T)
abline(h=0, col="gray70", lwd=2, lty=1)
#add vertical lines at x-axis proteomics mean +/- standard deviation
par(new=T)
abline(v=PO972PO1063PO1075_AvgFcMean+PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=1, lty=2)
par(new=T)
abline(v=PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=1, lty=2)
#add horizontal lines at y-axis RNAseq mean +/- standard deviation
par(new=T)
abline(h=RNAseq_FcMean+RNAseq_FcSd, col="gray70", lwd=1, lty=2)
par(new=T)
abline(h=RNAseq_FcMean-RNAseq_FcSd, col="gray70", lwd=1, lty=2)
#plot points
par(new=T)
plot(SubData$Stim4hr_AvgFc,
     SubData$log2FoldChangeStimUnstim,
     pch=19,
	 cex=1,
	 #cex=4*(abs(SubData$UbFcAvg_4hr_NormWcpFc)),
	 col="grey90",
	 #col=ifelse(SubData$GeneId.x=="Pi4k2a", "salmon", ifelse((SubData$GeneId.x=="Zap70" | SubData$GeneId.x=="Lat"), "gold", "dodgerblue")),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
)
#par(new=T)
#abline(fit <- lm(SubData$log2FoldChangeStimUnstim ~ SubData$Stim4hr_AvgFc), col="dodgerblue", lwd=2, lty=3) # regression line (y~x)
#legend("topright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)), text.col="grey50")
set <- list(ProtUnchRnaUnchUbUpreg, ProtUnchRnaDownregUbUpreg, ProtDownregRnaUpregUbUpreg, ProtDownregRnaUnchUbUpreg, ProtDownregRnaDownregUbUpreg, ProtUpregRnaUpregUbUpreg, ProtUpregRnaUnchUbUpreg, ProtUpregRnaDownregUbUpreg, ProtUnchRnaUpregUbUpreg)
col <- c("grey50", "dodgerblue", "yellow", "blue", "purple", "yellow", "red", "yellow", "darkcyan")
for ( i in 1:length(set) ){
	currset<-set[[i]]
	currcol<-col[i]
	TmpData<-SubData[row.names(SubData)%in%currset,]
	if ( dim(TmpData)[1] > 0 ){
		par(new=T)
		plot(TmpData$Stim4hr_AvgFc,
             TmpData$log2FoldChangeStimUnstim,
			 pch=1,
			 cex=3*(abs(TmpData$UbFcAvg_4hr_NormWcpFc)),
			 col=currcol,
			 lwd=2.5,
			 xlab="",ylab="",
			 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
			 xaxt='n', yaxt='n'
			)
		text(TmpData$Stim4hr_AvgFc,
		     TmpData$log2FoldChangeStimUnstim,
			 labels=TmpData$GeneId.x,
			 cex=0.5, pos=2, offset=0.2, col=currcol, font=2
			)
	}
}
#par(new=T)
#for ( i in 1:length(PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd) )
#{
#	currgene<-PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd[i]
#	currgeneid<-SubData[currgene, "GeneId.x"]
#	text(SubData[currgene, "Stim4hr_AvgFc"],
#	     SubData[currgene, "log2FoldChangeStimUnstim"],
#		 labels=currgeneid,
#		 cex=0.25, pos=2, offset=0.2, col="blue", font=2
#		)
#}

#plot(SubDataUbDown$Stim4hr_AvgFc,
#     SubDataUbDown$log2FoldChangeStimUnstim,
#	 pch=19,
#	 cex=1,
#	 cex=4*(abs(SubDataUbDown$UbFcAvg_4hr_NormWcpFc)),
#	 col="salmon",
##	 #col=ifelse(SubData$GeneId.x=="Pi4k2a", "salmon", ifelse((SubData$GeneId.x=="Zap70" | SubData$GeneId.x=="Lat"), "gold", "dodgerblue")),
#	 xlab="",ylab="",
#	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
#	 xaxt='n', yaxt='n'
#)
#par(new=T)
#abline(fit <- lm(SubDataII$log2FoldChangeStimUnstim ~ SubDataII$Stim4hr_AvgFc), col="olivedrab3", lwd=2, lty=3) # regression line (y~x)
#legend("topright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)), text.col="grey50")
#add points for genes in vector
#par(new=T)
#for ( i in 1:length(SlcTcellProts) ){
#for ( i in 1:length(sub) ){
#	currgene<-sub[i]
#	#currgene<-SlcTcellProts[i]
#	#currgene<-'Grap'
#	points(SubData[which(SubData$GeneId.x==currgene), "Stim4hr_AvgFc"],	
#	       SubData[which(SubData$GeneId.x==currgene), "log2FoldChangeStimUnstim"],
#	       pch=19,
#		   cex=1,
#		   #cex=4*(abs(SubData$UbFcAvg_4hr_NormWcpFc)),
#		   col="yellow"
#	      )
#	 text(SubData[which(SubData$GeneId.x==currgene), "Stim4hr_AvgFc"],
#	      SubData[which(SubData$GeneId.x==currgene), "log2FoldChangeStimUnstim"],
#		  labels=currgene,
#		  cex=0.25, pos=2, offset=0.2, col="yellow", font=2
#		 )
#	text(SubData[currgene, "Stim4hr_AvgFc"],
#	     SubData[currgene, "log2FoldChangeStimUnstim"],
#		 labels=currgene,
#		 cex=0.25, pos=2, offset=0.2, col="yellow", font=2
#		)
#
#}
#par(new=T)
#UbDiffProts<-union(PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd, PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd)
#UbIntProts<-intersect(UbDiffProts,PO972kggPO1063HLMkgg_UbProts_Int)
#SubDataIII<-WCPstimQuant[row.names(WCPstimQuant)%in%UbIntProts,]
#plot(SubDataIII$Stim4hr_AvgFc,
 #    SubDataIII$log2FoldChangeStimUnstim,
#     pch=1,
#	 cex=1,
#	 #cex=4*(abs(SubData$UbFcAvg_4hr_NormWcpFc)),
#	 col="red",
#	 #col=ifelse(SubData$GeneId.x=="Pi4k2a", "salmon", ifelse((SubData$GeneId.x=="Zap70" | SubData$GeneId.x=="Lat"), "gold", "dodgerblue")),
#	 xlab="",ylab="",
#	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
#	 xaxt='n', yaxt='n'
#)
#par(new=T)
#abline(fit <- lm(SubDataIII$log2FoldChangeStimUnstim ~ SubDataIII$Stim4hr_AvgFc), col="red", lwd=2, lty=3) # regression line (y~x)
#legend("topright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)), text.col="grey50")
#par(new=T)
#points(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Upreg_FcGt0,"Stim4hr_AvgFc"],
#       WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Upreg_FcGt0,"log2FoldChangeStimUnstim"],
#points(SubData$Stim4hr_AvgFc,
#       SubData$log2FoldChangeStimUnstim,
#	   pch=19,
#	   cex=4*(abs(SubData$UbFcAvg_4hr_NormWcpFc)),
#	   #cex=1,
#	   col="olivedrab3"
#	  )
#par(new=T)
#abline(fit <- lm(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Upreg_FcGt0,"log2FoldChangeStimUnstim"] ~ WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Upreg_FcGt0,"Stim4hr_AvgFc"]), col="olivedrab3", lwd=2, lty=3) # regression line (y~x)
#legend("topright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)), text.col="olivedrab3")
#par(new=T)
#points(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd,"Stim4hr_AvgFc"],
#       WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd,"log2FoldChangeStimUnstim"],
#	   pch=19,
#	   cex=1,
#	   col="salmon"
#	  )
#par(new=T)
#abline(fit <- lm(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd,"log2FoldChangeStimUnstim"] ~ WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd,"Stim4hr_AvgFc"]), col="salmon", lwd=2, lty=3) # regression line (y~x)
#legend("topright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)), text.col="salmon")
#par(new=T)
#text(WCPstimQuant[which(WCPstimQuant$GeneId.x=='Grap'), "Stim4hr_AvgFc"], 
#     WCPstimQuant[which(WCPstimQuant$GeneId.x=='Grap'), "log2FoldChangeStimUnstim"],
#	 labels="Grap",
#	 cex=0.25, pos=2, offset=0.2, col="black", font=2
#	)
#points(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanSd, "Stim4hr_AvgFc"],
#       WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanSd, "log2FoldChangeStimUnstim"],
#	   pch=1,
#	   cex=0.75,
#	   col="dodgerblue"
#)
#points(WCPstimQuant[ProtDiffRnaDiffUb, "Stim4hr_AvgFc"],
#       WCPstimQuant[ProtDiffRnaDiffUb, "log2FoldChangeStimUnstim"],
#	   pch=19,
#	   cex=0.75,
#	   col="salmon"
#)
axis(1,c(seq(xmin,xmax,by=0.25)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=0.5)),col="black",cex.axis=0.80)
mtext("wcp log2 fold change",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)),col="black",cex.axis=0.80)
mtext("RNAseq log2 fold change",2,line=2.75,cex=1.5)
dev.off()
}
##end proteomics vs RNAseq results plot

if ( FALSE )
{
##plot proteomics results vs RNAseq results
sub <- c(ProtDownregRnaUpreg,ProtDownregRnaUnch,ProtDownregRnaDownreg,
         ProtUpregRnaUpreg, ProtUpregRnaUnch, ProtUpregRnaDownreg,
		 ProtUnchRnaUpreg, ProtUnchRnaUnch, ProtUnchRnaDownreg)
SubData<-WCPstimQuant[row.names(WCPstimQuant)%in%sub,]
pdf(file.path(".", analysisdir, "ProteomicsRNAseq_Log2FoldChangeComp_w4hrUbProts.pdf"))
par(bty="l")
#xmin <- floor(min(WCPstimQuant$Stim4hr_AvgFc, na.rm=TRUE))
#xmax <- ceiling(max(WCPstimQuant$Stim4hr_AvgFc, na.rm=TRUE))
#ymin <- floor(min(WCPstimQuant$log2FoldChangeStimUnstim, na.rm=TRUE))
#ymax <- ceiling(max(WCPstimQuant$log2FoldChangeStimUnstim, na.rm=TRUE))
xmin <- -1
xmax <- 1
ymin <- -3
ymax <- 3
#make axes
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
#add horizontal and vertical lines at 0
par(new=T)
abline(v=0, col="gray70", lwd=2, lty=1)
par(new=T)
abline(h=0, col="gray70", lwd=2, lty=1)
#add vertical lines at x-axis proteomics mean +/- standard deviation
par(new=T)
abline(v=PO972PO1063PO1075_AvgFcMean+PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=1, lty=2)
par(new=T)
abline(v=PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=1, lty=2)
#add horizontal lines at y-axis RNAseq mean +/- standard deviation
par(new=T)
abline(h=RNAseq_FcMean+RNAseq_FcSd, col="gray70", lwd=1, lty=2)
par(new=T)
abline(h=RNAseq_FcMean-RNAseq_FcSd, col="gray70", lwd=1, lty=2)
#plot points
par(new=T)
plot(SubData$Stim4hr_AvgFc,
     SubData$log2FoldChangeStimUnstim,
     pch=19,
	 cex=1,
	 #cex=4*(abs(SubData$UbFcAvg_4hr_NormWcpFc)),
	 col="grey90",
	 #col=ifelse(SubData$GeneId.x=="Pi4k2a", "salmon", ifelse((SubData$GeneId.x=="Zap70" | SubData$GeneId.x=="Lat"), "gold", "dodgerblue")),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
)
par(new=T)
Ub4hr<-c("A0A0R4J1L0", "P43275", "Q9D0M0", "Q9JKF1", "D3Z0M9", "E9PW15", "F6SMY7", "O08914", "O70201", "O70400", "P13020", "P47968", "P70353", "Q4FJQ0", "Q4FJQ4", "Q4FK35", "Q545A2", "Q5SVP0", "Q5XJE5", "Q62120", "Q6QD59", "Q6ZWQ9", "Q6ZWV3", "Q80Y84", "Q8BK67", "Q8BR63", "Q91W67", "Q99LF4", "Q9CZN7", "Q9QUR7", "Q9WTV7")
for ( i in 1:length(Ub4hr) ){
	currprot<-Ub4hr[i]
	currgene<-SubData[currprot, "GeneId.x"]
	points(SubData[currprot, "Stim4hr_AvgFc"],
           SubData[currprot, "log2FoldChangeStimUnstim"],
		   pch=19,
		   #cex=3*(abs(TmpData$UbFcAvg_4hr_NormWcpFc)),
		   cex=1.5,
		   col="red",
		   #lwd=2.5,
		   xlab="",ylab="",
		   xlim=c(xmin,xmax), ylim=c(ymin,ymax),
		   xaxt='n', yaxt='n'
		  )
	text(SubData[currprot, "Stim4hr_AvgFc"],
	     SubData[currprot, "log2FoldChangeStimUnstim"],
		 labels=currgene,
		 cex=0.5, pos=2, offset=0.2, col="blue", font=2
	    )
}
axis(1,c(seq(xmin,xmax,by=0.25)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=0.5)),col="black",cex.axis=0.80)
mtext("wcp log2 fold change",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)),col="black",cex.axis=0.80)
mtext("RNAseq log2 fold change",2,line=2.75,cex=1.5)
dev.off()
}
##end proteomics vs RNAseq results plot


##plot proteomics results vs RNAseq results
##plot protein abundances for RNA counts
##determine if there is an RNA count "threshold"
##to see protein expression
pdf(file.path(".", analysisdir, "ProteomicsRNAseq_RnaMeanCountProtZscoreAbundance.pdf"))
par(bty="l")
xmin <- floor(min(RnaData$log10baseMean.x, na.rm=TRUE))
xmax <- ceiling(max(RnaData$log10baseMean.x, na.rm=TRUE))
ymin <- floor(min(WCPstimQuant$AbundanceInt_zscore, na.rm=TRUE))
ymax <- ceiling(max(WCPstimQuant$AbundanceInt_zscore, na.rm=TRUE))
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
box(bty="l", lwd=3)
par(new=T)
plot(WCPstimQuant$log10baseMean, WCPstimQuant$AbundanceInt_zscore,
     pch=19,
	 cex=0.5,
	 col="grey50",
	 #col=ifelse((!is.na(WCPstimQuant$Ub_0hrUnique) | !is.na(WCPstimQuant$Ub_4hrUnique) | !is.na(WCPstimQuant$UbFcAvg_4hr_NormWcpFc)), "cyan", "grey50"),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
)
par(new=T)
plot(UbProtsData$log10baseMean, UbProtsData$AbundanceInt_zscore,
     pch=19,
	 cex=0.5,
	 col="cyan",
	 #col=ifelse((!is.na(WCPstimQuant$Ub_0hrUnique) | !is.na(WCPstimQuant$Ub_4hrUnique) | !is.na(WCPstimQuant$UbFcAvg_4hr_NormWcpFc)), "cyan", "grey50"),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
)


#for ( i in 1:length(PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts) )
#{
#	currgene<-PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts[i]
#	points(WCPstimQuant[currgene, "log10baseMean"], 
#	       WCPstimQuant[currgene, "AbundanceInt_zscore"],
#		   pch=19,
#		   cex=1,
#		   col="cyan3"
#		  )
#}
axis(1,c(seq(xmin,xmax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=1)),col="black",lwd=2,cex.axis=1.5)
mtext("log10 RNAseq count",1,line=2.75,cex=1.75)
axis(2,c(seq(ymin,ymax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)),col="black",lwd=2,cex.axis=1.5)
mtext("wcp protein abundance",2,line=2.75,cex=1.75)
dev.off()


##plot proteomics results vs RNAseq results
tcrubupreg <- c("Cd3e", "Grap", "Zap70", "Cd247", "Rhoa", "Cd3g", "Prkcq", "Lat", "Cdc42", "Rac1", "Ube2n", "Nfkb1")
TCRUbId <- WCPstimQuant[WCPstimQuant$GeneId.x%in%tcrubupreg,]
TCRUbId_ordered <- TCRUbId[rev(order(TCRUbId$UbFcAvg_4hr_NormWcpFc)),] 
head(TCRUbId_ordered)
pdf(file.path(".", analysisdir, "ProteomicsRNAseq_Log2FoldChangeComp_UbNormToWcp_TCRProtUpregUb.pdf"))
par(bty="l")
xmin <- -0.5
xmax <- 1
ymin <- -2
ymax <- 2.5
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
par(new=T)
abline(v=0, col="gray70", lwd=2, lty=1)
par(new=T)
abline(h=0, col="gray70", lwd=2, lty=1)
#par(new=T)
#abline(v=PO972PO1063PO1075_AvgFcMean+PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=1, lty=2)
#par(new=T)
#abline(v=PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=1, lty=2)
#par(new=T)
#abline(h=RNAseq_FcMean+RNAseq_FcSd, col="gray70", lwd=1, lty=2)
#par(new=T)
#abline(h=RNAseq_FcMean-RNAseq_FcSd, col="gray70", lwd=1, lty=2)
par(new=T)
plot(TCRUbId_ordered$Stim4hr_AvgFc, TCRUbId_ordered$log2FoldChangeStimUnstim,
     pch=21,
	 cex=4*(abs(TCRUbId_ordered$UbFcAvg_4hr_NormWcpFc)),
	 #col="red",
	 col="grey40",
	 bg=ifelse((TCRUbId_ordered$GeneId.x=="Zap70" | TCRUbId_ordered$GeneId.x=="Lat" | TCRUbId_ordered$GeneId.x=="Nfkb1" | TCRUbId_ordered$GeneId.x=="Ube2n"), "grey50", "red"),
	 lwd=2.5,
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
)
#par(new=T)
#points(1,-2,
#       pch=19,
#	   cex=4*(abs(1)),,
#	   col="grey50"
#	  )
axis(1,c(seq(xmin,xmax,by=0.25)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=0.5)),col="black",lwd=2,cex.axis=1.5)
mtext("WCP log2 fold change",1,line=2.75,cex=1.75)
axis(2,c(seq(ymin,ymax,by=0.25)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=0.5)),col="black",lwd=2,cex.axis=1.5)
mtext("RNAseq log2 fold change",2,line=2.75,cex=1.75)
dev.off()


##plot proteomics results vs RNAseq results and Ubiquitination for specific pathways
#pathwayub<-c("Sqstm1", "Traf1", "Icam1", "Mcl1", "Slc2a3", "Eif1", "Tnip2", "Smad3", "Ddx58", "Pfkfb3", "Atp2b1", "Nfkb1", "Cd44", "Rela", "Tap1")
#pdf(file.path(".", analysisdir, "ProteomicsRNAseq_Log2FoldChangeComp_UbNormToWcp_PathwayUb_DNARepair.pdf"))
#pathwayub<-c("Ddb1", "Impdh2", "Pcna", "Ncbp2", "Pola1", "Gtf2b", "Hprt1", "Ssrp1", "Pom121", "Srsf6", "Polr2a", "Gtf2h1", "Arl6ip1", "Dut", "Nudt21", "Aprt", "Cetn2", "Itpa", "Lig1", "Nme1", "Hcls1", "Dctn4", "Rbx1")
#pdf(file.path(".", analysisdir, "ProteomicsRNAseq_Log2FoldChangeComp_UbNormToWcp_PathwayUb_DNARepair.pdf"))
#pathwayub<-c("Sfn", "Pcna", "Bax", "Pom121", "Cyfip2", "Hist1h1c", "S100a10", "Rb1", "Gnb2l1", "Rps12", "Txnip", "Rpl18", "Vamp8", "Dgka", "Slc3a2", "Tap1")
#pdf(file.path(".", analysisdir, "ProteomicsRNAseq_Log2FoldChangeComp_UbNormToWcp_PathwayUb_p53.pdf"))
#pathwayub<-c("Icam1", "Ccr7", "Sri", "Il4r", "Tapbp", "Gnai3", "Rhog", "Itgb3", "Atp2b1", "Nfkb1", "Nmi", "Slc7a1", "Lck", "Rela", "Sema4d")
#pdf(file.path(".", analysisdir, "ProteomicsRNAseq_Log2FoldChangeComp_UbNormToWcp_PathwayUb_InflammatoryResponse.pdf"))
#pathway name
#pathname<-'HALLMARK_TNFA_SIGNALING_VIA_NFKB'
#pathname<-'HALLMARK_E2F_TARGETS'
#pathname<-'HALLMARK_G2M_CHECKPOINT'
#pathname<-'HALLMARK_MTORC1_SIGNALING'
pathname<-'HALLMARK_E2FG2Moverlap_TARGETS'
#read in file
pathfilename<-paste(pathname,'_EnrichedUbiquitinatedGenes.txt',sep="")
pathvec<-scan(pathfilename,what='character')
#adjust gene id strings
pathwayub<-vector()
for ( i in 1:length(pathvec) ){
	ubgene<-pathvec[i]
	ubgene_lc<-tolower(ubgene)
	ubgene_formatted<-paste(toupper(substr(ubgene_lc, 1, 1)), substr(ubgene_lc, 2, nchar(ubgene_lc)), sep="")
	if(ubgene_formatted=='Gapdh'){ ubgene_formatted<-'GAPDH' }
	pathwayub<-c(pathwayub,ubgene_formatted)
}
print(pathwayub)
#assign file name based on path name
pathwayfigfile<-paste(pathname,'ProteomicsRNAseq_Log2FoldChangeComp_UbNormToWcp_PathwayUb.pdf',sep="")
pdf(file.path(".", analysisdir, pathwayfigfile))
PathUbId<-WCPstimQuant[WCPstimQuant$GeneId.x%in%pathwayub,]
par(bty="l")
#xmin <- -10
#xmax <- 10
#ymin <- -10
#ymax <- 10
xmin <- floor(min(PathUbId$Stim4hr_AvgFc, na.rm=TRUE))
xmax <- ceiling(max(PathUbId$Stim4hr_AvgFc, na.rm=TRUE))
ymin <- floor(min(PathUbId$log2FoldChangeStimUnstim, na.rm=TRUE))
ymax <- ceiling(max(PathUbId$log2FoldChangeStimUnstim, na.rm=TRUE))
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
#plot(WCPstimQuant[which(WCPstimQuant$GeneId.x=tcrub), "Stim4hr_AvgFc"], WCPstimQuant[which(WCPstimQuant$GeneId.x=tcrub), "log2FoldChangeStimUnstim"],
#plot(TCRUbId$Stim4hr_AvgFc, TCRUbId$log2FoldChangeStimUnstim,
#     pch=19,
#	 #cex=1.5,
#	 cex=2*(2^(abs(TCRUbId$UbFcAvg))),
#	 col="red",
#	 xlab="",ylab="",
#	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
#	 xaxt='n', yaxt='n'
#)
par(new=T)
abline(v=0, col="gray70", lwd=2, lty=1)
par(new=T)
abline(h=0, col="gray70", lwd=2, lty=1)
par(new=T)
abline(v=PO972PO1063PO1075_AvgFcMean+PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=1, lty=2)
par(new=T)
abline(v=PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=1, lty=2)
par(new=T)
abline(h=RNAseq_FcMean+RNAseq_FcSd, col="gray70", lwd=1, lty=2)
par(new=T)
abline(h=RNAseq_FcMean-RNAseq_FcSd, col="gray70", lwd=1, lty=2)
par(new=T)
plot(PathUbId$Stim4hr_AvgFc, PathUbId$log2FoldChangeStimUnstim,
	 
	 ######protein Ub downregulation by KeGG avg fold change +/- threshold
	 #####PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanThreshold<-row.names(subset(WCPstimQuant, WCPstimQuant$UbFcAvg_4hr_NormWcpFc <= -1*ubdeltathreshold))
	 

	 pch=ifelse((is.na(PathUbId$Ub_0hrUnique) & is.na(PathUbId$Ub_4hrUnique)), 16, 22),
	 cex=ifelse((is.na(PathUbId$Ub_0hrUnique) & is.na(PathUbId$Ub_4hrUnique)), ((5*(abs(PathUbId$UbFcAvg_4hr_NormWcpFc)))+1), 3),
	 lwd=2.5,
	 #col=ifelse(is.na(PathUbId$UbFcAvg_4hr_NormWcpFc), ifelse(!is.na(PathUbId$Ub_4hrUnique), "red", "dodgerblue"), ifelse(PathUbId$UbFcAvg_4hr_NormWcpFc>0, "red", "dodgerblue")),
	 col=ifelse(is.na(PathUbId$UbFcAvg_4hr_NormWcpFc), ifelse(!is.na(PathUbId$Ub_4hrUnique), "red", "dodgerblue"), ifelse(PathUbId$UbFcAvg_4hr_NormWcpFc > ubdeltathreshold, "red", ifelse(PathUbId$UbFcAvg_4hr_NormWcpFc <= -1*ubdeltathreshold, "dodgerblue","grey50"))),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
)
#par(new=T)
#points(1,-2,
#       pch=19,
#	   cex=4*(abs(1)),,
#	   col="grey50"
#	  )


for ( i in 1:length(pathwayub) )
{
	currgene<-pathwayub[i]
	x<-PathUbId[which(PathUbId$GeneId.x==currgene), "Stim4hr_AvgFc"]
	y<-PathUbId[which(PathUbId$GeneId.x==currgene), "log2FoldChangeStimUnstim"]
	print(currgene)
	print("x")
	print(x)
	print("y")
	print(y)
	text(PathUbId[which(PathUbId$GeneId.x==currgene), "Stim4hr_AvgFc"],
	     PathUbId[which(PathUbId$GeneId.x==currgene), "log2FoldChangeStimUnstim"],
		 labels=currgene,
		 cex=0.5, pos=2, offset=0.2, col="grey30", font=2
		)
	#text(x, y, labels=currgene, cex=0.25, pos=2, offset=0.2, col="black", font=2)
}
axis(1,c(seq(xmin,xmax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=1)),col="black",cex.axis=1.0)
mtext("WCP log2 fold change",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)),col="black",cex.axis=1.0)
mtext("RNAseq log2 fold change",2,line=2.75,cex=1.5)
dev.off()


##plot proteomics results vs RNAseq results
##plot all identified proteins for wcp and RNA
##include all Ub proteins identified by color
#sub <- c(ProtDownregRnaUpreg,ProtDownregRnaUnch)
SubDataUbUp<-WCPstimQuant[row.names(WCPstimQuant)%in%PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd,]
SubDataUbDown<-WCPstimQuant[row.names(WCPstimQuant)%in%PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd,]
SubDataUbAll<-WCPstimQuant[row.names(WCPstimQuant)%in%PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts,]
pdf(file.path(".", analysisdir, "ProteomicsRNAseq_Log2FoldChangeComp_AllUbProts.pdf"))
par(bty="l")
#xmin <- floor(min(WCPstimQuant$Stim4hr_AvgFc, na.rm=TRUE))
#xmax <- ceiling(max(WCPstimQuant$Stim4hr_AvgFc, na.rm=TRUE))
#ymin <- floor(min(WCPstimQuant$log2FoldChangeStimUnstim, na.rm=TRUE))
#ymax <- ceiling(max(WCPstimQuant$log2FoldChangeStimUnstim, na.rm=TRUE))
xmin <- -3
xmax <- 5
ymin <- -5
ymax <- 7
#make axes
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
#add horizontal and vertical lines at 0
par(new=T)
abline(v=0, col="gray70", lwd=2, lty=1)
par(new=T)
abline(h=0, col="gray70", lwd=2, lty=1)
#add vertical lines at x-axis proteomics mean +/- standard deviation
par(new=T)
abline(v=PO972PO1063PO1075_AvgFcMean+PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=1, lty=2)
par(new=T)
abline(v=PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=1, lty=2)
#add horizontal lines at y-axis RNAseq mean +/- standard deviation
par(new=T)
abline(h=RNAseq_FcMean+RNAseq_FcSd, col="gray70", lwd=1, lty=2)
par(new=T)
abline(h=RNAseq_FcMean-RNAseq_FcSd, col="gray70", lwd=1, lty=2)
#plot points
par(new=T)
plot(WCPstimQuant$Stim4hr_AvgFc,
     WCPstimQuant$log2FoldChangeStimUnstim,
	 pch=19,
	 cex=0.75,
	 col="grey80",
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
)
par(new=T)
plot(SubDataUbAll$Stim4hr_AvgFc,
     SubDataUbAll$log2FoldChangeStimUnstim,
     pch=19,
	 cex=1,
	 #cex=4*(abs(SubData$UbFcAvg_4hr_NormWcpFc)),
	 col="dodgerblue",
	 #col=ifelse(SubData$GeneId.x=="Pi4k2a", "salmon", ifelse((SubData$GeneId.x=="Zap70" | SubData$GeneId.x=="Lat"), "gold", "dodgerblue")),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
)
par(new=T)
#abline(fit <- lm(SubDataUbAll$log2FoldChangeStimUnstim ~ SubDataUbAll$Stim4hr_AvgFc), col="dodgerblue", lwd=2, lty=3) # regression line (y~x)
#legend("topright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)), text.col="grey50")
par(new=T)
plot(SubDataUbUp$Stim4hr_AvgFc,
     SubDataUbUp$log2FoldChangeStimUnstim,
     pch=19,
	 cex=1,
#	 #cex=4*(abs(SubData$UbFcAvg_4hr_NormWcpFc)),
	 col="olivedrab3",
#	 #col=ifelse(SubData$GeneId.x=="Pi4k2a", "salmon", ifelse((SubData$GeneId.x=="Zap70" | SubData$GeneId.x=="Lat"), "gold", "dodgerblue")),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
)
par(new=T)
#abline(fit <- lm(SubDataUbUp$log2FoldChangeStimUnstim ~ SubDataUbUp$Stim4hr_AvgFc), col="olivedrab3", lwd=2, lty=3) # regression line (y~x)
#legend("topright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)), text.col="grey50")
par(new=T)
plot(SubDataUbDown$Stim4hr_AvgFc,
     SubDataUbDown$log2FoldChangeStimUnstim,
     pch=19,
	 cex=1,
#	 #cex=4*(abs(SubData$UbFcAvg_4hr_NormWcpFc)),
	 col="salmon",
#	 #col=ifelse(SubData$GeneId.x=="Pi4k2a", "salmon", ifelse((SubData$GeneId.x=="Zap70" | SubData$GeneId.x=="Lat"), "gold", "dodgerblue")),
	 xlab="",ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
)
par(new=T)
#abline(fit <- lm(SubDataUbDown$log2FoldChangeStimUnstim ~ SubDataUbDown$Stim4hr_AvgFc), col="salmon", lwd=2, lty=3) # regression line (y~x)
#legend("topright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)), text.col="grey50")
#add points for genes in vector
#par(new=T)
#for ( i in 1:length(SlcTcellProts) ){
#for ( i in 1:length(sub) ){
#	currgene<-sub[i]
#	#currgene<-SlcTcellProts[i]
#	#currgene<-'Grap'
#	points(SubData[which(SubData$GeneId.x==currgene), "Stim4hr_AvgFc"],	
#	       SubData[which(SubData$GeneId.x==currgene), "log2FoldChangeStimUnstim"],
#	       pch=19,
#		   cex=1,
#		   #cex=4*(abs(SubData$UbFcAvg_4hr_NormWcpFc)),
#		   col="yellow"
#	      )
#	 text(SubData[which(SubData$GeneId.x==currgene), "Stim4hr_AvgFc"],
#	      SubData[which(SubData$GeneId.x==currgene), "log2FoldChangeStimUnstim"],
#		  labels=currgene,
#		  cex=0.25, pos=2, offset=0.2, col="yellow", font=2
#		 )
#	text(SubData[currgene, "Stim4hr_AvgFc"],
#	     SubData[currgene, "log2FoldChangeStimUnstim"],
#		 labels=currgene,
#		 cex=0.25, pos=2, offset=0.2, col="yellow", font=2
#		)
#
#}
#par(new=T)
#UbDiffProts<-union(PO972kggPO1063HLMkgg_UbNormFc_Upreg_AvgFcMeanSd, PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd)
#UbIntProts<-intersect(UbDiffProts,PO972kggPO1063HLMkgg_UbProts_Int)
#SubDataIII<-WCPstimQuant[row.names(WCPstimQuant)%in%UbIntProts,]
#plot(SubDataIII$Stim4hr_AvgFc,
 #    SubDataIII$log2FoldChangeStimUnstim,
#     pch=1,
#	 cex=1,
#	 #cex=4*(abs(SubData$UbFcAvg_4hr_NormWcpFc)),
#	 col="red",
#	 #col=ifelse(SubData$GeneId.x=="Pi4k2a", "salmon", ifelse((SubData$GeneId.x=="Zap70" | SubData$GeneId.x=="Lat"), "gold", "dodgerblue")),
#	 xlab="",ylab="",
#	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
#	 xaxt='n', yaxt='n'
#)
#par(new=T)
#abline(fit <- lm(SubDataIII$log2FoldChangeStimUnstim ~ SubDataIII$Stim4hr_AvgFc), col="red", lwd=2, lty=3) # regression line (y~x)
#legend("topright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)), text.col="grey50")
#par(new=T)
#points(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Upreg_FcGt0,"Stim4hr_AvgFc"],
#       WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Upreg_FcGt0,"log2FoldChangeStimUnstim"],
#points(SubData$Stim4hr_AvgFc,
#       SubData$log2FoldChangeStimUnstim,
#	   pch=19,
#	   cex=4*(abs(SubData$UbFcAvg_4hr_NormWcpFc)),
#	   #cex=1,
#	   col="olivedrab3"
#	  )
#par(new=T)
#abline(fit <- lm(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Upreg_FcGt0,"log2FoldChangeStimUnstim"] ~ WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Upreg_FcGt0,"Stim4hr_AvgFc"]), col="olivedrab3", lwd=2, lty=3) # regression line (y~x)
#legend("topright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)), text.col="olivedrab3")
#par(new=T)
#points(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd,"Stim4hr_AvgFc"],
#       WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd,"log2FoldChangeStimUnstim"],
#	   pch=19,
#	   cex=1,
#	   col="salmon"
#	  )
#par(new=T)
#abline(fit <- lm(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd,"log2FoldChangeStimUnstim"] ~ WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Downreg_AvgFcMeanSd,"Stim4hr_AvgFc"]), col="salmon", lwd=2, lty=3) # regression line (y~x)
#legend("topright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)), text.col="salmon")
#par(new=T)
#text(WCPstimQuant[which(WCPstimQuant$GeneId.x=='Grap'), "Stim4hr_AvgFc"], 
#     WCPstimQuant[which(WCPstimQuant$GeneId.x=='Grap'), "log2FoldChangeStimUnstim"],
#	 labels="Grap",
#	 cex=0.25, pos=2, offset=0.2, col="black", font=2
#	)
#points(WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanSd, "Stim4hr_AvgFc"],
#       WCPstimQuant[PO972kggPO1063HLMkgg_UbNormFc_Unch_AvgFcMeanSd, "log2FoldChangeStimUnstim"],
#	   pch=1,
#	   cex=0.75,
#	   col="dodgerblue"
#)
#points(WCPstimQuant[ProtDiffRnaDiffUb, "Stim4hr_AvgFc"],
#       WCPstimQuant[ProtDiffRnaDiffUb, "log2FoldChangeStimUnstim"],
#	   pch=19,
#	   cex=0.75,
#	   col="salmon"
#)
axis(1,c(seq(xmin,xmax,by=0.25)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=0.5)),col="black",cex.axis=0.80)
mtext("wcp log2 fold change",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)),col="black",cex.axis=0.80)
mtext("RNAseq log2 fold change",2,line=2.75,cex=1.5)
dev.off()
##end proteomics vs RNAseq results plot


##plot proteomics results vs RNAseq results (only for proteins identified with a KeGG)
pdf(file.path(".", analysisdir, "ProteomicsRNAseq_Log2FoldChangeComp_OnlyUbProts.pdf"))
par(bty="l")
xmin <- -6
xmax <- 4
ymin <- -5
ymax <- 7
plot(WCPstimQuant[PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts, "Stim4hr_AvgFc"], WCPstimQuant[PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts, "log2FoldChangeStimUnstim"],
     pch=19,
	 cex=0.75,
	 #####col="gray70",
	 col=ifelse((WCPstimQuant[PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts, "Stim4hr_AvgFc"]>(PO972PO1063PO1075_AvgFcMean+PO972PO1063PO1075_AvgFcSd) |
	             WCPstimQuant[PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts, "Stim4hr_AvgFc"]<(PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd) |
				 WCPstimQuant[PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts, "log2FoldChangeStimUnstim"]>(RNAseq_FcMean+RNAseq_FcSd) |
				 WCPstimQuant[PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts, "log2FoldChangeStimUnstim"]<(RNAseq_FcMean-RNAseq_FcSd)),
				 "darkorange", "gold"
				),
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xlab="", ylab="",
	 xaxt='n', yaxt='n'
)    
par(new=T)
abline(v=0, col="gray70", lwd=2, lty=1)
par(new=T)
abline(h=0, col="gray70", lwd=2, lty=1)
par(new=T)
abline(v=PO972PO1063PO1075_AvgFcMean+PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=2, lty=2)
par(new=T)
abline(v=PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=2, lty=2)
par(new=T)
abline(h=RNAseq_FcMean+RNAseq_FcSd, col="gray70", lwd=2, lty=2)
par(new=T)
abline(h=RNAseq_FcMean-RNAseq_FcSd, col="gray70", lwd=2, lty=2)
axis(1,c(seq(xmin,xmax,by=0.5)),label=F,col="black",cex.axis=1,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=1)),col="black",cex.axis=1)
mtext("proteomics 0-4hr TCR stim log2 fold change",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)),label=F,col="black",cex.axis=1,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)),col="black",cex.axis=1)
mtext("RNAseq 0-4hr TCR stim log2 fold change",2,line=2.75,cex=1.5)
dev.off()


##plot proteomics results vs KeGG quantification results (only for proteins identified with a KeGG)
alltcr<-c("Cd247","Ptprc","Cd3g","Lat","Lck","Itk","Grb2","Sos1","Zap70","Prkcq","Ikbkg","Nfkb1","Cblb","Cd3e","Rhoa","Nras","Kras","Ppp3r1")
upregtcr<-c("Lat","Zap70","Prkcq","Cd3e","Rhoa")
otherupreg<-c("Grap")
pdf(file.path(".", analysisdir, "ProteomicsKeGGUb_TCRproteinUbUpreg_Log2FoldChangeComp.pdf"))
par(bty="l")
#xmin <- -2
#xmax <- 2
#ymin <- -2
#ymax <- 3
xmin <- -6
xmax <- 4
ymin <- -4
ymax <- 4
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
#par(new=T)
#abline(v=0, col="gray70", lwd=2, lty=1)
#par(new=T)
#abline(h=0, col="gray70", lwd=2, lty=1)
#par(new=T)
#abline(a = 0, b = 1, col="gray70", lty=2, lwd=2)
par(new=T)
abline(h=PO972PO1063PO1075_AvgFcMean+PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=2, lty=3)
par(new=T)
abline(h=PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=2, lty=3)
par(new=T)
#abline(v=PO972kggPO1063HLMkgg_UbFc_Mean+PO972kggPO1063HLMkgg_UbFc_Sd, col="gray70", lwd=2, lty=3)
abline(v=PO972kggPO1063HLMkgg_UbNormFc_Mean+PO972kggPO1063HLMkgg_UbNormFc_Sd, col="gray70", lwd=2, lty=3)
par(new=T)
#abline(v=PO972kggPO1063HLMkgg_UbFc_Mean-PO972kggPO1063HLMkgg_UbFc_Sd, col="gray70", lwd=2, lty=3)
abline(v=PO972kggPO1063HLMkgg_UbNormFc_Mean-PO972kggPO1063HLMkgg_UbNormFc_Sd, col="gray70", lwd=2, lty=3)
par(new=T)
plot(WCPstimQuant[PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts, "Stim4hr_AvgFc"], WCPstimQuant[PO972kggPO1063HLMkgg_All0hr4hrFc_Union_UbProts, "UbFcAvg_4hr_NormWcpFc"],
     pch=19,
	 cex=0.75,
	 col="gray50",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xlab="", ylab="",
	 xaxt='n', yaxt='n'
)    
#for ( i in 1:length(alltcr) ){
#	currgene<-alltcr[i]
#	points(WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "Stim4hr_AvgFc"],
#	       WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "UbFcAvg"],
#		   pch=19,
#		   cex=1.5,
#		   col="grey50"
#		  )
#}
for ( i in 1:length(upregtcr) ){
	currgene<-upregtcr[i]
	points(WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "Stim4hr_AvgFc"],
	       WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "UbFcAvg"],
		   pch=19,
		   cex=1.5,
		   col="red"
		  )
}
for ( i in 1:length(otherupreg) ){
	currgene<-otherupreg[i]
	points(WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "Stim4hr_AvgFc"],
	       WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "UbFcAvg"],
		   pch=19,
		   cex=1.5,
		   col="red"
		  )
}
axis(1,c(seq(xmin,xmax,by=0.25)), labels=F, col="black",cex.axis=0.75, tck=-0.01)
axis(1,c(seq(xmin,xmax,by=0.5)), labels=T, col="black",cex.axis=0.8)
mtext("wcp log2 fold change",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.25)), labels=F, col="black",cex.axis=0.75, tck=-0.01)
axis(2,c(seq(ymin,ymax,by=0.5)), labels=T, col="black",cex.axis=0.8)
mtext("Ub log2 fold change",2,line=2.75,cex=1.5)
dev.off()


##plot proteomics results vs KeGG quantification results (only for proteins identified in PO1063 HLM KeGG experiment)
HLMWcpKeGGid<-row.names(subset(WCPstimQuant, (!is.na(WCPstimQuant$HLMStim4hr_fc) & !is.na(WCPstimQuant$intensityweightedmeanratio.PO1063kgg_HLM))))
HLMWcpKeGGid_Num<-length(HLMWcpKeGGid)
print(HLMWcpKeGGid_Num)
#modified scatter plot
alltcr<-c("Cd247","Ptprc","Cd3g","Lat","Lck","Itk","Grb2","Sos1","Zap70","Prkcq","Ikbkg","Nfkb1","Cblb","Cd3e","Rhoa","Nras","Kras","Ppp3r1")
upregtcr<-c("Lat","Zap70","Prkcq","Cd3e","Rhoa")
pdf(file.path(".", analysisdir, "ProteomicsKeGGUb_Log2FoldChangeComp_PO1063HLMidentifiedProteins.pdf"))
par(bty="l")
xmin <- -1
xmax <- 2
ymin <- -3
ymax <- 3
plot(WCPstimQuant[HLMWcpKeGGid,"HLMStim4hr_fc"], WCPstimQuant[HLMWcpKeGGid,"intensityweightedmeanratio.PO1063kgg_HLM"],
	 pch=19,
	 cex=0.75,
	 col="gray70",
	 xlab="", ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
)
par(new=T)
abline(v=0, col="gray70", lwd=2, lty=2)
par(new=T)
abline(h=0, col="gray70", lwd=2, lty=2)
for ( i in 1:length(alltcr) ){
	currgene<-alltcr[i]
	points(WCPstimQuant[which(WCPstimQuant$GeneId==currgene), "HLMStim4hr_fc"],
	       WCPstimQuant[which(WCPstimQuant$GeneId==currgene), "intensityweightedmeanratio.PO1063kgg_HLM"],
		   pch=19,
		   cex=1.25,
		   col="salmon"
		  )
}
for ( i in 1:length(upregtcr) ){
	currgene<-upregtcr[i]
	points(WCPstimQuant[which(WCPstimQuant$GeneId==currgene), "HLMStim4hr_fc"],
	       WCPstimQuant[which(WCPstimQuant$GeneId==currgene), "intensityweightedmeanratio.PO1063kgg_HLM"],
		   pch=19,
		   cex=1.25,
		   col="red"
		  )
}
axis(1,c(seq(xmin,xmax,by=0.25)), labels=F, col="black",cex.axis=0.75,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=1)), labels=T, col="black",cex.axis=1)
mtext("wcp log2 fold change",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.25)), labels=F, col="black",cex.axis=0.75,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)), labels=T, col="black",cex.axis=1)
mtext("Ub log2 fold change",2,line=2.75,cex=1.5)
dev.off()

##plot proteomics results vs KeGG quantification results (only for proteins identified in PO1063 HLM KeGG experiment)
#modified scatter plot
alltcr<-c("Cd247","Ptprc","Cd3g","Lat","Lck","Itk","Grb2","Sos1","Zap70","Prkcq","Ikbkg","Nfkb1","Cblb","Cd3e","Rhoa","Nras","Kras","Ppp3r1")
upregtcr<-c("Lat","Zap70","Prkcq","Cd3e","Rhoa")
pdf(file.path(".", analysisdir, "ProteomicsKeGGUb_Log2FoldChangeComp_PO1063HLMidentifiedProteins.pdf"))
par(bty="l")
xmin <- -1
xmax <- 2
ymin <- -3
ymax <- 3
plot(WCPstimQuant[HLMWcpKeGGid,"HLMStim4hr_fc"], WCPstimQuant[HLMWcpKeGGid,"intensityweightedmeanratio.PO1063kgg_HLM"],
	 pch=19,
	 cex=0.75,
	 col="gray50",
	 xlab="", ylab="",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
)
par(new=T)
abline(v=0, col="gray70", lwd=2, lty=2)
par(new=T)
abline(h=0, col="gray70", lwd=2, lty=2)
for ( i in 1:length(alltcr) ){
	currgene<-alltcr[i]
	points(WCPstimQuant[which(WCPstimQuant$GeneId==currgene), "HLMStim4hr_fc"],
	       WCPstimQuant[which(WCPstimQuant$GeneId==currgene), "intensityweightedmeanratio.PO1063kgg_HLM"],
		   pch=19,
		   cex=1.25,
		   col="salmon"
		  )
}
for ( i in 1:length(upregtcr) ){
	currgene<-upregtcr[i]
	points(WCPstimQuant[which(WCPstimQuant$GeneId==currgene), "HLMStim4hr_fc"],
	       WCPstimQuant[which(WCPstimQuant$GeneId==currgene), "intensityweightedmeanratio.PO1063kgg_HLM"],
		   pch=19,
		   cex=1.25,
		   col="red"
		  )
}
axis(1,c(seq(xmin,xmax,by=0.25)), labels=F, col="black",cex.axis=0.75,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=1)), labels=T, col="black",cex.axis=1)
mtext("protein expression: 0-4hr TCR stim log2 fold change",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.25)), labels=F, col="black",cex.axis=0.75,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)), labels=T, col="black",cex.axis=1)
mtext("protein Ub: 0-4hr TCR stim log2 fold change",2,line=2.75,cex=1.5)
dev.off()

if ( FALSE )
{
##plot proteomics results vs KeGG quantification results (only for proteins identified with a KeGG)
protrnacorr<-c(ProtUpregRnaUpregUb, ProtDownregRnaDownregUb)
pdf(file.path(".", analysisdir, "ProteomicsKeGGUb_Log2FoldChangeComp_ProtRnaCorrelated.pdf"))
par(bty="l")
xmin <- -3
xmax <- 3
ymin <- -3
ymax <- 3
plot(WCPstimQuant[protrnacorr, "Stim4hr_AvgFc"], KeGGnormLHdata[protrnacorr, "UbFcAvg"],
     pch=19,
	 cex=1,
	 #col="gray50",
	 col=ifelse(WCPstimQuant[protrnacorr,"Stim4hr_AvgFc"]<(PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd), "red", "steelblue1"),
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xlab="", ylab="",
	 xaxt='n', yaxt='n'
)    
par(new=T)
abline(v=0, col="gray70", lwd=2, lty=3)
par(new=T)
abline(h=0, col="gray70", lwd=2, lty=3)
#par(new=T)
#abline(a = 0, b = 1, col="gray70", lwd=2, lty=2)
#par(new=T)
#abline(fit <- lm(KeGGnormLHdata[protrnacorr, "UbFcAvg"] ~ WCPstimQuant[protrnacorr, "Stim4hr_AvgFc"]), col="gray70", lwd=2, lty=5) # regression line (y~x)
#legend("topright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)), text.col="gray")
axis(1,c(seq(xmin,xmax,by=0.5)), labels=F, col="black",cex.axis=0.75)
axis(1,c(seq(xmin,xmax,by=1)), labels=T, col="black",cex.axis=1)
mtext("protein expression: 0-4hr TCR stim log2 fold change",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)), labels=F, col="black",cex.axis=0.75)
axis(2,c(seq(ymin,ymax,by=1)), labels=T, col="black",cex.axis=1)
mtext("protein Ub: 0-4hr TCR stim log2 fold change",2,line=2.75,cex=1.5)
dev.off()

#####PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd

##plot proteomics results vs KeGG quantification results (only for proteins identified with a KeGG)
protrnacorr<-c(ProtUpregRnaUnchUb, ProtUpregRnaDownregUb, ProtUnchRnaUpregUb, ProtUnchRnaDownregUb)
pdf(file.path(".", analysisdir, "ProteomicsKeGGUb_Log2FoldChangeComp_ProtRnaOther.pdf"))
par(bty="l")
xmin <- -6
xmax <- 4
ymin <- -4
ymax <- 4
plot(WCPstimQuant[protrnacorr, "Stim4hr_AvgFc"], KeGGnormLHdata[protrnacorr, "UbFcAvg"],
     pch=19,
	 cex=0.75,
	 col="gray50",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xlab="", ylab="",
	 xaxt='n', yaxt='n'
)    
par(new=T)
abline(v=0, col="gray70", lty=1)
par(new=T)
abline(h=0, col="gray70", lty=1)
par(new=T)
abline(a = 0, b = 1, col="gray70", lty=2)
par(new=T)
points(WCPstimQuant[ProtUpregRnaUnchUb, "Stim4hr_AvgFc"],
       KeGGnormLHdata[ProtUpregRnaUnchUb, "UbFcAvg"],
	pch=19,
	cex=0.75,
	col="purple"
)
#####text(WCPstimQuant[ProtUpregRnaUnchUb, "Stim4hr_AvgFc"], KeGGnormLHdata[ProtUpregRnaUnchUb, "UbFcAvg"], labels=(WCPstimQuant[ProtUpregRnaUnchUb, "GeneId"]), cex=0.25, pos=1, offset=0.2, col="purple", font=2)
points(WCPstimQuant[ProtUpregRnaDownregUb, "Stim4hr_AvgFc"],
       KeGGnormLHdata[ProtUpregRnaDownregUb, "UbFcAvg"],
	pch=19,
	cex=0.75,
	col="cyan"
)
#####text(WCPstimQuant[ProtUpregRnaDownregUb, "Stim4hr_AvgFc"], KeGGnormLHdata[ProtUpregRnaDownregUb, "UbFcAvg"], labels=(WCPstimQuant[ProtUpregRnaDownregUb, "GeneId"]), cex=0.25, pos=1, offset=0.2, col="cyan", font=2)
points(WCPstimQuant[ProtUnchRnaUpregUb, "Stim4hr_AvgFc"],
       KeGGnormLHdata[ProtUnchRnaUpregUb, "UbFcAvg"],
	pch=19,
	cex=0.75,
	col="darkorange"
)
#####for ( i in 1:length(ProtUnchRnaUpregUb) ){
#####   prot<-ProtUnchRnaUpregUb[i]
#####   UbFc<-KeGGnormLHdata[prot, "UbFcAvg"]
#####if ( abs(UbFc)>0.5 ){
#####   text(WCPstimQuant[prot, "Stim4hr_AvgFc"], KeGGnormLHdata[prot, "UbFcAvg"], labels=(WCPstimQuant[prot, "GeneId"]), cex=0.25, pos=1, offset=0.2, col="darkorange", font=2)
#####}
#####}   
points(WCPstimQuant[ProtUnchRnaDownregUb, "Stim4hr_AvgFc"],
       KeGGnormLHdata[ProtUnchRnaDownregUb, "UbFcAvg"],
	pch=19,
	cex=0.75,
	col="gold"
)
axis(1,c(seq(xmin,xmax,by=0.5)), labels=T, col="black",cex.axis=0.75)
#axis(1,c(seq(xmin,xmax,by=1)), labels=T, col="black",cex.axis=1)
mtext("protein expression: 0-4hr TCR stim log2 fold change",1,line=2.75,cex=1.25)
axis(2,c(seq(ymin,ymax,by=0.5)), labels=T, col="black",cex.axis=0.75)
#axis(2,c(seq(ymin,ymax,by=1)), labels=T, col="black",cex.axis=1)
mtext("protein ubiquitylation: 0-4hr TCR stim log2 fold change",2,line=2.75,cex=1.25)
dev.off()

##plot proteomics results vs KeGG quantification results (only for proteins identified with a KeGG)
protrnacorr<-c(ProtDownregRnaUpregUb, ProtDownregRnaUnchUb)
pdf(file.path(".", analysisdir, "ProteomicsKeGGUb_Log2FoldChangeComp_ProtDownRnaUpUnch.pdf"))
par(bty="l")
xmin <- -6
xmax <- 4
ymin <- -4
ymax <- 4
plot(WCPstimQuant[protrnacorr, "Stim4hr_AvgFc"], KeGGnormLHdata[protrnacorr, "UbFcAvg"],
     pch=19,
	 cex=0.75,
	 #col="gray50",
	 col=ifelse(WCPstimQuant[protrnacorr,"log2FoldChangeStimUnstim"]>(RNAseq_FcMean+RNAseq_FcSd), "salmon", "blue"),
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xlab="", ylab="",
	 xaxt='n', yaxt='n'
)    
par(new=T)
abline(v=0, col="gray70", lty=1)
par(new=T)
abline(h=0, col="gray70", lty=1)
par(new=T)
abline(a = 0, b = 1, col="gray70", lty=2)
axis(1,c(seq(xmin,xmax,by=0.5)), labels=T, col="black",cex.axis=0.75)
#axis(1,c(seq(xmin,xmax,by=1)), labels=T, col="black",cex.axis=1)
mtext("protein expression: 0-4hr TCR stim log2 fold change",1,line=2.75,cex=1.25)
axis(2,c(seq(ymin,ymax,by=0.5)), labels=T, col="black",cex.axis=0.75)
#axis(2,c(seq(ymin,ymax,by=1)), labels=T, col="black",cex.axis=1)
mtext("protein ubiquitylation: 0-4hr TCR stim log2 fold change",2,line=2.75,cex=1.25)
dev.off()
}

print("nedd fig")
##plot fold change in KeGG ubiquitination abundance normalized to wcp fold change for + Nedd inhibitor vs - Nedd inhibitor
allcullin<-c("Cul1","Cul2","Cul3","Cul4a","Cul4b","Cul5")
#alltcr<-c("Cd247","Ptprc","Cd3g","Lat","Lck","Itk","Grb2","Sos1","Zap70","Prkcq","Ikbkg","Nfkb1","Cblb","Cd3e","Rhoa","Nras","Kras","Ppp3r1")
#set x and y lims
#xmin <- -2
#xmax <- 3
#ymin <- -4
#ymax <- 3
xmin <- floor(min(WCPstimQuant$PO1063HLMkgg_NormWcpFc, na.rm=TRUE))
xmax <- ceiling(max(WCPstimQuant$PO1063HLMkgg_NormWcpFc, na.rm=TRUE))
ymin <- floor(min(WCPstimQuant$PO1063HLPkgg_NormWcpFc, na.rm=TRUE))
ymax <- ceiling(max(WCPstimQuant$PO1063HLPkgg_NormWcpFc, na.rm=TRUE))
#open file
PlusMinusNeddInhFig<-file.path(".", analysisdir, "ProteinKeGGUb_Cullins_PlusMinusNeddInhibitor_UbFcComp.pdf")
#PlusMinusNeddInhFig<-file.path(".", analysisdir, "ProteinKeGGUb_TCRproteins_PlusMinusNeddInhibitor_UbFcComp.pdf")
#PlusMinusNeddInhFig<-file.path(".", analysisdir, "ProteinKeGGUb_PlusMinusNeddInhibitor_UbFcComp.pdf")
pdf(PlusMinusNeddInhFig)
#make axes
par(bty='l',lwd=2)
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
#make horizontal and vertical lines intersecting origin
par(new=T)
abline(h=0, v=0, col="gray50", lwd=1, lty=1)
#make x=y line
par(new=T)
abline(a = 0, b = 1, col="gray50", lwd=2, lty=2)
#make y=x+1 line
par(new=T)
abline(a = 1, b = 1, col="gray50", lwd=2.25, lty=3)
#make y=x-1 line
par(new=T)
abline(a = -1, b = 1, col="gray50", lwd=2.25, lty=3)
#plot points
par(new=T)
plot(WCPstimQuant$PO1063HLMkgg_NormWcpFc,
     WCPstimQuant$PO1063HLPkgg_NormWcpFc,
     pch=19,
     cex=1,
	 #col="grey70",
	 col=ifelse((WCPstimQuant$PO1063HLPkgg_NormWcpFc<WCPstimQuant$PO1063HLMkgg_NormWcpFc+1 & WCPstimQuant$PO1063HLPkgg_NormWcpFc>WCPstimQuant$PO1063HLMkgg_NormWcpFc-1), "gold", "grey80"),
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     xlab="", ylab="",
     xaxt='n', yaxt='n'
)
#plot diff reg proteins in different points
#par(new=T)
#for ( i in 1:length(PO1063HLMHLPkgg_UbDiffReg) ){
#	currprot<-PO1063HLMHLPkgg_UbDiffReg[i]
#	currgene<-WCPstimQuant[currprot,"GeneId.x"]
#	points(WCPstimQuant[currprot, "PO1063HLMkgg_NormWcpFc"],
#	       WCPstimQuant[currprot, "PO1063HLPkgg_NormWcpFc"],
#	       pch=19,
#		   cex=1.5,
#	       col="grey50"
#	)
#}
#plot cullins
par(new=T)
for ( i in 1:length(allcullin) ){
	currgene<-allcullin[i]
	points(WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "PO1063HLMkgg_NormWcpFc"],
	       WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "PO1063HLPkgg_NormWcpFc"],
		   pch=19,
		   cex=1.5,
		   col="purple"
	)
	#text(WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "PO1063HLMkgg_NormWcpFc"],
	#     WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "PO1063HLPkgg_NormWcpFc"],
	#	 labels=currgene,
	#	 cex=0.4, pos=2, offset=0.2, col="red", font=2
	#)

}
#plot TCR pathway proteins
#par(new=T)
#for ( i in 1:length(alltcr) ){
#	currgene<-alltcr[i]
#	points(WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "PO1063HLMkgg_NormWcpFc"],
#	       WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "PO1063HLPkgg_NormWcpFc"],
#		   pch=19,
#		   cex=2,
#		   col="red"
#	)
#	#text(WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "PO1063HLMkgg_NormWcpFc"],
#	#     WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "PO1063HLPkgg_NormWcpFc"],
#	#	 labels=currgene,
#	#	 cex=1, pos=2, offset=0.2, col="red", font=2
#	#)
#}
#plot proteins with decreased Ub and increased expression in different points
#par(new=T)
#for ( i in 1:length(HLMHLP_NormUbDownReg_WcpUpReg) ){
#	currprot<- HLMHLP_NormUbDownReg_WcpUpReg[i]
#	currgene<-WCPstimQuant[currprot,"GeneId.x"]
#	points(WCPstimQuant[currprot, "PO1063HLMkgg_NormWcpFc"],
#	       WCPstimQuant[currprot, "PO1063HLPkgg_NormWcpFc"],
#		   pch=19,
#		   cex=1.5,
#		   col="dodgerblue",
#		  )
#	text(WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "PO1063HLMkgg_NormWcpFc"],
#	     WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "PO1063HLPkgg_NormWcpFc"],
#		 labels=currgene,
#		 cex=0.4, pos=2, offset=0.2, col="dodgerblue", font=2
#		)
#}
#label axes
axis(1,c(seq(xmin,xmax,by=0.5)), labels=F, col="black",cex.axis=1.25, tck=-0.01)
axis(1,c(seq(xmin,xmax,by=1)), labels=T, col="black",lwd=2,cex.axis=1.5)
mtext("Ub log2 fold change (- Nedd inhibitor)",1,line=2.75,cex=1.75)
axis(2,c(seq(ymin,ymax,by=0.5)), labels=F, col="black",cex.axis=1.25, tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)), labels=T, col="black",lwd=2,cex.axis=1.5)
mtext("Ub log2 fold change (+ Nedd inhibitor)",2,line=2.75,cex=1.75)
dev.off()

stop("CULSTOP")

##plot fold change in wcp protein abundance for + Nedd inhibitor vs - Nedd inhibitor
#set x and y lims
xmin <- -2
xmax <- 3
ymin <- -3
ymax <- 3
#xmin <- floor(min(WCPstimQuant$HLMStim4hr_fc, na.rm=TRUE))
#xmax <- ceiling(max(WCPstimQuant$HLMStim4hr_fc, na.rm=TRUE))
#ymin <- floor(min(WCPstimQuant$HLPStim4hr_fc, na.rm=TRUE))
#ymax <- ceiling(max(WCPstimQuant$HLPStim4hr_fc, na.rm=TRUE))
#WCPNeddFig<-file.path(".", analysisdir, "ProteinWCPabundance_PlusMinusNeddInhibitor_FcComp.pdf")
WCPNeddFig<-file.path(".", analysisdir, "ProteinWCPabundance_WithUbAndExDiff_PlusMinusNeddInhibitor_FcComp.pdf")
pdf(WCPNeddFig)
#make axes
par(bty='l')
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
#make horizontal and vertical lines intersecting origin
par(new=T)
abline(h=0, v=0, col="gray70", lwd=1, lty=3)
#make x=y line
par(new=T)
abline(a = 0, b = 1, col="gray70", lwd=2, lty=2)
#plot points
par(new=T)
plot(WCPstimQuant$HLMStim4hr_fc,
     WCPstimQuant$HLPStim4hr_fc,
	 #pch=21,
	 pch=19,
	 cex=1.0,
	 #bg="grey70",
	 col="grey70",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
     xlab="", ylab="",
     xaxt='n', yaxt='n'
)
#plot diff reg proteins in different points
par(new=T)
for ( i in 1:length(PO1063HLPPO1063HLM_Upreg_AvgFcMeanSd) ){
	currprot<-PO1063HLPPO1063HLM_Upreg_AvgFcMeanSd[i]
	currgene<-WCPstimQuant[currprot,"GeneId.x"]
	points(WCPstimQuant[currprot,"HLMStim4hr_fc"],
	       WCPstimQuant[currprot,"HLPStim4hr_fc"],
		   #pch=21,
		   pch=19,
		   cex=1.5,
		   #bg="grey50",
		   col="grey50"
	)
}
#plot cullins in different points
par(new=T)
for ( i in 1:length(cullins) ){
	currgene<-cullins[i]
	points(WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "HLMStim4hr_fc"],
	       WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "HLPStim4hr_fc"],
		   #pch=21,
		   pch=19,
		   cex=1.5,
		   #bg="red",
		   col="red"
		  )
	text(WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "HLMStim4hr_fc"],
	     WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "HLPStim4hr_fc"],
		 labels=currgene,
		 cex=0.4, pos=2, offset=0.2, col="red", font=2
		)
}
#plot proteins with decreased Ub and increased expression in different points
par(new=T)
for ( i in 1:length(HLMHLP_NormUbDownReg_WcpUpReg) ){
	currprot<- HLMHLP_NormUbDownReg_WcpUpReg[i]
	currgene<-WCPstimQuant[currprot,"GeneId.x"]
	points(WCPstimQuant[currprot, "HLMStim4hr_fc"],
	       WCPstimQuant[currprot, "HLPStim4hr_fc"],
		   #pch=1,
		   pch=19,
		   cex=1.5,
		   col="dodgerblue",
		  )
	text(WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "HLMStim4hr_fc"],
	     WCPstimQuant[which(WCPstimQuant$GeneId.x==currgene), "HLPStim4hr_fc"],
		 labels=currgene,
		 cex=0.4, pos=2, offset=0.2, col="dodgerblue", font=2
		)
}
#label axes
axis(1,c(seq(xmin,xmax,by=0.5)), labels=F, col="black",cex.axis=1.0, tck=-0.01)
axis(1,c(seq(xmin,xmax,by=1)), labels=T, col="black",cex.axis=1.0)
mtext("wcp log2 fold change (- Nedd inhibitor)",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)), labels=F, col="black",cex.axis=1.0, tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)), labels=T, col="black",cex.axis=1.0)
mtext("wcp log2 fold change (+ Nedd inhibitor)",2,line=2.75,cex=1.5)
dev.off()

print("HDAC analysis...")

##analysis of HDAC proteins
#list of HDAC proteins
HDACgenes <- c("Hdac1", "Hdac2", "Hdac3", "Hdac8", "Hdac4", "Hdac5", "Hdac7", "Hdac9", "Hdac6", "Hdac10", "Sirt1", "Sirt2", "Sirt3", "Sirt4", "Sirt5", "Sirt6", "Sirt7", "Hdac11")
#protein abundance barplot
pdf(file.path(".", analysisdir, "HDACproteins_AbundanceZscores.pdf"))
quants <- vector()
genelabels <- vector()
for ( i in 1:length(HDACgenes) ){
	currgene <- HDACgenes[i]
	currzscore <- WCPstimQuant[which(WCPstimQuant$GeneId.x == currgene),"AbundanceInt_zscore0hr"]
	if ( identical(currzscore, numeric(0)) ){
		currzscore<-NA
	}
	quants <- c(quants, currzscore)
}
print(quants)
barplot(quants, beside=TRUE, names.arg=HDACgenes, xlab="", ylab="", ylim=c(-1.5,1.5), space=0.5, cex.axis=1, cex.names=1, las=2, col=c("darkblue", "darkblue", "darkblue", "darkblue", "darkcyan", "darkcyan", "darkcyan", "darkcyan", "darkcyan", "darkcyan", "dodgerblue", "dodgerblue", "dodgerblue", "dodgerblue", "dodgerblue", "dodgerblue", "dodgerblue", "purple"))
#width=c(2)
#space=c(1,10)
mtext("abundance Z-score", side=2, line=2.75, cex=1)
dev.off()
#protein and RNAseq fold change scatterplot 
SubData<-WCPstimQuant[WCPstimQuant$GeneId.x%in%HDACgenes,]
pdf(file.path(".", analysisdir, "ProteomicsRNAseq_Log2FoldChangeComp_HDACprots.pdf"))
par(bty="l")
xmin <- floor(min(SubData$Stim4hr_AvgFc, na.rm=TRUE))
xmax <- ceiling(max(SubData$Stim4hr_AvgFc, na.rm=TRUE))
ymin <- floor(min(SubData$log2FoldChangeStimUnstim, na.rm=TRUE))
ymax <- ceiling(max(SubData$log2FoldChangeStimUnstim, na.rm=TRUE))
#xmin <- -2
#xmax <- 1.5
#ymin <- -4
#ymax <- 4
#make axes
plot(NULL, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab="" ,ylab="", xaxt='n', yaxt='n')
#add horizontal and vertical lines at 0
par(new=T)
abline(v=0, col="gray70", lwd=2, lty=1)
par(new=T)
abline(h=0, col="gray70", lwd=2, lty=1)
#add vertical lines at x-axis proteomics mean +/- standard deviation
par(new=T)
abline(v=PO972PO1063PO1075_AvgFcMean+PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=1, lty=2)
par(new=T)
abline(v=PO972PO1063PO1075_AvgFcMean-PO972PO1063PO1075_AvgFcSd, col="gray70", lwd=1, lty=2)
#add horizontal lines at y-axis RNAseq mean +/- standard deviation
par(new=T)
abline(h=RNAseq_FcMean+RNAseq_FcSd, col="gray70", lwd=1, lty=2)
par(new=T)
abline(h=RNAseq_FcMean-RNAseq_FcSd, col="gray70", lwd=1, lty=2)
#plot points
par(new=T)
plot(SubData$Stim4hr_AvgFc,
     SubData$log2FoldChangeStimUnstim,
     pch=19,
	 cex=2,
	 #cex=4*(abs(SubData$UbFcAvg_4hr_NormWcpFc)),
	 col="grey50",
	 xlab="",ylab="",
	 #col=ifelse(is.na(SubData$UbFcAvg_4hr_NormWcpFc), "grey90", ifelse(SubData$UbFcAvg_4hr_NormWcpFc>0, "salmon", ifelse(SubData$UbFcAvg_4hr_NormWcpFc<0, "blue", "grey50"))), 
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xaxt='n', yaxt='n'
)
par(new=T)
for ( i in 1:length(HDACgenes) )
{
	currgene<-HDACgenes[i]
	print(currgene)
	currzscore <- WCPstimQuant[which(WCPstimQuant$GeneId.x == currgene),"AbundanceInt_zscore"]
	if ( identical(SubData[which(SubData$GeneId.x==currgene), "Stim4hr_AvgFc"], numeric(0)) | identical(SubData[which(SubData$GeneId.x==currgene), "log2FoldChangeStimUnstim"], numeric(0)) ){
		next
	}
	text(SubData[which(SubData$GeneId.x==currgene), "Stim4hr_AvgFc"],
	     SubData[which(SubData$GeneId.x==currgene), "log2FoldChangeStimUnstim"],
		 labels=currgene,
		 cex=0.25, pos=2, offset=0.2, col="blue", font=2
 	    )
}
axis(1,c(seq(xmin,xmax,by=0.25)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(1,c(seq(xmin,xmax,by=0.5)),col="black",cex.axis=0.80)
mtext("wcp log2 fold change",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)),label=F,col="black",cex.axis=1.0,tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)),col="black",cex.axis=0.80)
mtext("RNAseq log2 fold change",2,line=2.75,cex=1.5)
dev.off()
##end proteomics vs RNAseq results plot
stop("EndAll")


par(new=T)
points(WCPstimQuant[ProtDownregRnaUpregUb, "Stim4hr_AvgFc"],
       KeGGnormLHdata[ProtDownregRnaUpregUb, "UbFcAvg"],
	   pch=19,
	   cex=1,
	   col="salmon"
)
#####text(WCPstimQuant[ProtDownregRnaUpregUb, "Stim4hr_AvgFc"], KeGGnormLHdata[ProtDownregRnaUpregUb, "UbFcAvg"], labels=(WCPstimQuant[ProtDownregRnaUpregUb, "GeneId"]), cex=0.25, pos=2, offset=0.2, col="salmon", font=2)
points(WCPstimQuant[ProtDownregRnaUnchUb, "Stim4hr_AvgFc"],
       KeGGnormLHdata[ProtDownregRnaUnchUb, "UbFcAvg"],
	   pch=19,
	   cex=1,
	   col="blue"
)
#####text(WCPstimQuant[ProtDownregRnaUnchUb, "Stim4hr_AvgFc"], KeGGnormLHdata[ProtDownregRnaUnchUb, "UbFcAvg"], labels=(WCPstimQuant[ProtDownregRnaUnchUb, "GeneId"]), cex=0.25, pos=2, offset=0.2, col="blue", font=2)
points(WCPstimQuant[ProtDownregRnaDownregUb, "Stim4hr_AvgFc"],
       KeGGnormLHdata[ProtDownregRnaDownregUb, "UbFcAvg"],
	   pch=19,
	   cex=1,
	   col="red"
)
#####text(WCPstimQuant[ProtDownregRnaDownregUb, "Stim4hr_AvgFc"], KeGGnormLHdata[ProtDownregRnaDownregUb, "UbFcAvg"], labels=(WCPstimQuant[ProtDownregRnaDownregUb, "GeneId"]), cex=0.25, pos=2, offset=0.2, col="red", font=2)
points(WCPstimQuant[ProtUpregRnaUpregUb, "Stim4hr_AvgFc"],
       KeGGnormLHdata[ProtUpregRnaUpregUb, "UbFcAvg"],
	   pch=19,
	   cex=1,
	   col="steelblue1"
)
#####text(WCPstimQuant[ProtUpregRnaUpregUb, "Stim4hr_AvgFc"], KeGGnormLHdata[ProtUpregRnaUpregUb, "UbFcAvg"], labels=(WCPstimQuant[ProtUpregRnaUpregUb, "GeneId"]), cex=0.25, pos=1, offset=0.2, col="steelblue1", font=2)
points(WCPstimQuant[ProtUpregRnaUnchUb, "Stim4hr_AvgFc"],
       KeGGnormLHdata[ProtUpregRnaUnchUb, "UbFcAvg"],
	   pch=19,
	   cex=1,
	   col="purple"
)
#####text(WCPstimQuant[ProtUpregRnaUnchUb, "Stim4hr_AvgFc"], KeGGnormLHdata[ProtUpregRnaUnchUb, "UbFcAvg"], labels=(WCPstimQuant[ProtUpregRnaUnchUb, "GeneId"]), cex=0.25, pos=1, offset=0.2, col="purple", font=2)
points(WCPstimQuant[ProtUpregRnaDownregUb, "Stim4hr_AvgFc"],
       KeGGnormLHdata[ProtUpregRnaDownregUb, "UbFcAvg"],
	   pch=19,
	   cex=1,
	   col="cyan"
)
#####text(WCPstimQuant[ProtUpregRnaDownregUb, "Stim4hr_AvgFc"], KeGGnormLHdata[ProtUpregRnaDownregUb, "UbFcAvg"], labels=(WCPstimQuant[ProtUpregRnaDownregUb, "GeneId"]), cex=0.25, pos=1, offset=0.2, col="cyan", font=2)
points(WCPstimQuant[ProtUnchRnaUpregUb, "Stim4hr_AvgFc"],
       KeGGnormLHdata[ProtUnchRnaUpregUb, "UbFcAvg"],
	   pch=19,
	   cex=1,
	   col="darkorange"
)
#####for ( i in 1:length(ProtUnchRnaUpregUb) ){
#####	prot<-ProtUnchRnaUpregUb[i]
#####	UbFc<-KeGGnormLHdata[prot, "UbFcAvg"]
	#####if ( abs(UbFc)>0.5 ){
	#####	text(WCPstimQuant[prot, "Stim4hr_AvgFc"], KeGGnormLHdata[prot, "UbFcAvg"], labels=(WCPstimQuant[prot, "GeneId"]), cex=0.25, pos=1, offset=0.2, col="darkorange", font=2)
	#####}
#####}
points(WCPstimQuant[ProtUnchRnaDownregUb, "Stim4hr_AvgFc"],
       KeGGnormLHdata[ProtUnchRnaDownregUb, "UbFcAvg"],
	   pch=19,
	   cex=1,
	   col="gold"
)
#####for ( i in 1:length(ProtUnchRnaDownregUb) ){
	#####prot<-ProtUnchRnaDownregUb[i]
	#####UbFc<-KeGGnormLHdata[prot, "UbFcAvg"]
	#####if ( abs(UbFc)>0.5 ){
	#####	text(WCPstimQuant[prot, "Stim4hr_AvgFc"], KeGGnormLHdata[prot, "UbFcAvg"], labels=(WCPstimQuant[prot, "GeneId"]), cex=0.25, pos=1, offset=0.2, col="gold", font=2)
	#####}
#####}
axis(1,c(seq(xmin,xmax,by=0.5)), labels=T, col="black",cex.axis=0.75)
#axis(1,c(seq(xmin,xmax,by=1)), labels=T, col="black",cex.axis=1)
mtext("protein expression: 0-4hr TCR stim log2 fold change",1,line=2.75,cex=1.25)
axis(2,c(seq(ymin,ymax,by=0.5)), labels=T, col="black",cex.axis=0.75)
#axis(2,c(seq(ymin,ymax,by=1)), labels=T, col="black",cex.axis=1)
mtext("protein ubiquitylation: 0-4hr TCR stim log2 fold change",2,line=2.75,cex=1.25)
dev.off()


stop("Stopped")


if ( FALSE )
{
#KeGG unique up and down regulated
KeGG_Downreg_Unique<-setdiff(PO972_Downreg_MeanSd, Reduce(union, list(PO1075PO972_MeanSd_Downreg,PO1063PO972_MeanSd_Downreg,PO1075PO1063PO972_MeanSd_Downreg)))
KeGG_Downreg_Unique_Num<-length(KeGG_Downreg_Unique)
print(KeGG_Downreg_Unique_Num)
KeGG_Upreg_Unique<-setdiff(PO972_Upreg_MeanSd, Reduce(union, list(PO1075PO972_MeanSd_Upreg,PO1063PO972_MeanSd_Upreg,PO1075PO1063PO972_MeanSd_Upreg)))
KeGG_Upreg_Unique_Num<-length(KeGG_Upreg_Unique)
print(KeGG_Upreg_Unique_Num)
KeGG_DiffReg_Unique<-c(KeGG_Downreg_Unique,KeGG_Upreg_Unique)
KeGG_DiffReg_Unique_Num<-length(KeGG_DiffReg_Unique)
print(KeGG_DiffReg_Unique_Num)
#count dist of unique up and down regulated genes
KeGG_DiffReg_Unique_KeGGCounts <- WCPstimQuant[KeGG_DiffReg_Unique,"IGD0hr4hr_Maxzscore"]
KeGG_DiffReg_Unique_KeGGCounts_Num <- length(KeGG_DiffReg_Unique_KeGGCounts)
print(KeGG_DiffReg_Unique_KeGGCounts_Num)
plotnormdistribution(KeGG_DiffReg_Unique_KeGGCounts, "NA", "NA", 0.1, "NA", "NA", "NA", "NA", "KeGG 0-4hr max intensity zscore", "norm freq", file.path(".", analysisdir, "KeGG_DiffReg_Unique_IntensityDistribution.pdf"))
#KeGG common up and donw regulated
KeGG_Downreg_Common<-Reduce(union, list(PO1075PO972_MeanSd_Downreg,PO1063PO972_MeanSd_Downreg,PO1075PO1063PO972_MeanSd_Downreg))
KeGG_Downreg_Common_Num<-length(KeGG_Downreg_Common)
print(KeGG_Downreg_Common_Num)
KeGG_Upreg_Common<-Reduce(union, list(PO1075PO972_MeanSd_Upreg,PO1063PO972_MeanSd_Upreg,PO1075PO1063PO972_MeanSd_Upreg))
KeGG_Upreg_Common_Num<-length(KeGG_Upreg_Common)
print(KeGG_Upreg_Common_Num)
KeGG_DiffReg_Common<-c(KeGG_Downreg_Common,KeGG_Upreg_Common)
KeGG_DiffReg_Common_Num<-length(KeGG_DiffReg_Common)
print(KeGG_DiffReg_Common_Num)
#count dist of common up and down regulated genes
KeGG_DiffReg_Common_KeGGCounts <- WCPstimQuant[KeGG_DiffReg_Common,"IGD0hr4hr_Maxzscore"]
KeGG_DiffReg_Common_KeGGCounts_Num <- length(KeGG_DiffReg_Common_KeGGCounts)
print(KeGG_DiffReg_Common_KeGGCounts_Num)
plotnormdistribution(KeGG_DiffReg_Common_KeGGCounts, "NA", "NA", 0.1, "NA", "NA", "NA", "NA", "KeGG 0-4hr max intensity zscore", "norm freq", file.path(".", analysisdir, "KeGG_DiffReg_Common_IntensityDistribution.pdf"))
#sig test
t.test(KeGG_DiffReg_Unique_KeGGCounts,KeGG_DiffReg_Common_KeGGCounts,paired=FALSE,var.equal=FALSE,na.action="na.omit")
#print(xyTest)
distlist<-list(KeGG_DiffReg_Unique_KeGGCounts,KeGG_DiffReg_Common_KeGGCounts)
plotnormdistributionS(distlist,"NA", "NA", 0.1, "NA", "NA", "NA", "NA", "KeGG 0-4hr max intensity zscore", "norm freq", file.path(".", analysisdir, "KeGG_DiffReg_UniqueAndCommon_IntensityDistribution.pdf"))


#TUBE unique up and down regulated
Tube_Downreg_Unique<-setdiff(PO1075_Downreg_MeanSd, Reduce(union, list(PO1075PO972_MeanSd_Downreg,PO1075PO1063_MeanSd_Downreg,PO1075PO1063PO972_MeanSd_Downreg)))
Tube_Downreg_Unique_Num<-length(Tube_Downreg_Unique)
print(Tube_Downreg_Unique_Num)
Tube_Upreg_Unique<-setdiff(PO1075_Upreg_MeanSd, Reduce(union, list(PO1075PO972_MeanSd_Upreg,PO1075PO1063_MeanSd_Upreg,PO1075PO1063PO972_MeanSd_Upreg)))
Tube_Upreg_Unique_Num<-length(Tube_Upreg_Unique)
print(Tube_Upreg_Unique_Num)
#count dist of unique up and down regulated genes
Tube_DiffReg_Unique<-c(Tube_Downreg_Unique,Tube_Upreg_Unique)
Tube_DiffReg_Unique_Num<-length(Tube_DiffReg_Unique)
print(Tube_DiffReg_Unique_Num)
Tube_DiffReg_Unique_TubeCounts <- WCPstimQuant[Tube_DiffReg_Unique,"WCPMint_zscore"]
Tube_DiffReg_Unique_TubeCounts_Num <- length(Tube_DiffReg_Unique_TubeCounts)
print(Tube_DiffReg_Unique_TubeCounts_Num)
plotnormdistribution(Tube_DiffReg_Unique_TubeCounts, "NA", "NA", 1, "NA", "NA", "NA", "NA", "Tube LHratio Intensity Zscore", "norm freq", file.path(".", analysisdir, "Tube_DiffReg_Unique_InentsityDistribution.pdf"))
#TUBE common up and donw regulated
Tube_Downreg_Common<-Reduce(union, list(PO1075PO972_MeanSd_Downreg,PO1075PO1063_MeanSd_Downreg,PO1075PO1063PO972_MeanSd_Downreg))
Tube_Downreg_Common_Num<-length(Tube_Downreg_Common)
print(Tube_Downreg_Common_Num)
Tube_Upreg_Common<-Reduce(union, list(PO1075PO972_MeanSd_Upreg,PO1075PO1063_MeanSd_Upreg,PO1075PO1063PO972_MeanSd_Upreg))
Tube_Upreg_Common_Num<-length(Tube_Upreg_Common)
print(Tube_Upreg_Common_Num)
Tube_DiffReg_Common<-c(Tube_Downreg_Common,Tube_Upreg_Common)
Tube_DiffReg_Common_Num<-length(Tube_DiffReg_Common)
print(Tube_DiffReg_Common_Num)
#count dist of common up and down regulated genes
Tube_DiffReg_Common_TubeCounts <- WCPstimQuant[Tube_DiffReg_Common,"WCPMint_zscore"]
Tube_DiffReg_Common_TubeCounts_Num <- length(Tube_DiffReg_Common_TubeCounts)
print(Tube_DiffReg_Common_TubeCounts_Num)
plotnormdistribution(Tube_DiffReg_Common_TubeCounts, "NA", "NA", 0.1, "NA", "NA", "NA", "NA", "Tube LHratio Intensity Zscore", "norm freq", file.path(".", analysisdir, "Tube_DiffReg_Common_IntensityDistribution.pdf"))
#sig test
t.test(Tube_DiffReg_Unique_TubeCounts,Tube_DiffReg_Common_TubeCounts,paired=FALSE,var.equal=FALSE,na.action="na.omit")
#print(xyTest)
distlist<-list(Tube_DiffReg_Unique_TubeCounts,Tube_DiffReg_Common_TubeCounts)
plotnormdistributionS(distlist,"NA", "NA", 0.1, "NA", "NA", "NA", "NA", "Tube LHratio Intensity Zscore", "norm freq", file.path(".", analysisdir, "Tube_DiffReg_UniqueAndCommon_IntensityDistribution.pdf"))
					

#Nedd unique up and down regulated
Nedd_Downreg_Unique<-setdiff(PO1063_Downreg_MeanSd, Reduce(union, list(PO1063PO972_MeanSd_Downreg,PO1075PO1063_MeanSd_Downreg,PO1075PO1063PO972_MeanSd_Downreg)))
Nedd_Downreg_Unique_Num<-length(Nedd_Downreg_Unique)
print(Nedd_Downreg_Unique_Num)
Nedd_Upreg_Unique<-setdiff(PO1063_Upreg_MeanSd, Reduce(union, list(PO1063PO972_MeanSd_Upreg,PO1075PO1063_MeanSd_Upreg,PO1075PO1063PO972_MeanSd_Upreg)))
Nedd_Upreg_Unique_Num<-length(Nedd_Upreg_Unique)
print(Nedd_Upreg_Unique_Num)
#count dist of unique up and down regulated genes
Nedd_DiffReg_Unique<-c(Nedd_Downreg_Unique,Nedd_Upreg_Unique)
Nedd_DiffReg_Unique_Num<-length(Nedd_DiffReg_Unique)
print(Nedd_DiffReg_Unique_Num)
Nedd_DiffReg_Unique_NeddCounts <- WCPstimQuant[Nedd_DiffReg_Unique,"HLMint_zscore"]
Nedd_DiffReg_Unique_NeddCounts_Num <- length(Nedd_DiffReg_Unique_NeddCounts)
print(Nedd_DiffReg_Unique_NeddCounts_Num)
plotnormdistribution(Nedd_DiffReg_Unique_NeddCounts, "NA", "NA", 1, "NA", "NA", "NA", "NA", "Nedd LHratio Intensity Zscore", "norm freq", file.path(".", analysisdir, "Nedd_DiffReg_Unique_IntensityDistribution.pdf"))
#Nedd common up and donw regulated
Nedd_Downreg_Common<-Reduce(union, list(PO1063PO972_MeanSd_Downreg,PO1075PO1063_MeanSd_Downreg,PO1075PO1063PO972_MeanSd_Downreg))
Nedd_Downreg_Common_Num<-length(Nedd_Downreg_Common)
print(Tube_Downreg_Common_Num)
Nedd_Upreg_Common<-Reduce(union, list(PO1063PO972_MeanSd_Upreg,PO1075PO1063_MeanSd_Upreg,PO1075PO1063PO972_MeanSd_Upreg))
Nedd_Upreg_Common_Num<-length(Nedd_Upreg_Common)
print(Nedd_Upreg_Common_Num)
Nedd_DiffReg_Common<-c(Nedd_Downreg_Common,Nedd_Upreg_Common)
Nedd_DiffReg_Common_Num<-length(Nedd_DiffReg_Common)
print(Nedd_DiffReg_Common_Num)
#count dist of common up and down regulated genes
Nedd_DiffReg_Common_NeddCounts <- WCPstimQuant[Nedd_DiffReg_Common,"HLMint_zscore"]
Nedd_DiffReg_Common_NeddCounts_Num <- length(Nedd_DiffReg_Common_NeddCounts)
print(Nedd_DiffReg_Common_NeddCounts_Num)
plotnormdistribution(Nedd_DiffReg_Common_NeddCounts, "NA", "NA", 0.1, "NA", "NA", "NA", "NA", "Nedd LHratio Intensity Zscore", "norm freq", file.path(".", analysisdir, "Nedd_DiffReg_Common_IntensityDistribution.pdf"))
#sig test
t.test(Nedd_DiffReg_Unique_NeddCounts,Nedd_DiffReg_Common_NeddCounts,paired=FALSE,var.equal=FALSE,na.action="na.omit")
#print(xyTest)
distlist<-list(Nedd_DiffReg_Unique_NeddCounts,Nedd_DiffReg_Common_NeddCounts)
plotnormdistributionS(distlist,"NA", "NA", 0.1, "NA", "NA", "NA", "NA", "Nedd LHratio Intensity Zscore", "norm freq", file.path(".", analysisdir, "Nedd_DiffReg_UniqueAndCommon_IntensityDistribution.pdf"))
}


stop("stopped")




########## old stuff from a different project ##########

plotnormdistributionS <- function(distlist, breaksmin, breaksmax, binsize, filename, xlabel, ylabel) {

	####mean and sd
	###dist_mean<-mean(dist, na.rm=TRUE)
	###dist_sd<-sd(dist, na.rm=TRUE)

	#break points
	if ( breaksmin=="NA"| breaksmax=="NA" ){
		breaksmin<-NA
		breaksmax<-NA
		for ( distnum in 1:length(distlist)){
			currentdist <- distlist[[distnum]]
			currentmin<-floor(min(currentdist,na.rm=TRUE))
			currentmax<-ceiling(max(currentdist,na.rm=TRUE))
			if ( is.na(breaksmin) ){
				breaksmin<-currentmin
			}else if ( currentmin < breaksmin ){
				breaksmin<-currentmin
			}
			if ( is.na(breaksmax) ){
				breaksmax<-currentmax
			}else if ( currentmax > breaksmax ){
				breaksmax<-currentmax
			}
		}
	}

	#bin size
	if ( binsize=="NA" ){
		binsize<-1
	}

	#print(breaksmin)
	#print(breaksmax)
	#print(binsize)

	#need $mids, $counts

	pdf(filename)
	par(bty="l")
	color<-c("red", "blue")
	xmin<-breaksmin
	xmax<-breaksmax
	for ( distnum in 1:length(distlist)){
	currentcol<-color[distnum]
	dist <- distlist[[distnum]]
	#print(dist)
	#histogram and norm frequency
	hist_calc<-hist(dist, breaks=seq(breaksmin,breaksmax,by=binsize), plot=FALSE, right=FALSE)
	hist_sum<-sum(hist_calc$counts)
	hist_normfreq<-hist_calc$counts/hist_sum
	hist_normfreq<-c(hist_normfreq, 0)

	ymin<-0
	hist_normfreq_max<-max(hist_normfreq)
	hist_normfreq_max_round<-round(hist_normfreq_max,1)
	if ( hist_normfreq_max_round<hist_normfreq_max ){
		ymax<-hist_normfreq_max_round+0.05
	}else{
		ymax<-signif(hist_normfreq_max_round,2)
	}

	#intensityHeavy fold change distribution
	###pdf(filename)
	###par(bty="l")
	plot(hist_calc$breaks, hist_normfreq,     
		 #type="l", col="gray50", lwd=3,
		 type="l", col=currentcol, lwd=3,
		 xlab="", ylab="",
		 #xlim=c(breaksmin,breaksmax),
		 #ylim=c(ymin,ymax),
		 xlim=c(xmin,xmax),
		 xaxt='n', yaxt='n'
		)
	par(new=T)
	}
	###par(new=T)
	###abline(v=dist_mean, col="gray70")
	###par(new=T)
	###abline(v=dist_mean+dist_sd, col="gray70", lty=2)
	###par(new=T)
	###abline(v=dist_mean-dist_sd, col="gray70", lty=2)
	###par(new=T)
	###abline(v=dist_mean+(2*dist_sd), col="gray70", lty=2)
	###par(new=T)
	###abline(v=dist_mean-(2*dist_sd), col="gray70", lty=2)
	axis(1,at=c(seq(xmin,xmax,by=0.50)),labels=F,col="black",cex.axis=1)
	axis(1,at=c(seq(xmin,xmax,by=1)),col="black",cex.axis=1.1)
	mtext(xlabel,1,line=2.75,cex=1.5)
	axis(2,at=c(seq(ymin,ymax,by=0.01)),labels=F,col="black",cex.axis=1)
	axis(2,at=c(seq(ymin,ymax,,by=0.05)),col="black",cex.axis=1.1)
	mtext(ylabel,2,line=2.75,cex=1.5)
	dev.off()
	
	return

}

pdf("Nedd-TUBE_LHratiocomp.pdf")
par(bty='l')
plot(WCPstimQuant[HLMWCPMLHratioId,"normLHratio.PO1063_HLM"], WCPstimQuant[HLMWCPMLHratioId,"normLHratio.PO1075_WCPM"],
     pch=19, cex=1.25, col="darkgray",
	 #xlim=c(-3,3), ylim=c(-3,3),
	 #xaxt='n', yaxt='n',
	 xlab="", ylab="",
	)
par(new=T)
abline(a = 0, b = 1, col="gray70", lty=2)
par(new=T)
abline(h=0, col="gray70", lty=2)
par(new=T)
abline(v=0, col="gray70", lty=2)
#par(new=T)
#abline(lm(filtereddata$intensityweightedmean.PO972kgg_4diff0 ~ filtereddata$intensityweightedmean.PO1063kgg_HLM), col="gray70", lwd=2) # regression line (y~x)
#axis(1,at=c(seq(-3,3,by=1)),col="black",cex.axis=1.0)
mtext("Nedd 0-4 hr stim fold change",1,line=2.75,cex=1.25)
#axis(2,at=c(seq(-3,3,by=1)),col="black",cex.axis=1.0)
mtext("TUBE 0-4 hr stim fold change",2,line=2.75,cex=1.25)
#add points
points(WCPstimQuant["Q2TB02","normLHratio.PO1063_HLM"],WCPstimQuant["Q2TB02","normLHratio.PO1075_WCPM"], pch=19, col="red")
text(WCPstimQuant["Q2TB02","normLHratio.PO1063_HLM"],WCPstimQuant["Q2TB02","normLHratio.PO1075_WCPM"], labels="Nfkbid", cex=0.75, pos=2, col="red", font=2)
points(WCPstimQuant["Q3UV17","normLHratio.PO1063_HLM"],WCPstimQuant["Q3UV17","normLHratio.PO1075_WCPM"], pch=19, col="red")
text(WCPstimQuant["Q3UV17","normLHratio.PO1063_HLM"],WCPstimQuant["Q3UV17","normLHratio.PO1075_WCPM"], labels="Krt76", cex=0.75, pos=2, col="red", font=2)
points(WCPstimQuant["Q544J7","normLHratio.PO1063_HLM"],WCPstimQuant["Q544J7","normLHratio.PO1075_WCPM"], pch=19, col="red")
text(WCPstimQuant["Q544J7","normLHratio.PO1063_HLM"],WCPstimQuant["Q544J7","normLHratio.PO1075_WCPM"], labels="Irf8", cex=0.75, pos=2, col="red", font=2)
points(WCPstimQuant["Q569U6","normLHratio.PO1063_HLM"],WCPstimQuant["Q569U6","normLHratio.PO1075_WCPM"], pch=19, col="red")
text(WCPstimQuant["Q569U6","normLHratio.PO1063_HLM"],WCPstimQuant["Q569U6","normLHratio.PO1075_WCPM"], labels="Junb", cex=0.75, pos=2, col="red", font=2)
points(WCPstimQuant["Q64287.2","normLHratio.PO1063_HLM"],WCPstimQuant["Q64287.2","normLHratio.PO1075_WCPM"], pch=19, col="red")
text(WCPstimQuant["Q64287.2","normLHratio.PO1063_HLM"],WCPstimQuant["Q64287.2","normLHratio.PO1075_WCPM"], labels="Irf4", cex=0.75, pos=2, col="red", font=2)
points(WCPstimQuant["Q9JL25.2","normLHratio.PO1063_HLM"],WCPstimQuant["Q9JL25.2","normLHratio.PO1075_WCPM"], pch=19, col="red")
text(WCPstimQuant["Q9JL25.2","normLHratio.PO1063_HLM"],WCPstimQuant["Q9JL25.2","normLHratio.PO1075_WCPM"], labels="Rgs1", cex=0.75, pos=2, col="red", font=2)
points(WCPstimQuant["Q5BL07.2","normLHratio.PO1063_HLM"],WCPstimQuant["Q5BL07.2","normLHratio.PO1075_WCPM"], pch=19, col="red")
text(WCPstimQuant["Q5BL07.2","normLHratio.PO1063_HLM"],WCPstimQuant["Q5BL07.2","normLHratio.PO1075_WCPM"], labels="Pex1", cex=0.75, pos=2, col="red", font=2)
points(WCPstimQuant["P07750","normLHratio.PO1063_HLM"],WCPstimQuant["P07750","normLHratio.PO1075_WCPM"], pch=19, col="red")
text(WCPstimQuant["P07750","normLHratio.PO1063_HLM"],WCPstimQuant["P07750","normLHratio.PO1075_WCPM"], labels="Il4", cex=0.75, pos=2, col="red", font=2)
dev.off()


##########

##subset of proteins with iBAQ quantification in 0-4 hours whole cell proteome
UnstimStim4hrQuantifiedWCP<-row.names(subset(WCPiBAQdata, !is.na(WCPiBAQdata$Unstim_Stim4hr_fc)))
UnstimStim4hrQuantifiedWCPNum<-length(UnstimStim4hrQuantifiedWCP)
#head(UnstimStim4hrQuantifiedWCP)
print(UnstimStim4hrQuantifiedWCPNum)

##subset of proteins with KeGG-based Ub quantification in 0-4 hours KeGG enrichment
UnstimStim4hrQuantifiedUb<-row.names(subset(KeGGUbnormLHdata, !is.na(intensityweightedmean.PO972kgg_4diff0)))
UnstimStim4hrQuantifiedUbNum<-length(UnstimStim4hrQuantifiedUb)
#head(UnstimStim4hrQuantifiedUb)
print(UnstimStim4hrQuantifiedUbNum)

#intersection of WCP and Ub quantified proteins for 0-4 hours
UnstimStim4hrQuantifiedWCPUb<-intersect(UnstimStim4hrQuantifiedWCP,UnstimStim4hrQuantifiedUb)
UnstimStim4hrQuantifiedWCPUbNum<-length(UnstimStim4hrQuantifiedWCPUb)
#head(UnstimStim4hrQuantifiedWCPUbNum)
print(UnstimStim4hrQuantifiedWCPUbNum)

##subset of proteins with iBAQ quantification in 0-1 hours whole cell proteome
UnstimStim1hrQuantifiedWCP<-row.names(subset(WCPiBAQdata, !is.na(WCPiBAQdata$Unstim_Stim1hr_fc)))
UnstimStim1hrQuantifiedWCPNum<-length(UnstimStim1hrQuantifiedWCP)
#head(UnstimStim1hrQuantifiedWCPNum)
print(UnstimStim1hrQuantifiedWCPNum)

##subset of proteins with KeGG-based Ub quantification in 0-1 hours KeGG enrichment
UnstimStim1hrQuantifiedUb<-row.names(subset(KeGGUbnormLHdata, !is.na(intensityweightedmean.PO972kgg_1diff0)))
UnstimStim1hrQuantifiedUbNum<-length(UnstimStim1hrQuantifiedUb)
#head(UnstimStim1hrQuantifiedUb)
print(UnstimStim1hrQuantifiedUbNum)

#intersection of WCP and Ub quantified proteins for 0-1 hours
UnstimStim1hrQuantifiedWCPUb<-intersect(UnstimStim1hrQuantifiedWCP,UnstimStim1hrQuantifiedUb)
UnstimStim1hrQuantifiedWCPUbNum<-length(UnstimStim1hrQuantifiedWCPUb)
#head(UnstimStim1hrQuantifiedWCPUbNum)
print(UnstimStim1hrQuantifiedWCPUbNum)

#union of proteins with quantificaiton in WCP and Ub for 0-4 hours and 0-1 hours
UnstimStim1hr4hrQuantifiedWCPUb<-union(UnstimStim4hrQuantifiedWCPUb,UnstimStim1hrQuantifiedWCPUb)
UnstimStim1hr4hrQuantifiedWCPUbNum<-length(UnstimStim1hr4hrQuantifiedWCPUb)
#head(UnstimStim1hr4hrQuantifiedWCPUb)
print(UnstimStim1hr4hrQuantifiedWCPUbNum)

##create data frame with wcp and Ub quantification 0-4 and 0-1 hour stim comparison
WCPKeGGCombined<-data.frame()
for( i in 1:UnstimStim1hr4hrQuantifiedWCPUbNum ){
	prot<-UnstimStim1hr4hrQuantifiedWCPUb[i]
	WCPKeGGCombined[prot,"UnstimStim4hrWCP"]<-WCPiBAQdata[prot,"Unstim_Stim4hr_fc"]
	WCPKeGGCombined[prot,"UnstimStim4hrUb"]<-KeGGUbnormLHdata[prot, "intensityweightedmean.PO972kgg_4diff0"]
	WCPKeGGCombined[prot,"UnstimStim1hrWCP"]<-WCPiBAQdata[prot,"Unstim_Stim1hr_fc"]
	WCPKeGGCombined[prot,"UnstimStim1hrUb"]<-KeGGUbnormLHdata[prot, "intensityweightedmean.PO972kgg_1diff0"]
}
#print(WCPKeGGCombined)
#head(WCPKeGGCombined)
dim(WCPKeGGCombined)

##make heatmap
#convert combined dataframe to matrix
WCPKeGGCombined_Mat <- as.matrix(WCPKeGGCombined[,1:ncol(WCPKeGGCombined)])  # transform column 2-5 into a matrix
rownames(WCPKeGGCombined_Mat) <- row.names(WCPKeGGCombined)	# assign row names
#summary(WCPKeGGCombined_Mat)
boxplot(WCPKeGGCombined_Mat, main="unscaled (raw) data", ylab="LH ratio | fold change", cex.axis=0.75, cex.lab=0.75)

#filter to remove outliers
wcp4hr_mean<-mean(WCPKeGGCombined_Mat[,1], na.rm=TRUE)
wcp4hr_sd<-sd(WCPKeGGCombined_Mat[,1], na.rm=TRUE)
ub4hr_mean<-mean(WCPKeGGCombined_Mat[,2], na.rm=TRUE)
ub4hr_sd<-sd(WCPKeGGCombined_Mat[,2], na.rm=TRUE)
wcp1hr_mean<-mean(WCPKeGGCombined_Mat[,3], na.rm=TRUE)
wcp1hr_sd<-sd(WCPKeGGCombined_Mat[,3], na.rm=TRUE)
ub1hr_mean<-mean(WCPKeGGCombined_Mat[,4], na.rm=TRUE)
ub1hr_sd<-sd(WCPKeGGCombined_Mat[,4], na.rm=TRUE)
#
#summary(WCPKeGGCombined_Mat)
for ( rownum in 1:dim(WCPKeGGCombined_Mat)[1] ){
     if ( !is.na(WCPKeGGCombined_Mat[rownum,1]) & (WCPKeGGCombined_Mat[rownum,1] < (wcp4hr_mean-3*wcp4hr_sd) | WCPKeGGCombined_Mat[rownum,1] > (wcp4hr_mean+3*wcp4hr_sd)) ){
	      WCPKeGGCombined_Mat[rownum,1]<-NA
	 }
     if ( !is.na(WCPKeGGCombined_Mat[rownum,2]) & (WCPKeGGCombined_Mat[rownum,2] < (ub4hr_mean-3*ub4hr_sd) | WCPKeGGCombined_Mat[rownum,2] > (ub4hr_mean+3*ub4hr_sd)) ){
	      WCPKeGGCombined_Mat[rownum,2]<-NA
	 }
     if ( !is.na(WCPKeGGCombined_Mat[rownum,3]) & (WCPKeGGCombined_Mat[rownum,3] < (wcp1hr_mean-3*wcp1hr_sd) | WCPKeGGCombined_Mat[rownum,3] > (wcp1hr_mean+3*wcp1hr_sd)) ){
	      WCPKeGGCombined_Mat[rownum,3]<-NA
	 }
     if ( !is.na(WCPKeGGCombined_Mat[rownum,4]) & (WCPKeGGCombined_Mat[rownum,4] < (ub1hr_mean-3*ub1hr_sd) | WCPKeGGCombined_Mat[rownum,4] > (ub1hr_mean+3*ub1hr_sd)) ){
	      WCPKeGGCombined_Mat[rownum,4]<-NA
	 }
}
#summary(WCPKeGGCombined_Mat)
boxplot(WCPKeGGCombined_Mat, main="unscaled (raw) data, outliers removed", ylab="LH ratio | fold change", cex.axis=0.75, cex.lab=0.75)

#proteins with WCP and Ub quantified for 0-4 hour
UnstimStim4hrWCPUbNonNa<-row.names(subset(WCPKeGGCombined_Mat, !is.na(WCPKeGGCombined_Mat[,1]) & !is.na(WCPKeGGCombined_Mat[,2])))
UnstimStim4hrWCPUbNonNaNum<-length(UnstimStim4hrWCPUbNonNa)
print(UnstimStim4hrWCPUbNonNaNum)

#proteins with WCP and Ub quantified for 0-1 hour
UnstimStim1hrWCPUbNonNa<-row.names(subset(WCPKeGGCombined_Mat, !is.na(WCPKeGGCombined_Mat[,3]) & !is.na(WCPKeGGCombined_Mat[,4])))
UnstimStim1hrWCPUbNonNaNum<-length(UnstimStim1hrWCPUbNonNa)
print(UnstimStim1hrWCPUbNonNaNum)

#heatmap(mat_data_scaled_filtered_trans_subset, Colv=NA, scale="none")
quantile_range <- quantile(WCPKeGGCombined_Mat, probs = seq(0, 1, 0.01), na.rm=TRUE)
col_breaks <- seq(quantile_range["5%"], quantile_range["95%"], 0.01)
col_palette  <- colorRampPalette(c("blue", "gray", "red"))(length(col_breaks) - 1)
heatmap.2(WCPKeGGCombined_Mat[UnstimStim4hrWCPUbNonNa,1:2],
		  main = "unscaled (raw) data, outliers removed, 0-4hr Stim",
	      notecol="black",      # change font color of cell labels to black
	      density.info="histogram",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  margins =c(9,12),     # widens margins around plot
		  dendrogram="row",     # only draw a row dendrogram
		  Colv="NA",            # turn off column clustering
		  col=col_palette,
		  scale="none",
		  cexRow=0.1,
		  cexCol=0.8,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1.3,
		  sepwidth=c(0.0, 0.0)
		 )
heatmap.2(WCPKeGGCombined_Mat[UnstimStim1hrWCPUbNonNa,3:4],
		  main = "unscaled (raw) data, outliers removed, 0-1hr Stim",
	      notecol="black",      # change font color of cell labels to black
	      density.info="histogram",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  margins =c(9,12),     # widens margins around plot
		  dendrogram="row",     # only draw a row dendrogram
		  Colv="NA",            # turn off column clustering
		  col=col_palette,
		  scale="none",
		  cexRow=0.1,
		  cexCol=0.8,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1.3,
		  sepwidth=c(0.0, 0.0)
		 )

#rescale to between -1 and 1 to correct for different SILAC and iBAQ fold change ranges
WCPKeGGCombined_Mat_RescaledRange <- matrix(nrow=length(row.names(WCPKeGGCombined_Mat)), ncol = dim(WCPKeGGCombined_Mat)[2])
rownames(WCPKeGGCombined_Mat_RescaledRange) <- row.names(WCPKeGGCombined_Mat) # assign row names
colnames(WCPKeGGCombined_Mat_RescaledRange) <- colnames(WCPKeGGCombined_Mat) #assign col names
WCPKeGGCombined_Mat_RescaledRange[,1] <- rescale(WCPKeGGCombined_Mat[,1], to=c(-1,1))
WCPKeGGCombined_Mat_RescaledRange[,2] <- rescale(WCPKeGGCombined_Mat[,2], to=c(-1,1))
WCPKeGGCombined_Mat_RescaledRange[,3] <- rescale(WCPKeGGCombined_Mat[,3], to=c(-1,1))
WCPKeGGCombined_Mat_RescaledRange[,4] <- rescale(WCPKeGGCombined_Mat[,4], to=c(-1,1))
#summary(WCPKeGGCombined_Mat_RescaledRange)
boxplot(WCPKeGGCombined_Mat_RescaledRange, main="scaled: range -1 to 1", ylab="LH ratio | fold change", cex.axis=0.75, cex.lab=0.75)

#heatmap(mat_data_scaled_filtered_trans_subset, Colv=NA, scale="none")
quantile_range <- quantile(WCPKeGGCombined_Mat_RescaledRange, probs = seq(0, 1, 0.01), na.rm=TRUE)
col_breaks <- seq(quantile_range["5%"], quantile_range["95%"], 0.01)
col_palette  <- colorRampPalette(c("blue", "gray", "red"))(length(col_breaks) - 1)
heatmap.2(WCPKeGGCombined_Mat_RescaledRange[UnstimStim4hrWCPUbNonNa,1:2],
		  main = "scaled: range -1 to 1, 0-4hr Stim",
	      notecol="black",      # change font color of cell labels to black
	      density.info="histogram",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  margins =c(9,12),     # widens margins around plot
		  dendrogram="row",     # only draw a row dendrogram
		  Colv="NA",            # turn off column clustering
		  col=col_palette,
		  scale="none",
		  cexRow=0.1,
		  cexCol=0.8,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1.3
		 )
heatmap.2(WCPKeGGCombined_Mat_RescaledRange[UnstimStim1hrWCPUbNonNa,3:4],
		  main = "scaled: range -1 to 1, 0-1hr Stim",
	      notecol="black",      # change font color of cell labels to black
	      density.info="histogram",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  margins =c(9,12),     # widens margins around plot
		  dendrogram="row",     # only draw a row dendrogram
		  Colv="NA",            # turn off column clustering
		  col=col_palette,
		  scale="none",
		  cexRow=0.1,
		  cexCol=0.8,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1.3
		 )

#rescale to mean=0 and sd=1 to correct for different SILAC and iBAQ fold change ranges
WCPKeGGCombined_Mat_RescaledMeanSd <- matrix(nrow=length(row.names(WCPKeGGCombined_Mat)), ncol = dim(WCPKeGGCombined_Mat)[2])
rownames(WCPKeGGCombined_Mat_RescaledMeanSd) <- row.names(WCPKeGGCombined_Mat) # assign row names
colnames(WCPKeGGCombined_Mat_RescaledMeanSd) <- colnames(WCPKeGGCombined_Mat) #assign col names
WCPKeGGCombined_Mat_RescaledMeanSd <- scale(WCPKeGGCombined_Mat)
#summary(WCPKeGGCombined_Mat_RescaledMeanSd)
boxplot(WCPKeGGCombined_Mat_RescaledMeanSd, main="scaled: mean=0, sd=1", ylab="LH ratio | fold change", cex.axis=0.75, cex.lab=0.75)

#heatmap(mat_data_scaled_filtered_trans_subset, Colv=NA, scale="none")
quantile_range <- quantile(WCPKeGGCombined_Mat_RescaledMeanSd, probs = seq(0, 1, 0.01), na.rm=TRUE)
col_breaks <- seq(quantile_range["5%"], quantile_range["95%"], 0.01)
col_palette  <- colorRampPalette(c("blue", "gray", "red"))(length(col_breaks) - 1)
heatmap.2(WCPKeGGCombined_Mat_RescaledMeanSd[UnstimStim4hrWCPUbNonNa,1:2],
	      main = "scaled: mean=0 sd=1, 0-4hr Stim",
		  notecol="black",      # change font color of cell labels to black
	      density.info="histogram",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  margins =c(9,12),     # widens margins around plot
		  dendrogram="row",     # only draw a row dendrogram
		  Colv="NA",            # turn off column clustering
		  col=col_palette,
		  scale="none",
		  cexRow=0.1,
		  cexCol=0.8,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1.3
		 )
heatmap.2(WCPKeGGCombined_Mat_RescaledMeanSd[UnstimStim1hrWCPUbNonNa,3:4],
	      main = "scaled: mean=0 sd=1, 0-1hr Stim",
		  notecol="black",      # change font color of cell labels to black
	      density.info="histogram",  # turns off density plot inside color legend
		  trace="none",         # turns off trace lines inside the heat map
		  margins =c(9,12),     # widens margins around plot
		  dendrogram="row",     # only draw a row dendrogram
		  Colv="NA",            # turn off column clustering
		  col=col_palette,
		  scale="none",
		  cexRow=0.1,
		  cexCol=0.8,
		  vline=NULL,
		  hline=NULL,
		  key=TRUE,
		  keysize=1.3
		 )

###clustering
#nc<-50
#seed<-1234
#wss <- (nrow(WCPKeGGCombined_Mat)-1)*sum(apply(WCPKeGGCombined_Mat,2,var))
#for (i in 2:nc){
#set.seed(seed)
#wss[i] <- sum(kmeans(WCPKeGGCombined_Mat, centers=i)$withinss)
#}
#plot(1:nc, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")


#kmeans(WCPKeGGCombined_Mat, centers=i, 


#set seed to get same clustering result each time
seed<-1234
set.seed(seed)


#clusters...
#
#####WCPKeGGCombined_4hrClus<-WCPKeGGCombined_Mat[,1:2]
#####head(WCPKeGGCombined_4hrClus)
#####WCPKeGGCombined_4hrClus<-WCPKeGGCombined_4hrClus[WCPKeGGCombined_4hrClus[,1] is.na(]
#####stop("stopped")

numclus<-18
ClusMatUnstim4hr<-WCPKeGGCombined_Mat[UnstimStim4hrWCPUbNonNa,1:2]
dim(ClusMatUnstim4hr)

ClusMatUnstim1hr<-WCPKeGGCombined_Mat[UnstimStim1hrWCPUbNonNa,3:4]
dim(ClusMatUnstim1hr)

##########cl04hr <- kmeans(WCPKeGGCombined_Mat[UnstimStim4hrWCPUbNonNa,1:2], centers=numclus, nstart = 25)
cl04hr <- kmeans(ClusMatUnstim4hr, centers=numclus, nstart = 25)
##########plot(WCPKeGGCombined_Mat[UnstimStim4hrWCPUbNonNa,1],WCPKeGGCombined_Mat[UnstimStim4hrWCPUbNonNa,2],
plot(ClusMatUnstim4hr[,1], ClusMatUnstim4hr[,2],
     col = cl04hr$cluster, pch = cl04hr$cluster,
	 xlim=c(-1.75,1.75),
	 main="kmeans 18 clusters, 0-4hr Stim", xlab="unstim-stim4hr WCP", ylab="unstim-stim4hr Ub")
#points(cl$centers, col = 1:numclus, pch = 8)
abline(h=0,v=0)

##########cl01hr <- kmeans(WCPKeGGCombined_Mat[UnstimStim1hrWCPUbNonNa,3:4], centers=numclus, nstart = 25)
cl01hr <- kmeans(ClusMatUnstim1hr, centers=numclus, nstart = 25)
##########plot(WCPKeGGCombined_Mat[UnstimStim1hrWCPUbNonNa,3],WCPKeGGCombined_Mat[UnstimStim1hrWCPUbNonNa,4],
plot(ClusMatUnstim1hr[,1], ClusMatUnstim1hr[,2],
     col = cl01hr$cluster, pch = cl01hr$cluster,
	 xlim=c(-1.75,1.75),
	 main="kmeans 18 clusters, 0-1hr Stim", xlab="unstim-stim1hr WCP", ylab="unstim-stim1hr Ub")
#points(cl$centers, col = 1:numclus, pch = 8)
abline(h=0,v=0)

stop("stopped")

cl04hrclus<-cl04hr$cluster
Clus04Stim<-data.frame(row.names=row.names(WCPKeGGCombined))
Clus04Stim[,"clusterNum4hrStim"]<-NA
for ( clusnum in 1:18 ){
	clusmems<-row.names(ClusMatUnstim4hr[cl04hr$cluster==clusnum,])
	#print(clusmems)
	for ( i in 1:length(clusmems) ){
		prot<-clusmems[i]
		Clus04Stim[prot,"clusterNum4hrStim"]<-clusnum
	}
}
print(Clus04Stim)

cl01hrclus<-cl01hr$cluster
Clus01Stim<-data.frame(row.names=row.names(WCPKeGGCombined))
Clus01Stim[,"clusterNum1hrStim"]<-NA
for ( clusnum in 1:18 ){
	clusmems<-row.names(ClusMatUnstim1hr[cl01hr$cluster==clusnum,])
	#print(clusmems)
	for ( i in 1:length(clusmems) ){
		prot<-clusmems[i]
		Clus01Stim[prot,"clusterNum1hrStim"]<-clusnum
	}
}
print(Clus01Stim)



#attach cluster numbers to WCPKeGGCombined dataframe and write
WCPKeGGCombined <- cbind(WCPKeGGCombined, Clus04Stim)
WCPKeGGCombined <- cbind(WCPKeGGCombined, Clus01Stim)
dim(WCPKeGGCombined)
print(WCPKeGGCombined)
write.table(WCPKeGGCombined, file="WCP0-4hr_Ub0-4hr_WCP0-1hr_Ub0-1hr_QuantificationAndClustering_TEST.txt", sep="\t", quote=FALSE, na="NA", dec=".")
#print(WCPKeGGCombined_Mat)
#

stop("stopped")

#
numclus<-15
cl <- kmeans(WCPKeGGCombined_Mat[,1:2], centers=numclus, nstart = 25)
plot(WCPKeGGCombined_Mat[,1],WCPKeGGCombined_Mat[,2],
     col = cl$cluster, pch = cl$cluster,
	 main="kmeans 15 clusters", xlab="unstim-stim4hr WCP", ylab="unstim-stim4hr Ub")
points(cl$centers, col = 1:numclus, pch = 8)
abline(h=0,v=0)
#
numclus<-12
cl <- kmeans(WCPKeGGCombined_Mat[,1:2], centers=numclus, nstart = 25)
plot(WCPKeGGCombined_Mat[,1],WCPKeGGCombined_Mat[,2],
     col = cl$cluster, pch = cl$cluster,
	 main="kmeans 12 clusters", xlab="unstim-stim4hr WCP", ylab="unstim-stim4hr Ub")
points(cl$centers, col = 1:numclus, pch = 8)
abline(h=0,v=0)
#
numclus<-9
cl <- kmeans(WCPKeGGCombined_Mat[,1:2], centers=numclus, nstart = 25)
plot(WCPKeGGCombined_Mat[,1],WCPKeGGCombined_Mat[,2],
     col = cl$cluster, pch = cl$cluster,
	 main="kmeans 9 clusters", xlab="unstim-stim4hr WCP", ylab="unstim-stim4hr Ub")
points(cl$centers, col = 1:numclus, pch = 8)
abline(h=0,v=0)

plot(WCPKeGGCombined_Mat[,3],WCPKeGGCombined_Mat[,2],
     col = WCPKeGGCombined_Mat[,4], pch = WCPKeGGCombined_Mat[,4],
	 main="kmeans 18 clusters", xlab="unstim-stim1hr WCP", ylab="unstim-stim4hr Ub")
abline(h=0,v=0)

