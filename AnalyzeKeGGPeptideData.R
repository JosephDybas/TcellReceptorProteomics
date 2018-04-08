

##disable scientific notation
options(scipen = 999)


source("/Users/dybasj/JoeRLib/ProteomicsAnalysisFunctions.R")
source("/Users/dybasj/JoeRLib/PlottingFunctions.R")


##calculate 4hr-0hr ratios
diffcomp <- function(x,y)
{

	if ( is.na(x) | is.na(y) ){
		result<-NA
	}else{
		result<-as.numeric(x)-as.numeric(y)
	}
	
	return (result)

}
##end calculate 4hr-0hr ratios

#####read data and combine into "peptidedata" dataframe
##read Log2 transformed, normalized H/L ratio quantification data
##formatted with experiments as columns and samples as rows
HLratiodata <- read.table("./datafiles/KeGGData_RatioHLnorm_SILAC_Log2Norm.txt", sep="\t", row.names=1, header=TRUE)
#since the MaxQuant ratios are H/L and L/H is informative for fold changes, multiply log2 scores by -1 (inverse)
HLratiodata <- -(HLratiodata)
##change column names to specify HLratio
HLratiodatacolnames <- colnames(HLratiodata)
HLratiodatacolnamesadj <- paste("HLratioLog2Norm", HLratiodatacolnames, sep=".")
colnames(HLratiodata) <- HLratiodatacolnamesadj
#dim(HLratiodata)
#head(HLratiodata)


##read norm H/L ratio count data
##formatted with experiments as columns and samples as rows
HLcountdata <- read.table("./datafiles/KeGGData_RatioHLcount_SILAC_Raw.txt", sep="\t", row.names=1, header=TRUE)
##change column names to specify HLcount
HLcountdatacolnames <- colnames(HLcountdata)
HLcountdatacolnamesadj <- paste("HLcount", HLcountdatacolnames, sep=".")
colnames(HLcountdata) <- HLcountdatacolnamesadj
#dim(HLcountdata)
#head(HLcountdata)


##read Log2 transformed intensity quantification data
##formatted with experiments as columns and samples as rows
IntensityLog2data <- read.table("./datafiles/KeGGData_Intensity_SILAC_Log2.txt", sep="\t", row.names=1, header=TRUE)
##change column names to specify Intensity
intensityLog2datacolnames <- colnames(IntensityLog2data)
intensityLog2datacolnamesadj <- paste("IntensityLog2", intensityLog2datacolnames, sep=".")
colnames(IntensityLog2data) <- intensityLog2datacolnamesadj
#dim(IntensityLog2data)
#head(IntensityLog2data)


##read Log2 transformed, normalized intensity quantification data
##formatted with experiments as columns and samples as rows
IntensityLog2Normdata <- read.table("./datafiles/KeGGData_Intensity_SILAC_Log2Norm.txt", sep="\t", row.names=1, header=TRUE)
##change column names to specify Intensity
intensitylog2normdatacolnames <- colnames(IntensityLog2Normdata)
intensitylog2normdatacolnamesadj <- paste("IntensityLog2Norm", intensitylog2normdatacolnames, sep=".")
colnames(IntensityLog2Normdata) <- intensitylog2normdatacolnamesadj
#dim(IntensityLog2Normdata)
#head(IntensityLog2Normdata)


##read Log2 transformed, normalized intensity L quantification data
##formatted with experiments as columns and samples as rows
IntensityLLog2Normdata <- read.table("./datafiles/KeGGData_IntensityL_SILAC_Log2Norm.txt", sep="\t", row.names=1, header=TRUE)
##change column names to specify Intensity
intensityLlog2normdatacolnames <- colnames(IntensityLLog2Normdata)
intensityLlog2normdatacolnamesadj <- paste("IntensityLLog2Norm", intensityLlog2normdatacolnames, sep=".")
colnames(IntensityLLog2Normdata) <- intensityLlog2normdatacolnamesadj
#dim(IntensityLLog2Normdata)
#head(IntensityLLog2Normdata)



##read Log2 transformed, normalized intensity H quantification data
##formatted with experiments as columns and samples as rows
IntensityHLog2Normdata <- read.table("./datafiles/KeGGData_IntensityH_SILAC_Log2Norm.txt", sep="\t", row.names=1, header=TRUE)
##change column names to specify Intensity
intensityHlog2normdatacolnames <- colnames(IntensityHLog2Normdata)
intensityHlog2normdatacolnamesadj <- paste("IntensityHLog2Norm", intensityHlog2normdatacolnames, sep=".")
colnames(IntensityHLog2Normdata) <- intensityHlog2normdatacolnamesadj
#dim(IntensityHLog2Normdata)
#head(IntensityHLog2Normdata)


##combine dataframes
KeGGdata <- merge(HLratiodata, HLcountdata, by="row.names", all=TRUE)
rownames(KeGGdata)=KeGGdata$Row.names
KeGGdata$Row.names <- NULL
KeGGdata <- merge(KeGGdata, IntensityLog2data, by="row.names", all=TRUE)
rownames(KeGGdata)=KeGGdata$Row.names
KeGGdata$Row.names <- NULL
KeGGdata <- merge(KeGGdata, IntensityLog2Normdata, by="row.names", all=TRUE)
rownames(KeGGdata)=KeGGdata$Row.names
KeGGdata$Row.names <- NULL
KeGGdata <- merge(KeGGdata, IntensityLLog2Normdata, by="row.names", all=TRUE)
rownames(KeGGdata)=KeGGdata$Row.names
KeGGdata$Row.names <- NULL
KeGGdata <- merge(KeGGdata, IntensityHLog2Normdata, by="row.names", all=TRUE)
rownames(KeGGdata)=KeGGdata$Row.names
KeGGdata$Row.names <- NULL
#dim(KeGGdata)
#head(KeGGdata)


##split unique proteinID_position row names into protein ID column and site position column
proteinIDsitenames<-rownames(KeGGdata)
for ( i in 1:length(proteinIDsitenames) ){
	currentID <- proteinIDsitenames[i]
	tmpID <- strsplit(currentID, "_")
	protID <- tmpID[[1]][1]
	pos <- tmpID[[1]][2]
	KeGGdata[currentID,"protID"]=protID
	KeGGdata[currentID,"position"]=pos
}
dim(KeGGdata)
#head(KeGGdata)

##read iBAQ quantification protein IDs
ProtIDdata <- read.table("./datafiles/ProteinGroupData_ProteinIdGeneIdProteinName_LFQ.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors=FALSE)
ProtIDdata$protID<-row.names(ProtIDdata)
dim(ProtIDdata)
#head(ProtIDdata)


##combine with KeGGdata dataframe
#####KeGGdata <- merge(KeGGdata, ProtIDdata, by=c("protID"), all.x=TRUE)
#####rownames(KeGGdata)=KeGGdata$Row.names
#####KeGGdata$Row.names <- NULL
#####dim(KeGGdata)
#####head(KeGGdata)


KeGGdata$GeneId<-apply(KeGGdata, MARGIN=1, FUN=function(x) ProtIDdata[x[["protID"]],"GeneId"])
KeGGdata$ProtName<-apply(KeGGdata, MARGIN=1, FUN=function(x) ProtIDdata[x[["protID"]],"ProtName"])
#dim(KeGGdata)


#######merging causes some added lines to have NA where there should be 0
######change NA to 0 where necessary
#####changezero<-c("HLcount.PO1063kgg_HLM", "HLcount.PO1063kgg_HLP", "HLcount.PO972kgg_0hr", "HLcount.PO972kgg_1hr", "HLcount.PO972kgg_4hr")
#####KeGGdata[changezero][is.na(KeGGdata[changezero])] <- 0
######head(KeGGdata)
####### end change NA to 0 where necessary


##calculate the peptide based HL ratio change for 4hr compared to 0hr
KeGGdata$HLratioLog2Norm.PO972kgg_4diff0<-apply(KeGGdata, MARGIN=1, FUN=function(x) diffcomp(x[["HLratioLog2Norm.PO972kgg_4hr"]],x[["HLratioLog2Norm.PO972kgg_0hr"]]))

##calculate the peptide based HL ratio average for PO1063kggHLM and PO972kgg4diff0
KeGGdata$HLratioLog2Norm.PO0163kggHLMPO972kgg4diff0_AvgHLratio<-apply(KeGGdata, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["HLratioLog2Norm.PO1063kgg_HLM"]],x[["HLratioLog2Norm.PO972kgg_4diff0"]])), na.rm=TRUE))


###peptide set analysis
#PO1063HLM Ub peptides identified at 0hr and 4hr
PO1063HLMkgg_0hr4hr_UbPep <- row.names(subset(KeGGdata, !is.na(KeGGdata$HLratioLog2Norm.PO1063kgg_HLM)))
PO1063HLMkgg_0hr4hr_UbPep_Num <- length(PO1063HLMkgg_0hr4hr_UbPep)
print("PO1063HLMkgg 0hr-4hr fold change peptides identified")
print(PO1063HLMkgg_0hr4hr_UbPep_Num)
#
#PO1063HLM Ub peptides identified at 0hr unique
PO1063HLMkgg_0hr_Unique_UbPep <- row.names(subset(KeGGdata, (!is.na(KeGGdata$IntensityHLog2Norm.PO1063kgg_HLM) & is.na(KeGGdata$HLratioLog2Norm.PO1063kgg_HLM))))
PO1063HLMkgg_0hr_Unique_UbPep_Num <- length(PO1063HLMkgg_0hr_Unique_UbPep)
print("PO1063HLMkgg 0hr unique peptides identified")
print(PO1063HLMkgg_0hr_Unique_UbPep_Num)
#
#PO1063HLM Ub peptides identified at 4hr unique
PO1063HLMkgg_4hr_Unique_UbPep <- row.names(subset(KeGGdata, (!is.na(KeGGdata$IntensityLLog2Norm.PO1063kgg_HLM) & is.na(KeGGdata$HLratioLog2Norm.PO1063kgg_HLM))))
PO1063HLMkgg_4hr_Unique_UbPep_Num <- length(PO1063HLMkgg_4hr_Unique_UbPep)
print("PO1063HLMkgg 4hr unique peptides identified")
print(PO1063HLMkgg_4hr_Unique_UbPep_Num)
#
#PO1063HLM all identified Ub peptides
PO1063HLMkgg_All_UbPep <- Reduce(union, list(PO1063HLMkgg_0hr_Unique_UbPep, PO1063HLMkgg_4hr_Unique_UbPep, PO1063HLMkgg_0hr4hr_UbPep))
PO1063HLMkgg_All_UbPep_Num <- length(PO1063HLMkgg_All_UbPep)
print("PO1063HLMkgg all Ub peptides identified")
print(PO1063HLMkgg_All_UbPep_Num)
#
#PO972 Ub peptides identified at 0hr and 4hr
PO972kgg_0hr4hr_UbPep <- row.names(subset(KeGGdata, !is.na(KeGGdata$HLratioLog2Norm.PO972kgg_4diff0)))
PO972kgg_0hr4hr_UbPep_Num <- length(PO972kgg_0hr4hr_UbPep)
print("PO972kgg 0hr-4hr fold change peptides identified")
print(PO972kgg_0hr4hr_UbPep_Num)
#
#PO972 Ub peptides identified at 0hr unique
PO972kgg_0hr_Unique_UbPep <- row.names(subset(KeGGdata, (!is.na(KeGGdata$IntensityLLog2Norm.PO972kgg_0hr) & is.na(KeGGdata$HLratioLog2Norm.PO972kgg_4diff0))))
PO972kgg_0hr_Unique_UbPep_Num <- length(PO972kgg_0hr_Unique_UbPep)
print("PO972kgg 0hr unique peptides identified")
print(PO972kgg_0hr_Unique_UbPep_Num)
#
#PO972 Ub peptides identified at 4hr unique
PO972kgg_4hr_Unique_UbPep <- row.names(subset(KeGGdata, (!is.na(KeGGdata$IntensityLLog2Norm.PO972kgg_4hr) & is.na(KeGGdata$HLratioLog2Norm.PO972kgg_4diff0))))
PO972kgg_4hr_Unique_UbPep_Num <- length(PO972kgg_4hr_Unique_UbPep)
print("PO972kgg 4hr unique peptides identified")
print(PO972kgg_4hr_Unique_UbPep_Num)
#
#PO972 all identified Ub peptides
PO972kgg_All_UbPep <- Reduce(union, list(PO972kgg_0hr_Unique_UbPep, PO972kgg_4hr_Unique_UbPep, PO972kgg_0hr4hr_UbPep))
PO972kgg_All_UbPep_Num <- length(PO972kgg_All_UbPep)
print("PO972kgg all Ub peptides identified")
print(PO972kgg_All_UbPep_Num)
#
#PO1063HLM-PO972 all identified peptides union
PO972kggPO1063HLMkgg_All_Union_UbPep <- union(PO1063HLMkgg_All_UbPep, PO972kgg_All_UbPep)
PO972kggPO1063HLMkgg_All_Union_UbPep_Num <- length(PO972kggPO1063HLMkgg_All_Union_UbPep)
print("PO972kgg-PO1063HLMkgg all Ub peptides identified (union)")
print(PO972kggPO1063HLMkgg_All_Union_UbPep_Num)
#
#PO1063HLM-PO972 all identified peptides union
PO972kggPO1063HLMkgg_All_Int_UbPep <- intersect(PO1063HLMkgg_All_UbPep, PO972kgg_All_UbPep)
PO972kggPO1063HLMkgg_All_Int_UbPep_Num <- length(PO972kggPO1063HLMkgg_All_Int_UbPep)
print("PO972kgg-PO1063HLMkgg all common Ub peptides identified (intersection)")
print(PO972kggPO1063HLMkgg_All_Int_UbPep_Num)
#
#PO1063HLM-PO972 all quantified peptide with 0-4hr fold change union
PO972kggPO1063HLMkgg_0hr4hr_Union_UbPep <- union(PO1063HLMkgg_0hr4hr_UbPep, PO972kgg_0hr4hr_UbPep)
PO972kggPO1063HLMkgg_0hr4hr_Union_UbPep_Num <- length(PO972kggPO1063HLMkgg_0hr4hr_Union_UbPep)
print("PO972kgg-PO1063HLMkgg all Ub peptides quantified with 0hr4hr fold change (union)")
print(PO972kggPO1063HLMkgg_0hr4hr_Union_UbPep_Num)
#

#PO1063HLM-PO972 all quantified peptide with 0-4hr fold change intersection
PO972kggPO1063HLMkgg_0hr4hr_Int_UbPep <- intersect(PO1063HLMkgg_0hr4hr_UbPep, PO972kgg_0hr4hr_UbPep)
PO972kggPO1063HLMkgg_0hr4hr_Int_UbPep_Num <- length(PO972kggPO1063HLMkgg_0hr4hr_Int_UbPep)
print("PO972kgg-PO1063HLMkgg all Ub peptides quantified with 0hr4hr fold change (intersection)")
print(PO972kggPO1063HLMkgg_0hr4hr_Int_UbPep_Num)
#
###

##
PO972kgg_UbPep_Data <- KeGGdata[row.names(KeGGdata)%in%PO972kgg_All_UbPep,]
PO972kgg_redundantprots <- PO972kgg_UbPep_Data$protID
PO972kgg_Prots <- unique(PO972kgg_redundantprots)
PO972kgg_Prots_Num <- length(PO972kgg_Prots)
print(PO972kgg_Prots_Num)
#
PO1063HLMkgg_UbPep_Data <- KeGGdata[row.names(KeGGdata)%in%PO1063HLMkgg_All_UbPep,]
PO1063HLMkgg_redundantprots <- PO1063HLMkgg_UbPep_Data$protID
PO1063HLMkgg_Prots <- unique(PO1063HLMkgg_redundantprots)
PO1063HLMkgg_Prots_Num <- length(PO1063HLMkgg_Prots)
print(PO1063HLMkgg_Prots_Num)
#
PO972kggPO1063HLMkgg_Prot_Int <- intersect(PO972kgg_Prots,PO1063HLMkgg_Prots)
PO972kggPO1063HLMkgg_Prot_Int_Num <- length(PO972kggPO1063HLMkgg_Prot_Int)
print("PO972kggPO1063HLMkgg_Prot_Int_Num")
print(PO972kggPO1063HLMkgg_Prot_Int_Num)
#
PO972kggPO1063HLMkgg_Prot_Union <- union(PO972kgg_Prots,PO1063HLMkgg_Prots)
PO972kggPO1063HLMkgg_Prot_Union_Num <- length(PO972kggPO1063HLMkgg_Prot_Union)
print("PO972kggPO1063HLMkgg_Prot_Union_Num")
print(PO972kggPO1063HLMkgg_Prot_Union_Num)
#
PO972kgg_Prot_Unique <- setdiff(PO972kgg_Prots, PO1063HLMkgg_Prots)
PO972kgg_Prot_Unique_Num <- length(PO972kgg_Prot_Unique)
print("PO972kgg_Prot_Unique_Num")
print(PO972kgg_Prot_Unique_Num)
#
PO1063HLMkgg_Prot_Unique <- setdiff(PO1063HLMkgg_Prots, PO972kgg_Prots)
PO1063HLMkgg_Prot_Unique_Num <- length(PO1063HLMkgg_Prot_Unique)
print("PO1063HLMkgg_Prot_Unique_Num")
print(PO1063HLMkgg_Prot_Unique_Num)


##number of proteins and number of peptides per protein
PO972kggPO1063HLMkgg_UbPep_Data <- KeGGdata[row.names(KeGGdata)%in%PO972kggPO1063HLMkgg_All_Union_UbPep,]
dim(PO972kggPO1063HLMkgg_UbPep_Data)
allproteins <- PO972kggPO1063HLMkgg_UbPep_Data$protID
allproteins_num <- length(allproteins)
print("allproteins_num")
print(allproteins_num)
uniqueproteins <- unique(allproteins)
uniqueproteins_num <- length(uniqueproteins)
print("uniqueproteins_num")
print(uniqueproteins_num)


write.table(PO972kggPO1063HLMkgg_UbPep_Data, file="KeGGpeptide_SupplementalTable_2.txt", sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)


##venn diagram of number of Ub peptides for PO972, PO1063HLM and common
plotpairwisevenn(PO972kgg_All_UbPep_Num,
                 PO1063HLMkgg_All_UbPep_Num,
				 PO972kggPO1063HLMkgg_All_Int_UbPep_Num,
				 "", "", "mediumpurple", "mediumpurple1", "mediumpurple", "mediumpurple1", TRUE,
				 "PO972kgg-PO1063HLMkgg_NumUbPeptides_VennDiagram.pdf")
##end venn diagram of number of Ub peptides for PO972, PO1063HLM and common


##zscores
KeGGdataSub<-subset(KeGGdata, select=c(IntensityLLog2Norm.PO1063kgg_HLM, IntensityLLog2Norm.PO1063kgg_HLP, IntensityLLog2Norm.PO972kgg_0hr, IntensityLLog2Norm.PO972kgg_1hr, IntensityLLog2Norm.PO972kgg_4hr, IntensityHLog2Norm.PO1063kgg_HLM, IntensityHLog2Norm.PO1063kgg_HLP, IntensityHLog2Norm.PO972kgg_0hr, IntensityHLog2Norm.PO972kgg_1hr, IntensityHLog2Norm.PO972kgg_4hr))
KeGGdataSub_colmeans <- apply(KeGGdataSub, MARGIN=2, FUN=function(x) mean(x, na.rm=TRUE))
KeGGdataSub_colstdev <- apply(KeGGdataSub, MARGIN=2, FUN=function(x) sd(x, na.rm=TRUE))
KeGGdata$IntensityLLog2Norm.PO1063kgg_HLM.zscore<-apply(KeGGdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["IntensityLLog2Norm.PO1063kgg_HLM"]]),KeGGdataSub_colmeans[["IntensityLLog2Norm.PO1063kgg_HLM"]],KeGGdataSub_colstdev[["IntensityLLog2Norm.PO1063kgg_HLM"]]))
KeGGdata$IntensityLLog2Norm.PO1063kgg_HLP.zscore<-apply(KeGGdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["IntensityLLog2Norm.PO1063kgg_HLP"]]),KeGGdataSub_colmeans[["IntensityLLog2Norm.PO1063kgg_HLP"]],KeGGdataSub_colstdev[["IntensityLLog2Norm.PO1063kgg_HLP"]]))
KeGGdata$IntensityLLog2Norm.PO972kgg_0hr.zscore<-apply(KeGGdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["IntensityLLog2Norm.PO972kgg_0hr"]]),KeGGdataSub_colmeans[["IntensityLLog2Norm.PO972kgg_0hr"]],KeGGdataSub_colstdev[["IntensityLLog2Norm.PO972kgg_0hr"]]))
KeGGdata$IntensityLLog2Norm.PO972kgg_1hr.zscore<-apply(KeGGdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["IntensityLLog2Norm.PO972kgg_1hr"]]),KeGGdataSub_colmeans[["IntensityLLog2Norm.PO972kgg_1hr"]],KeGGdataSub_colstdev[["IntensityLLog2Norm.PO972kgg_1hr"]]))
KeGGdata$IntensityLLog2Norm.PO972kgg_4hr.zscore<-apply(KeGGdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["IntensityLLog2Norm.PO972kgg_4hr"]]),KeGGdataSub_colmeans[["IntensityLLog2Norm.PO972kgg_4hr"]],KeGGdataSub_colstdev[["IntensityLLog2Norm.PO972kgg_4hr"]]))
KeGGdata$IntensityHLog2Norm.PO1063kgg_HLM.zscore<-apply(KeGGdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["IntensityHLog2Norm.PO1063kgg_HLM"]]),KeGGdataSub_colmeans[["IntensityHLog2Norm.PO1063kgg_HLM"]],KeGGdataSub_colstdev[["IntensityHLog2Norm.PO1063kgg_HLM"]]))
KeGGdata$IntensityHLog2Norm.PO1063kgg_HLP.zscore<-apply(KeGGdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["IntensityHLog2Norm.PO1063kgg_HLP"]]),KeGGdataSub_colmeans[["IntensityHLog2Norm.PO1063kgg_HLP"]],KeGGdataSub_colstdev[["IntensityHLog2Norm.PO1063kgg_HLP"]]))
KeGGdata$IntensityHLog2Norm.PO972kgg_0hr.zscore<-apply(KeGGdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["IntensityHLog2Norm.PO972kgg_0hr"]]),KeGGdataSub_colmeans[["IntensityHLog2Norm.PO972kgg_0hr"]],KeGGdataSub_colstdev[["IntensityHLog2Norm.PO972kgg_0hr"]]))
KeGGdata$IntensityHLog2Norm.PO972kgg_1hr.zscore<-apply(KeGGdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["IntensityHLog2Norm.PO972kgg_1hr"]]),KeGGdataSub_colmeans[["IntensityHLog2Norm.PO972kgg_1hr"]],KeGGdataSub_colstdev[["IntensityHLog2Norm.PO972kgg_1hr"]]))
KeGGdata$IntensityHLog2Norm.PO972kgg_4hr.zscore<-apply(KeGGdata, MARGIN=1, FUN=function(x) zscores(as.numeric(x[["IntensityHLog2Norm.PO972kgg_4hr"]]),KeGGdataSub_colmeans[["IntensityHLog2Norm.PO972kgg_4hr"]],KeGGdataSub_colstdev[["IntensityHLog2Norm.PO972kgg_4hr"]]))
KeGGdata$IntensityHLog2Norm.PO972kggPO1063kggHLM.avgzscore<-apply(KeGGdata, MARGIN=1, FUN=function(x) mean(as.numeric(c(x[["IntensityHLog2Norm.PO1063kgg_HLM.zscore"]],x[["IntensityHLog2Norm.PO972kgg_0hr.zscore"]])),na.rm=TRUE))
##write output of peptide level KeGG data analysis
#directory to which output is written
analysisdir<-"./"
#write output
write.table(KeGGdata, file=file.path(".", analysisdir, "KeGGpeptideData_PositionBasedHLRatioAndIntensityQuantifications_PO1063-PO972_TEST.txt"), sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)

PO1063kggHLMPO972kgg0hr_UbKeGGLysines_AvgZscore<-c(KeGGdata["P62983_6","IntensityHLog2Norm.PO972kggPO1063kggHLM.avgzscore"],
                                                   KeGGdata["P62983_11","IntensityHLog2Norm.PO972kggPO1063kggHLM.avgzscore"],
												   KeGGdata["P62983_27","IntensityHLog2Norm.PO972kggPO1063kggHLM.avgzscore"],
												   KeGGdata["P62983_29","IntensityHLog2Norm.PO972kggPO1063kggHLM.avgzscore"],
												   KeGGdata["P62983_33","IntensityHLog2Norm.PO972kggPO1063kggHLM.avgzscore"],
												   KeGGdata["P62983_48","IntensityHLog2Norm.PO972kggPO1063kggHLM.avgzscore"],
												   KeGGdata["P62983_63","IntensityHLog2Norm.PO972kggPO1063kggHLM.avgzscore"])
UbKeGGLysinesNames<-c("K6", "K11", "K27", "K29", "K33", "K48", "K63")
pdf(file.path(".", analysisdir, "UbiquitinKeGGLysines_PO972kgg0hrPO1063kggHLM_IntensityHzscore_barplot.pdf"), height=5, width=7)
par(mar=c(5,6,4,2)+0.1)
ymin<-0.0
ymax<-ceiling(max(PO1063kggHLMPO972kgg0hr_UbKeGGLysines_AvgZscore, na.rm=TRUE))
barplot(PO1063kggHLMPO972kgg0hr_UbKeGGLysines_AvgZscore, beside=T, space=0.5, ylab="", yaxt='n', names=UbKeGGLysinesNames, cex.names=1.5, las=1, cex.axis=1.5, width=5, ylim=c(ymin,ymax), col="grey50")
box(bty="l", lwd=3)
axis(2,at=c(seq(ymin,ymax,by=0.5)),labels=F,col="black",cex.axis=1,lwd=3,tck=-0.01)
axis(2,at=c(seq(ymin,ymax,by=1.0)),labels=T,col="black",cex.axis=1.5,lwd=3,tck=-0.01,las=1)
#mtext("relative abundance", side=2, line=2.75, cex=1.5)
dev.off()
##end plot MSMS counts at 0 and 4 hours for TCR activation markers



PO1063kggHLMPO972kgg0hr_UbKeGGLysines_AvgHLratio<-c(KeGGdata["P62983_6","HLratioLog2Norm.PO0163kggHLMPO972kgg4diff0_AvgHLratio"],
                                                    KeGGdata["P62983_11","HLratioLog2Norm.PO0163kggHLMPO972kgg4diff0_AvgHLratio"],
												    KeGGdata["P62983_27","HLratioLog2Norm.PO0163kggHLMPO972kgg4diff0_AvgHLratio"],
												    KeGGdata["P62983_29","HLratioLog2Norm.PO0163kggHLMPO972kgg4diff0_AvgHLratio"],
												    KeGGdata["P62983_33","HLratioLog2Norm.PO0163kggHLMPO972kgg4diff0_AvgHLratio"],
												    KeGGdata["P62983_48","HLratioLog2Norm.PO0163kggHLMPO972kgg4diff0_AvgHLratio"],
												    KeGGdata["P62983_63","HLratioLog2Norm.PO0163kggHLMPO972kgg4diff0_AvgHLratio"])
UbKeGGLysinesNames<-c("K6", "K11", "K27", "K29", "K33", "K48", "K63")
pdf(file.path(".", analysisdir, "UbiquitinKeGGLysines_PO972kgg0hrPO1063kggHLM_HLratio_barplot.pdf"), height=5, width=7)
par(mar=c(5,6,4,2)+0.1)
#ymin<-floor(min(PO1063kggHLMPO972kgg0hr_UbKeGGLysines_AvgHLratio, na.rm=TRUE))
#ymax<-ceiling(max(PO1063kggHLMPO972kgg0hr_UbKeGGLysines_AvgHLratio, na.rm=TRUE))
ymin<- -0.2
ymax<- 0.5
barplot(PO1063kggHLMPO972kgg0hr_UbKeGGLysines_AvgHLratio, beside=T, space=0.5, ylab="", yaxt='n', names=UbKeGGLysinesNames, cex.names=1.5, las=1, cex.axis=1.5, width=5, ylim=c(ymin,ymax), col=c("grey50", "grey50", "grey50", "salmon", "salmon", "grey50", "salmon"))
abline(h=0, lwd=3)
axis(2,at=c(seq(ymin,ymax,by=0.05)),labels=F,col="black",cex.axis=1,lwd=3,tck=-0.01)
axis(2,at=c(seq(ymin,ymax,by=0.1)),labels=T,col="black",cex.axis=1.5,lwd=3,tck=-0.01,las=1)
#mtext("log2 fold change",2,line=2.75,cex=1.5)
dev.off()

stop("Z")

#####quantify protein-based ubiquitylation

#function to weight means by either intensity or count
weight <- function(x,y)
{

	if ( is.na(x) | is.na(y) ){
		wtval<-NA
	}else{
		wtval<-x*y
	}

	return(wtval)

}
#end weight function

#function to weight means by either intensity or count, normalized to total for that protein
normweight <- function(x,y,z)
{

	if ( is.na(x) | is.na(y) ){
		wtval<-NA
	}else{
		wtval<-x*(y/z)
	}

	return(wtval)

}
#end weight function

##unique protein ids in dataset
proteinIDs<-unique(KeGGdata$protID, incompariables=FALSE)
proteinIDsNum<-length(proteinIDs)
#print(proteinIDsNum)
#head(proteinIDs)


##combine peptide quantification in to protein-based quantification
#initialize data frame to store protein-level information
KeGGProteinData<-data.frame()
#each experiment in dataset
experiments<-c("PO1063kgg_HLM", "PO1063kgg_HLP", "PO972kgg_0hr", "PO972kgg_1hr", "PO972kgg_4hr")
#parse each experiment
for ( i in 1:length(experiments) )
{

	#current expriment
	expID<-experiments[i]
	#data column names
	LHcountid<-paste("HLcount.", expID, sep="")
	normLHid<-paste("HLratioLog2Norm.", expID, sep="")
	intensitylog2normid<-paste("IntensityLog2Norm.", expID, sep="")
	intensityLlog2normid<-paste("IntensityLLog2Norm.", expID, sep="")
	intensityHlog2normid<-paste("IntensityHLog2Norm.", expID, sep="")
	intensitylog2id<-paste("IntensityLog2.", expID, sep="")
	intensitywtnormLHid<-paste("intensityweightednormHLratio.", expID, sep="")



	#evaluate all proteins identified in the combined dataset
	#parse each unique protein id
	for ( i in 1:length(proteinIDs) )
	{
	
		#current protein
		tempprot<-proteinIDs[i]
		#print(tempprot)

		#current gene
		tempgene<-ProtIDdata[tempprot,"GeneId"]
		#print(tempgene)
		
		#current protein id
		tempid<-ProtIDdata[tempprot,"ProtName"]
		#print(tempid)

		KeGGProteinData[tempprot,"GeneId"]<-tempgene
		KeGGProteinData[tempprot,"ProtName"]<-tempid

		#peptide data for current protein
		tempdata<-subset(KeGGdata, KeGGdata$protID==tempprot)
		#print(tempdata)
		
		#num of identified sites for current protein
		numsites<-sum(!is.na(tempdata[,LHcountid]))
		#print(numsites)

		#populate protein data with number of sites for current protein
		colid<-paste("numsites.", expID, sep="")
		KeGGProteinData[tempprot,colid]<-numsites
		#print(numsites)

		#max intensity peptide value
		maxinten<-ifelse(!all(is.na(tempdata[,intensitylog2normid])), max(tempdata[,intensitylog2normid], na.rm=TRUE), NA)
		#populate protein data frame
		colid<-paste("maxinten.", expID, sep="")
		KeGGProteinData[tempprot,colid]<-maxinten
		#print(maxinten)

		#avg intensityL peptide value
		avgintenL<-ifelse(!all(is.na(tempdata[,intensityLlog2normid])), mean(tempdata[,intensityLlog2normid], na.rm=TRUE), NA)
		#populate protein data frame
		colid<-paste("avgintenL.", expID, sep="")
		KeGGProteinData[tempprot,colid]<-avgintenL
		#print(avgintenL)

		#avg intensityH peptide value
		avgintenH<-ifelse(!all(is.na(tempdata[,intensityHlog2normid])), mean(tempdata[,intensityHlog2normid], na.rm=TRUE), NA)
		#populate protein data frame
		colid<-paste("avgintenH.", expID, sep="")
		KeGGProteinData[tempprot,colid]<-avgintenH
		#print(avgintenH)

		#skip "0" sites, these proteins were identified in a different study
		#the number of sites are counted by those that have an HLcount
		if ( numsites>0 )
		{
		
			###populate protein data with number of sites for current protein
			##colid<-paste("numsites.", expID, sep="")  #################could move this to above the loop i think
			##KeGGProteinData[tempprot,colid]<-numsites
			###print(numsites)
			
			###max intensity peptide value
			##maxinten<-max(tempdata[,intensitylog2normid], na.rm=TRUE)
			###populate protein data frame
			##colid<-paste("maxinten.", expID, sep="") #################could move this to above the loop i think
			##KeGGProteinData[tempprot,colid]<-maxinten
			###print(maxinten)
			
			#max intensity peptide #################this is incorrect
			##########maxintenpep<-row.names(tempdata[which(tempdata[,intensitylog2normid] == max(tempdata[,intensitylog2normid], na.rm = TRUE)), ])[1]
			maxintenpep<-row.names(tempdata[which(tempdata[,intensitylog2normid]==maxinten), ])[1]
			#populate protein data frame
			colid<-paste("maxintenpep.", expID, sep="")
			KeGGProteinData[tempprot,colid]<-maxintenpep
			#print(maxintenpep)

			#max intensity peptide normLHratio
			maxintenratio<-tempdata[maxintenpep,normLHid]
			#populate protein data frame
			colid<-paste("maxintenpepratio.", expID, sep="")
			KeGGProteinData[tempprot,colid]<-maxintenratio

			#intensity weighted mean
			#calc total intensity
			totintensity<-sum(tempdata[,intensitylog2id], na.rm=TRUE)
			#weight ratios by by raw intensity
			intensityweighted<-apply(tempdata, MARGIN=1, FUN=function(x) weight(as.numeric(x[[normLHid]]),as.numeric(x[[intensitylog2id]])))
			#populate temp data frame with weighted intensity
			tempdata[,intensitywtnormLHid]<-intensityweighted
			#calc intensity weighted mean
			intensityweightedmean<-sum(tempdata[,intensitywtnormLHid], na.rm=TRUE)/totintensity

			###########weight by intensity norm to total intensity
			##########totintensity<-sum(tempdata[,intensitylog2id], na.rm=TRUE)
			##########intensityweighted<-apply(tempdata, MARGIN=1, FUN=function(x) normweight(as.numeric(x[[normLHid]]),as.numeric(x[[intensitylog2id]]),as.numeric(totintensity)))
			###########populate temp data frame with weighted intensity
			##########tempdata[,intensitywtnormLHid]<-intensityweighted
			###########calc weighted mean
			##########intensityweightedmean<-mean(tempdata[,intensitywtnormLHid], na.rm=TRUE)
			
			#populate protein data frame
			colid<-paste("intensityweightedmeanratio.", expID, sep="")
			KeGGProteinData[tempprot,colid]<-intensityweightedmean
			#print(intensityweightedmean)

		}else{
		
			#colid<-paste("numsites.", expID, sep="")
			#KeGGProteinData[tempprot,colid]<-numsites
			#colid<-paste("maxinten.", expID, sep="")
			#KeGGProteinData[tempprot,colid]<-NA
			colid<-paste("maxintenpep.", expID, sep="")
			KeGGProteinData[tempprot,colid]<-NA
			colid<-paste("maxintenpepratio.", expID, sep="")
			KeGGProteinData[tempprot,colid]<-NA
			colid<-paste("intensityweightedmeanratio.", expID, sep="")
			KeGGProteinData[tempprot,colid]<-NA

		} #end numsites>0 loop

	} #end protein site loop
	
} #end experiment loop

#dim(KeGGProteinData)
#head(KeGGProteinData)


#intensity weighted
#0hr-4hr
###intwtLHmeancomp4diff0<-apply(KeGGProteinData, MARGIN=1, FUN=function(x) diffcomp(x[["intensityweightedmeanratio.PO972kgg_4hr"]],x[["intensityweightedmeanratio.PO972kgg_0hr"]]))
###KeGGProteinData[,"intensityweightedmeanratio.PO972kgg_4diff0"]=intwtLHmeancomp4diff0
KeGGProteinData$intensityweightedmeanratio.PO972kgg_4diff0<-apply(KeGGProteinData, MARGIN=1, FUN=function(x) diffcomp(x[["intensityweightedmeanratio.PO972kgg_4hr"]],x[["intensityweightedmeanratio.PO972kgg_0hr"]]))
#0hr-1hr
###intwtLHmeancomp1diff0<-apply(KeGGProteinData, MARGIN=1, FUN=function(x) diffcomp(x[["intensityweightedmeanratio.PO972kgg_1hr"]],x[["intensityweightedmeanratio.PO972kgg_0hr"]]))
###KeGGProteinData[,"intensityweightedmeanratio.PO972kgg_1diff0"]=intwtLHmeancomp1diff0
KeGGProteinData$intensityweightedmeanratio.PO972kgg_1diff0<-apply(KeGGProteinData, MARGIN=1, FUN=function(x) diffcomp(x[["intensityweightedmeanratio.PO972kgg_1hr"]],x[["intensityweightedmeanratio.PO972kgg_0hr"]]))
#1hr-4hr
###intwtLHmeancomp4diff1<-apply(KeGGProteinData, MARGIN=1, FUN=function(x) diffcomp(x[["intensityweightedmeanratio.PO972kgg_4hr"]],x[["intensityweightedmeanratio.PO972kgg_1hr"]]))
###KeGGProteinData[,"intensityweightedmeanratio.PO972kgg_4diff1"]=intwtLHmeancomp4diff1
KeGGProteinData$intensityweightedmeanratio.PO972kgg_4diff1<-apply(KeGGProteinData, MARGIN=1, FUN=function(x) diffcomp(x[["intensityweightedmeanratio.PO972kgg_4hr"]],x[["intensityweightedmeanratio.PO972kgg_1hr"]]))


##write table of protein Ub quantification data
write.table(KeGGProteinData, file=file.path(".", analysisdir, "KeGGproteinData_PeptideBasedTotalProteinQuantification_PO1063-PO972_0810_NewWeighting_1203.txt"), sep="\t", quote=FALSE, na="NA", dec=".", row.names=TRUE, col.names=NA)


print("PO972kgg vs PO1063HLMkgg correlations for 0-4hr fold change for ubiquitinated peptides")
ubpepcortest_UbFc_pearson<-cor.test(KeGGdata$HLratioLog2Norm.PO972kgg_4diff0, KeGGdata$HLratioLog2Norm.PO1063kgg_HLM, alternative="two.sided", method="pearson", conf.level = 0.95)
print(ubpepcortest_UbFc_pearson)
ubpepcortest_UbFc_spearman<-cor.test(KeGGdata$HLratioLog2Norm.PO972kgg_4diff0, KeGGdata$HLratioLog2Norm.PO1063kgg_HLM, alternative="two.sided", method="spearman", conf.level = 0.95)
print(ubpepcortest_UbFc_spearman)
print("PO972kgg vs PO1063HLMkgg correlations for 0-4hr fold change for ubiquitinated proteins")
ubprotcortest_UbFc_pearson<-cor.test(KeGGProteinData$intensityweightedmeanratio.PO972kgg_4diff0, KeGGProteinData$intensityweightedmeanratio.PO1063kgg_HLM, alternative="two.sided", method="pearson", conf.level = 0.95)
print(ubprotcortest_UbFc_pearson)
ubprotcortest_UbFc_spearman<-cor.test(KeGGProteinData$intensityweightedmeanratio.PO972kgg_4diff0, KeGGProteinData$intensityweightedmeanratio.PO1063kgg_HLM, alternative="two.sided", method="spearman", conf.level= 0.95)
print(ubprotcortest_UbFc_spearman)

##density scatter plot of comparison of PO1063 vs PO972
#peptide-based Ub fold change
makedensitycoloredscatterplot(KeGGdata$HLratioLog2Norm.PO972kgg_4diff0, KeGGdata$HLratioLog2Norm.PO1063kgg_HLM, "NA", "NA", "NA", "NA", "peptide Ub (exp 1) log2 fold change", "peptide Ub (exp 2) log2 fold change", file.path(".", analysisdir, "PO972PO1063_0hr4hrStim_KeGGPeptideQuantFoldChangeComp_DensityPlot.pdf"))
#protein-based Ub fold change
makedensitycoloredscatterplot(KeGGProteinData$intensityweightedmeanratio.PO972kgg_4diff0, KeGGProteinData$intensityweightedmeanratio.PO1063kgg_HLM, "NA", "NA", "NA", "NA", "protein Ub (exp 1) log2 fold change", "protein Ub (exp 2) log2 fold change", file.path(".", analysisdir, "PO972PO1063_0hr4hrStim_KeGGProteinQuantFoldChangeComp_DensityPlot.pdf"))
##end density scatter plot of comparison of PO1063 vs PO972 peptide-based Ub fold change



if ( FALSE )
{

## the following are not used.... should delete?

##plot TCR pathway data
TCRpath<-c("Lat","Zap70","Prkcq","Cd3e","Rhoa","Grap", "Pi4k2a", "Pip5k1a")
xmin <- -3
xmax <- 3
ymin <- -4
ymax <- 3
TCRpathwayNeddfig<-file.path(".", analysisdir, "ProteinKeGGUbiquitylationVsUbWithNeddInhibitor.pdf")
pdf(TCRpathwayNeddfig)
par(bty='l')
plot(KeGGProteinData$intensityweightedmeanratio.PO1063kgg_HLM,
     KeGGProteinData$intensityweightedmeanratio.PO1063kgg_HLP,
	 pch=20,
	 cex=0.75,
	 col="grey50",
	 xlim=c(xmin,xmax), ylim=c(ymin,ymax),
	 xlab="", ylab="",
	 xaxt='n', yaxt='n'
)
par(new=T)
abline(a = 0, b = 1, col="grey", lwd=2, lty=2)
par(new=T)
abline(a = 1, b = 1, col="grey", lwd=2, lty=3)
par(new=T)
abline(a = -1, b = 1, col="grey", lwd=2, lty=3)
for ( i in 1:length(TCRpath) ){
	currgene<-TCRpath[i]
	points(KeGGProteinData[which(KeGGProteinData$GeneId==currgene), "intensityweightedmeanratio.PO1063kgg_HLM"],
	       KeGGProteinData[which(KeGGProteinData$GeneId==currgene), "intensityweightedmeanratio.PO1063kgg_HLP"],,
	       pch=20,
	       cex=1.75,
	       col="goldenrod2"
	)
}
axis(1,c(seq(xmin,xmax,by=0.5)), labels=F, col="black",cex.axis=1.0, tck=-0.01)
axis(1,c(seq(xmin,xmax,by=1)), labels=T, col="black",cex.axis=1.0)
mtext("Ub 0-4hr TCR stim log2 fold change (control)",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)), labels=F, col="black",cex.axis=1.0, tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)), labels=T, col="black",cex.axis=1.0)
mtext("Ub 0-4hr TCR stim log2 fold change (+ Nedd inhibitor)",2,line=2.75,cex=1.5)
dev.off()
}

if ( FALSE )
{
##plot Cullins data
cullins<-c("Cul1", "Cul3", "Cul4a", "Cul4b", "Cul5")
neddinhibdelta<-rownames(subset(KeGGProteinData, KeGGProteinData$intensityweightedmeanratio.PO1063kgg_HLP <= -2))
#print(neddinhibdelta)
#abline(fit <- lm(KeGGProteinData$intensityweightedmeanratio.PO1063kgg_HLP ~ KeGGProteinData$intensityweightedmeanratio.PO1063kgg_HLM), col="grey70", lwd=2) # regression line (y~x)
#par(new=T)
#legend("topright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)), text.col="grey30")
for ( i in 1:length(neddinhibdelta) ){
  currgene<-neddinhibdelta[i]
  points(KeGGProteinData[currgene, "intensityweightedmeanratio.PO1063kgg_HLM"],
       KeGGProteinData[currgene, "intensityweightedmeanratio.PO1063kgg_HLP"],
	   pch=20,
	   cex=1,
	   col="dodgerblue"
  )
  text(KeGGProteinData[currgene, "intensityweightedmeanratio.PO1063kgg_HLM"],
     KeGGProteinData[currgene, "intensityweightedmeanratio.PO1063kgg_HLP"], 
	 labels=(KeGGProteinData[currgene, "GeneId"]),
	 cex=0.75, pos=2, offset=0.2, col="dodgerblue", font=2
  )
}
par(new=T)
for ( i in 1:length(cullins) ){
  currcul<-cullins[i]
  points(KeGGProteinData[rownames(KeGGProteinData[which(KeGGProteinData$GeneId == currcul),]), "intensityweightedmeanratio.PO1063kgg_HLM"],
       KeGGProteinData[rownames(KeGGProteinData[which(KeGGProteinData$GeneId == currcul),]), "intensityweightedmeanratio.PO1063kgg_HLP"],
	   pch=20,
	   cex=1,
	   col="red"
  )
  text(KeGGProteinData[rownames(KeGGProteinData[which(KeGGProteinData$GeneId == currcul),]), "intensityweightedmeanratio.PO1063kgg_HLM"],
     KeGGProteinData[rownames(KeGGProteinData[which(KeGGProteinData$GeneId == currcul),]), "intensityweightedmeanratio.PO1063kgg_HLP"], 
	 labels=(KeGGProteinData[rownames(KeGGProteinData[which(KeGGProteinData$GeneId == currcul),]), "GeneId"]),
	 cex=0.75, pos=2, offset=0.2, col="red", font=2
  )
}
  points(KeGGProteinData[rownames(KeGGProteinData[which(KeGGProteinData$GeneId == "Cish"),]), "intensityweightedmeanratio.PO1063kgg_HLM"],
       KeGGProteinData[rownames(KeGGProteinData[which(KeGGProteinData$GeneId == "Cish"),]), "intensityweightedmeanratio.PO1063kgg_HLP"],
	   pch=20,
	   cex=1,
	   col="red"
  )
  text(KeGGProteinData[rownames(KeGGProteinData[which(KeGGProteinData$GeneId == "Cish"),]), "intensityweightedmeanratio.PO1063kgg_HLM"],
     KeGGProteinData[rownames(KeGGProteinData[which(KeGGProteinData$GeneId == "Cish"),]), "intensityweightedmeanratio.PO1063kgg_HLP"], 
	 labels=(KeGGProteinData[rownames(KeGGProteinData[which(KeGGProteinData$GeneId == "Cish"),]), "GeneId"]),
	 cex=0.75, pos=2, offset=0.2, col="red", font=2
  )


  points(KeGGProteinData[rownames(KeGGProteinData[which(KeGGProteinData$GeneId == "Slfn2"),]), "intensityweightedmeanratio.PO1063kgg_HLM"],
       KeGGProteinData[rownames(KeGGProteinData[which(KeGGProteinData$GeneId == "Slfn2"),]), "intensityweightedmeanratio.PO1063kgg_HLP"],
	   pch=20,
	   cex=1,
	   col="red"
  )
  text(KeGGProteinData[rownames(KeGGProteinData[which(KeGGProteinData$GeneId == "Slfn2"),]), "intensityweightedmeanratio.PO1063kgg_HLM"],
     KeGGProteinData[rownames(KeGGProteinData[which(KeGGProteinData$GeneId == "Slfn2"),]), "intensityweightedmeanratio.PO1063kgg_HLP"], 
	 labels=(KeGGProteinData[rownames(KeGGProteinData[which(KeGGProteinData$GeneId == "Slfn2"),]), "GeneId"]),
	 cex=0.75, pos=2, offset=0.2, col="red", font=2
  )
axis(1,c(seq(xmin,xmax,by=0.5)), labels=F, col="black",cex.axis=1.0, tck=-0.01)
axis(1,c(seq(xmin,xmax,by=1)), labels=T, col="black",cex.axis=1.0)
mtext("Ub 0-4hr TCR stim log2 fold change (control)",1,line=2.75,cex=1.5)
axis(2,c(seq(ymin,ymax,by=0.5)), labels=F, col="black",cex.axis=1.0, tck=-0.01)
axis(2,c(seq(ymin,ymax,by=1)), labels=T, col="black",cex.axis=1.0)
mtext("Ub 0-4hr TCR stim log2 fold change (+ Nedd inhibitor)",2,line=2.75,cex=1.5)
}
