rm(list=ls())
setwd("/Users/malcolm/Dropbox/Research/CFseqAnalysis")
dataT <- read.table("annotatedOTUsAsColumnsJoshUnLoggedWithSampleDays.txt", header=TRUE, sep="\t")
conT = dataT[,23:172]  # Hack/Magic number defining the columns of the excel sheet related to the OTU sequence counts
otuT  <- read.table("Cons.txt", header=FALSE, sep="\t")
conSum = colSums(conT,na.rm=TRUE)
percentTotalSeqs=conSum/sum(conSum)*100

treatments = factor( dataT$treatmentString ) 
treatments = relevel( treatments, ref ="Stable")

#----use column V5 (genus)----
allTaxa= otuT$V5
rawTotal  = rowSums(conT)
rawMean   = rowMeans(conT)
rawAll = array(0,dim=c(130,length(allTaxa)))
for(col in 1:length(allTaxa)){
  rawAll[,col] = conT[,col] 
}
colnames(rawAll) <- allTaxa


logTaxa = matrix(data=NA,dim(rawAll)[1],dim(rawAll)[2]) #rep(c(1),ncol(rawAll))
pVals.Lmoriginal = rep(c(1),ncol(rawAll))
pVals.Lslope = rep(c(1),ncol(rawAll))
pVals.Lintercept = rep(c(1),ncol(rawAll))
pVals.R2val      = rep(c(1),ncol(rawAll))
pVals.1original = rep(c(1),ncol(rawAll))
for( col in 1:ncol(rawAll)){
  logTaxa[,col] = log10(rawAll[,col]/rawTotal*rawMean + 1)
  lm1 = lm( logTaxa[,col] ~ dataT$sampleDays + treatments, x=TRUE, y=TRUE)
  pVals.Lslope[col]      = lm1$coefficients[2]
  pVals.Lintercept[col]  = lm1$coefficients[1]
  pVals.R2val[col]       = summary(lm1)$r.squared
  myAnova = anova( lm1 )
  pVals.Lmoriginal[col]  = myAnova$"Pr(>F)"[1]
  #myLm    = lm( logTaxa[,col] ~ dataT$sampleDays * treatments, x=TRUE, y=TRUE)
  #pVals.1original[col] = myLm
}
print("P-10% FDR, Original: ")
print(sum(pVals.Lmoriginal<0.10))
pVals.lmadj = p.adjust(pVals.Lmoriginal, "BH")
print("P-10% FDR, BH Corrected: ")
print(sum(pVals.lmadj<0.10))

###########################################################################
hist(pVals.Lmoriginal)
hist(pVals.lmadj,breaks=30)

dummy150 = 1:150
sigIdxs =  (dummy150<99) #  hack/Magic number, OTU's above 100 had a very low number of sequences
                         # & pVals.lmadj<0.10
sigTaxa = data.frame(name = as.character(allTaxa[sigIdxs]), p_value_adj = pVals.lmadj[sigIdxs], p_value_orig=pVals.Lmoriginal[sigIdxs], slope = pVals.Lslope[sigIdxs], intercept = pVals.Lintercept[sigIdxs], R2 = pVals.R2val[sigIdxs], OTUnum=dummy150[sigIdxs], percentTotalSeqs=percentTotalSeqs[sigIdxs])
sigSorted = sigTaxa[with(sigTaxa, order(p_value_adj)), ]
print(sigSorted)
write.csv(sigSorted, "significantTaxa10percentFDR.csv", row.names=FALSE)


sigValues = logTaxa[,sigIdxs]

for(plotNum in 1:9){
  plot(dataT$sampleDays,logTaxa[,sigSorted$OTUnum[plotNum]])#sigValues[,plotNum])
  y1=sigSorted$intercept[plotNum] + sigSorted$slope[plotNum]*dataT$sampleDays
  lines(dataT$sampleDays,y1,col="RED")
}


#-------------------------------RICHNESS PLOTS--------------------------------------------------------
# dataT$richness76622 takes into account normalization based on total number of sequences
subset = !is.na(rawAll[,1]) # use this to create a subset excluding NA values 
myLmR = lm( dataT$richness76622 ~ dataT$sampleDays + treatments, x=TRUE, y=TRUE )
coefs = myLmR$coefficients
plot( dataT$sampleDays, dataT$richness76622)
y1=coefs[1]+coefs[2]*myLmR$x[,2]
y2=coefs[1]+coefs[2]*myLmR$x[,2] + coefs[3]
y3=coefs[1]+coefs[2]*myLmR$x[,2] + coefs[4]
y4=coefs[1]+coefs[2]*myLmR$x[,2] + coefs[5]
lines(dataT$sampleDays[subset],y1,col="RED")
lines(dataT$sampleDays[subset],y2,col="MAGENTA")
lines(dataT$sampleDays[subset],y3,col="BLUE")
lines(dataT$sampleDays[subset],y4,col="GREEN")








