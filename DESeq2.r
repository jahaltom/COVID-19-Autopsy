library(DESeq2)
library(dplyr)

#Read in count information. 
countData = read.table("SRP.tsv",header=TRUE,sep = '\t',quote="")

#X an Y gene names can be the same. This makes them unique. Row is set to this unique value. 
rownames(countData)= paste(countData$Gene_name,countData$Gene_stable_ID, countData$chr,sep="$")
#Extract count columns
countData=select(countData,contains("SRR"))
countData=round(countData,0) 


##Read in expermental design
metadata = read.table("AutopsyStudiesCurated.tsv",header=TRUE,row.names=1,sep = '\t')
metadata=subset(metadata,study_accession == "SRP")


#Should return TRUE
#all(rownames(metadata) == colnames(countData))

##Make DEseq2 object 
dds = DESeqDataSetFromMatrix(countData = countData,colData = metadata,design = ~ disease)   
dds = DESeq(dds)
#Contrast case vs control
result = results(dds, contrast=c("disease","COVID-19","Non-Covid-19"))
## Remove rows with NA
result = result[complete.cases(result),]
#Put GeneID as column 
result = cbind(ID = rownames(result), result)


write.table(result,"SRP_Covid19VSNonCovid19_DGE.tsv" ,sep = '\t',row.names = FALSE)
