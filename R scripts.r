#### Integration of DNA Methylation and Gene Expression in Lung Adenocarcinoma ####

####################################################################################
################## emQTL Analysis #########################################
####################################################################################
library(msigdbr) # use msigdbr() to annotate transcripts as gene symbols
library(dplyr) # subset(), group_by()
library(TCGAutils) #package required to identify origin of biopsy
library(matrixStats) # used function rowIQRs

#Import TCGA data
# RNA-Seq TCGA-LUAD
seq_raw <- readRDS("/open/archive/TCGA/DATA/LUAD/RNA/mRNA/TCGA-LUAD.htseq_fpkm-uq.rds")
# DNA methylation TCGA-LUAD
meth_raw <- readRDS("/open/archive/TCGA/DATA/LUAD/DNA/Methylation/TCGA-LUAD.methylation450_processed_imp.rds")

################# RNA-Seq preprocessing ##############
#Remove decimals in Ensembl IDs to harmonize transcripts names
seq_raw$Ensembl_ID <- sapply(strsplit(seq_raw$Ensembl_ID,".",fixed = T), getElement, 1)

# Transcript annotation with MsigDB 
#Extracting Ensembl IDs based on msigDB
genesets <- msigdbr(species = "Homo sapiens")
#Convert ensembl ID to gene symbol
gene_convert <- subset(genesets, select=c("human_ensembl_gene","human_gene_symbol"))
#Grouping by gene symbol. Filtering out duplicates.
gene_convert <- gene_convert %>% group_by(human_ensembl_gene,human_gene_symbol) %>% filter(row_number() == 1)
seq <- merge(seq_raw,gene_convert, by.x = "Ensembl_ID", by.y = "human_ensembl_gene") # Merging RNA-Seq data and conversion table
seq$Ensembl_ID = NULL #Removing Ensembl_ID column
seq <- seq %>% group_by(human_gene_symbol) %>% dplyr::summarise_all(median) # Grouping by median
rownames(seq) <- seq$human_gene_symbol # Set symbols as rownames
seq <- data.matrix(seq, rownames.force = TRUE)
seq <- seq[,-1] # Remove symbol column
seq <- seq[rowSums(seq) != 0,] #remove genes with zero expression level

#Matching column names with columns names in methylation array data
colnames(seq) <- gsub(".","-",colnames(seq),fixed = TRUE)
my_samples <- colnames(seq) #Extract vector with sample names
my_samples_01 <- my_samples[TCGAsampleSelect(my_samples,01)] #Extract primary samples
primary_seq <- seq[,colnames(seq) %in% my_samples_01]
primary_seq <- primary_seq[rowSums(primary_seq) != 0,] #Filter out zero expression level

# Match RNA-Seq and DNA methylation array data across tumors 
seq_final <- primary_seq[,colnames(primary_seq) %in% colnames(meth_raw)]
## DNA methylation in TCGA-LUAD preprocessing
meth <- meth_raw[,colnames(meth_raw) %in% colnames(seq_final)] # extract tumors common for RNA-Seq and DNA mthhylation
gene_indx <- match(colnames(seq_final),colnames(meth))  #Match the colnames in meth with colnames in seq#############
														#The second input will change while the first one stays the same.
meth_ordered <- meth[,gene_indx]
meth_samples <- colnames(meth_raw) #Extract normal tissue biopsy Meth

#### BEFORE ANY FURTHER ANALYSIS, CHECK IF ALL SAMPLES ARE UNIQUE and comes from unique patients###
seq <- seq_final
colnames(seq) <- gsub("(-01).*","\\1",colnames(seq))
seq_df <- as.data.frame(seq)
seq_no_db <- sapply(unique(colnames(seq_df)), function(x) rowMeans(seq_df[grepl(x,colnames(seq_df))]))

#############################################################################################
#################### Split the data into Training and Validation cohorts/subsets ##################

set.seed(42)

sample_name <- colnames(seq_no_db)
discovery <- sample(sample_name,length(sample_name)/2)
seq_discovery <- seq_no_db[,colnames(seq_no_db) %in% discovery]
seq_no_zero <- seq_discovery[rowSums(seq_discovery) != 0,] #Remove rows with zero expression
seq_discovery <- seq_no_zero

############ Working with DISCOVERY subset###################
####################################################
meth <- meth_ordered
colnames(meth) <- gsub("(-01).*","\\1",colnames(meth))
meth_df <- as.data.frame(meth)
meth_no_db <- sapply(unique(colnames(meth_df)), function(x) rowMeans(meth_df[grepl(x,colnames(meth_df))]))
meth_discovery <- meth_no_db[,colnames(meth_no_db) %in% colnames(seq_discovery)]

#Harmonizing sample ID in RAN-seq and DNA methylation microarray.
samp_indx <- match(colnames(seq_discovery),colnames(meth_discovery)) 
meth_ordered <- meth_discovery[,samp_indx]
# Check the order
identical(colnames(seq_discovery),colnames(meth_ordered))

# Extract rows with IQR > 0.1
meth_discovery_final <- meth_ordered[rowIQRs(meth_ordered)>0.1,]

############Filtering Gender_related CpGs#######################
CpG_info <- read.table("/data2/thomas/breast450k/NewProbeinfo.txt", header = TRUE, sep ="\t")
only_cpg <- as.data.frame(rownames(meth_discovery_final))
colnames(only_cpg) = "CpG"
complete <- merge(CpG_info,only_cpg, by.x = "Probe", by.y = "CpG" )
Y <- subset(complete, Chr == "Y")
X <- subset(complete, Chr == "X")

meth_discovery_final_filt <- meth_discovery_final[setdiff(rownames(meth_discovery_final), X$Probe),]
meth_discovery_final_filt_final <- meth_discovery_final_filt[setdiff(rownames(meth_discovery_final_filt), Y$Probe),]

###################################################################
########## emQTL analysis in DISCOVERY subset ##########################
r <- cor(t(seq_discovery),t(meth_discovery_final_filt_final))
saveRDS(r, paste0(my_directory,"/corr_r_SexFiltered_discovery_",dim(r)[1],"x",dim(r)[2],".rds"))
n <- ncol(seq_discovery)
df <- n-2
t <- (r*sqrt(df))/(sqrt(1-r^2))
p <- 2*pt(abs(t),df,lower.tail=FALSE)# pt gives the distribution function
saveRDS(p,paste0(my_directory,"/p_value_SexFiltered_discovery.rds"))
pcut_discovery <- 0.05/(as.numeric(nrow(p))*as.numeric(ncol(p))) #8.612673e-12
rowsums <- rowSums(p < pcut_discovery) > 0  
colsums <- colSums(p < pcut_discovery) > 0 
p_filt <- p[rowsums,colsums]    
num_signif <- sum(p < pcut_discovery) 			
r_filt <- r[rowsums,colsums]
saveRDS(r_filt,paste0(my_directory,"/open/tmp/Anr_filt_SexFiltered_discovery_",dim(r_filt)[1],"x",dim(r_filt)[2],".rds"))
 
#########################################################################################
#######################Working with VALIDATION cohort/subset############################
#VALIDATION RNA-Seq ####################################################################
seq_validation <- seq_no_db[,!colnames(seq_no_db) %in% discovery]
seq_val_no_zero <- seq_validation[rowSums(seq_validation) != 0,]
seq_validation <- seq_val_no_zero
meth_validation <- meth_no_db[,colnames(meth_no_db) %in% colnames(seq_validation)]

#Sorting data in meth
gene_indx <- match(colnames(seq_validation),colnames(meth_validation)) 
meth_valid_ordered <- meth_validation[,gene_indx]

#Check if the colnames are ordered with respect to samples
identical(colnames(seq_validation),colnames(meth_valid_ordered))

##############Filter out CpGs related to gender##########
only_cpg <- as.data.frame(rownames(meth_valid_ordered))
colnames(only_cpg) = "CpG"
complete <- merge(CpG_info,only_cpg, by.x = "Probe", by.y = "CpG" )
Y <- subset(complete, Chr == "Y")
X <- subset(complete, Chr == "X")

meth_validation_final_filt <- meth_valid_ordered[setdiff(rownames(meth_valid_ordered), X$Probe),]
meth_validation_final_filt_final <- meth_validation_final_filt[setdiff(rownames(meth_validation_final_filt), Y$Probe),]

#Check if the colnames are ordered with respect to samples
identical(colnames(seq_validation),colnames(meth_validation_final_filt_final))

# Filter data based on significant CpGs and Genes
meth_valid <- meth_validation_final_filt_final
meth_filt_validation <- meth_valid[rownames(meth_valid) %in% colnames(p_filt),]

saveRDS(meth_filt_validation,paste0(my_directory,"/meth_filtByPval_validation_SexFiltered_discovery_",dim(meth_filt_validation)[1],"x",dim(meth_filt_validation)[2],".rds"))
seq_filt_validation <- seq_validation[rownames(seq_validation) %in% rownames(p_filt),]
seq_no_zero_validation <- seq_filt_validation[rowSums(seq_filt_validation) != 0,]
saveRDS(seq_no_zero_validation,paste0(my_directory,"/seq_filtByPval_validation_SexFiltered_",dim(seq_no_zero_validation)[1],"x",dim(seq_no_zero_validation)[2],".rds"))

##############################################################################
################ CORRELATION IN VALIDATION DATA ##################################
########################################################################

r <- cor(t(seq_no_zero_validation),t(meth_filt_validation))
saveRDS(r,paste0(my_directory,"/r_complete_validation_SexFiltered_",dim(r)[1],"x",dim(r)[2],".rds"))
n <- ncol(seq_no_zero_validation)
df <- n-2
t <- (r*sqrt(df))/(sqrt(1-r^2))
p <- 2*pt(abs(t),df,lower.tail=FALSE)# pt gives the distribution function

saveRDS(p,paste0("/open/tmp/Anastasia/LUAD/LUAD_postpilot/LUAD_latest/raw_data_unique_samp/EMQTL/SexFiltered/p_complete_validation_SexFiltered_",dim(p)[1],"x",dim(p)[2],".rds"))

##################################################
##Filtering discovery data with respect to significant values both in discovery and validation set ####

# Convert pval matrix to boolean based on pcut
#import filtered p-val data
p_discovery <- p_filt
#import complete p-val data for pcut calculation
p_discovery_complete <- readRDS(my_directory,"/p_value_SexFiltered_discovery.rds")
#Import validation p-val data
p_validation <- p
# Adjust dimensions of matrices
p_discovery_matched <- p_discovery[rownames(p_discovery) %in% rownames(p_validation),]
# Check the order of col and rownames in p_discovery_matched and p_validation
identical(colnames(p_discovery_matched),colnames(p_validation))

# Convert pval matrix to boolean based on pcut
pcut_discovery <- 0.05/(as.numeric(nrow(p_discovery_complete))*as.numeric(ncol(p_discovery_complete)))
num_signif <- sum(p_discovery < pcut_discovery) #10 402 154
p_bool_discovery <- p_discovery_matched < pcut_discovery
pcut_validation <- 0.05/num_signif #  4.806697e-09
p_bool_validation <- p_validation < pcut_validation
# filtering based on signif values both in discovery and validation data.
bool_signif_p <- p_bool_discovery & p_bool_validation
 
#####Filtering based on number of associations 
rowsums <- rowSums(bool_signif_p == TRUE) > 5 
colsums <- colSums(bool_signif_p == TRUE) > 5 

p_final <- p_discovery[rowsums,colsums]    
r_discovery <- readRDS(my_directory,"/corr_r_SexFiltered_discovery_37431x155096.rds")
r_final <- r_discovery[rownames(r_discovery) %in% rownames(p_final),colnames(r_discovery) %in% colnames(p_final)]
r <- -r_final
r <- t(r)

saveRDS(r,paste0(my_directory,"/r_final_discovery_forClustering_SexFiltered_",dim(r)[1],"x",dim(r)[2],".rds"))
#########################################################################

#########################################################################################################################
############################ Gene Set Enrichment Analysis (GSEA) in R ##################################################################
##################################################################################################################
library(dplyr)
library(ggplot2)
library(EnsDb.Hsapiens.v79)

#create collection name convertor
	num_clusters = 5
#### Define N for the enrichment analysis #####
	msig <- readRDS(paste0(my_directory, "/MsigDB_anotated_seq_38048x585.rds")
	msig_genes <- rownames(msig)
	all_gensets <- readRDS(my_directory,"/MSigDB_GSEA_data.rds")
	all_genes <- unique(all_gensets$Gene) #
	GenesInUni <- length(intersect(all_genes,msig_genes))# number of overlapping my genes with all total unique genes in defined gene sets

### Create dictionary for gene collection ###
	name_col <- unique(all_gensets$GeneSetCollection)  #"H"  "C1" "C2" "C3" "C4" "C5" "C6" "C7" "C8"
	name_col <- name_col[-(2)]
	keys <- c("H: hallmark gene sets",
			#"C1: positional gene sets", # not adequate data
			"C2: curated gene sets",
			"C3: regulatory target gene sets",
			"C4: computational gene sets",
			"C5: ontology gene sets",
			"C6: oncogenic signature gene sets",
			"C7: immunologic signature gene sets",
			"C8: cell type signature gene sets")			

	dict <- list()
	for (i in 1:length(name_col)){
	dict[i] = keys[i]
	}	
	names(dict) <- name_col
	print(dict)
	
for (c in 1:length(name_col)){
geneset_collection <- name_col[c]
### Choose collection to use for analysis ###
	genesets <- all_gensets[all_gensets$GeneSetCollection %in% geneset_collection,]
	if(geneset_collection=="C8"){genesets <- genesets[grep("LUNG",genesets$GeneSetName),]} #from scRNA include only data for lung
### exclude HPO (Phenotypes)
	include_hpo = FALSE # CHANGE if needed
	if(include_hpo==FALSE){genesets <- genesets[!genesets$GeneSetName%in%unique(genesets[grep("HP_",genesets$GeneSetName),]$GeneSetName),]}
#######
	res <- data.frame(matrix(nrow=length(unique(genesets$GeneSetName)),ncol=9))
	colnames(res) <- c("GeneSetName","GenesInGeneset","GenesInOverlap","p_value","FDR_q_value","FE","RatioAndFE","CollectionName","Only_geneset_names")
	res$GeneSetName <- unique(genesets$GeneSetName)
	res$CollectionName <- rep(geneset_collection, nrow(res))
	Gesets_no_dup <- genesets[!duplicated(genesets$GeneSetName),]
	res$GenesInGeneset <- Gesets_no_dup$GenesInGeneset # to get the column with number of genes per gene set

######## My genes ###########
	num_clusters <- 5 #CHANGE if required
	read_biclusters <- readLines(paste0(my_directory,"/List_all_Genes_",num_clusters,"_biclusters.txt")) #read the txt file with biclusters
	m_bicluster <- as.matrix(read_biclusters) #convert list to a matrix
	sep_matrix <- strsplit(m_bicluster," ") # separate genes inside the list
	counter = 1
	name_clust <- c("Cell Cycle", "Mixed I", "Immune Ig", "Mixed II", "Immune") # Change	

	for (line in 1:length(sep_matrix)){ #

		if  (length(sep_matrix[[line]]) == 0){
			next
			}

	bicluster <- unlist(sep_matrix[line])
	out <- split(genesets,f=genesets$GeneSetName); out <- out[match(res$GeneSetName,names(out))]
	empty_vector <- vector()
	for(i in 1:nrow(res)){empty_vector[i] <- length(intersect(out[[i]]$Gene,bicluster))} # we only receive a number of overlapping genes.
	res$GenesInOverlap <- empty_vector
	genes_length <- length(intersect(unique(genesets$Gene),bicluster))# number of overlapping my genes with all total unique genes in defined gene sets

	N <- GenesInUni # total number of genes in define set collection
	k <- genes_length # intersect
	pvals <- vector()
	fe <- vector()
	ratioFE <- vector()
	
	for(i in 1:nrow(res)){ # for every gene set
	q <- res$GenesInOverlap[i]
	m <- res$GenesInGeneset[i]
	n <- N-m
	pvals[i] <- phyper(q-1,m,n,k,lower.tail=FALSE)
	fe[i] <- formatC(q*N/(m*k),digits = 1, format = "f")
	ratioFE[i] <- paste0(fe[i]," ","(",q,"/",m,")")
	}
	res$p_value <- pvals
	res$FE <- fe
	res$FE <- formatC(res$FE, digits = 1, format = "f")
	res$RatioAndFE <- ratioFE
	res$FDR_q_value <- p.adjust(pvals,method="fdr")
	geneset_threshold <- 50
	pcut <- 0.05

	res <- res[res$GenesInGeneset>=geneset_threshold,]
	res <- res[res$FDR_q_value<=pcut & res$GenesInOverlap>=0,] #FDR less or equal to 0.05, only rows with more than 0 genes in overlap are included
	res <- res[order(res$p_value,decreasing=FALSE),]
	res <- res[order(res$p_value,decreasing=FALSE),]

############## Bar PLOTS ###############

	res[1:10,] %>% #top 10 results
	  ggplot( aes(x=reorder(GeneSetName[1:10],-FDR_q_value[1:10]), y=-log(FDR_q_value[1:10],10), fill=as.numeric(FE))) +
	  geom_bar(stat="identity", color="black", width=1) + 
	  coord_flip() +
	  xlab("") +
	  ylab("-log10(Benjamini-Hochberg corrected p-value)") +
	  scale_y_continuous(expand = c(0, 0)) +
	  scale_x_discrete(expand = c(0, 0)) +
	  labs(fill = "FE") +
	  ggtitle(paste0(name_clust[counter],"   ",dict[[geneset_collection]])) + #" (",length(bicluster)," genes",")")
	  theme_classic() +  
	  theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1), plot.margin = grid::unit(c(1,0,1,0), "mm")) 
	ggsave(paste0(my_directory,"/",geneset_collection,"/",name_clust[counter],"_",counter,"of",num_clusters,".png"), height = 2.5 , width = 7,dpi = 600)
	counter = counter + 1
	}
}	
################################################################################3

#########################################################################################################################
############################ Inter-Bicluster correlation in R ##################################################################
##################################################################################################################

library(ggplot2)
library(ggpubr) # for stat_cor
library(robustbase) # for median calculations

SimilarityPlot <- function(data,main,x_clust1_name,y_clust2_name,outputfilename) { #legend_name_b1,legend_name_b2

ggplot(data = data,aes(x = bicluster1, y = bicluster2)) +										#CHANGE
	geom_point()+
	labs(title = main) +
	#stat_smooth(method = "lm",col = "red") +
	scale_color_grey()+
	theme_classic() +
	xlab(paste0(x_clust1_name)) +
    ylab(paste0(y_clust2_name)) +
	stat_cor(label.y.npc="top",method = "pearson") + #default
	stat_smooth(method = "lm", color = "black") 
ggsave(paste0(outputfilename,"_",n,".png"),height = 4,width = 4,dpi = 600)
 
}

num_clusters <- 5 										 
n <- 453 # 226 or 453
clust1 <- 3 #Change 
clust2 <- 5 #Change	

#### Expression #####
gene1 <- readRDS(paste0(my_directory,num_clusters,"_Biclusters/Genes_Bicluster_",clust1,"of",num_clusters,"_",n,".rds"))
gene2 <- readRDS(paste0(my_directory,num_clusters,"_Biclusters/Genes_Bicluster_",clust2,"of",num_clusters,"_",n,".rds"))
gene1_mean <- as.data.frame(colMeans(gene1))
gene2_mean <- as.data.frame(colMeans(gene2))
average <- cbind(gene1_mean,gene2_mean)
names(average) <- c("bicluster1","bicluster2")
data = average
main = "Average Expression"
outputfilename = paste0(my_directory, "_Ig_Immune_",main)
x_clust1_name <- "Immune Ig/Hormone-like"	#Change								
y_clust2_name <- "Immune"					#Change
SimilarityPlot(data,main,x_clust1_name,y_clust2_name,outputfilename) 

#### Meth ###
meth1 <- readRDS(paste0(my_directory,num_clusters,"_Biclusters/Meth_Bicluster_",clust1,"of",num_clusters,"_",n,".rds"))
meth2 <- readRDS(paste0(my_directory,num_clusters,"_Biclusters/Meth_Bicluster_",clust2,"of",num_clusters,"_",n,".rds"))
meth1_mean <- as.data.frame(colMeans(meth1))
meth2_mean <- as.data.frame(colMeans(meth2))

############ Create data frame with data for analysis ##########
average <- cbind(meth1_mean,meth2_mean)
names(average) <- c("bicluster1","bicluster2")
data = average
main = "Average Methylation"
outputfilename = paste0(my_directory,"Ig_Immune_",main)
x_clust1_name <- "Immune Ig/Hormone-like"			#Change								
y_clust2_name <- "Immune"							#Change
SimilarityPlot(data,main,x_clust1_name,y_clust2_name,outputfilename) 
######################################################################################################################


#########################################################################################################################
############################ Intra-Bicluster correlation in R ##################################################################
##################################################################################################################
library(ggplot2)
library(ggpubr) # for stat_cor
library(robustbase) # for median calculations

setwd(my_directory)
num_clusters <- 5
n <- 453 # nymber of samples	
clust_name = c("Cell Cycle","Hormone-like I","Immune Ig/Hormone-like","Hormone-like II","Immune")

for (i in 1:num_clusters){
meth <- readRDS(paste0(my_directory,num_clusters,"_Biclusters/Meth_Bicluster_",i,"of",num_clusters,"_",n,".rds"))
express <- readRDS(paste0(my_directory,num_clusters,"_Biclusters/Genes_Bicluster_",i,"of",num_clusters,"_",n,".rds"))

meth_mean <- as.data.frame(colMeans(meth))
express_mean <- as.data.frame(colMeans(express))

average <- cbind(meth_mean,express_mean)
names(average) <- c("meth","expression")


ggplot(data = average,aes(x = meth, y = expression)) +									
	geom_point()+
	labs(title = paste0(clust_name[i])) +
	theme_classic() +
	xlab("Average Methylation") +
    ylab("Average Expression") +
	stat_cor(label.y.npc="bottom",method = "pearson") +
	stat_smooth(method = "lm", color = "black") 
ggsave(paste0("Bicluster_",i,"of",num_clusters,"_n=",n,".png"),height = 4,width = 4,dpi = 600)
}
###########################################################

#########################################################################################################################
############################ Enrichment of emQTL-CpGs at cis-regulatory elements in R ##################################################################
##################################################################################################################
###### Creating intersect regions ######

setwd(my_directory)
options(scipen=100) # scientific penalty for writing number in scientific format. Here we "turn off" the scientific annotation.

refseq<-read.table(paste0(my_directory, "/refseq_genes_hg19_.txt"),header=T,sep="\t",na.strings="NA",row.names=NULL,stringsAsFactors=F)
 refseq$TSS<-rep(NA,nrow(refseq))
 refseq$TSS[refseq$strand=="+"]<-refseq$txStart[refseq$strand=="+"]
 refseq$TSS[refseq$strand=="-"]<-refseq$txEnd[refseq$strand=="-"]
 refseq<-refseq[,c(9,2,10,3)]
 refseq$strand[refseq$strand=="+"] <- "1"
 refseq$strand[refseq$strand=="-"] <- "-1"
 refseq$strand <- as.integer(refseq$strand)
 refseq$Chr <- gsub("chr","",refseq$chrom,fixed=T)
 refseq$Chr <- gsub("X","23",refseq$Chr,fixed=T)
 refseq$Chr <- gsub("Y","24",refseq$Chr,fixed=T)
refseq <- refseq[refseq$Chr %in% as.character(1:24) , ]
 refseq$Chr <- as.integer(refseq$Chr)
refseq <- refseq[ , -2]
colnames(refseq)[1] <- "Gene"
colnames(refseq)[3] <- "Strand"

emptyChr <- c()
for (i in 1:nrow(refseq)){
	emptyChr[i] = paste0("chr",refseq$Chr[i])
}
refseq$chrom <- emptyChr

################## Prosessing gene data into BED file ###############
windowsizes_string <- c("1000000","500000","100000","10000","1000","500","200")
window_num = c(1000000,500000,100000,10000,1000,500,200)

for (size in 1:(length(window_num)-1)){
emptyStart <- c()
emptyEnd <- c()

for (i in 1:nrow(refseq)){
	if  (refseq$Strand[i] == 1){
		emptyEnd[i] = as.numeric(refseq$TSS[i]+(window_num[size]-window_num[size+1]))
		emptyStart[i] = as.numeric(refseq$TSS[i])
	}
	else if (refseq$Strand[i] == -1){
		emptyEnd[i] = as.numeric(refseq$TSS[i])
		emptyStart[i] = as.numeric(refseq$TSS[i]-(window_num[size]-window_num[size+1]))
	}
}
refseq_new = refseq
refseq_new$chromStart <- emptyStart
refseq_new$chromEnd <- emptyEnd
refseq_new <- refseq_new[c("chrom","chromStart","chromEnd", "Gene")]
sum(refseq_new$chromStart < 0) #80 # måtte filtrere pga cutoff på minst avstand, derfor fikk noen negative verdier.
sum(refseq_new$chromEnd < 0)
refseq_new <- refseq_new[refseq_new$chromStart > 0,]
names(refseq_new) <- NULL
print(sum(is.na(refseq_new)))
write.table(refseq_new,paste0("refseq_genes_hg19_",window_num[size],"bp.bed"),sep="\t",row.names=F,col.names=F, quote=F)
}

########################## Intersection  in Linux #################
cd my_directory # Define directory where files to be intersected are found
bedtools intersect -wa -wb -a Probeinfo450K_hg19.bed -b refseq_genes_hg19_1000000bp.bed > Intersect_Probeinfo450K_x_refseq_genes_hg19_1000000bp.bed # Change for every region.
#######################################################################################################

##### Enruchemtn of CpGs in the defined regions #######
########## CpG Enrichment #########
probeinfo <- read.table(paste0(my_directory, "Probeinfo2017.txt"),header=T,sep="\t",na.strings="NA",row.names=NULL,stringsAsFactors=FALSE)
probeinfo <- probeinfo[!probeinfo$Chr=="Y",]
probeinfo <- probeinfo[!probeinfo$Chr=="X",]

N <- nrow(probeinfo); print(paste("N =",N))
bicluster_CpG <- readRDS(paste0(my_directory, "/CpG_list_5_biclusters.rds")
bicluster_Genes <- readRDS(my_directory, "/Genes_list_5_biclusters.rds")
num_clusters = 5
name_clust <- c("Cell Cycle","Hormone-like I", "Immune Ig/Hormone-like", "Hormone-like II", "Immune") # Change
windowsizes_string <- c("-1000000","-500000","-100000","-10000","-1000","-500","-200","+200","+500","+1000","+10000","+100000","+500000","+1000000")
windowsizes <- as.numeric(windowsizes_string)
res <- matrix(NA,nrow=num_clusters,ncol=length(windowsizes))
rownames(res) <- name_clust
colnames(res) <- windowsizes_string
library(pheatmap)

for(size in 1:length(windowsizes_string)){
  print(size)             
                int <- read.table(paste0("Intersect_Probeinfo450K_x_refseq_genes_hg19_",windowsizes_string[size],"bp.bed"),sep="\t", header=F, row.names=NULL)
               
                for(i in 1:num_clusters){
                              
                               bicluster <- name_clust[i]
                               print(bicluster)
                               genes <- bicluster_Genes[[i]]
                               probes <- bicluster_CpG[[i]]
                               int_genes <- int[int$V8 %in% genes,]
                               int_genes_probes <- int[int$V8 %in% genes & int$V4 %in% probes,]
 
                               # k --- number of balls drawn from the urn
                               k <- as.numeric(length(probes));print(paste("k =",k))
                              
                               # m --- number of white balls in the urn
                               m <- length(unique(int_genes$V4)); print(paste("m =",m))
                              
                               n <- N-m; print(paste("n =",n))
                              
                               # q --- number of white balls (success) drawn without replacement from an urn
                               q <- length(unique(int_genes_probes$V4)); print(paste("q =",q))
 
 
                               # enrichment
                               fe <- q*N/(as.numeric(m)*as.numeric(k))
                               p <- phyper(q=q, m=m, n=n, k=k, lower.tail=F)
                              
                               print(paste("FE =",fe))
                               print(paste("p =",p))
                               res[i,size] <- fe
 
                }
}
 
m <- res
annot_col <- NA
annot_row <- NA
k_hc <- NA
k_hr <- NA
dist <- NA
link <- NA
tag <- paste("emQTL_biclusters_distanseEnrichment")
paletteLength <- 100
cols <- colorRampPalette(c("white","red"))(n=paletteLength)
mybreaks <- seq(from=1, to=3, length.out=paletteLength+1) 
png(paste("test_UpandDownstream_Heatmap_",tag,".png",sep=""),height=400,width=700)
pheatmap(m, color=cols, show_rownames=T, show_colnames=T,
 clustering_method=link, cluster_rows=F, cluster_cols=F, clustering_distance_rows=dist, clustering_distance_cols=dist,
 annotation_col=annot_col, annotation_row=annot_row, annotation_colors=annot_colors,
cutree_rows=k_hr,cutree_cols=k_hc,breaks=mybreaks)
dev.off()
##################################################################

#########################################################################################################################
############################ SEGWAY in R ##################################################################
##################################################################################################################

#CHECK OUT DATA used for SEGWAY
#probeinfo <- readRDS("/open/work/Jorgen/Data/SEGWAY/Processed data/probeinfo_segway.rds")
#cell.line <- "A549" #MCF7 or HCT-116							#CHANGE
#temp <- temp[,colnames(temp)%in%cell.line,drop=FALSE]
#probeinfo_SEGWAY <- readRDS("/open/work/Jorgen/Data/SEGWAY/Processed data/A549_probeinfo.rds")

library(dplyr)

#probeinfo <- readRDS("/open/work/Jorgen/Data/SEGWAY/Processed data/A549_probeinfo.rds") #This is for segway. 
#chromatin.states <- colnames(table(unique(probeinfo)))
probeinfo <- readRDS("/open/work/Jorgen/Data/SEGWAY/Processed data/probeinfo_segway.rds")
chromatin.states <- colnames(table(unique(probeinfo)))
cell.line <- "A549" #MCF7 or HCT-116
							#CHANGE
temp <- probeinfo[,colnames(probeinfo)%in%cell.line,drop=FALSE]
temp$Probe <- rownames(temp)

probeinfo_illumina <- read.table("/data2/thomas/breast450k/Probeinfo2017.txt",header=T,sep="\t",na.strings="NA",row.names=NULL,stringsAsFactors=FALSE)
probeinfo_illumina <- probeinfo_illumina[,c("Probe","Chr")]

CpG_X <- probeinfo_illumina[probeinfo_illumina$Chr == "X",]
CpG_Y <- probeinfo_illumina[probeinfo_illumina$Chr == "Y",]
probeinfo_illumina <- probeinfo_illumina[!probeinfo_illumina$Probe %in% CpG_X$Probe,]
probeinfo_illumina <- probeinfo_illumina[!probeinfo_illumina$Probe %in% CpG_Y$Probe,]

temp <- temp[intersect(temp$Probe,probeinfo_illumina$Probe),]
temp <- select(temp,-"Probe")

num_clusters <- 5							#CHANGE
read_CpGs <- readLines(paste0("/open/tmp/Anastasia/LUAD/Gender_Filtered_LUAD/All_Biclusters/List_all_CpGs_",num_clusters,"_biclusters.txt")) #read the txt file with biclusters
m_bicluster <- as.matrix(read_CpGs) #convert list to a matrix
sep_matrix <- strsplit(m_bicluster," ") # separate genes inside the list


counter = 1 
for (line in 1:length(sep_matrix)){

	if  (length(sep_matrix[[line]]) == 0){
		next
		}
  probes <- unlist(sep_matrix[line])
  #print(my.probes)
 # -- Freq.y all probes
Freq.y <- data.frame(matrix(nrow=length(chromatin.states),ncol=ncol(temp))); colnames(Freq.y) <- colnames(temp); rownames(Freq.y) <- chromatin.states
te <- lapply(temp,table)
Freq.y <- data.frame(unlist(te))

# -- Freq.x probes
Freq.x <- data.frame(matrix(nrow=length(chromatin.states),ncol=ncol(temp))); colnames(Freq.x) <- colnames(temp); rownames(Freq.x) <- chromatin.states
selection <- temp[rownames(temp)%in%probes,,drop=FALSE]
te <- lapply(selection,table)
Freq.x <- data.frame(unlist(te))

pvalue.cutoff <- 0.05
p.adjust.method <- "BH"
fold.enrichment.cutoff <- 0

freq <- merge(Freq.x,Freq.y,by="row.names",all=TRUE)
colnames(freq) <- c("names","Freq.x","Freq.y")

freq$FE <- NA
freq$p.value <- NA

for(i in 1:nrow(freq)){
q <- as.numeric(freq$Freq.x[i])
m <- as.numeric(freq$Freq.y[i])
N <- as.numeric(nrow(temp))
n <- as.numeric(N-m)
k <- as.numeric(length(probes))
freq$p.value[i] <- phyper(q-1,m,n,k,lower.tail=FALSE)
freq$FE[i] <- q*N/(m*k)}

freq$BH <- p.adjust(freq$p.value,method=p.adjust.method)

freq$Freq.x[is.na(freq$Freq.x)] <- 0
freq$p.value[is.na(freq$p.value)] <- 1
freq$BH[is.na(freq$BH)] <- 1
freq$FE[is.na(freq$FE)] <- 0

freq <- freq[freq$BH<=pvalue.cutoff & freq$FE>=fold.enrichment.cutoff,]
 
freq <- freq[order(freq$p.value,decreasing=FALSE),]

freq$Region <- freq$names

freq$Region <- unlist(strsplit(freq$Region,".",fixed=TRUE))[seq(from=2,to=length(unlist(strsplit(freq$Region,".",fixed=TRUE))),by=2)]
freq$names <- unlist(strsplit(freq$names,".",fixed=TRUE))[seq(from=1,to=length(unlist(strsplit(freq$names,".",fixed=TRUE))),by=2)]

freq$p.value <- signif(freq$p.value,digits=4)
freq$BH <- signif(freq$BH,digits=4)
freq$FE <- signif(freq$FE,digits=4)

result <- freq

print("--------------") 
print(counter)
print(result) 
print("--------------") 
 
  colnames(result)[which(names(result)=="names")] <- "Cell line"
  write.table(result, paste0(my_directory, "/Segway_bicluster_",counter,"of",num_clusters,".txt"))
  counter = counter + 1
}

### create a barplot from SEGWAY results ###
setwd(my_directory)
library(dplyr)
library(ggplot2)

name_clust <- c("Cell Cycle","Hormone I", "Immune Ig", "Hormone II", "Immune")
for (i in 1:num_clusters){
res = readRDS(paste0(my_directory, "Segway_bicluster_",i,"of",num_clusters,".txt"))
	res %>%
	  ggplot() +
	  geom_bar(aes(x=-log(BH,10), y=reorder(Region,-(as.numeric(BH))) , fill=FE), stat="identity", colour="black", width=1,) + 
	  scale_fill_gradient(low="white", high="red") +
	  ylab("") +
	  xlab("-log10(Benjamini-Hochberg corrected p-value)") +
	  scale_x_continuous(expand = c(0, 0)) +
	  scale_y_discrete(expand = c(0, 0)) +
	  labs(fill = "FE") +
	  ggtitle(paste0(name_clust[i])) + 
	  theme_classic() +  
	  theme(axis.line = element_line(size = 0.2, colour = "black", linetype=1), plot.margin = grid::unit(c(1,0,1,0), "mm")) 

	ggsave(paste0("Barplot_SEGWAY_LUAD_",name_clust[i],".png"), height = 1.5 , width = 4,dpi = 600)
}
##############################################################################################


#########################################################################################################################
############################ enrichment of emQTL-CpGs in TENET-defined enhancer regions in R ##################################################################
##################################################################################################################
# Use positions only with all TRUE (TENET-annotated enhancers and promoters simultaniously for every cell line/tissue sample 
# Probeinfo lifted over to hg38
a <- read.csv(paste0(my_directory, "/open/tmp/Anastasia/LUAD/Gender_Filtered_LUAD/TENET/Probeinfo450K_hg38.bed"), sep = "\t", header = F)  
write.table(a,paste0("/open/tmp/Anastasia/LUAD/Gender_Filtered_LUAD/TENET/Probeinfo450K_hg38_integers.bed"),sep="\t",row.names=F,col.names=F, quote=F)

# enhancer and open chromatin data direectly from TENET
b1 <- read.csv(paste0(my_directory, "TENET/tenet_enhancers_bed_allTRUE_7270x3.bed"), sep = "\t", header = F) # 7 270    3
write.table(b1,paste0(my_directory, "/TENET/tenet_enhancers_bed_allTRUE_7270x3_integers.bed"),sep="\t",row.names=F,col.names=F, quote=F)
b2 <- read.csv(paste0(my_directory, "/TENET/tenet_openChr_bed_allTRUE_13837x3.bed"), sep = "\t", header = F) # 13 837     3
write.table(b2,paste0(my_directory, "/TENET/tenet_openChr_bed_allTRUE_13837x3_integers.bed"),sep="\t",row.names=F,col.names=F, quote=F)

##### intersect in Linux
# Intersect with enhancer data
cd "/open/tmp/Anastasia/LUAD/Gender_Filtered_LUAD/TENET" # Define directory where files to be intersected are found
bedtools intersect -wa -wb -a Probeinfo450K_hg38_integers.bed -b tenet_enhancers_bed_allTRUE_7270x3_integers.bed  > Intersect_Probeinfo450K_hg38_x_tenet_enhancers_allTRUE_122488_integers.bed 

### R ####
#Filter out non-overlapping CpGs.
filt1 <- a[a$V4 %in% unique(axb1$V4),]
write.table(filt1,paste0(my_directory, "/Probeinfo450K_after_1st_inersect_allTRUE_integers",dim(filt1)[1],"x",dim(filt1)[2],"_integers.bed"),sep="\t",row.names=F,col.names=F, quote=F)
# Intersect with open chromatin regions data
bedtools intersect -wa -wb -a Probeinfo450K_after_1st_inersect_allTRUE_integers80207x4_integers.bed -b tenet_openChr_bed_allTRUE_13837x3.bed  > Intersect_2nd_Probeinfo450K_hg38_x_tenet_openchr_allTRUE_385429_integers.bed 

# Enrichement in TENET active regulatory regions ############
library(dplyr)
library(ggplot2)

setwd(my_directory)

probeinfo <- read.table(paste0(my_directory, "/Probeinfo2017.txt"),header=T,sep="\t",na.strings="NA",row.names=NULL,stringsAsFactors=FALSE)
probeinfo <- probeinfo[!probeinfo$Chr=="Y",]
probeinfo <- probeinfo[!probeinfo$Chr=="X",]
probeinfo_cpg <- unique(probeinfo$Probe) 

# Import bed file with enhancer in open chromatin data
tenet <- read.csv(paste0(my_directory, "/TENET/Intersect_2nd_Probeinfo450K_hg38_x_tenet_openchr_allTRUE_385429_integers.bed"), sep = "\t", header = F) # 75492     7
tenet_cpg <- unique(tenet$V4)

num_clusters = 5
name_clust <- c("Cell Cycle","Hormone I", "Immune Ig", "Hormone II", "Immune") # Change
biclusters_empty <- c()
pval_empty <- c()
df <- data.frame(matrix(ncol = 4, nrow = num_clusters))
name <- c("Bicluster","pval","BH","FE")
names(df) <- name
for(i in 1:num_clusters){
	bicluster <- readRDS(paste0(my_directory, "/CpGs_Bicluster_",i,"of",num_clusters,".rds"))
	my_probes <- rownames(bicluster)							        # 
	N = as.numeric(length(probeinfo_cpg))		        
	m = as.numeric(length(intersect(probeinfo_cpg,tenet_cpg)))  #Freq.y 
	n = as.numeric(N-m)													
	q = as.numeric(length(intersect(my_probes,tenet_cpg)))		#Freq.x		
	k = as.numeric(length(my_probes))
	pval <- phyper(q-1,m,n,k,lower.tail=FALSE)
	FE <- q*N/(m*k)
	df$Bicluster[i] <- name_clust[i]
	df$pval[i] <- pval
	df$FE[i] <- FE
}
df$BH <- p.adjust(df$pval, method="BH")
	df %>%
	  ggplot() + 
	  geom_bar(aes(x=-log(BH,10), y=Bicluster, fill=FE), stat="identity", colour="black", width=1,) + 
	  scale_fill_gradient(low="white", high="red") +
	  ylab("") +
	  xlab("-log10(Benjamini-Hochberg corrected p-value)") +
	  scale_x_continuous(expand = c(0, 0)) +
	  scale_y_discrete(expand = c(0, 0)) +
	  labs(fill = "FE") +
	  ggtitle("Enrichment in enhancers within open chromatin") + 
	  theme_classic() +  
	  theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1), plot.margin = grid::unit(c(1,0,1,0), "mm")) 

	ggsave(paste0("TENET_enrichment.png"), height = 2.5 , width = 4,dpi = 600)
###########################################################################################

#########################################################################################################################
############################ enrichment of emQTL-CpGs and genes in IM-PET loops in R ##################################################################
#########################################################################################################################
setwd(my_directory)
options(scipen=100) # scientific penalty for writing number in scientific format. Here we "turn off" the scientific annotation.

#https://4dgenome.research.chop.edu/Download.html link for IM-PET and 3C data
loopinfo <- read.table("4DGenome_HomoSapiens_hg19.txt", header = TRUE, sep = "\t")
loopinfo <- loopinfo[loopinfo$Cell.Tissue=="A549",]

# Separate footA and FootB
loopinfo$ID <- 1:nrow(loopinfo)
InteractorA <- loopinfo[c("InteractorAChr","InteractorAStart","InteractorAEnd","Agene","ID")]
rownames(InteractorA) = NULL
names(InteractorA) = NULL
write.table(InteractorA,paste0("InteractorA.bed"),sep="\t",row.names=F,col.names=F, quote=F)

InteractorB <- loopinfo[c("InteractorBChr","InteractorBStart","InteractorBEnd","Bgene","ID")]
rownames(InteractorB) = NULL
names(InteractorB) = NULL
write.table(InteractorB,paste0("InteractorB.bed"),sep="\t",row.names=F,col.names=F, quote=F)

#### Run in LINUX ####
 cd "my_directory" # Define directory where files to be intersected are found
 bedtools intersect -wa -wb -a InteractorB.bed -b Probeinfo450K_hg19.bed > Intersected_InteractorB_X_Probeinfo450K_hg19.bed 
 bedtools intersect -wa -wb -a InteractorB.bed -b refseq_genes_hg19_200bp_24971x4.bed > Intersected_InteractorB_X_refseq_genes_hg19_200bp_24971x4.bed 
 bedtools intersect -wa -wb -a InteractorA.bed -b Probeinfo450K_hg19.bed > Intersected_InteractorA_X_Probeinfo450K_hg19.bed 
 bedtools intersect -wa -wb -a InteractorA.bed -b refseq_genes_hg19_200bp_24971x4.bed > Intersected_InteractorA_X_refseq_genes_hg19_200bp_24971x4.bed 

##### Run in R ######
library(dplyr)
library(ggplot2)

cpgA <- read.csv("Intersected_InteractorA_X_Probeinfo450K_hg19.bed", sep = "\t", header = F) 
geneA <- read.csv("Intersected_InteractorA_X_refseq_genes_hg19_200bp_24971x4.bed", sep = "\t", header = F)

cpgB <- read.csv("Intersected_InteractorB_X_Probeinfo450K_hg19.bed", sep = "\t", header = F) 
geneB <- read.csv("Intersected_InteractorB_X_refseq_genes_hg19_200bp_24971x4.bed", sep = "\t", header = F)

cpgAxgeneB <- merge(cpgA,geneB,by ="V5")
cpgBxgeneA <- merge(cpgB,geneA,by ="V5")
CpgGene = rbind(cpgAxgeneB,cpgBxgeneA) # 267,028
final = CpgGene[c("V1.x","V5","V9.x","V9.y")]
names(final) <- c("chr","LoopID","CpG","Gene")

## There are duplicates despite the unique IDs
noID <- final[,c("CpG","Gene")]
nodp <- noID[!duplicated(noID),]

probeinfo <- read.table("Probeinfo450K_hg19.bed", header = FALSE)
probeinfo <- probeinfo[,c("V1","V4")]
names(probeinfo) <- c("chr","CpG")

geneinfo <- read.table("refseq_genes_hg19_200bp_24971x4.bed", header = FALSE)
geneinfo <- geneinfo[,c("V1","V4")]
names(geneinfo) <- c("chr","Gene")


cis_cpg = as.data.frame(table(probeinfo$chr))
colnames(cis_cpg) = c("chr", "probecount")
cis_gene = as.data.frame(table(geneinfo$chr))
colnames(cis_gene) = c("chr", "genecount")

cis_all = merge(cis_cpg,cis_gene, by = "chr")
cis_all$cis_frec = cis_all$probecount*cis_all$genecount
cis_all$cis_frec 

name_clust <- c("Cell Cycle","Hormone I", "Immune Ig", "Hormone II", "Immune") # Change

N = as.numeric(sum(cis_all$cis_frec)) # all possible cis CpG-gene pair
m = as.numeric(nrow(nodp)) 
n = N - m 

num_clusters = 5
df <- data.frame(matrix(ncol = 8, nrow = num_clusters))
name <- c("Bicluster","pval","FE","N","n","m","q","k")
names(df) <- name
### import biclusters #####
bicluster_CpG <- readRDS(paste0(my_directory,"CpG_list_5_biclusters.rds"))
bicluster_Genes <- readRDS(paste0(my_directory, "/open/work/Anastasia/LUAD/Ind_biclusters/Genes_list_5_biclusters.rds"))
## Definitions for the enrichment analysis 
# N: number of possible CpG-gene pairs (450K and refseq) on same chromosome
# q: number of CpG-gene pairs in bicluster in loops
# m: number of CpG-gene pairs (450K and refseq) in loops 
# n=N-m: number of possible CpG-gene pairs (450K and refseq) on same chromosome not in loops
# k: number of CpG-gene pairs on same chromosome in bicluster

for (i in 1:5){
my_genes <- bicluster_Genes[[i]]
my_probes <- bicluster_CpG[[i]]
common_cis <- cis_all[cis_all$Gene %in% my_genes & cis_all$CpG %in% my_probes,]
		q = as.numeric(sum(nodp$Gene %in% my_genes & nodp$CpG %in% my_probes))# within a loop
		k = as.numeric(nrow(common_cis))
        pval <- phyper(q-1,m,n,k,lower.tail=FALSE)
        FE <- (q*N)/(m*k)
        df$Bicluster[i] <- name_clust[i]
        df$pval[i] <- pval
        df$FE[i] <- FE
        df$N[i] <- N
        df$n[i] <- n
        df$m[i] <- m
        df$k[i] <- k
        df$q[i] <- q
}
	df %>%
	  ggplot() + 
	  geom_bar(aes(x=-log(p.adjust(df$pval, method="BH"),10), y=Bicluster, fill=FE), stat="identity", colour="black", width=1,) + 
	  scale_fill_gradient(low="white", high="red") +
	  ylab("") +
	  xlab("-log10(Benjamini-Hochberg corrected p-value)") +
	  scale_x_continuous(expand = c(0, 0)) +
	  scale_y_discrete(expand = c(0, 0)) +
	  labs(fill = "FE") +
	  ggtitle("Enrichment at IM-PET loops ") + 
	  theme_classic() +  
	  theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1), plot.margin = grid::unit(c(1,0,1,0), "mm")) 

	ggsave(paste0("Looping_enrichment.png"), height = 2.5 , width = 4,dpi = 600)
#############################################################


#########################################################################################################################
############################ enrichment of emQTL-CpG at TFBSs defined by UniBind in R ##################################################################
#########################################################################################################################


######################################################


#########################################################################################################################
############################ correlation with ASCAT purity score (R) ##################################################################
#########################################################################################################################
library(dplyr)
library(ggplot2)
library(ggprism)
library(patchwork)
library(magrittr)
library(ggpubr)

num_clusters <- 5 				
n_samp <- 453
clust_name = c("Cell Cycle","Hormone-like I","Immune Ig/Hormone-like","Hormone-like II","Immune")

temp <- read.csv(paste0(my_directory, "ASCAT_v2.5.2/snp6.ascat2.metadata.tsv",sep = "\t"))
temp <- temp[,colnames(temp)%in%c("tumor.aliquot.submitter_id","project","purity")]
temp$sample_ID <- substring(temp$tumor.aliquot.submitter_id,1,16)
temp$sample_ID <- gsub("(-01).*","\\1",temp$sample_ID) #Matching sample IDs with IDs in methylation and RNA-seq data.

#Plotting data METHYLATION vs purity 
for (i in 1:num_clusters){
	bicluster = readRDS(paste0(my_directory"Meth_Bicluster_",i,"of",num_clusters,"_",n_samp,".rds"))
	temp_filt <- temp[temp$sample_ID %in% colnames(bicluster),]
	temp_final <- temp_filt[,c("sample_ID","purity")] #isolating only patient ID and purity columns
	temp_ndp <- temp_final %>% group_by(sample_ID) %>% dplyr::summarise_all(mean) #removing dublicates by averagig them
	rownames(temp_ndp) <- temp_ndp$sample_ID
	temp_df <- as.data.frame(temp_ndp)
	no_ones <- temp_df[temp_df$purity < 1,]
	rownames(no_ones) = no_ones$sample_ID
	bicluster_average <- as.data.frame(colMeans(bicluster))
	names(bicluster_average) <- c("Mean")
	merged_df <- merge(bicluster_average,no_ones, by = "row.names")
	a <- merged_df
	
	ggplot(data = a,aes(x = purity, y = Mean)) +								
	labs(title = paste0(clust_name[i])) +
	scale_color_grey()+
	theme_classic() +
	xlab("ASCAT score") +
    ylab("Average Methylation") +
	stat_smooth(method = "lm", color = "black") 
ggsave(paste0("ASCAT_Bicluster_",i,"of",num_clusters,".png"),height = 4,width = 4,dpi = 600)
}

for (i in 1:num_clusters){
	bicluster = readRDS(paste0(my_directory, "/Genes_Bicluster_",i,"of",num_clusters,"_",n,".rds"))
	temp_filt <- temp[temp$sample_ID %in% colnames(bicluster),]
	temp_final <- temp_filt[,c("sample_ID","purity")] #isolating only patient ID and purity columns
	temp_ndp <- temp_final %>% group_by(sample_ID) %>% dplyr::summarise_all(mean) #removing dublicates by averagig them
	rownames(temp_ndp) <- temp_ndp$sample_ID
	temp_df <- as.data.frame(temp_ndp)
	no_ones <- temp_df[temp_df$purity < 1,] 
	rownames(no_ones) = no_ones$sample_ID
	bicluster_average <- as.data.frame(colMeans(bicluster))
	names(bicluster_average) <- c("Mean")
	merged_df <- merge(bicluster_average,no_ones, by = "row.names")
	a <- merged_df
	
	ggplot(data = a,aes(x = purity, y = Mean)) +										
	geom_point()+
	labs(title = paste0(clust_name[i])) +
	scale_color_grey()+
	theme_classic() +
	xlab("ASCAT score") +
    ylab("Average Expression") +
	stat_smooth(method = "lm", color = "black") 
ggsave(paste0("ASCAT_Bicluster",i,"of",num_clusters,".png"),height = 4,width = 4,dpi = 600)
}
#######################################################


#########################################################################################################################
############################ single-cell RNA-Seq (R) ##################################################################
#########################################################################################################################

#Creating UMAP
# Loading required packages
   lapply(c("Matrix","readr","Seurat","ggplot2","dplyr","patchwork","umap"),require,character.only=TRUE)
# Importing seurat object 
	cancer.type = "Lung"
	data <- readRDS(paste("/open/work/Jorgen/Data/scRNA-seq/",cancer.type,"/",cancer.type,"_metadata.rds",sep=""))
ndata <- subset(data, PatientNumber == c(3,4,6)) #Extracting patients with LUAD		
# Set working directory
	setwd(my_directory)
	
png(paste0("umap_n3_LUAD.png"), width=500, height=400) 
plot <- DimPlot(ndata, reduction = "umap") #tsne  reduction = "tsne"
plot + labs(title = "Clustering of 13,670 tumor cells from 3 LUAD patients")
dev.off()

# Creating dotplots fro every bicluster
	num_clusters = 5
	a <- paste0(my_directory, "Gender_filtered_Top30_minimum_val_list_",num_clusters,"_Bicluster") # most negativly correlated Genes to a CpG. 
	mytop <- readRDS(paste0(a,".rds")) # Change
	name_clust <- c("Cell Cycle", "Mixed I", "Immune Ig", "Mixed II", "Immune") # Change	
	for(i in 1:num_clusters){
		vec <- vector()
		vec <- mytop[[i]]$Gene
		my.genes <- vec
	# Select only genes found in the scRNA-seq dataset
		genes <- my.genes[my.genes%in%data@assays$RNA@counts@Dimnames[[1]]]
		notfound <- my.genes[!my.genes%in%data@assays$RNA@counts@Dimnames[[1]]]
		features <- genes
	png(paste0(,"DotPlot_scRNA_",name_clust[i],".png"), width=650, height=480) 
	fig = DotPlot(object=data,features=c(features),scale=TRUE,dot.scale=5,scale.max=20) + #https://satijalab.org/seurat/reference/dotplot From Seurat
		  RotatedAxis() +
		  ggtitle(paste0(name_clust[i]," ",str_title)) + ##"Expression by cell type"
		  theme(axis.text.x = element_text(angle = 90,size=10),
		  plot.title = element_text(hjust = 0.5),
		  legend.title=element_text(size=10),
		  legend.text=element_text(size=8),
		  legend.position="bottom", legend.box = "horizontal")+
		  ylab("Cell types") +
		  xlab("Genes")
	print("Missing genes from scRNA data")
	print(name_clust[i])
	print(notfound)
	print(paste("__________________________"))
	print(fig)
	dev.off()
	}
	
	
#########################################################################################################################
############################ Unsupervised Hierarchical Clustering and Survival Analyses (R) ##################################################################
#########################################################################################################################
library(pheatmap)
library(survival)
library(survminer)
library(grid)
library(stringr)
# Example for patients from the Oslo cohort
setwd(my_directory)
num_clusters = 5
n_samp <- 453	
name_clust <- c("Cell Cycle","Hormone-like I","Immune IgHormone-like","Hormone-like II","Immune") # Change

clinic = read.table(paste0(my_directory, "/complete_clinical_Oslo_meth_164.txt"), header = TRUE)
smoking = read.table(my_directory, "/Oslo_methylation_135.txt", header = TRUE)
smoking = smoking[,c("SampleID","Smoking.status", "Stage")]
names(smoking) = c("sampleID","Smoking status", "Stage")
clinic_filt = merge(clinic,smoking,by = "sampleID", all = TRUE)
names(clinic_filt)[9] = "Gender"
clinic_filt$Gender = str_to_title(clinic_filt$Gender)
names(clinic_filt)[13] = "Expression subtype"

	for (sample in 1:length(clinic_filt$Stage)){
		if  (is.na(clinic_filt$Stage[sample]) == TRUE){
			next
		}
		else if  (clinic_filt$Stage[sample] == "") {
			clinic_filt$Stage[sample] = NA
		}
		else if  (clinic_filt$Stage[sample] == 1) {
			clinic_filt$Stage[sample] = "I"
		}
		else if  (clinic_filt$Stage[sample] == 2) {
			clinic_filt$Stage[sample] = "II"
		}
		else if  (clinic_filt$Stage[sample] == 3) {
			clinic_filt$Stage[sample] = "III"
		}
		else if  (clinic_filt$Stage[sample] == 4) {
			clinic_filt$Stage[sample] = "IV"
		}

	}

	
		for (sample in 1:length(clinic_filt$'Smoking status')){
		if  (is.na(clinic_filt$'Smoking status'[sample]) == TRUE){
			next
		}
		else if  (clinic_filt$'Smoking status'[sample] == "") {
			clinic_filt$'Smoking status'[sample] = NA
		}
		else if  (clinic_filt$'Smoking status'[sample] == "Never smoker") {
			clinic_filt$'Smoking status'[sample] = "Never-smoker"
		}

		else if  (clinic_filt$'Smoking status'[sample] == "Current smoker") {
			clinic_filt$'Smoking status'[sample] = "Current"
		}
		else if  (clinic_filt$'Smoking status'[sample] == "Ex-smoker" ) {
			clinic_filt$'Smoking status'[sample] = "Former"
		}
	}

	
array_oslo <- readRDS(my_directory, "/Validation_Oslo/Prosessed_data/Methylation_LUAD_456946x164.rds")

for (i in 1:num_clusters){
k = 2
bicluster = readRDS(paste0(my_directory,"/CpGs_Bicluster_",i,"of",num_clusters,"_",n_samp,".rds"))
	cpg <- intersect(rownames(bicluster),rownames(array_oslo))
	length(cpg)
	bicluster_oslo <- array_oslo[rownames(array_oslo) %in% cpg,]
	dim(bicluster_oslo)
	bicluster <- bicluster_oslo

fig <- pheatmap(bicluster,
          clustering_distance_rows = "euclidean",
          clustering_distance_cols = "euclidean",
          clustering_method = "ward.D2",
		  show_rownames = FALSE,
		  show_colnames = FALSE
		)

# SURV PLOTS #####################################
samp_sort <- sort(cutree(fig$tree_col,k = k)) # obtain two groups
cuttree <- samp_sort

#Based only on dendrogram
cuttree2 <- cuttree
df2 <- as.data.frame(cuttree2)

#Merge groups from clustering and clinical data by sample names
e <- merge(clinic_filt, df2, by.x = "sampleID", by.y = "row.names")
e$OS_time_surgery <- e$OS_time_surgery/365.25

# use Kaplan-Meier method
e_surv_object = Surv(e$OS_time_surgery,e$OS_surgery)
e_fit <- survfit(e_surv_object ~ cuttree2, data = e) 
e_fit
e_gr1_median <- format(summary(e_fit)$table[,"median"][1], digits = 3)
e_gr2_median <- format(summary(e_fit)$table[,"median"][2], digits = 3)

e_n_1 <- sum(e$cuttree2 == 1)
e_n_2 <- sum(e$cuttree2 == 2)
e$cuttree2 <- factor(e$cuttree2)

res_cox <- coxph(Surv(OS_time_surgery,OS_surgery) ~ cuttree2, data = e)
a <- summary(res_cox)
coef <- format(a$coefficients[2],digits = 3)

png(filename = paste0(name_clust[i],"/Oslo_Meth_Surv_",n_samp,"_",name_clust[i],"_k=",k,"_",i,"of",num_clusters,".png"),width = 9, height = 7, units = "in", res = 600)       #  
fig1 <- ggsurvplot(e_fit, 
           data = e,
           #title = paste0(name_clust[i],": DNA methylation"),
		   surv.median.line = "hv",
           pval = TRUE,
		   pval.method=TRUE,
           pval.size = 6,
           #pval.coord = c(2.5, 1),
           xlab = "Time (years)",
           #surv.scale = "percent",
           legend.title = paste0("Dendrogram groups (k=2)"),
           legend.labs = c(paste("\n",paste0("Group I (n=",e_n_1,")"), paste("Median: ",e_gr1_median), sep ="\n"), 
                           paste("\n",paste0("Group II (n=",e_n_2,")"),paste("Median: ",e_gr2_median, "  HR:",coef), sep ="\n") 
                            ),
           legend = c(0.8,0.8),
           font.legend = c(16),
		   font.main = c(20),
		   font.x = c(18),
		   font.y = c(18),
		   break.time.by = 2.5,
		   font.tickslab = c(16),
		   tables.theme = theme_cleantable()
           )
print(fig1)
dev.off()


################ k = 3 #############################################
k = 3
samp_sort_3 <- sort(cutree(fig$tree_col,k = k)) # obtain two groups
#################################################################################
################## SURV PLOTS #####################################
cuttree3 <- samp_sort_3
df3 <- as.data.frame(cuttree3)

f <- merge(e, df3, by.x = "sampleID", by.y = "row.names")
f_surv_object = Surv(f$OS_time_surgery,f$OS_surgery)
f_fit <- survfit(f_surv_object ~ cuttree3, data = f) 
f_fit
f_gr1_median <- format(summary(f_fit)$table[,"median"][1], digits = 3)
f_gr2_median <- format(summary(f_fit)$table[,"median"][2], digits = 3)
f_gr3_median <- format(summary(f_fit)$table[,"median"][3], digits = 3)

f_n_1 <- sum(f$cuttree3 == 1)
f_n_2 <- sum(f$cuttree3 == 2)
f_n_3 <- sum(f$cuttree3 == 3)

f$cuttree3 <- factor(f$cuttree3)

res_cox <- coxph(Surv(OS_time_surgery,OS_surgery) ~ cuttree3, data = f)
a <- summary(res_cox)
coef2 <- format( a$coefficients[1,2],digits = 3)
coef3 <- format( a$coefficients[2,2],digits = 3)

p_val_2 <- format(a$coefficients[1,5],digits = 3)
p_val_3 <- format(a$coefficients[2,5],digits = 3)

png(filename = paste0(name_clust[i],"/Oslo_Meth_Surv_",n_samp,"_",name_clust[i],"_k=",k,"_",i,"of",num_clusters,".png"),width = 9, height = 8, units = "in", res = 600)       #  
fig2 <- ggsurvplot(f_fit, 
           data = f,
           #title = paste0(name_clust[i],": DNA Methylation"),
		   surv.median.line = "hv",
           pval = TRUE,
		   pval.method=TRUE,
           pval.size = 6,
           xlab = "Time (years)",
           legend.title = paste0("Dendrogram groups (k=3)"),
           legend.labs = c(paste("\n",paste0("Group I (n=",f_n_1,")"), paste("Median: ",f_gr1_median), sep ="\n"), 
                           paste("\n",paste0("Group II (n=",f_n_2,")"),paste("Median: ",f_gr2_median),paste("HR: ",coef2,"  p = ",p_val_2), sep ="\n"),
                           paste("\n",paste0("Group III (n=",f_n_3,")"),paste("Median: ",f_gr3_median),paste("HR: ",coef3,"  p = ",p_val_3),sep ="\n")), 						   
           legend = c(0.8,0.8),
           font.legend = c(16),
		   font.main = c(20),
		   font.x = c(18),
		   font.y = c(18),
		   break.time.by = 2.5,
		   font.tickslab = c(16),
		   tables.theme = theme_cleantable()
           )
print(fig2)
dev.off()


################ HEATMAPS med Unsupervised Hierarchical clusterint ####################


h = f

	for (sample in 1:nrow(h)){
		if  (h$cuttree2[sample] == 1){
			h$'Groups (k=2)'[sample] = "I"
		} # for some reason dosn't work
		else if  (h$cuttree2[sample] == 2) {
			h$'Groups (k=2)'[sample] = "II"
		}
	}
	
unique(h$'Groups (k=2)')


	for (sample in 1:nrow(h)){
		if  (h$cuttree3[sample] == 1){
			h$'Groups (k=3)'[sample] = "I"
		} 
		else if  (h$cuttree3[sample] == 2) {
			h$'Groups (k=3)'[sample] = "II"
		}
		else if  (h$cuttree3[sample] == 3) {
			h$'Groups (k=3)'[sample] = "III"
		}
	}
unique(h$'Groups (k=3)')


annotation_col = h[ , colnames(h)%in%c("Gender","Smoking status","Stage","Expression subtype","Groups (k=2)","Groups (k=3)")] 

c <- c("Expression subtype","Gender","Smoking status","Stage","Groups (k=2)","Groups (k=3)") # denne vektoren bestemmer rekkefølgen. «PAM50» blir øverst og «pCR» nederst
annotation_col <- annotation_col[, match(c,colnames(annotation_col))]
rownames(annotation_col) <- h$sampleID
paletteLength <- 299
purity_col <- colorRampPalette(c("black","gray"))(n=paletteLength)				 
col = list(
'Expression subtype' = c("TRU" = "blue", "PP" = "yellow","PI" = "orange"),
Gender = c("Male" = "blue3", "Female" = "brown3"),
'Smoking status' = c( "Never-smoker" = "deepskyblue", "Former" = "deepskyblue4","Current" = "gray17"),
Stage = c("I" = "chartreuse3", "II" = "darkgoldenrod1", "III" = "brown2","IV" = "black"),
'Groups (k=2)' = c("I" = "red", "II" = "darkturquoise"),
'Groups (k=3)' = c("I" = "brown1", "II" = "green3", "III" = "cornflowerblue")
)

my_palette = colorRampPalette(c("blue","white","red"))(100) # palett for the gene expression/methylation data

	png(paste0(name_clust[i],"/Oslo_Meth_heatmap_",i,"of",num_clusters,".png"), width = 2500, height = 2500, res=300)                                      
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")

fig3 = pheatmap(bicluster,
          clustering_distance_rows = "euclidean",
          clustering_distance_cols = "euclidean",
          clustering_method = "ward.D2",
          annotation_col = annotation_col,
          annotation_colors = col,
		  color = my_palette,
          legend = TRUE,
		  show_rownames = FALSE,
		  show_colnames = FALSE,
		  fontsize = 10
         )
		 setHook("grid.newpage", NULL, "replace")
	grid.text(paste0(dim(bicluster)[2], " Oslo Tumors"), y=-0.025, x = 0.4, gp=gpar(fontsize=12.5))
	grid.text(paste0(dim(bicluster)[1]," CpGs"), x=-0.02, rot=90, gp=gpar(fontsize=12.5))
print(fig3)

dev.off()

}
#######################################################################################


#########################################################################################################################
############################ BoxPlots Expression/Methylation across Expression Subtypes (R) ##################################################################
#########################################################################################################################

setwd(my_directory)
library(ggplot2)
library(ggprism)
library(patchwork)
library(magrittr)
library(ggpubr)
library(reshape2)#use "melt" function to sort your data

## import clinical data
clinic <- read.table(my_directory,"/LUAD_clinicalMatrix.txt", header = TRUE, sep ="\t")
clinic_subtype <- subset(clinic, select = c(sampleID,Expression_Subtype))
clinic_subtype <- clinic_subtype[!clinic_subtype$Expression_Subtype == "",]

empty = c()
for (i in 1:nrow(clinic_subtype)){
if (clinic_subtype$Expression_Subtype[i] == "Bronchioid"){
empty[i] = "TRU"
}
else if (clinic_subtype$Expression_Subtype[i] == "Squamoid"){
empty[i] = "PI"
}
else if (clinic_subtype$Expression_Subtype[i] == "Magnoid"){
empty[i] = "PP"
}
}
clinic_subtype$Expression_subtype = empty
num_clusters = 5   #CHANGE
name_clust <- c("Cell Cycle","Hormone I","Immune Ig","Hormone II","Immune")

for(i in 1:num_clusters){
print(i)
	bicluster <- readRDS(paste0(my_directory, "/CpGs_Bicluster_",i,"of",num_clusters,"_453.rds")) # or expression data 
	bicluster_average <- as.data.frame(colMeans(bicluster))
	final <- merge(bicluster_average, clinic_subtype, by.x = "row.names", by.y = "sampleID")
	names(final) <- c("sampleID", "Mean","Exp_sub_old","Expression_subtype")
	final$Expression_subtype =  factor(final$Expression_subtype, levels=c("TRU","PP","PI"))
	
	png(filename = paste0("Wilcox_boxplot_meth_mean_vs_expsub_",name_clust[i],".png"), width = 16, height = 13, res = 600, units = "cm")       #  
		fig = ggplot(final, aes_string(x = "Expression_subtype", y = "Mean", fill = "Expression_subtype")) +
		geom_boxplot(aes_string(x = "Expression_subtype", y = "Mean", fill = "Expression_subtype"), colour = "black") + 
		labs(x = "Expression Subtype", y = "Average Methylation", fill = "Expression_subtype") + #title = name_clust[i],
		stat_compare_means(method = "wilcox.test",comparisons = my_comparisons, size = 5) +
		theme_classic(base_size = 15) + 
		theme(legend.position = "none") + 
		scale_fill_manual(values = c("blue","yellow","orange")) 
	print(fig)
	dev.off()
}

##################################################################################################################
############################ Distribution of bicluster CpGs across Illumina Annotated regions ######################
################################################################################################################
### Location of bicluster-CpGs ####

library(dplyr) 
library(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())
probeinfo <- read.table(my_directory, "/NewProbeinfo.txt",header=T,sep="\t",na.strings="NA",row.names=NULL,stringsAsFactors=F)

location = probeinfo[c("Probe","Relation_CGI")]

for (i in 1:nrow(location)){if (is.na(location$Relation_CGI[i])){location$Relation_CGI[i] = "Open Sea"}}

geneinfo <- read.table(paste0(my_directory, "/450KprobesWmatchedgenes.txt"),header=T,sep="\t",na.strings="NA",row.names=NULL,quote="",stringsAsFactors=F)
names(geneinfo) = c("Probe","Gene","ProbeID","Relation_gene")

a = merge(location,geneinfo, all.x = TRUE)
for (i in 1:nrow(a)){
print(i)
	if (a$Relation_CGI[i] == "N_Shore"){
	a$Relation_CGI[i] = "Shore"
	}
	else if (a$Relation_CGI[i] == "S_Shore"){
	a$Relation_CGI[i] = "Shore"
	}
	else if (a$Relation_CGI[i] == "N_Shelf"){
	a$Relation_CGI[i] = "Shelf"
	}
	else if (a$Relation_CGI[i] == "S_Shelf"){
	a$Relation_CGI[i] = "Shelf"
	}
}

for (i in 1:nrow(a)){
print(i)
if (is.na(a$Relation_gene[i])){
a$Relation_gene[i] = "Intergenic"
}
}

a$Relation_gene = factor(a$Relation_gene, levels = c("TSS1500", "TSS200", "5UTR","1stExon","Body","3UTR","Intergenic"))
a$Relation_CGI =factor(a$Relation_CGI, levels = c("Island","Shore","Shelf","Open Sea"))
df <- a %>%
  group_by(Relation_CGI, Relation_gene) %>%
  summarise(counts = n()) %>%
  mutate(cumulative = cumsum(counts))   

setwd(my_directory)
ggplot(df, aes(x = Relation_CGI, y = counts)) +
		geom_bar(
		aes(color = Relation_gene, fill = Relation_gene),
		stat = "identity", position = position_stack(),
    ) +
	coord_flip() +
	theme_classic() +
	scale_y_continuous(expand = c(0, 0)) +
  xlab("") +
  ylab("CpGs included in the 450K array") +
  theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1),axis.text=element_text(size=12), plot.margin = grid::unit(c(1,0,1,0), "mm")) 
ggsave(paste0("Illumina_Probes_Distribution_450K.png"), height = 4 , width = 7,dpi = 600)

### For every bicluster ###
bicluster_CpG <- readRDS(paste0(my_directory,"/CpG_list_5_biclusters.rds"))
biclucter_name = c("Cell Cycle","Hormone-like I","Immune Ig/Hormone-like","Hormone-like II","Immune")

for (i in 1:5){
my_probes <- bicluster_CpG[[i]]
get_location = a[a$Probe %in% my_probes,]

df <- get_location %>%
  group_by(Relation_CGI, Relation_gene) %>%
  summarise(counts = n()) 

ggplot(df, aes(x = Relation_CGI, y = counts)) +
		geom_bar(
		aes(color = Relation_gene, fill = Relation_gene),
		stat = "identity", position = position_stack(),
    ) +
	coord_flip() +
	theme_classic() +
	scale_y_continuous(expand = c(0, 0)) +
    xlab("") +
    ylab("CpGs included in the 450K array") +
    ggtitle(paste0(biclucter_name[i])) +
    theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1),axis.text=element_text(size=12), plot.margin = grid::unit(c(1,0,1,0), "mm")) 

ggsave(paste0("CpG_distribution_Bicluster_",i,".png"), height = 4 , width = 7,dpi = 600)
}
###########################################################

#UNIBIND
counter = 1 
 p.adjust.method="BH"

for (line in 1:length(sep_matrix)){

	if  (length(sep_matrix[[line]]) == 0){
		next
		}
  my.probes <- unlist(sep_matrix[line])
  m <- readRDS("/open/work/Jorgen/Data/UniBind/450k_annotation_UniBind_300_by_TF.rds")  
  res <- data.frame(matrix(nrow=ncol(m),ncol=5)); colnames(res) <- c("TF","Freq.x","Freq.y","FE","p_val")
  res$TF <- colnames(m)
  res$Freq.y <- colSums(m)
  temp <- m[rownames(m)%in%my.probes,]
  res$Freq.x <- colSums(temp)

	for(i in 1:nrow(res)){
	q <- res[i,"Freq.x"]
	m2 <- res[i,"Freq.y"]
	N <- nrow(m)
	n <- N - m2
	k <- sum(my.probes%in%rownames(m))
	res$p_val[i] <- phyper(q-1,m2,n,k,lower.tail=FALSE)
	res$FE[i] <- q*N/(m2*k)
}

res$p.adjust <- p.adjust(res$p_val,method=p.adjust.method)
res <- res[order(res$p_val,decreasing=FALSE),]}

 saveRDS(res, paste0(my_directory,"/UniBind_bicluster_",counter,"of",num_clusters,".rds"))
}

for (file_rds in 1:num_clusters){
  unibind_rds = readRDS(paste0(my_directory,"/UniBind_bicluster_",file_rds,"of",num_clusters,".rds"))

  unibind_rds[1:10,] %>%
  ggplot( aes(x=reorder(TF[1:10],-p.adjust[1:10]), y=-log(p.adjust[1:10],10), fill=enrichment)) +
  geom_bar(stat="identity", color="black", width=1) + 
  coord_flip() +
  xlab("") +
  ylab("-log10(Benjamini-Hochberg corrected p-value)") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(fill = "Fold enrichment") +
  ggtitle(paste0("Bicluster ",file_rds)) +
  theme_classic() 
  ggsave(paste0("/Barplot_UniBind_bicluster_",file_rds,"of",num_clusters,".png"), height = 2.5 , width = 10,dpi = 600)
  counter = counter + 1
}
########################################################
