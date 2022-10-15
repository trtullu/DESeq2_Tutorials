

## Load the required packages
library(tidyverse)
library(knitr)
library(DESeq2)


## Import the data
counts <- readRDS("~/Desktop/Toker_PD_Project/ForShreejoy_BigRNA_Data/countMatrix.Rds")
metadata <- readRDS("~/Desktop/Toker_PD_Project/ForShreejoy_BigRNA_Data/MetadataFull.Rds")

#-----------------------------------------------------------------------
## Metadata Pre-processing
# Step 1) Update the Cohorts names (you may do the same for other covars that will be in the ~Design)
metadata$Cohort[metadata$Cohort== "Barcelona Brain Bank"] <- "Barcelona"
metadata$Cohort[metadata$Cohort== "MS-PD Brain Bank"] <- "UK"
metadata$Cohort[metadata$Cohort== "Norway Brain Bank"] <- "Norway"
metadata$Cohort[metadata$Cohort== "KCL- London Neurodegenerative Diseases Brain Bank"] <- "UK"
metadata$Cohort[metadata$Cohort== "Netherlands Brain Bank"] <- "Netherlands"
metadata$Cohort[metadata$Cohort== "Oxford Brain Bank"] <- "UK"

# Step 2) Remove samples form metadata if they don't have RNAseq counts data 
metadata <- metadata%>% filter(!is.na(RNAseq_id_ParkOme2)) #595 samples

# you will need to remove NAs from all the covars that will be used 
# in the ~Design, such as Age, PMI, PH, etc. Use the following code:
# metadata <- metadata %>% filter(!is.na(PMI)) 


# Step 3) Set rownames
rownames(metadata) <- metadata$RNAseq_id_ParkOme2

# Step 4) Take a look at numerical variables
Features.summary <- metadata %>% group_by(Cohort) %>% 
  summarise(Agemean= mean(Age),
            Agemax= max(Age),
            RINmean = mean(RIN),
            RINmin = min(RIN),
            RINmax = max(RIN), 
            DV200mean = mean(DV200), 
            DV200min = min(DV200), 
            DV200max = max(DV200),
            PMImean=mean(PMI),
            PMImin= min(PMI),
            PMImax= max(PMI)) 

# Step 5) Check the distribution of numerical covars by histogram:
num.vars <- metadata[, c(4,5,6,11, 13,14, 18, 21,23)]
num.vars %>% gather() %>%
  ggplot(aes(x=value)) + 
  geom_histogram(fill="slategray2", alpha=.9) +
  theme_bw() +
  facet_wrap(~key, scales="free")

# Step 6) Correlation Matrix for numerical variables
#Let’s move on to bivariate statistics. 
#We are plotting a correlation matrix, in order to 
#a) check if we have features that are highly 
# correlated (which is problematic for some algorithms), 
# and b) get a first feeling about which features are correlated 
# with the target (MT.Grouping) and which are not:

CovarCor = cor(num.vars, method = "pearson", use = "complete.obs")
round(CovarCor, 2)

CovarCor %>% as.data.frame %>% mutate(var2=rownames(.)) %>%
  pivot_longer(!var2, values_to = "value") %>%
  ggplot(aes(x=name,y=var2,fill=abs(value),label=round(value,2))) +
  geom_tile() + geom_label() + xlab("") + ylab("") +
  ggtitle("Correlation matrix of our predictors") +
  labs(fill="Correlation\n(absolute):")


# Step 7) Boxplots of the associations between our continuous variables
# and the grouping variable (like Diagnosis)
metadata[, c(10, 12, 28)]%>%
  pivot_longer(!Diagnosis, values_to = "value") %>%
  ggplot(aes(x=factor(Diagnosis), y=value, color=factor(Diagnosis))) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(size=.7, width=.1, alpha=.5) +
  scale_color_manual(values=alpha(c("#E0AC57", "#F7D5A1", "#F5CAC3", "#84A59D", "#BB9590","#F28482"), 0.8)) +
  theme_bw() +
  facet_wrap(~name, scales="free")+
  theme(strip.text = element_text(size=12, face = "bold"),
        axis.text = element_text(size=8, face="bold"),
        legend.text = element_text(size=5))


# Step 8) Categorical variables
## we just use simple stacked barplots to show the 
## differences between groups:
metadataata [, c(7, 17, 19)]%>% 
  pivot_longer(!MT_Group, values_to = "value") %>%
  ggplot(aes(x=factor(value), fill=factor(MT_Group))) +
  scale_fill_manual(values=alpha(c("#67a9cf", "#ef8a62", "#fddbc7"), 0.7)) +
  geom_bar(position="fill", alpha=.7)+
  theme_minimal() +
  facet_wrap(~name, scales="free")

#------------------------------------------------------------------------------

## Prepare the Data for DESEq2 Analysis

# Step 1) Scale metadata for numerical variables
metadata$Age = scale(metadata$Age)
metadata$PMI = scale(metadata$PMI)


# Step 2) Factor and level the variables 
metadata$MT_Group <- factor(metadata$MT_Group, levels=c("NonMT_PD", "MT_PD"))
metadata$Cohort <- factor(metadata$Cohort, levels= c("Barcelona", "Norway"))
metadata$Sex <- factor(metadata$Sex, levels=c("M", "F"))


# Step 3) Counts matrix processing
# Subset counts for samples presented in metadata and make a matrix
cts <- counts[, metadata$RNAseq_id_ParkOme2] 

cts <- round(cts) %>%  # convert values to integers
  as.data.frame() %>%
  dplyr::filter(rowSums(.) > 0) # Only keep rows that have total counts above the cutoff 



# Step 4) Make sure that the columns of the count data are in the same order as rows names of the metadata
all(colnames(cts) == rownames(metadata)) # same order check
all(colnames(cts) %in% rownames(metadata)) 


#----------------------------------------------------------------------------

## Run the DESeq2 Analysis

# Step 1) Make the DESeq2 object
ddsFun <- function(x, y){
  dds <- DESeqDataSetFromMatrix(countData=x, 
                                colData=y,
                                design= ~ Age+PMI+Cohort+Sex+MT_Group) # change the design formula based on your project objective
  dds <- estimateSizeFactors(dds)
  ids <- rowSums(counts(dds, normalized=TRUE) >= 10 ) >= 3 
  dds <- dds[ids, ]
  return(dds)
}

dds <- ddsFun(x= cts.Train, y=metadata)

# Step 2) Define the reference levels
dds $MT_Group <- relevel(dds $MT_Group, "NonMT_PD")
dds $Cohort<- relevel(dds $Cohort, "Barcelona")
dds $Sex <- relevel(dds $Sex, "M" )

# Step 3) Run the DE analysis
dds  <- DESeq(dds )

#save
saveRDS(dds  , file =file.path("~/dds _with_more.than.10.Counts.Rds"))


# Step 4) Get the DEG results
resultsNames(dds)
res <- results(dds , contrast=c("MT_Group", "MT_PD", "NonMT_PD"))

res.sig <- res %>% as.data.frame() %>%
  filter(!is.na(padj)) %>% filter (padj<0.01)%>% #9703  DEG 
  filter (! between(log2FoldChange, -1, 1 )) %>% #859 DEG
  #filter (! log2FoldChange < -2.5) #%>% #852 DEG
  #filter (! log2FoldChange > 2.5) #848 DEG
# The stepwise filtering of genes with extreme log2FC values was done to narrow down the list of DEGs (as predictors) in the final Train data.

# Step 5) Add genes names/symbols column to the result file
library(org.Hs.eg.db)  
library(AnnotationDbi)

res.sig$symbol <- mapIds(org.Hs.eg.db, 
                         keys=rownames(res.sig), 
                         keytype= "ENSEMBL",
                         column= "SYMBOL" ) #This returns NA for 355 genes. 

# Save the results
write.csv(res.sig, file="DESEq2_results.csv")

#---------------------------------------------------------------------------
# Plot the results  
plotMA(res.sig, main= "MAplot_ Age related DEGs") 

### Volcano Plot: visualize the significant tags
library(ggrepel)

res.sig$diffexpressed <- "NO"
res.sig$diffexpressed[res.sig$log2FoldChange > 0.5] <- "UP"
res.sig$diffexpressed[res.sig$log2FoldChange < - 0.5] <- "DOWN"

table(res.sig$diffexpressed)
# DOWN    NO    UP 
# 49     13942    23 

# plot
ggplot(data=res.sig, 
       aes(x=log2FoldChange, y=-log10(padj), color=diffexpressed)) +
  geom_point() +
  theme_minimal() + 
  scale_color_manual(values=c("steelblue", "grey56","red3"))+
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -2.5, 0, 2.5, 5))+
  geom_text_repel(data = res.sig %>% filter(diffexpressed==c("UP", "DOWN")),
                  aes(label=symbol),
                  box.padding   = 0.1,
                  point.padding = 0.1,
                  force         = 100,
                  segment.size  = 0.1,
                  direction     = "x") +
  geom_vline(xintercept=c(-0.5, 0.5), col="grey70") +
  geom_hline(yintercept=-log10(0.05), col="grey50") + 
  ggtitle('Differentially expressed genes related to Age') + 
  ylab('-Log10(Adjusted P-Value)') + 
  xlab('-Log2 Fold Change') +
  theme(axis.line = element_line(size=1),
        axis.text= element_text(size=9, face="bold"),
        axis.title = element_text(size=10, face="bold"),
        plot.title = element_text(size=14, face="bold")
  )


#-----------------------------------------------------------------------------
## Data transformations and visualization
### *Count data transformations*
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE) #The running times are shorter 
                              # when using blind=FALSE
# this gives log2(n + 1)
ntd <- normTransform(dds)


## Effects of transformations on the variance
library("vsn")
meanSdPlot(assay(vsd))
meanSdPlot(assay(ntd))
meanSdPlot(assay(rld))

# Save
saveRDS(vsd, file =file.path("~/vsd.Allsamples_With_more.than.10.Counts.Rds"))





