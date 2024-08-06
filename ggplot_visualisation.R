#loading necessary packages
library(ggplot2)
library(ggfortify)

#Reading datafiles into R
ss = read.table("C:/Users/2938235s/Documents/Stats/Report/Sample_Sheet.csv", header=TRUE, row.names=1, sep="\t")
bg = read.table("C:/Users/2938235s/Documents/Stats/Report/Annotations.csv", header=TRUE, row.names=1, sep="\t")
em = read.table("C:/Users/2938235s/Documents/Stats/Report/Expression_Table.csv", header=TRUE, row.names=1, sep="\t")
de_gout = read.table("C:/Users/2938235s/Documents/Stats/Report/DE_Gout_vs_HC.csv", header=TRUE, row.names=1, sep="\t")
de_SA = read.table("C:/Users/2938235s/Documents/Stats/Report/DE_SA_vs_HC.csv", header=TRUE, row.names=1, sep="\t")

#Matching em column names to ss row names as one is in lower case and the other in upper case letters.
names(em)= row.names(ss)

#Looking at age and sex distribution
ss$GROUP=as.factor(ss$GROUP)
ss$SEX=as.factor(ss$SEX)
str(ss)
summary(ss)
ss_healthy = subset (ss, GROUP == "HEALTHY")
ss_gout = subset (ss, GROUP == "GOUT")
ss_SA = subset (ss, GROUP == "SA")
summary(ss_healthy)
summary(ss_gout)
summary(ss_SA)

#Checking for significant differences between clinical groups
shapiro.test(ss_healthy$AGE)
shapiro.test(ss_healthy$NEUTROPHILS)
shapiro.test(ss_healthy$MONOCYTES)

shapiro.test(ss_gout$AGE)
shapiro.test(ss_gout$NEUTROPHILS)
shapiro.test(ss_gout$MONOCYTES)

shapiro.test(ss_SA$AGE)
shapiro.test(ss_SA$NEUTROPHILS)
shapiro.test(ss_SA$MONOCYTES)

bartlett.test(AGE ~ GROUP, data=ss)
bartlett.test(NEUTROPHILS ~ GROUP, data=ss)
bartlett.test(MONOCYTES ~ GROUP, data=ss)

age_box = ggplot(ss, aes(fill=GROUP, x = GROUP, y=AGE)) + geom_boxplot(alpha=0.5)
age_box

kruskal.test(AGE ~ GROUP, data=ss)
kruskal.test(NEUTROPHILS ~ GROUP, data=ss)
kruskal.test(MONOCYTES ~ GROUP, data=ss)

#Aim 2, comparing significant genes with HC
#Creating annotated data
em_annotated = merge(bg, em, by.x=0, by.y=0)
de_gout_annotated = merge(bg, de_gout, by.x=0, by.y=0)
de_SA_annotated = merge(bg, de_SA, by.x=0, by.y=0)

#Creating variable of unannotated genes for background check
unannotated_ids = row.names(em)[!(row.names(em) %in% row.names(em_annotated))]

#Renaming "Row.names" column to gene id
names(em_annotated)[1] = "gene_id"
names(de_gout_annotated)[1] = "gene_id"
names(de_SA_annotated)[1] = "gene_id"

#Renaming the actual row names to gene id (tried to rename to the gene symbols, however this was complicated due to gene duplicates)
row.names(em_annotated) = em_annotated[,"gene_id"]
row.names(de_gout_annotated) = de_gout_annotated[,"gene_id"]
row.names(de_SA_annotated) = de_SA_annotated[,"gene_id"]

#Creating subsets of significant genes
de_gout_annotated_sig = subset(de_gout_annotated, p.adj<0.05)
de_SA_annotated_sig = subset(de_SA_annotated, p.adj<0.05)

#Expression table of significant genes for gout and SA
samples_gout = row.names(subset(ss, GROUP == c("GOUT")))
samples_SA = row.names(subset(ss, GROUP == c("SA")))
samples_healthy = row.names(subset(ss, GROUP == c("HEALTHY")))

sig_gout_id = row.names(de_gout_annotated_sig)
sig_SA_id = row.names(de_SA_annotated_sig)

em_gout_id_sig = em_annotated[row.names(de_gout_annotated_sig), samples_gout]
em_SA_id_sig = em_annotated[row.names(de_SA_annotated_sig), samples_SA]

#Creating column of absolute values of log2fold in de sig dataframes to easily isolate the genes with the highest log2fold change, whether that be negative or positive.
de_SA_annotated_sig$absolute = abs(de_SA_annotated_sig$log2fold)
de_gout_annotated_sig$absolute = abs(de_gout_annotated_sig$log2fold)

#making most differentiated subset of the de sig dataframes by ordering them in descending order by the "absolute" column and isolating top 2.
de_SA_annotated_sig_subset = head(de_SA_annotated_sig[order(-de_SA_annotated_sig$absolute),], 2)
de_gout_annotated_sig_subset = head(de_gout_annotated_sig[order(-de_gout_annotated_sig$absolute),], 2)
combined_subsets_id = c(row.names(de_SA_annotated_sig_subset), row.names(de_gout_annotated_sig_subset))
combined_subsets_id = unique(combined_subsets_id)
combined_subsets_id
                       
#Applying subset to em tables
em_gout_id_sig_subset = em_annotated[row.names(de_gout_annotated_sig_subset),names(em_gout_id_sig)]
em_SA_id_sig_subset = em_annotated[row.names(de_SA_annotated_sig_subset),names(em_SA_id_sig)]

#Making graphs of distribution of significant SA genes
gene_data = em_SA_id_sig_subset["ENSG00000137869",]
gene_data = data.frame(t(gene_data))
names(gene_data) = "gene"
CYP19A1SA = ggplot(gene_data, aes(x=gene)) + geom_density(colour="blue", fill="blue", alpha=1) +
  labs(x="CYP19A1  expression", y="Count", title="CYP19A1 expression distribution in SA patients")

gene_data = em_SA_id_sig_subset["ENSG00000170439",]
gene_data = data.frame(t(gene_data))
names(gene_data) = "gene"
METTL7BSA = ggplot(gene_data, aes(x=gene)) + geom_density(colour="blue", fill="blue", alpha=1) +
  labs(x="METTL7B  expression", y="Count", title="METTL7B expression distribution in SA patients")

#And of gout significant genes

gene_data = em_gout_id_sig_subset["ENSG00000229807",]
gene_data = data.frame(t(gene_data))
names(gene_data) = "gene"
XIST = ggplot(gene_data, aes(x=gene)) + geom_density(colour="blue", fill="blue", alpha=1) +
  labs(x="XIST  expression", y="Count", title="XIST expression distribution in gout patients")

gene_data = em_gout_id_sig_subset["ENSG00000170439",]
gene_data = data.frame(t(gene_data))
names(gene_data) = "gene"
METTL7Bgout = ggplot(gene_data, aes(x=gene)) + geom_density(colour="blue", fill="blue", alpha=1) +
  labs(x="METTL7B  expression", y="Count", title="METTL7B expression distribution in gout patients")

#Merging tables to be able to test effect of age, neutrophils and monocytes on gene expression
em_transposed = t(em)
ss_sig_id = ss

for (i in de_gout_annotated_sig_subset[,"gene_id"])
  ss_sig_id[,i]=as.numeric(em_transposed[,i])
  
for (i in de_SA_annotated_sig_subset[,"gene_id"])
  ss_sig_id[,i]=as.numeric(em_transposed[,i])

#Testing correlation with age
correlations = as.data.frame(matrix(0, ncol=3, nrow=nrow(em)))
names(correlations) = c(, "p", "p.adj")
row.names(de_goutSA) = row.names(em)

cor(ss$AGE,ss_sig_id$ENSG00000229807)
cor(ss$ENSG00000229807, ss$AGE, method = c("spearman"))

#ggp = ggplot(ss, aes(x=ENSG00000229807)) + geom_histogram()
#ggp
#ggp = ggplot(ss, aes(x=AGE, y=ENSG00000229807)) + geom_point()
#ggp

#Question 3 needs a lot more work

#Testing for overlap of significant genes in SA and gout by a hypergeometric test
#Merging de tables and changing column names
de_merged = merge(de_gout, de_SA, by=0)
names(de_merged) = c("gene_id", "log2fold_gout", "p_gout", "padj_gout", "log2fold_SA", "p_SA", "padj_SA")

#Sets, overlaps and differences
gout=nrow(subset(de_merged, padj_gout<0.05))
SA=nrow(subset(de_merged, padj_SA<0.05))
overlap=nrow(subset(de_merged, padj_gout<0.05 & padj_SA<0.05))
gout_only=gout-overlap
SA_only=SA-overlap
total=nrow(de_merged)

#Hypergeometric test
phyper(overlap-1, gout, total-gout, SA,lower.tail= FALSE)
#Result of hypergeometric test == 0

#Percentage of overlap in each disease calculated
#overlap/SA*100=52.6%
#overlap/gout*100=65.9%

#Finding log2fold change and p values between gout and SA
em_gout = em[samples_gout]
em_SA = em[samples_SA]

de_goutSA = as.data.frame(matrix(0, ncol=3, nrow=nrow(em)))
names(de_goutSA) = c("Log2fold", "p", "p.adj")
row.names(de_goutSA) = row.names(em)

for (row in 1:nrow(em))
{
  gene_data_gout = as.numeric(em_gout[row,])
  gene_data_SA = as.numeric(em_SA[row,])
  mean_gout = mean(gene_data_gout)
  mean_SA = mean(gene_data_SA)
  log2fold = log2(mean_SA) - log2(mean_gout) 
  p = t.test(gene_data_gout,gene_data_SA)
  p = p$p.value
  de_goutSA[row,"Log2fold"] = log2fold
  de_goutSA[row,"p"] = p
}

#Correcting p values with bonferroni method
pvalues=de_goutSA$p
p.adj = p.adjust(pvalues, method="fdr")
de_goutSA$p.adj = p.adj

#task 19
gene_data = em["ENSG00000163520",]
gene_data = data.frame(t(gene_data))
names(gene_data) = "gene"
gene_data$age = ss$AGE_halves
ggp = ggplot(gene_data, aes(fill=ss$AGE_halves, x = age, y=gene)) + geom_boxplot(alpha=1) +
  labs(x = "Age group", y = "Count")
ggp
