#load libraries ####
library(ggplot2)
library(amap)
library(reshape2)
library(ggrepel)
library(STRINGdb)
library(clusterProfiler)
library(org.Hs.eg.db)
library(svglite)


#fixed variables ####
p_thres = 0.001
fold_thres = 2

human_db="org.Hs.eg.db"

p_colour = "forestgreen"
s_colour = "orchid4"
m_colour = "skyblue3"

group_colours = c(p_colour, s_colour, m_colour)
gradient_colours = c("blue", "white", "red")

general_path = "C:/Users/2938235s/Documents/DataExploration/Assessment/Plots"

#Functions ####
#Theme to be used in all plots
my_theme = theme(
  panel.border= element_rect(colour="black", fill = NA),
  plot.title = element_text(size=20, hjust = 0.5),
  axis.text.x = element_text(size=12),
  axis.text.y = element_text(size=12),
  axis.title.x = element_text(size=15),
  axis.title.y = element_text(size=15),
  panel.grid = element_line(colour="lightgrey"),
  panel.background = element_rect(fill = NA),
  legend.key = element_blank()
)
#Theme to be used in facet plots
exp_facet_theme= theme(
  panel.grid = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=1),
  strip.text.x = element_text(size = 10, family="sans", face="bold", vjust=1),
  panel.spacing = unit(0.5, "lines"),
  legend.position = "none"
)
#Function to add sig column, -log10p.adj and correcting the infinite p.adj to the lowesrt finite p.adj
expand_de = function(path, p_threshold, fold_threshold)
{
  de = read.table(path, header=TRUE, row.names=1, sep= "\t")
  de$sig = as.factor(de$p.adj < p_threshold & abs(de$log2fold) > fold_threshold)
  smallest_p=min(de$p.adj[de$p.adj>0])
  de$p.adj[de$p.adj==0]= smallest_p
  de$mlog10p = -log10(de$p.adj)
  return(de)
}
#Function to make a master table with one de table
make_master=function(em, annotations, de, p_threshold, fold_threshold)
{
  master = merge(em, annotations, by.x=0, by.y=0)
  master = merge(master, de, by.x=1, by.y=0)
  row.names(master) = master[,"SYMBOL"]
  names(master)[1] = "ensemblID"
  master = na.omit(master)
  order(master[,"p.adj"], decreasing=FALSE)
  sorted_order = order(master[,"p"], decreasing=FALSE)
  master = master[sorted_order,]
  em_columns = names(em)
  master$mean = rowMeans(master[,em_columns])
  master$sig = as.factor(master$p.adj < p_threshold & abs(master$log2fold) > fold_threshold)
  master_sig = subset(master, master$p.adj < p_threshold & abs(master$log2fold) > fold_threshold)
  master_sig_up = subset(master_sig, log2fold>0)
  master_sig_down = subset(master_sig, log2fold<0)
  master_non_sig = subset(master, master$p.adj > p_threshold | abs(master$log2fold) < fold_threshold)
  #Including a column for direction of regulation
  master_non_sig$direction = "a"
  master_sig_up$direction = "c"
  master_sig_down$direction = "b"
  master = rbind(master_non_sig, master_sig_up, master_sig_down)
  return(master)
}
#Functin for makinga master table with two de tables
make_master_ABN=function(em, annotations, de1, de2, p_threshold, log2fold_threshold)
{
  master = merge(de1, de2, by.x=0, by.y=0, suffixes=c(".x",".y"))
  master = merge(annotations, master, by.x=0, by.y=1)
  master = merge(em, master, by.x=0, by.y=1)
  
  row.names(master) = master[,"SYMBOL"]
  names(master)[1] = "ensemblID"
  master = na.omit(master)
  em_columns = names(em)
  master$mean = rowMeans(master[,em_columns])
  return(master)
}
#Plotting PCA
plot_pca = function(em_scaled, groups, output_path)
{
  pca = prcomp(t(as.matrix(sapply(em_scaled, as.numeric))))
  pca_coordinates = data.frame(pca$x)
  
  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  x_axis_label = paste("PC1 ", " (",prop_x, "%)",sep="")
  y_axis_label = paste("PC2 ", " (",prop_y, "%)",sep="")
  
  ggp = ggplot(pca_coordinates, aes(x=PC1, y=PC2, colour = sample_groups)) + 
    geom_point() +
    my_theme +
    labs(title="PCA plot", x=x_axis_label, y=y_axis_label) +
    scale_color_manual(values=group_colours, name="Sample group")
  ggsave(output_path, height = 3, width = 4)
  return(ggp)
}
#Plotting expression density
plot_exp_density = function(em, ncolumns, output_path)
{
  em_symbols_melted = melt(em)
  ggp = ggplot(em_symbols_melted, aes(x=log10(value), colour=variable, fill=variable)) + 
    geom_density(alpha = 0.6)+
    facet_wrap(~variable, ncol=ncolumns)+
    labs(title="Expression density", x="log10(expression value)", y="Density")+
    scale_fill_manual(values=rep(group_colours, each=3))+
    scale_colour_manual(values=rep(group_colours, each=3))+
    my_theme +
    exp_facet_theme
  ggsave(output_path, height = 5, width = 5)
  return(ggp)
}
#Plotting MA plot
plot_MA = function(master, title, p_threshold, fold_threshold, ylimit, output_path)
{
  ggp = ggplot(master, aes(x=log10(mean), y=log2fold, colour= direction)) + 
    geom_point() +
    my_theme +
    labs(title=title, x="log10(mean expression)", y="log2(fold change)") +
    ylim(c(-ylimit, ylimit)) +  
    geom_hline(yintercept=fold_thres, linetype="dashed", color = "grey", linewidth=0.6) +
    geom_hline(yintercept=-fold_thres, linetype="dashed", color = "grey", linewidth=0.6) +
    geom_hline(yintercept=0, color = "lightgrey", linewidth=0.5) +
    scale_colour_manual(values = c("black", "blue", "red"), labels=c("Insignificant","Downregulated","Upregulated"), name="Change")
  ggsave(output_path, height = 4, width = 4)
  return(ggp)
}
#Plotting volcano plot
plot_volcano = function(master, de_groups, p_threshold, fold_threshold, xlimit, ylimit, output_path)
{
  master_sig = subset(master, master$p.adj < p_threshold & abs(master$log2fold) > fold_threshold)
  master_sig_up = subset(master_sig, log2fold>0)
  master_sig_down = subset(master_sig, log2fold<0)
  sorted_order_up = order(master_sig_up[,"p.adj"], decreasing=FALSE)
  sorted_order_down = order(master_sig_down[,"p.adj"], decreasing=FALSE)
  master_sig_up = master_sig_up[sorted_order_up,]
  master_sig_down = master_sig_down[sorted_order_down,]
  master_sig_up_top10 = master_sig_up[1:5,]
  master_sig_down_top10 = master_sig_down[1:5,]
  
  ggp = ggplot(master, aes(x=log2fold,y=mlog10p, colour=direction)) + 
    geom_point() +
    labs(title=sprintf("Volcano plot %s", de_groups), x="log2(fold change)", y="-log10(p-value)") +
    my_theme+
    geom_vline(xintercept= -fold_threshold, linetype="dashed", color = "grey", linewidth=0.5) +
    geom_vline(xintercept= fold_threshold, linetype="dashed", color = "grey", linewidth=0.5) +
    geom_hline(yintercept=-log10(p_threshold), linetype="dashed", color = "grey", linewidth=0.5) +
    xlim(c(-xlimit, xlimit)) +
    ylim(c(0, ylimit)) +
    scale_colour_manual(values = c("black", "blue", "red"), labels=c("Insignificant","Downregulated","Upregulated"), name="Change")+
    geom_text_repel(data = master_sig_up_top10, aes(label=SYMBOL), show.legend = FALSE) +
    geom_text_repel(data = master_sig_down_top10, aes(label=SYMBOL), show.legend = FALSE)
  ggsave(output_path, height = 5, width = 5.5)
  return(ggp)
}
#Plotting heatmap
plot_heatmap = function(em_scaled_variant, title, x, y, output_path)
{
  em_scaled_variant = na.omit(em_scaled_variant)
  colours = c("blue","red")
  
  hm.matrix = as.matrix(em_scaled_variant)
  y.dist = Dist(hm.matrix, method="spearman")
  x.dist = Dist(t(hm.matrix), method="spearman")
  
  y.cluster = hclust(y.dist, method="average")
  x.cluster = hclust(x.dist, method="average")
  
  y.dd = as.dendrogram(y.cluster)
  x.dd = as.dendrogram(x.cluster)
  
  y.dd.reorder = reorder(y.dd,0,FUN="average")
  x.dd.reorder = reorder(x.dd,0,FUN="average")
  
  y.order = order.dendrogram(y.dd.reorder)
  x.order = order.dendrogram(x.dd.reorder)
  #Giving choice of clustering when calling the function
  if (x ==TRUE & y ==TRUE)
  {
    hm.matrix_clustered = hm.matrix[y.order,x.order]
  }
  
  if (x ==FALSE & y ==TRUE)
  {
    hm.matrix_clustered = hm.matrix[y.order,]
  }
  
  if (x ==TRUE & y ==FALSE)
  {
    hm.matrix_clustered = hm.matrix[,x.order]
  }
  
  if (x ==FALSE & y ==FALSE)
  {
    hm.matrix_clustered = hm.matrix
  }
  
  hm.matrix_clustered = melt(hm.matrix_clustered)
  
  ggp = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile() + 
    scale_fill_gradientn(colours = colorRampPalette(gradient_colours)(100), name = "Expression\n(scaled)") + 
    labs(title= title, y = "Genes", x ="Samples") + 
    my_theme +
    theme(axis.text.y = element_blank(), axis.text.x =element_text(angle =30, vjust=0.6), axis.ticks=element_blank(), legend.spacing.x = unit(0.25, 'cm'), panel.grid = element_blank())
  ggsave(output_path, height = 5, width = 5.5)
  return(ggp)
}
#Plotting a simple rug to go under any heatmap plotted from above function.
plot_rug = function(groups_column, output_path)
{
  groups_data = as.matrix(as.numeric(as.factor(groups_column)))
  groups_data = melt(groups_data)
  
  ggp = ggplot(groups_data, aes(x = value, y = Var2, fill = value)) + geom_tile() + 
    scale_fill_gradientn(colours = group_colours)+ 
    geom_tile(linetype="blank")+ 
    labs(x = "", y = "") + 
    my_theme +
    theme(legend.position="none", legend.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks=element_blank(), aspect.ratio = 0.03, panel.grid = element_blank())
  ggsave(output_path)
  return(ggp)
}
#Plotting a foldvfold plot with significant genes of each de table show in automated different colours by adding a temporary new column
plot_foldvfold= function(master, fold1, fold2, de1name, de2name, output_path)
{
  master_shared=subset(master, sig.x==TRUE & sig.y==TRUE)
  master_x= subset(master, sig.x==TRUE & sig.y==FALSE)
  master_y= subset(master, sig.x==FALSE & sig.y==TRUE)
  master_non_sig= subset(master, sig.x==FALSE & sig.y==FALSE)
  
  master_non_sig$fold_sig = "a"
  master_y$fold_sig = "b"
  master_x$fold_sig = "c"
  master_shared$fold_sig ="d"
  master = rbind(master_non_sig, master_y, master_x, master_shared)

  cor1= cor.test(fold1, fold2, method = "spearman")
  cor2= cor.test(master_shared$log2fold.x, master_shared$log2fold.y, method = "spearman")
  ggp = ggplot(master, aes(x=log2fold.x,y=log2fold.y, colour=fold_sig)) + 
    geom_point() +
    scale_colour_manual(values = c("black", p_colour, m_colour, "red"), labels=c("Not sig", "SVP sig only", "MVS sig only", "Shared sig"), name="Significant group") +
    labs(x= de1name, y=de2name, subtitle= sprintf("All genes: R = %s\nShared significant genes: R = %s", cor1$estimate, cor2$estimate))+
    my_theme
  ggsave(output_path, height = 4.5, width = 4.5)
  return(ggp)
}
#Plotting multigene box plots
plot_box_multigene = function(em_scaled_variant, groups, title, output_path)
{
  gene_data = na.omit(em_scaled_variant)
  gene_data = data.frame(t(gene_data))
  gene_data$sample_group = groups
  
  gene_data.m = melt(gene_data, id.vars="sample_group")
  ggp=ggplot(gene_data.m, aes(x=variable, y=value, colour=sample_group, fill=sample_group))+
    geom_boxplot(width = 0.5,linewidth=0.4, outlier.size = 0, alpha = 0.7) +
    my_theme +
    labs(x="Gene", y="Expression (scaled)", title=title)+
    scale_fill_manual(values = group_colours, name="Sample group")+ 
    scale_colour_manual(values = group_colours, name="Sample group") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(output_path, height = 3, width = 10)
  return(ggp)
}
#Function for turning gene symbols into entrezIDs and find overrepresentation, to go within the next function
get_pathways=function(genes, ontology, p_adj_cutoff, organism_db)
{
  genes_entrez = bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = organism_db)
  ora_results = enrichGO(gene = genes_entrez$ENTREZID, OrgDb = organism_db, readable = T, ont = ontology, pvalueCutoff = p_adj_cutoff, qvalueCutoff = 0.10)
  return(ora_results)
}
#Major function to plot all ora and gsea plots
plot_ora_gsea = function(master, p_adj_cutoff, ontology, database, sample_groups, general_path)
{
  master_sig = subset(master, sig.x == TRUE | sig.y == TRUE)
  master_sig_shared = subset(master, sig.x == TRUE & sig.y == TRUE)
  #Defining each signature
  signature_1 = row.names(subset(master_sig_shared, (log2fold.x > 0) & (log2fold.y > 0)))
  signature_2 = row.names(subset(master_sig_shared, (log2fold.x < 0) & (log2fold.y < 0)))
  signature_3 = row.names(subset(master_sig_shared, (log2fold.x > 0) & (log2fold.y < 0)))
  signature_4 = row.names(subset(master_sig_shared, (log2fold.x < 0) & (log2fold.y > 0)))
  signature_5 = row.names(subset(master_sig, (sig.x == FALSE) & (sig.y ==TRUE & log2fold.y > 0)))
  signature_6 = row.names(subset(master_sig, (sig.x == FALSE) & (sig.y ==TRUE & log2fold.y < 0)))
  signature_7 = row.names(subset(master_sig, (sig.x == TRUE & log2fold.x > 0) & (sig.y ==FALSE)))
  signature_8 = row.names(subset(master_sig, (sig.x == TRUE & log2fold.x < 0) & (sig.y ==FALSE)))
  #Finding pathways
  ora_results1= get_pathways(signature_1, ontology, p_adj_cutoff, database)
  ora_results2= get_pathways(signature_2, ontology, p_adj_cutoff, database)
  ora_results3= get_pathways(signature_3, ontology, p_adj_cutoff, database)
  ora_results4= get_pathways(signature_4, ontology, p_adj_cutoff, database)
  ora_results5= get_pathways(signature_5, ontology, p_adj_cutoff, database)
  ora_results6= get_pathways(signature_6, ontology, p_adj_cutoff, database)
  ora_results7= get_pathways(signature_7, ontology, p_adj_cutoff, database)
  ora_results8= get_pathways(signature_8, ontology, p_adj_cutoff, database)
  #Making list for looping over
  list = c(ora_results1, ora_results2, ora_results3, ora_results4, ora_results5, ora_results6, ora_results7, ora_results8)
  #Organising plots in output folder
  dir.create(sprintf("%s/ORA", general_path))
  #Making a counter to keep track of which signature each plot belongs to
  resultcount = 0
  for (result in list)
  {
    resultcount=resultcount+1
    #Preventing errors caused by no significant overrepresentation in a signature
    if (length(result$geneID)!=0)
    {
      dir.create(sprintf("%s/ORA/signature%s", general_path, resultcount))
      output_path=(sprintf("%s/ORA/signature%s", general_path, resultcount))
      
      bar = barplot(result, showCategory=10)
      path= sprintf("%s/signature%s_%s_bar.svg", output_path, resultcount, ontology)
      ggsave(path, bar, height = 5, width = 5)
      
      dot = dotplot(result, showCategory=10)
      path= sprintf("%s/signature%s_%s_dot.svg", output_path, resultcount, ontology)
      ggsave(path, dot, height = 5, width = 5)
      
      go = goplot(result, showCategory = 10)
      path= sprintf("%s/signature%s_%s_go.svg", output_path, resultcount, ontology)
      ggsave(path, go, height = 5, width = 5)
      
      cnet = cnetplot(result, categorySize="pvalue")
      path= sprintf("%s/signature%s_%s_cnet.svg", output_path, resultcount, ontology)
      ggsave(path, cnet, height = 5, width = 5)
      
      #Finding the enriched genes of each pathway and creating new dataframe
      gene_sets = result$geneID
      description = result$Description
      p.adj = result$p.adjust
      
      ora_results_table = data.frame(cbind(gene_sets, p.adj))
      row.names(ora_results_table)=description
      
      dir.create(sprintf("%s/ORA/signature%s/enriched", general_path, resultcount))
      output_path=(sprintf("%s/ORA/signature%s/enriched", general_path, resultcount))
      
      #Catching errors for tables with less than 10 significantly overrepresented pathways
      if (nrow(ora_results_table)>9)
      {
        j=10
      }
      else
      {
        j= nrow(ora_results_table)
      }
      
      for (i in 1:j)
      {
        
        enriched_gene_set = as.character(ora_results_table[i,1])
        candidate_genes = unlist(strsplit(enriched_gene_set, "/")) 
        
        title= sprintf("Signature %s (%s, enriched pathway %s)", resultcount, ontology, i)
        path= sprintf("%s/signature%s_%s_enriched%s_box.svg", output_path, resultcount, ontology,i)
        ggp=plot_box_multigene(em_scaled[candidate_genes,], sample_groups, title, path)
  
        path= sprintf("%s/signature%s_%s_enriched%s_heat.svg", output_path, resultcount, ontology, i)
        ggp= plot_heatmap(em_scaled[candidate_genes,], title, FALSE, TRUE, path)
      }
    }
  }
  #New directory for gsea results
  dir.create(sprintf("%s/gsea", general_path))
  gseacounter=0
  #looping over each of the DE data sets
  for (foldchange in c("log2fold.x", "log2fold.y"))
  {
    gseacounter = gseacounter + 1   
    dir.create(sprintf("%s/gsea/foldchange%s", general_path, gseacounter))
    output_path=(sprintf("%s/gsea/foldchange%s", general_path, gseacounter))
    
    gsea_input = master[,foldchange]
    names(gsea_input) = row.names(master)
    gsea_input = na.omit(gsea_input)
    gsea_input = sort(gsea_input, decreasing = TRUE)
    
    gse_results = gseGO(geneList=gsea_input, 
                        ont =ontology, 
                        keyType = "SYMBOL", 
                        nPerm = 10000, 
                        minGSSize = 3, 
                        maxGSSize = 800, 
                        pvalueCutoff = p_adj_cutoff, 
                        verbose = TRUE, 
                        OrgDb = database, 
                        pAdjustMethod = "none")
    
    ridge = ridgeplot(gse_results, showCategory = 10)
    path= sprintf("%s/gse%s_%s_ridge.svg", output_path, gseacounter, ontology)
    ggsave(path,  height=7, width=5)
    
    #Making dataframe of the enriched pathways
    gene_sets = gse_results$core_enrichment
    description=gse_results$Description
    p.adj = gse_results$p.adjust
    
    gse_results_table = data.frame(cbind(gene_sets, p.adj))
    row.names(gse_results_table)=description
    #Preventing errors from less than 10 significantly overrepresented pathways
    if (nrow(gse_results_table)>9)
    {
      j=10
    }
    else
    {
      j= nrow(gse_results_table)
    }
    for (i in 1:j)
    {
      enriched_gene_set = as.character(gse_results_table [i,1])
      candidate_genes = unlist(strsplit(enriched_gene_set, "/")) 
      
      title= sprintf("Foldchange %s (%s, enriched pathway %s)", gseacounter, ontology, i)
      path= sprintf("%s/gse_%s_enriched%s_box.svg", output_path, ontology, i)
      ggp=plot_box_multigene(em_scaled[candidate_genes,], sample_groups, title, path)

      path= sprintf("%s/gse%s_%s_enriched%s_heat.svg", output_path, gseacounter, ontology, i)
      ggp= plot_heatmap(em_scaled[candidate_genes,], title,  FALSE, TRUE, path)
    }
  }
}


#load tables####
em = read.table("C:/Users/2938235s/Documents/DataExploration/Assessment/EM.csv", header=TRUE, row.names=1, sep= "\t")
annotations = read.table("C:/Users/2938235s/Documents/DataExploration/Assessment/Human_Background_GRCh38.p13.csv", header=TRUE, row.names=1, sep= "\t")
ss = read.table("C:/Users/2938235s/Documents/DataExploration/Assessment/sample_sheet.csv ", header=TRUE, row.names=1, sep= "\t")
#Only loading two de tables as I am not interested in directly comparing Prolif with Senes MTD
de_mvs = expand_de("C:/Users/2938235s/Documents/DataExploration/Assessment/DE_Senes_MtD_vs_Senes.csv", p_thres, fold_thres)
de_svp = expand_de("C:/Users/2938235s/Documents/DataExploration/Assessment/DE_Senes_vs_Prolif.csv", p_thres, fold_thres)

#parsing####
sample_names = row.names(ss)
ss$sample = sample_names
sample_groups = ss$SAMPLE_GROUP
#A master for each of the DE tables
master_mvs= make_master(em, annotations, de_mvs, p_thres, fold_thres)
master_svp= make_master(em, annotations, de_svp, p_thres, fold_thres)
#A master containing both DE tables
master=make_master_ABN(em, annotations, de_mvs, de_svp, p_thres, fold_thres)

em_symbols = master[, as.vector(row.names(ss))]
em_scaled= data.frame(t(scale(data.frame(t(em_symbols)))))

#Making master subsets with either fold change being significant
master_sig = subset(master, sig.x==TRUE | sig.y ==TRUE)
#Making master subsets with both fold changes being significant
master_sig_shared = subset(master, sig.x==TRUE & sig.y ==TRUE)

#Selecting 150 random significant genes for a heatmap
sig_genes_random150 = sample(row.names(master_sig), 150)

#Plotting####
pca=plot_pca(em_scaled, sample_groups, sprintf("%s/PCA.svg", general_path))

density= plot_exp_density(em_symbols, 3, sprintf("%s/density.svg", general_path))

ma_mvs=plot_MA(master_mvs, "Senes MTD vs Senes", p_thres, fold_thres, 11.5, sprintf("%s/ma_mvs.svg", general_path))
ma_svp=plot_MA(master_svp, "Senes vs Prolif", p_thres, fold_thres, 10, sprintf("%s/ma_svp.svg", general_path))

vol_mvs=plot_volcano(master_mvs, "Senes MtD vs Senes", p_thres, fold_thres, 12, 310, sprintf("%s/vol_mvs.svg", general_path))
vol_svp=plot_volcano(master_svp, "Senes vs Prolif", p_thres, fold_thres, 10, 300, sprintf("%s/vol_svp.svg", general_path))

heat=plot_heatmap(em_scaled[sig_genes_random150,], "150 random significant genes", FALSE, TRUE, sprintf("%s/heat_random150.svg", general_path))

rug=plot_rug(sample_groups, sprintf("%s/rug.svg", general_path))

#hypergeometric test
x.n = nrow(subset(master, sig.x == TRUE))
y.n = nrow(subset(master, sig.y == TRUE))
overlap.n = nrow(subset(master, sig.x == TRUE & sig.y == TRUE))
x_only = x.n - overlap.n
y_only = y.n - overlap.n
total.n = nrow(master)
hyper = phyper(overlap.n-1, y.n, total.n-y.n, x.n, lower.tail= FALSE)
hyper_factor=as.factor(hyper<p_thres)

fold_mvs_svp=plot_foldvfold(master, master$log2fold.x, master$log2fold.y, "Senes MtD vs Senes log2(fold change)", "Senes vs Prolif log2(fold change)", sprintf("%s/foldvfold.svg", general_path))

#Defining each signature. This could be moved into the ora/gsea function, but I'm  very short on time
signature_1 = subset(master_sig_shared, (log2fold.x > 0) & (log2fold.y > 0))
signature_2 = subset(master_sig_shared, (log2fold.x < 0) & (log2fold.y < 0))
signature_3 = subset(master_sig_shared, (log2fold.x > 0) & (log2fold.y < 0))
signature_4 = subset(master_sig_shared, (log2fold.x < 0) & (log2fold.y > 0))
signature_5 = subset(master_sig, (sig.x == FALSE) & (sig.y ==TRUE & log2fold.y > 0))
signature_6 = subset(master_sig, (sig.x == FALSE) & (sig.y ==TRUE & log2fold.y < 0))
signature_7 = subset(master_sig, (sig.x == TRUE & log2fold.x > 0) & (sig.y ==FALSE))
signature_8 = subset(master_sig, (sig.x == TRUE & log2fold.x < 0) & (sig.y ==FALSE))

heat=plot_heatmap(em_scaled[row.names(signature_1),], "Signature 1", FALSE, TRUE, sprintf("%s/heat_sig1.svg", general_path))
heat=plot_heatmap(em_scaled[row.names(signature_2),], "Signature 2", FALSE, TRUE, sprintf("%s/heat_sig2.svg", general_path))
heat=plot_heatmap(em_scaled[row.names(signature_3),], "Signature 3", FALSE, TRUE, sprintf("%s/heat_sig3.svg", general_path))
heat=plot_heatmap(em_scaled[row.names(signature_4),], "Signature 4", FALSE, TRUE, sprintf("%s/heat_sig4.svg", general_path))
heat=plot_heatmap(em_scaled[row.names(signature_5),], "Signature 5", FALSE, TRUE, sprintf("%s/heat_sig5.svg", general_path))
heat=plot_heatmap(em_scaled[row.names(signature_6),], "Signature 6", FALSE, TRUE, sprintf("%s/heat_sig6.svg", general_path))
heat=plot_heatmap(em_scaled[row.names(signature_7),], "Signature 7", FALSE, TRUE, sprintf("%s/heat_sig7.svg", general_path))
heat=plot_heatmap(em_scaled[row.names(signature_8),], "Signature 8", FALSE, TRUE, sprintf("%s/heat_sig8.svg", general_path))

#Making an em table with a column for signature
metagene1 = data.frame(colMeans(em_scaled[row.names(signature_1),], na.rm=TRUE))
names(metagene1)="mean"
metagene1$group = sample_groups
metagene1$signature="Signature 1"
metagene2 = data.frame(colMeans(em_scaled[row.names(signature_2),], na.rm=TRUE))
names(metagene2)="mean"
metagene2$group = sample_groups
metagene2$signature="Signature 2"
metagene3 = data.frame(colMeans(em_scaled[row.names(signature_3),], na.rm=TRUE))
names(metagene3)="mean"
metagene3$group = sample_groups
metagene3$signature="Signature 3"
metagene4 = data.frame(colMeans(em_scaled[row.names(signature_4),], na.rm=TRUE))
names(metagene4)="mean"
metagene4$group = sample_groups
metagene4$signature="Signature 4"
metagene5 = data.frame(colMeans(em_scaled[row.names(signature_5),], na.rm=TRUE))
names(metagene5)="mean"
metagene5$group = sample_groups
metagene5$signature="Signature 5"
metagene6 = data.frame(colMeans(em_scaled[row.names(signature_6),], na.rm=TRUE))
names(metagene6)="mean"
metagene6$group = sample_groups
metagene6$signature="Signature 6"
metagene7 = data.frame(colMeans(em_scaled[row.names(signature_7),], na.rm=TRUE))
names(metagene7)="mean"
metagene7$group = sample_groups
metagene7$signature="Signature 7"
metagene8 = data.frame(colMeans(em_scaled[row.names(signature_8),], na.rm=TRUE))
names(metagene8)="mean"
metagene8$group = sample_groups
metagene8$signature="Signature 8"

metagenes=rbind(metagene1, metagene2, metagene3, metagene4, metagene5, metagene6, metagene7, metagene8)

#Plotting a faceted boxplot of gene expression in each signature. The input is not in the format that my boxplot function takes.
ggp= ggplot(metagenes, aes(x=group, y=mean, colour=group))+
  geom_boxplot(width = 0.5,linewidth=0.4, outlier.size = 0, alpha = 0.7) +
  facet_wrap(~signature, ncol=2) +
  my_theme+
  labs(x="", y="Expression (scaled)")+
  scale_fill_manual(values = group_colours, name="Sample group")+ 
  scale_colour_manual(values = group_colours, name="Sample group") +
  theme(panel.grid.minor=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(sprintf("%s/metagenes_box.svg", general_path), height=5.5, width=4)

pathways=plot_ora_gsea(master, p_thres, "BP", human_db, sample_groups, general_path)
