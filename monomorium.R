library(ggplot2)
library(RMySQL)
library(edgeR)
library(DESeq)
library(GOstats)
library(GSEABase)
library(gplots)
library(gtools)
library(pgirmess)
library(class)
library(gridExtra)
library(reshape)
library(WGCNA)

mydb = dbConnect(MySQL(), user='monomorium', password='pharaonis', dbname='monomorium', host='ECOEVO.UNIT.OIST.JP')
contrasts <- dbGetQuery(mydb,"SELECT * FROM contrasts") # read a priori contrasts
contrasts <- contrasts[,sort(names(contrasts))]  # this makes the contrasts matrix compatible with design matrix, which is sorted by level
factors <- dbGetQuery(mydb,"SELECT * FROM factors")  #read treatments

#Note: MP12 (age15), MP14 (age18), MP19 (prot) have relatively low mapping rates (<70%)

counts <- dbGetQuery(mydb,"SELECT gene_id,MP1,MP2,MP3,MP4,MP5,MP6,MP7,MP8,MP9,MP10,MP11,MP12,MP13,MP14,MP15,MP16,MP17,MP18,MP19,MP20,MP21,MP22,MP23,MP24 FROM expected_counts_genes")
row.names(counts) <- counts$gene_id # change first column to row names
counts <- subset(counts,select=-c(gene_id))  
fpkm <- dbGetQuery(mydb,"SELECT gene_id,MP1,MP2,MP3,MP4,MP5,MP6,MP7,MP8,MP9,MP10,MP11,MP12,MP13,MP14,MP15,MP16,MP17,MP18,MP19,MP20,MP21,MP22,MP23,MP24 FROM fpkm_genes")
row.names(fpkm) <- fpkm$isoform
fpkm <- subset(fpkm,select=-c(gene_id))

keep=rowSums(fpkm>= 1)>= ncol(counts)/2  #select isoforms where at least half of the libraries have FPKM >1
#create design matrix, and label rows and columns according to our contrasts
design <- model.matrix(~0+factors$factor,data=counts)
rownames(design) <- colnames(counts)
colnames(design) <- names(contrasts)

pdf('/Users/sasha/Dropbox/Manuscripts/monomorium age polyethism/plots/coverage cutoff.pdf')
ggplot(melt(fpkm),aes(x=value,color=variable))+geom_density()+scale_x_log10()+theme_bw()+geom_vline(xintercept=1, color="red")+xlab('expected counts')+ylab('density')+ theme(legend.position="none")
dev.off()

dge <- DGEList(counts=round(counts[keep,]),group=factors$factor)   #apply filtering!
dge <- calcNormFactors(dge)
dge <- estimateGLMCommonDisp(dge,design,verbose=TRUE)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge,design)
lrt <- glmLRT(fit)
# plot NMDS as a sanity check

mdsplot <- plotMDS(dge,top=500)
plot(mdsplot, main="",xlab="Dimension 1",ylab="Dimension 2")
text(mdsplot$cmdscale.out[,1],mdsplot$cmdscale.out[,2],factors$factor)

#plot heatmap for ages
ages <- grep("age",factors$factor)
cdsTop <- newCountDataSet(round(counts[rownames(topTags(lrt,1000)),ages]),factors$factor[ages])
cdsTop <- estimateSizeFactors( cdsTop )
cdsTop <- estimateDispersions( cdsTop)
vsd <- getVarianceStabilizedData(cdsTop)[1:100,]
col_dendro <- reorder(as.dendrogram(hclust(dist(t(vsd)), method="average")),1/match(factors$factor[ages],mixedsort(factors$factor)))
mycols <- colorpanel(n=19,low="white",mid="grey46",high="black")
#libraries <- which(apply(selected[grepl(i,rownames(selected)),],1,max)!=0)
color_codes <- data.frame(factor = c("age0", "age3", "age6", "age9", "age12", "age15", "age18", "groom", "nurse", "troph", "carb", "prot"), color = c("#CCFFFF","#66CCFF","#3399FF","#0066FF","#0000FF","#0000CC","#000066","#d53e4f","#fc8d59","#fee08b","#e6f598","#99d594"))
colors <-  as.vector(color_codes$color[match(factors$factor[ages], color_codes$factor)])
heatmap.2(atan(vsd),key=FALSE,trace="none",density.info="none",col=mycols,Colv=col_dendro,ColSideColors=as.vector(colors),labCol=factors$factor[ages])

#    __             __       
#   / /_____ ______/ /_______
#  / __/ __ `/ ___/ //_/ ___/
# / /_/ /_/ (__  ) ,< (__  ) 
# \__/\__,_/____/_/|_/____/  
                           

#examine task-specific gene expression
et_forager <- glmLRT(fit, contrast=makeContrasts((carb+prot)/2 - (nurse+groom+troph)/3  , levels=design))  #logFC is computed as X - Y, so positive values are upregulated in focal contast
pdf('/Users/sasha/Dropbox/Manuscripts/monomorium age polyethism/plots/diffrential expression.pdf')
par(mfrow=c(2,2))
summary(de_forager <- decideTestsDGE(et_forager, p=0.05, adjust="BH"))
plotSmear(et_forager, de.tags=rownames(fit)[as.logical(de_forager)])
title("noragers vs. others")
et_nurse <- glmLRT(fit, contrast=makeContrasts(nurse - (carb+prot+groom+troph)/4, levels=design)) 
summary(de_nurse <- decideTestsDGE(et_nurse, p=0.05, adjust="BH"))
plotSmear(et_nurse, de.tags=rownames(fit)[as.logical(de_nurse)])
title("nurses vs. others")
et_groom <- glmLRT(fit, contrast=makeContrasts(groom - (carb+prot+nurse+troph)/4, levels=design)) 
summary(de_groom <- decideTestsDGE(et_groom, p=0.05, adjust="BH"))
plotSmear(et_groom, de.tags=rownames(fit)[as.logical(de_groom)])
title("grooming vs. others")
et_troph <- glmLRT(fit, contrast=makeContrasts(troph - (carb+prot+nurse+groom)/4, levels=design)) 
summary(de_troph <- decideTestsDGE(et_troph, p=0.05, adjust="BH"))
plotSmear(et_troph, de.tags=rownames(fit)[as.logical(de_troph)])
title("trophallaxis vs. others")
dev.off()
par(mfrow=c(1,1))
et_nurse_forager <- glmLRT(fit, contrast=makeContrasts(nurse - (carb + prot)/2, levels=design)) 
summary(de_nurse_forager <- decideTestsDGE(et_nurse_forager, p=0.05, adjust="BH"))


#     _____                          __      
#    / __(_)_______     ____ _____  / /______
#   / /_/ / ___/ _ \   / __ `/ __ \/ __/ ___/
#  / __/ / /  /  __/  / /_/ / / / / /_(__  ) 
# /_/ /_/_/   \___/   \__,_/_/ /_/\__/____/  
                                           

#test to see if evolutionary rates differ between forager and nurse genes
dnds <- dbGetQuery(mydb,"SELECT * FROM dNdS") # read evolutionary rates vs S. invicta
nurse_names <- rownames(et_nurse_forager$table[de_nurse_forager == 1,])
forager_names <- rownames(et_nurse_forager$table[de_nurse_forager == -1,])
nde_names <- rownames(et_nurse_forager$table[de_nurse_forager == 0 ,])
rates <- data.frame(rate = c(dnds[na.omit(match(nurse_names,dnds$gene)),2],
	dnds[na.omit(match(forager_names,dnds$gene)),2],
	dnds[na.omit(match(nde_names,dnds$gene)),2]), 
	task = c(rep("nurse",length(na.omit(match(nurse_names,dnds$gene)))),
		rep("forager",length(na.omit(match(forager_names,dnds$gene)))),
		rep("NDE",length(na.omit(match(nde_names,dnds$gene))))))
ggplot(data=rates, aes(x=factor(task),y=rate))+geom_boxplot(notch=TRUE)+scale_y_log10()
kruskal.test(rate ~ task, data = rates)
kruskalmc(rate ~ task, data = rates)
kruskal.test(rate ~ task, data = rates[rates$task != "NDE",])

# rates$task_upreg <- 1
# rates$task_upreg[rates$task == "NDE"] <- 0
# ggplot(data=rates, aes(x=factor(task_upreg),y=rate))+geom_boxplot(notch=TRUE)+scale_y_log10()
# kruskal.test(rate ~ task_upreg, data = rates)

#are more highly expressed genes evolving faster?
cor.test(lrt$table[dnds$gene,"logCPM"],dnds$rate, method="spearman")
plot(lrt$table[dnds$gene,"logCPM"],log(dnds$rate))

#are forager genes or nurse, genes more highly expressed?
#caveat more power to detect differential expression in highly expressed genes
cpm_task <- data.frame(cpm = lrt$table$logCPM, task = de_nurse_forager)
rownames(cpm_task) <- rownames(counts[keep,])
cpm_task[forager_names,"task"] <- "forager"
cpm_task[nurse_names,"task"] <- "nurse"
cpm_task[cpm_task[,"task"] == 0,"task"] <- "nde"
kruskalmc(cpm~ task, cpm_task)
p0 <- ggplot(cpm_task,aes(x=factor(task),y=cpm))+geom_boxplot()+theme_bw()+xlab("")+ylab("log counts per million reads")

#comparing results with fire ant foraging study by Manfredini et al 2013

manfredini <- dbGetQuery(mydb, '
SELECT DISTINCT mp_id, microarray.gene_id, expression,rate FROM `mp vs sinv`
JOIN dNdS
ON `mp vs sinv`.mp_id = dNdS.gene_id
LEFT JOIN (SELECT DISTINCT gene_id, IF (logFC>0,"forager","nest") AS expression FROM  manfredini_expression
JOIN manfredini_names 
ON manfredini_expression.manfredini_id = manfredini_names.manfredini_id) AS microarray
ON sinv_id = microarray.gene_id 
GROUP BY mp_id
') #a few probes target the same gene, so just choose one of them
rownames(manfredini) <- manfredini$mp_id
fisher.test(matrix(c(
	nrow(subset(manfredini, expression=="forager")), # number of sinv forager genes with mp homologs
	nrow(manfredini[is.na(manfredini$expression),]), 
	table(forager_names %in% subset(manfredini,expression=="forager")$mp_id)["TRUE"],
	table(forager_names %in% subset(manfredini,expression=="forager")$mp_id)["FALSE"]),ncol=2,byrow=T),alternative="less")

fisher.test(matrix(c(
	nrow(subset(manfredini, expression=="nest")), # number of sinv nest genes with mp homologs
	nrow(manfredini[is.na(manfredini$expression),]), 
	table(nurse_names %in% subset(manfredini,expression=="nest")$mp_id)["TRUE"],
	table(nurse_names %in% subset(manfredini,expression=="nest")$mp_id)["FALSE"]),ncol=2,byrow=T),,alternative="less")

#comparing patterns of gene expression with honeybees
hbee <- dbGetQuery(mydb, '
SELECT `blast amel vs mp`.amel, `blast amel vs mp`.mp,`alaux dge`.task FROM `blast amel vs mp` JOIN `blast mp vs amel`
ON `blast amel vs mp`.amel = `blast mp vs amel`.amel AND `blast amel vs mp`.mp = `blast mp vs amel`.mp
LEFT JOIN `alaux dge` 
ON `alaux dge`.gene = SUBSTR(`blast amel vs mp`.amel,1,7)')

# test to see whether there is an enrichment in mp forager genes that have bee homologs
fisher.test(matrix(c(
	nrow(subset(hbee, task=="forager")), # number of honeybee forager genes with mp homologs
	nrow(hbee[is.na(hbee$task),]), 
	table(forager_names %in% subset(hbee,task=="forager")$mp)["TRUE"],
	table(forager_names %in% subset(hbee,task=="forager")$mp)["FALSE"]),ncol=2,byrow=T),alternative="less")

fisher.test(matrix(c(
	nrow(subset(hbee, task=="nurse")),  # total hb nurse gene count
	nrow(hbee[is.na(hbee$task),]), 		# number of NDE  
	table(nurse_names %in% subset(hbee,task=="nurse")$mp)["TRUE"],  # nurse genes matching mp nurse genes
	table(nurse_names %in% subset(hbee,task=="nurse")$mp)["FALSE"]),ncol=2,byrow=T),alternative="less") # genes not matching mp nurse genes

#plot heatmap for foragers, nurse and ages
polyethism <- grepl("age|nurse|forager",gsub("carb|prot","forager",factors$factor))
forager_nurse_et_names <- rownames(et_nurse_forager$table[de_nurse_forager !=0,]) #names of DE genes
#plot the 50 genes most significantly differentiated between nurses and workers
strong_dge <- order(et_nurse_forager$table[forager_nurse_et_names,"PValue"],decreasing=FALSE)[0:200]
strong_names <- rownames(et_nurse_forager$table[forager_nurse_et_names,][strong_dge,])
col_order <- mixedorder(gsub("nurse","0",factors$factor[polyethism]))
cdsTop <- newCountDataSet(round(counts[strong_names,polyethism][,col_order]),factors$factor[polyethism])
cdsTop <- estimateSizeFactors( cdsTop )
cdsTop <- estimateDispersions( cdsTop, method="pooled",fitType="local")
vsd <- getVarianceStabilizedData(cdsTop)
mycols <- colorpanel(n=19,low="black",high="white")
mycols <- colorpanel(n=19,low="white",high="darkgreen")
color_codes <- data.frame(factor = c("age0", "age3", "age6", "age9", "age12", "age15", "age18", "groom", "nurse", "troph", "prot","carb"), color = c("#CCFFFF","#66CCFF","#3399FF","#0066FF","#0000FF","#0000CC","#000066","#d53e4f","#FFFF00","#fee08b","#660066","#660066"))
colors <-  as.vector(color_codes$color[match(factors$factor[polyethism][col_order], color_codes$factor)])
pdf('/Users/sasha/Dropbox/Manuscripts/monomorium age polyethism/plots/expression heatmap.pdf')
heatmap.2(log(1+as.matrix(counts[strong_names,polyethism][1:20,col_order])),key=FALSE,trace="none",col=mycols,labRow=FALSE,density.info="none",Colv=FALSE, ColSideColors=as.vector(colors),labCol=factors$factor[polyethism][col_order],hclustfun=function(c){hclust(c, method="ward")},margins = c(10, 2))
dev.off()

#     __ __      _   ___   __
#    / //_/     / | / / | / /
#   / ,< ______/  |/ /  |/ / 
#  / /| /_____/ /|  / /|  /  
# /_/ |_|    /_/ |_/_/ |_/                              

#supervised learning of task-specific genes. 
forager_nurse <- grepl("nurse|carb|prot",factors$factor)
tasks <- ! grepl("age",factors$factor)
train <- as.data.frame(t(counts[forager_nurse_et_names,forager_nurse]))
test <- as.data.frame(t(counts[forager_nurse_et_names,! tasks]))
cl <- factor(gsub("carb|prot","forager",factors$factor[forager_nurse]))
# k==2 and k==3 give close to the same result
# http://saravananthirumuruganathan.wordpress.com/2010/05/17/a-detailed-introduction-to-k-nearest-neighbor-knn-algorithm/
#knn performs well
#http://www.stat.cmu.edu/~jiashun/Research/software/GenomicsData/papers/dudoit.pdf
model <- knn(train, test, cl, k = 3, prob=TRUE)
cbind(model,factors$factor[! tasks],attr(model,"prob"))[mixedorder(factors$factor[! tasks]),]
# of the behavioral categories
                                                                         
#    __________     __                          
#   / ____/ __ \   / /____  _________ ___  _____
#  / / __/ / / /  / __/ _ \/ ___/ __ `__ \/ ___/
# / /_/ / /_/ /  / /_/  __/ /  / / / / / (__  ) 
# \____/\____/   \__/\___/_/  /_/ /_/ /_/____/  
                                              

#load go terms, and append 
go <- dbGetQuery(mydb,'SELECT gene, GO FROM blast2go JOIN genes_isoforms ON genes_isoforms.isoform = blast2go.isoform WHERE go != ""')
go$evidence <- "ISS"
go <- go[,c("GO","evidence","gene")]
dbDisconnect(mydb)
universe <- unique(go$gene)
goFrame=GOFrame(go[go$gene %in% universe,],organism="Monomorium pharaonis")
goAllFrame=GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

# GO terms upregulated in foragers
forager_upreg_go <- hyperGTest(GSEAGOHyperGParams(name = "worker upregulated",
	geneSetCollection=gsc,geneIds = intersect(forager_names,universe),
	universeGeneIds=universe,ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over"))

# GO terms downregulated in foragers
forager_downreg_go <- hyperGTest(GSEAGOHyperGParams(name = "worker upregulated",
	geneSetCollection=gsc,geneIds = intersect(forager_names,universe),
	universeGeneIds=universe,ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "under"))
summary(forager_downreg_go)

# GO terms upregulated in nurses
nurse_upreg_go <- hyperGTest(GSEAGOHyperGParams(name = "worker upregulated",
	geneSetCollection=gsc,geneIds = intersect(nurse_names,universe),
	universeGeneIds=universe,ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over"))                                      

#  _       _______________   _____ 
# | |     / / ____/ ____/ | / /   |
# | | /| / / / __/ /   /  |/ / /| |
# | |/ |/ / /_/ / /___/ /|  / ___ |
# |__/|__/\____/\____/_/ |_/_/  |_|
                                 

# Note, this analysis takes hours...
datt <- t(counts[keep,])
adjacency = adjacency(datt, power = 12,type="signed");
TOM = TOMsimilarity(adjacency,TOMType="signed");
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(1-TOM), method = "average");

# set the minimum module size to something relatively large
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = 1-TOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)

# Calculate eigengenes
MEList = moduleEigengenes(datt, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
METree = flashClust(as.dist(1-cor(MEs)), method = "average");

plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(datt, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# plotting the fabulous ridiculogram
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# how many genes in each module?
table(moduleColors)

# Define numbers of genes and samples
nGenes = ncol(datt);
nSamples = nrow(datt);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datt)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);


#go-term enrichment in modules
# for (i in table(moduleColors)) {
for (i in c("darkturquoise")) {
	module_genes <- rownames(counts[keep,])[moduleColors == i]
	module_upreg <- hyperGTest(GSEAGOHyperGParams(name = paste(i,"upregulated"),
	geneSetCollection=gsc,geneIds = intersect(module_genes,universe),
	universeGeneIds=universe,ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over"))                                      	
	print(paste(i,"upregulated"))
	print(summary(module_upreg))
}


# module-trait correlations
plot_order <-  c(8,1,5,6,7,2,3,4,9)
traits = cbind(subset(design, select = -c(troph,groom, carb,prot)),forager=design[,"prot"]+ design[,"carb"])[,plot_order]

moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

pdf('/Users/sasha/Dropbox/Manuscripts/monomorium age polyethism/plots/module trait.pdf')
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = colnames(traits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 1,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()

moduleTraitCor_stack <- melt(moduleTraitCor)
colnames(moduleTraitCor_stack) <- c("module","stage","corr")
moduleTraitCor_stack$stage <- factor(moduleTraitCor_stack$stage,unique(moduleTraitCor_stack[,2]))
ggplot(subset(moduleTraitCor_stack, stage !="nurse" & stage != "forager"),aes(x=stage,y=corr,color=module,group=module)) + geom_line()+theme_bw()

pvals <-  c()
for (i in 1:5) 
	pvals <- c(pvals , cor.test(t(moduleTraitCor[,2:8])[,i],seq(0,18,3))$p.value)
p.adjust(pvals,method="fdr")

#############
# centrality vs selection

geneModuleMembership = as.data.frame(signedKME(datt, MEs));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

#mean_centrality <- rowSums(t(t(abs(geneModuleMembership))*as.vector(table(moduleColors))))
mean_centrality <-rowMeans(abs(geneModuleMembership))
genes_with_rates <- intersect(names(mean_centrality),manfredini$mp_id)
centrality_rate <- data.frame(rate=manfredini[genes_with_rates,"rate"],centrality=mean_centrality[names(mean_centrality) %in% genes_with_rates])
cor.test(centrality_rate$rate,centrality_rate$centrality,method="spearman")
#expression and connectivity
cor.test(lrt$table$logFC,mean_centrality,method="spearman")
plot(lrt$table$logFC,mean_centrality)

#effect of connectivity on presence of fa homolog
mean_centrality <- data.frame(centrality = mean_centrality, homolog = rep(1,length(mean_centrality)))
rownames(mean_centrality) <- rownames(counts[keep,])
mean_centrality[rownames(manfredini)[which(is.na(manfredini$gene_id))],"homolog"] <- 0
kruskal.test(centrality ~ factor(lost), data=mean_centrality)
ggplot(mean_centrality,aes(x=factor(homolog),y=centrality))+geom_boxplot()

scientific_10 <- function(x)  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
p1 <- ggplot(centrality_rate,aes(centrality,rate))+geom_point(alpha=.2)+scale_y_log10(label=scientific_10)+ stat_smooth(method="lm", se=FALSE,size=1)+theme_bw()+xlab("connectivity")

# are proteins conserved in certain roles more central?

colorGenes <- data.frame(colors=moduleColors)
rownames(colorGenes) <- rownames(counts[keep,])

#forager genes
centrality_rate$forager <- "all others"
centrality_rate[rownames(subset(manfredini, expression=="forager")),"forager"] <- "conserved"
kruskal.test(centrality ~ factor(forager), data=centrality_rate)

forager_table <- table(colorGenes[rownames(centrality_rate[centrality_rate$forager=="conserved",]),])
chisq.test(forager_table,table(moduleColors)/length(moduleColors)*sum(forager_table))

centrality_rate$nest <- "all others"
centrality_rate[rownames(subset(manfredini, expression=="nest")),"nest"] <- "conserved"
kruskal.test(centrality ~ factor(nest), data=centrality_rate)
nest_table <- table(colorGenes[rownames(centrality_rate[centrality_rate$nest=="conserved",]),])

chisq.test(nest_table,table(moduleColors)/length(moduleColors)*sum(nest_table))

#figure out centrality of monomorium polyethism biased genes
centrality <- data.frame(centrality=mean_centrality)
rownames(centrality) <- rownames(counts[keep,])
centrality$conserved <- "all others"
centrality[forager_names,"conserved"] <- "forager"
centrality[nurse_names,"conserved"] <- "nurse"
kruskalmc(centrality ~ conserved, data = centrality_rate)
p2 <- ggplot(centrality,aes(x=factor(conserved),y=centrality))+geom_boxplot(notch=TRUE)+theme_bw()+ylab("connectivity")+xlab("")


# p2 <- ggplot(na.omit(melt(centrality_rate,id=c("rate","centrality"))),aes(x=factor(value),y=centrality))+geom_boxplot(notch=TRUE)+facet_grid(~variable)+theme_bw()+xlab("")+theme(strip.background = element_blank())
pdf('/Users/sasha/Dropbox/Manuscripts/monomorium age polyethism/plots/centrality.pdf',width=8,height=4)
#grid.arrange(arrangeGrob(p1,widths = unit(0.6, "npc")),p2,widths=c(1,5),nrow=1)
grid.arrange(p0,p1,p2,nrow=1)
dev.off()

#               _          
#    ____ ___  (_)_________
#   / __ `__ \/ / ___/ ___/
#  / / / / / / (__  ) /__  
# /_/ /_/ /_/_/____/\___/  
                         

#how many nurse vs forager genes have homologs in FA?
table(forager_names %in% go$gene)
table(nurse_names %in% go$gene)
genes_withGO <- data.frame(rate=manfredini[intersect(rownames(mean_centrality), rownames(manfredini)),"rate"], forager=rep(NA,length(intersect(rownames(mean_centrality), rownames(manfredini)))),nurse=rep(NA,length(intersect(rownames(mean_centrality), rownames(manfredini)))))
rownames(genes_withGO) <- rownames(manfredini[intersect(rownames(mean_centrality), rownames(manfredini)),])
genes_withGO$conn <- mean_centrality[intersect(rownames(mean_centrality), rownames(manfredini)), "centrality"]
genes_withGO[forager_names[forager_names %in% go$gene],"forager"] <- 1
genes_withGO[forager_names[! forager_names %in% go$gene],"forager"] <- 0
genes_withGO[nurse_names[nurse_names %in% go$gene],"nurse"] <- 1
genes_withGO[nurse_names[! nurse_names %in% go$gene],"nurse"] <- 0
kruskal.test(rate ~ factor(forager),genes_withGO)
kruskal.test(rate ~ factor(nurse),genes_withGO)
ggplot(genes_withGO,aes(x=factor(nurse),y=log(rate)))+geom_boxplot()
ggplot(genes_withGO,aes(x=factor(forager),y=log(rate)))+geom_boxplot()
ggplot(genes_withGO,aes(x=factor(nurse),y=conn))+geom_boxplot()
ggplot(genes_withGO,aes(x=factor(forager),y=conn))+geom_boxplot()
kruskalmc(conn~factor(forager),genes_withGO)

