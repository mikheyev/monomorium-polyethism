library(ggplot2)
library(RMySQL)
library(edgeR)
library(DESeq)
library(GOstats)
library(GSEABase)
library(gplots)
library(gtools)
library(pgirmess)

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

#mods for 
mdsplot <- plotMDS(dge$Table[])
plot(mdsplot, main="",xlab="Dimension 1",ylab="Dimension 2")
text(mdsplot$cmdscale.out[,1],mdsplot$cmdscale.out[,2],factors$factor)


#examine task-specific gene expression
et_forager <- glmLRT(fit, contrast=makeContrasts((carb+prot)/2 - (nurse+groom+troph)/3  , levels=design))  #logFC is computed as X - Y, so positive values are upregulated in focal contast
summary(de_forager <- decideTestsDGE(et_forager, p=0.05, adjust="BH"))
et_nurse <- glmLRT(fit, contrast=makeContrasts(nurse - (carb+prot+groom+troph)/4, levels=design)) 
summary(de_nurse <- decideTestsDGE(et_nurse, p=0.05, adjust="BH"))
et_groom <- glmLRT(fit, contrast=makeContrasts(groom - (carb+prot+nurse+troph)/4, levels=design)) 
summary(de_groom <- decideTestsDGE(et_groom, p=0.05, adjust="BH"))
et_troph <- glmLRT(fit, contrast=makeContrasts(troph - (carb+prot+nurse+groom)/4, levels=design)) 
summary(de_groom <- decideTestsDGE(et_troph, p=0.05, adjust="BH"))
et_nurse_forager <- glmLRT(fit, contrast=makeContrasts(nurse - (carb + prot)/2, levels=design)) 
summary(de_nurse_forager <- decideTestsDGE(et_nurse_forager, p=0.05, adjust="BH"))

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

rates$task_upreg <- 1
rates$task_upreg[rates$task == "NDE"] <- 0
ggplot(data=rates, aes(x=factor(task_upreg),y=rate))+geom_boxplot(notch=TRUE)+scale_y_log10()
kruskal.test(rate ~ task_upreg, data = rates)

#are more highly expressed genes evolving faster?
cor.test(lrt$table[dnds$gene,"logCPM"],dnds$rate, method="spearman")
plot(lrt$table[dnds$gene,"logCPM"],log(dnds$rate))

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
	table(forager_names %in% subset(hbee,task=="forager")$mp)["FALSE"]),ncol=2,byrow=T))

fisher.test(matrix(c(
	nrow(subset(hbee, task=="nurse")),  # total hb nurse gene count
	nrow(hbee[is.na(hbee$task),]), 		# number of NDE  
	table(nurse_names %in% subset(hbee,task=="nurse")$mp)["TRUE"],  # nurse genes matching mp nurse genes
	table(nurse_names %in% subset(hbee,task=="nurse")$mp)["FALSE"]),ncol=2,byrow=T)) # genes not matching mp nurse genes

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
mycols <- colorpanel(n=19,low="black",mid="#C0C0C0",high="white")
color_codes <- data.frame(factor = c("age0", "age3", "age6", "age9", "age12", "age15", "age18", "groom", "nurse", "troph", "prot","carb"), color = c("#CCFFFF","#66CCFF","#3399FF","#0066FF","#0000FF","#0000CC","#000066","#d53e4f","#FFFF00","#fee08b","#660066","#660066"))
colors <-  as.vector(color_codes$color[match(factors$factor[polyethism][col_order], color_codes$factor)])
heatmap.2(log(1+as.matrix(counts[strong_names,polyethism][0:50,col_order])),key=FALSE,trace="none",col=mycols,labRow=FALSE,density.info="none",Colv=FALSE, ColSideColors=as.vector(colors),labCol=factors$factor[polyethism][col_order],hclustfun=function(c){hclust(c, method="ward")},margins = c(10, 2))

#supervised learning of task-specific genes. 
forager_nurse <- grepl("nurse|carb|prot",factors$factor) 
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

# GO terms upregulated in nurses
nurse_upreg_go <- hyperGTest(GSEAGOHyperGParams(name = "worker upregulated",
	geneSetCollection=gsc,geneIds = intersect(nurse_names,universe),
	universeGeneIds=universe,ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over"))

