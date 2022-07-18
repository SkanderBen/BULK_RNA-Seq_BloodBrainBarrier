library(edgeR)
library(fgsea)
library(qusage)
library(DESeq2)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichR)
library(VennDiagram)
library(ggpubr)
library(ggrepel)

## Functions ##

## Plot a gene (boxplot)
# plot.gene("SAT1",cols=con_colors)
plot.gene <-function(symbol,x='group',cols){
  
  geneid = rownames(subset(conv,SYMBOL==symbol))
  meta.tmp = meta
  meta.tmp$expr = E$E[geneid,rownames(meta)]
  if(x=='condition'){
    pl = ggplot(meta.tmp,aes(condition,expr,fill=group))+geom_boxplot(width=.75,size=.75)+theme_bw()+scale_fill_manual(values=cols)+ scale_x_discrete(labels= c("BBB", "BMB"))+
      labs(title = symbol, x = "Group", y = "Expr") +scale_color_manual(labels = c("Untreated", "ACM", "Infammation"))
  }else if(x=='group'){
    pl = ggplot(meta.tmp,aes(group,expr,fill=condition))+geom_boxplot(width=.75,size=.75)+theme_bw()+scale_fill_manual(values=cols)+ scale_x_discrete(labels= c("BBB", "BMB"))+
      labs(title = symbol, x = "Group", y = "Expr") +scale_color_manual(labels = c("Untreated", "ACM", "Infammation")) + xlab("") + ylab("")
  }
  return(pl)
}

## Method#2 pour les GSEA par geneset
fisher_enrichment<-function(cluster_markers,pathway,universe,pathway_name){
enrichment_tables=list() 
  for(cluster in 1:length(cluster_markers)){
    cluster_name=names(cluster_markers)[cluster]
    output <- lapply(pathway, function(x) {
    freq.table <- table(factor(universe %in% as.character(cluster_markers[[cluster]]), 
                          levels = c(TRUE,FALSE)), 
                    	  factor(universe %in% x, 
                          levels = c(TRUE,FALSE)))

      fit <- fisher.test(freq.table, alternative = "greater")
      interSection <- intersect(cluster_markers[[cluster]], x)
      interSection <- paste(interSection, collapse = ",")
      return(value = c(NOM_pval = fit$p.value, INTERSECT = interSection,"OR"=fit$estimate))})
    
    term_names=character()
    for (pathway.term in 1:length(output)){
      term_names[pathway.term]=pathway[[pathway.term]][1]
    }
  
    results=data.frame(do.call("rbind",output))
    results$fdr=p.adjust(as.numeric(as.character(results$NOM_pval)),method = "BH")
    results=results[order(results$fdr),]
    enrichment_tables[[cluster_name]]=results
  }

  return(enrichment_tables)
}

goBP = read.gmt('/Users/macbookpro/OneDrive/Documents/UDEM/Master/BulkData/c5.go.bp.v7.5.1.symbols.gmt')
goMF = read.gmt('/Users/macbookpro/OneDrive/Documents/UDEM/Master/BulkData/c5.go.mf.v7.5.1.symbols.gmt')
goCC = read.gmt('/Users/macbookpro/OneDrive/Documents/UDEM/Master/BulkData/c5.go.cc.v7.5.1.symbols.gmt')

counts = read.csv('/Users/macbookpro/OneDrive/Documents/UDEM/Master/BulkData/raw_counts.csv',row.names=1)
ref = select(org.Hs.eg.db, keys=rownames(counts), keytype='ENSEMBL',columns=c('ENSEMBL','SYMBOL'))
ref=ref[!duplicated(ref[,1]),]
rownames(ref)=ref[,1]
meta = read.table('/Users/macbookpro/OneDrive/Documents/UDEM/Master/BulkData/meta_data.txt')

con_colors = c('deepskyblue','olivedrab2','darkorange')
names(con_colors)=c('untreated','astrocytes','inflammation')

grp_colors = c('tomato','mediumpurple1')
names(grp_colors) = c('H','M')

## limma DEA
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

fit <- cmdscale(dist(t(d$counts)),eig=TRUE, k=2)
meta$MDS1 = fit$points[,1]
meta$MDS2 = fit$points[,2]
mds.plot = ggplot(meta,aes(MDS1,MDS2,fill=condition,shape=group))+geom_point(size=4)+scale_fill_manual(values=con_colors)+scale_shape_manual(values=c(21,23))+theme_bw()

mm <- model.matrix(~0 + group:condition,data=meta)
E <- voom(d, mm, plot = T)
meta$condition = factor(meta$condition,levels=c('untreated','astrocytes','inflammation'))
#faire pca sur E$E normalisé
library(PCAtools)

groups = c(rep("H", 3), rep("HA", 3), rep("HI", 3),rep("M", 3), rep("MA", 3), rep("MI", 3))
cols = c('deepskyblue','olivedrab2','darkorange','deepskyblue','olivedrab2','darkorange')[factor(groups)]
cond = c("H","H","H","F","F","F")[factor(groups)]

#groups = c(rep("BBB_untreated", 3), rep("BBB_astrocyte", 3), rep("BBB_inflammation", 3),rep("BMB_untreated", 3), rep("BMB_astrocyte", 3), rep("BMB_inflammation", 3))

metadata <- data.frame(groups = groups,cond=cond)
rownames(metadata) <- colnames(E$E)
#rownames(metadata) <- paste0('BBB_', 1:nrow(metadata))
p <- pca(E$E, metadata = metadata)

plot1<-biplot(p, colby = 'groups', shape = 'groups'  ,hline = 0, vline = 0,legendPosition = 'right')
plot2<- plot1  +  scale_color_manual( name= "Groups",labels = c("H:BBB_untreated", "HA:BBB_astrocyte", "HI:BBB_inflammation","M:BMB_untreated","MA:BMB_astrocyte","MI:BMB_inflammation"), 
                                     values= c('deepskyblue','olivedrab2','darkorange','deepskyblue','olivedrab2','darkorange')) + 
  scale_shape_manual(name= "Groups", labels = c("H:BBB_untreated", "HA:BBB_astrocyte", "HI:BBB_inflammation","M:BMB_untreated","MA:BMB_astrocyte","MI:BMB_inflammation"), 
                     values= c(16,16,16,18,18,18) )+
  theme(legend.title = element_text(size=7), legend.text = element_text(size=6))
       #showLoadings = TRUE,labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
plot2

#plot3 <- biplot(p, colby = 'groups', showLoadings = TRUE,labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
#plot3


## Inflammation
# The chemokine (C-X-C motif) ligand (CXCL) family plays an important role in inflammation.
##ENSG00000138755 -> CXCL9
##ENSG00000169245 -> CXCL10
##ENSG00000169248 -> CXCL11
# https://pubmed.ncbi.nlm.nih.gov/26110930/
##ENSG00000131203 -> IDO1
# Ubiquitin D
##ENSG00000213886 -> UBD

## BMB_untreated & BMB_astrocyte
##ENSG00000145864 -> GABRB2
##ENSG00000165474 -> GJB2
##ENSG00000064655 -> EYA2
##ENSG00000147676 -> MAL2
##ENSG00000164764 -> SBSPON

##FIG : La stimulation des cytokines a induit des changements transcriptionnels significatifs dans les CE associées 
#à la BBB/BMB. Ces changements comprenaient des gènes tels que ENSG00000138755 (CXCL9), ENSG00000169245 (CXCL10), 
#ENSG00000169248 (CXCL11) et ENSG00000131203 (IDO1). 


fit <- cmdscale(dist(t(E$E)),eig=TRUE, k=2)
meta$voomMDS1 = fit$points[,1]
meta$voomMDS2 = fit$points[,2]
mds.plot.voom = ggplot(meta,aes(voomMDS1,voomMDS2,fill=condition,shape=group))+geom_point(size=4)+scale_fill_manual(values=con_colors)+scale_shape_manual(values=c(21,23))+theme_bw()

colnames(mm)=gsub(':','.',colnames(mm))
fit <- lmFit(E, mm)


## Comparaison entre les barrieres
contr <- makeContrasts(
            A = groupH.conditionUntreated - groupM.conditionUntreated,
            B = groupH.conditionACM - groupM.conditionACM,
            C = groupH.conditionInflammation - groupM.conditionInflammation, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

res.A <- topTable(tmp, sort.by = "P", n = Inf, coef='A')
res.B <- topTable(tmp, sort.by = "P", n = Inf, coef='B')
res.C <- topTable(tmp, sort.by = "P", n = Inf, coef='C')

res.A$condition = 'untreated'
res.B$condition = 'astrocytes'
res.C$condition = 'inflammation'

res.A$geneid = rownames(res.A)
res.B$geneid = rownames(res.B)
res.C$geneid = rownames(res.C)

res.bbb_vs_bmb = rbind(res.A,res.B,res.C)
res.bbb_vs_bmb$symbol = ref[res.bbb_vs_bmb$geneid,2]

venn.diagram(
    x=list(na.omit(subset(res.bbb_vs_bmb,condition=='untreated' & adj.P.Val<0.01 & abs(logFC)>1)$symbol), 
           na.omit(subset(res.bbb_vs_bmb,condition=='astrocytes' & adj.P.Val<0.01 & abs(logFC)>1)$symbol), 
           na.omit(subset(res.bbb_vs_bmb,condition=='inflammation' & adj.P.Val<0.01 & abs(logFC)>1)$symbol)),
    filename='/Users/macbookpro/OneDrive/Documents/UDEM/Master/BulkData/Venn.BBB_vs_BMB.png',
    category.names=c('untreated','astrocytes','inflammation'),
    fill=con_colors,
    col='black',
    output=TRUE
    )

degs.A = subset(res.bbb_vs_bmb, condition == 'untreated' & adj.P.Val<0.01 & abs(logFC)>1)$symbol
degs.B = subset(res.bbb_vs_bmb, condition == 'astrocytes' & adj.P.Val<0.01 & abs(logFC)>1)$symbol
degs.C = subset(res.bbb_vs_bmb, condition == 'inflammation' & adj.P.Val<0.01 & abs(logFC)>1)$symbol

res.bbb_vs_bmb$sign = 'NS'
res.bbb_vs_bmb[which(res.bbb_vs_bmb$adj.P.Val<0.01 & res.bbb_vs_bmb$logFC>1),'sign']='UP'
res.bbb_vs_bmb[which(res.bbb_vs_bmb$adj.P.Val<0.01 & res.bbb_vs_bmb$logFC<(-1)),'sign']='DN'

volcano.plot.bbb_vs_bmb <- ggplot(res.bbb_vs_bmb,aes(x=logFC,y=-log10(P.Value),col=sign))+geom_point()+theme_bw()+scale_color_manual(values=c('deepskyblue','grey75','tomato'))+facet_wrap(~condition) + ggtitle("BBB vs. BMB")
volcano.plot.bbb_vs_bmb


## GSEA geneset approach: method #1
inf.specific.degs = degs.C[!degs.C%in%degs.A & !degs.C%in%degs.B]
listEnrichrDbs()[,3]
enr.go.bp.sp = enrichr(inf.specific.degs, 'GO_Biological_Process_2021')

common.core = na.omit(degs.C[degs.C%in%degs.A & degs.C%in%degs.B])
enr.go.bp.cc = enrichr(as.character(common.core), 'GO_Biological_Process_2021')

## GSEA geneset approach: method #2
genes = list(common = common.core, specific = inf.specific.degs)
fish = fisher_enrichment(genes,go,unique(res.bbb_vs_bmb$symbol),'GO_BP')


### Comparaisons entre les conditions (groupH & M)
contr <- makeContrasts(
  A = groupH.conditionUntreated - groupH.conditionACM,
  B = groupH.conditionInflammation - groupH.conditionUntreated,
  C = groupH.conditionInflammation - groupH.conditionACM,
  D = groupM.conditionUntreated - groupM.conditionACM,
  E = groupM.conditionInflammation - groupM.conditionUntreated,
  F = groupM.conditionInflammation - groupM.conditionACM, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

res.A <- topTable(tmp, sort.by = "P", n = Inf, coef='A')
res.B <- topTable(tmp, sort.by = "P", n = Inf, coef='B')
res.C <- topTable(tmp, sort.by = "P", n = Inf, coef='C')
res.D <- topTable(tmp, sort.by = "P", n = Inf, coef='D')
res.E <- topTable(tmp, sort.by = "P", n = Inf, coef='E')
res.F <- topTable(tmp, sort.by = "P", n = Inf, coef='F')

res.A$condition = 'BBB-UvA'
res.B$condition = 'BBB-UvI'
res.C$condition = 'BBB-AvI'
res.D$condition = 'BMB-UvA'
res.E$condition = 'BMB-UvI'
res.F$condition = 'BMB-AvI'

res.A$geneid = rownames(res.A)
res.B$geneid = rownames(res.B)
res.C$geneid = rownames(res.C)
res.D$geneid = rownames(res.D)
res.E$geneid = rownames(res.E)
res.F$geneid = rownames(res.F)

res.conditionH = rbind(res.A,res.B,res.C,res.D,res.E,res.F)
res.A$symbol = ref[res.A$geneid,2]
res.B$symbol = ref[res.B$geneid,2]
res.C$symbol = ref[res.C$geneid,2]
res.D$symbol = ref[res.D$geneid,2]
res.E$symbol = ref[res.E$geneid,2]
res.F$symbol = ref[res.F$geneid,2]



venn.diagram(
    x=list(na.omit(subset(res.conditionH,condition=='UvA' & adj.P.Val<0.01 & abs(logFC)>2)$symbol), 
           na.omit(subset(res.conditionH,condition=='UvI' & adj.P.Val<0.01 & abs(logFC)>2)$symbol), 
           na.omit(subset(res.conditionH,condition=='AvI' & adj.P.Val<0.01 & abs(logFC)>2)$symbol)),
    filename='/Users/macbookpro/OneDrive/Documents/UDEM/Master/BulkData/Venn.ConditionH.png',
    category.names=c('UvA','UvI','AvI'),
    fill=con_colors,
    col='black')

res.conditionH$sign = 'NS'
res.conditionH[which(res.conditionH$adj.P.Val<0.01 & res.conditionH$logFC>1),'sign']='UP'
res.conditionH[which(res.conditionH$adj.P.Val<0.01 & res.conditionH$logFC<(-1)),'sign']='DN'

res.C$sign = 'NS'
res.C[which(res.C$adj.P.Val<0.01 & res.C$logFC>1),'sign']='UP'
res.C[which(res.C$adj.P.Val<0.01 & res.C$logFC<(-1)),'sign']='DN'
res.F$sign = 'NS'
res.F[which(res.F$adj.P.Val<0.01 & res.F$logFC>1),'sign']='UP'
res.F[which(res.F$adj.P.Val<0.01 & res.F$logFC<(-1)),'sign']='DN'


res.conditionH$symbol <- as.character(res.conditionH$symbol)


#res.C$genelabels <- ifelse(res.C$symbol %in% c("ICAM1", "VCAM1","CLMP", "JAM2","JAM3","F11R","CXADR"), as.character(res.C$symbol), NA)

#volcano.plot.conditionH <- ggplot(res.conditionH,aes(x=logFC,y=-log10(P.Value),label=genelabels,col=sign ))+geom_point()+theme_bw()+scale_color_manual(values=c('deepskyblue','grey75','darkorange')) + ggtitle("Astrocyte Conditioned Medium VS Inflammatory Condition") + facet_wrap(~condition)
#volcano.plot.conditionH


res.conditionH<-res.conditionH[order(res.conditionH$P.Value),]
res.C<-res.C[order(res.C$P.Value),]
res.F<-res.F[order(res.F$P.Value),]

LabelAll <- as_labeller(c(
  'BBB-UvA'="BBB-Unt Vs ACM",
  'BBB-UvI'="BBB-Unt Vs Inf",
  'BBB-AvI'="BBB-ACM Vs Inf",
  'BMB-UvA'="BMB-Unt Vs ACM",
  'BMB-UvI'="BMB-Unt Vs Inf",
  'BMB-AvI'="BMB-ACM Vs Inf"
))


res.conditionH$condition = factor(res.conditionH$condition, levels=c("BBB-UvA", "BBB-UvI", "BBB-AvI", "BMB-UvA", "BMB-UvI", "BMB-AvI"))
volcano.plot.all <- ggplot(res.conditionH,aes(x=logFC,y=-log10(P.Value)))+geom_point(aes(col=t))+theme_bw()+scale_color_gradient2(low='deepskyblue',mid='grey85',high='olivedrab2') +
  theme(legend.position = "none",strip.background = element_rect(fill=c("grey95")))+ facet_wrap(~condition, labeller=LabelAll)
volcano.plot.all


volcano.plot.H <- ggplot(res.conditionH,aes(x=logFC,y=-log10(P.Value)))+geom_point(aes(col=t))+theme_bw()+scale_color_gradient2(low='deepskyblue',mid='grey85',high='darkorange') +
  geom_point(data=res.A[1:5,], col="red") +geom_label_repel(data=res.A[1:5,],aes(label=symbol),force=8,fill = alpha(c("white"),0.5)) + 
  geom_point(data=res.B[1:5,], col="red") +geom_label_repel(data=res.B[1:5,],aes(label=symbol),force=8,fill = alpha(c("white"),0.5)) + 
  geom_point(data=res.C[1:5,], col="red") +geom_label_repel(data=res.C[1:5,],aes(label=symbol),force=8,fill = alpha(c("white"),0.5)) + 
  geom_point(data=res.D[1:5,], col="red") +geom_label_repel(data=res.D[1:5,],aes(label=symbol),force=50,fill = alpha(c("white"),0.5)) + 
  geom_point(data=res.E[1:5,], col="red") +geom_label_repel(data=res.E[1:5,],aes(label=symbol),force=50,fill = alpha(c("white"),0.5)) + 
  geom_point(data=res.F[1:5,], col="red") +geom_label_repel(data=res.F[1:5,],aes(label=symbol),force=50,fill = alpha(c("white"),0.5)) + 
  theme(legend.position = "none",strip.background = element_rect(fill=c("grey95")))+ facet_wrap(~condition , labeller=LabelAll)
volcano.plot.H

genelabels <- c("ICAM1", "VCAM1","CLMP", "JAM2","JAM3","F11R","CXADR")
volcano.plot.C <- ggplot(res.C,aes(x=logFC,y=-log10(P.Value)))+geom_point(aes(col=t))+theme_bw()+scale_color_gradient2(low='deepskyblue',mid='grey85',high='darkorange') + ggtitle("HBECs") +
  geom_point(data=res.C[match(genelabels, res.C$symbol),], col="red") +geom_label_repel(data=res.C[match(genelabels, res.C$symbol),],aes(label=genelabels),force=5) + theme(legend.position = "none") 
volcano.plot.C

volcano.plot.F <- ggplot(res.F,aes(x=logFC,y=-log10(P.Value)))+geom_point(aes(col=t))+theme_bw()+scale_color_gradient2(low='deepskyblue',mid='grey85',high='darkorange') + ggtitle("HMECs") +
  geom_point(data=res.F[match(genelabels, res.F$symbol),], col="red") +geom_label_repel(data=res.F[match(genelabels, res.F$symbol),],aes(label=genelabels),force=25)+ theme(legend.position = "none") 
volcano.plot.F


degs.A = subset(res.conditionH, condition == 'UvA' & adj.P.Val<0.01 & abs(logFC)>1)$symbol
degs.B = subset(res.conditionH, condition == 'UvI' & adj.P.Val<0.01 & abs(logFC)>1)$symbol
degs.C = subset(res.conditionH, condition == 'AvI' & adj.P.Val<0.01 & abs(logFC)>1)$symbol
degs.F = subset(res.conditionH, condition == 'AvI' & adj.P.Val<0.01 & abs(logFC)>1)$symbol

common = as.character(na.omit(degs.A[degs.A%in%degs.B & degs.A%in%degs.C]))
common.inf = as.character(na.omit(degs.B[!degs.B%in%degs.A & degs.B%in%degs.C]))
specificA.inf = as.character(na.omit(degs.C[!degs.C%in%degs.A & !degs.C%in%degs.B]))
specificU.inf = as.character(na.omit(degs.B[!degs.B%in%degs.A & !degs.B%in%degs.C]))

#degsC
conv = select(org.Hs.eg.db,keys=rownames(res.C),keytype='ENSEMBL',columns='SYMBOL')
conv=conv[!duplicated(conv[,1]),]
rownames(conv)=conv[,1]
res.C$symbol = conv[rownames(res.C),'SYMBOL']
degsC = subset(res.C,abs(logFC)>2&adj.P.Val <0.01)$symbol
#degsF
conv1 = select(org.Hs.eg.db,keys=rownames(res.F),keytype='ENSEMBL',columns='SYMBOL')
conv1=conv[!duplicated(conv[,1]),]
rownames(conv1)=conv1[,1]
res.F$symbol = conv1[rownames(res.F),'SYMBOL']
degsF = subset(res.F,abs(logFC)>2&adj.P.Val <0.01)$symbol
#
commonCF = as.character(na.omit(degsC[degsC%in%degsF]))
specificC = as.character(na.omit(degsC[!degsC%in%degsF]))
specificF = as.character(na.omit(degsF[!degsF%in%degsC]))

#
write.table(na.omit(commonCF), "commonCF.txt", quote= F, row.names= F, col.names = F )
write.table(na.omit(specificC), "specificC1.txt", quote= F, row.names= F, col.names = F )
write.table(na.omit(specificF), "specificF1.txt", quote= F, row.names= F, col.names = F )
write.table(na.omit(res.F$symbol), "universe.txt", quote= F, row.names= F, col.names = F )

#upset
degsUpset = list(
  BBB = rownames(na.omit(subset(res.C, condition=='BBB-AvI' &adj.P.Val<0.01& abs(logFC)>0.25))),
  BMB = rownames(na.omit(subset(res.F, condition=='BMB-AvI' &adj.P.Val<0.01& abs(logFC)>0.25)))
)
library(UpSetR)
upset(fromList(degsUpset),sets = c("BBB", "BMB"),order.by='freq',decreasing = c(F),  text.scale = c(1.8, 2,1,1.4,3,3))



venn.diagram(
  x=list(rownames(na.omit(subset(res.C,condition=='BBB' & adj.P.Val<0.01 & abs(logFC)>0.25))), 
         rownames(na.omit(subset(res.F,condition=='BMB' & adj.P.Val<0.01 & abs(logFC)>0.25)))),
  filename='/Users/macbookpro/OneDrive/Documents/UDEM/Master/BulkData/Venn.BBBvsBMB_IvA.png',
  category.names=c('BBB','BMB'),
  main = "Astrocyte Conditioned Medium VS Inflammatory Condition",
  main.pos = c(0.5, 0.95),
  fill=c('deepskyblue', 'olivedrab2'),
  col='black')



###### PanglaoDB_Augmented_2021
# Read in the data
x <- scan("/Users/macbookpro/OneDrive/Documents/UDEM/Master/BulkData/PanglaoDB_Augmented_2021.txt", what="", sep="\n")
# Separate elements by one or more whitepace
y <- strsplit(x, "[[:space:]]+")
# Extract the first vector element and set it as the list element name
names(y) <- sapply(y, `[[`, 1)
#names(y) <- sapply(y, function(x) x[[1]]) # same as above
# Remove the first vector element from each list element
goPangloaDB <- lapply(y, `[`, -1)
#y <- lapply(y, function(x) x[-1]) # same as above

##Gene set enrichment analysis
list = list(common=commonCF, specificC=specificC,specificF=specificF)
fishBP = fisher_enrichment(list,goBP,unique(res.C$symbol))
fishCC = fisher_enrichment(list,goCC,unique(res.C$symbol))
fishMF = fisher_enrichment(list,goMF,unique(res.C$symbol))
#fishPang = fisher_enrichment(list,goPangloaDB,unique(res.C$symbol))

##top5 pathway
topPath<- function(fish ,specific, p ,n){
  l=n*10
  result <- matrix(NA, l, l)
  
  colnames(result) <- rownames(result) <- rownames(getElement(fish, specific))[1:l]
  for(i in 1:l){
    for(j in 1:l){
      minlist = 0
      if(length(unlist(strsplit(fish[[specific]][["INTERSECT"]][i], split = ",")))>=length(unlist(strsplit(fish[[specific]][["INTERSECT"]][j], split = ",")))) minlist=length(unlist(strsplit(fish[[specific]][["INTERSECT"]][j], split = ","))) else minlist=length(unlist(strsplit(fish[[specific]][["INTERSECT"]][i], split = ",")))
      result[i, j] <- length(intersect(unlist(strsplit(fish[[specific]][["INTERSECT"]][i], split = ",")), unlist(strsplit(fish[[specific]][["INTERSECT"]][j], split = ","))))/(minlist)      
      result[j, i] <- result[i, j]     
    }
  }
  listBPspecific = list(c(rownames(result)[1]))
  i=1
  while(i < l){
    j=i
    while(j < l){
      if(result[i, j] < p){ 
        i=j
        listBPspecific <- c(listBPspecific, rownames(result)[j])
        j=l
      }else {
        j=j+1
      }
    }
    if(length(listBPspecific)>n-1) 
      i=l
  }
  return(listBPspecific)
}
##top5 pathway BP(1 14 17 21 22)
tabCBP <- fishBP$specificC 
tabCBP$fdr <- as.numeric(tabCBP$fdr) 
tabCBP$OR.odds.ratio <- as.numeric(tabCBP$OR.odds.ratio) 
tabCBP$NOM_pval <- as.numeric(tabCBP$NOM_pval) 
topCBP <- topPath(fishBP ,"specificC", 0.5, 5)
tabCBP <- tabCBP[match(topCBP, rownames(tabCBP)),]
rownames(tabCBP) <- gsub("GOBP_", "",rownames(tabCBP)) 

#goLen<-lapply(goBP,length)
#goLen<-data.frame(t(data.frame(goLen)))
#tabCBP$Length <-goLen$t.data.frame.goLen..[match(rownames(tabCBP), rownames(goLen))]
#tabCBP <- subset(tabCBP, tabCBP$Length<1000)
#tabCBP <- head(tabCBP[order(tabCBP$NOM_pval,  decreasing= F),], n = 5)

tabFBP <- fishBP$specificF
tabFBP$fdr <- as.numeric(tabFBP$fdr)
tabFBP$OR.odds.ratio <- as.numeric(tabFBP$OR.odds.ratio) 
tabFBP$NOM_pval <- as.numeric(tabFBP$NOM_pval)
topFBP <- topPath(fishBP ,"specificF", 0.5, 5)
tabFBP <- tabFBP[match(topFBP, rownames(tabFBP)),]
rownames(tabFBP) <- gsub("GOBP_", "",rownames(tabFBP)) 

##top5 pathway BP common
tabComBP <- fishBP$common
tabComBP$fdr <- as.numeric(tabComBP$fdr) 
tabComBP$OR.odds.ratio <- as.numeric(tabComBP$OR.odds.ratio) 
tabComBP$NOM_pval <- as.numeric(tabComBP$NOM_pval) 
topComBP <- topPath(fishBP ,"common", 0.7, 5)
tabComBP <- tabComBP[match(topComBP, rownames(tabComBP)),]

#tabFBP$Length <-goLen$t.data.frame.goLen..[match(rownames(tabFBP), rownames(goLen))]
#tabFBP <- subset(tabFBP, tabFBP$Length<1000)
#tabFBP <- head(tabFBP[order(tabFBP$NOM_pval,  decreasing= F),], n = 5)


##top5 pathway MF
tabCMF <- fishMF$specificC
tabCMF$fdr <- as.numeric(tabCMF$fdr)
tabCMF$OR.odds.ratio <- as.numeric(tabCMF$OR.odds.ratio) 
tabCMF$NOM_pval <- as.numeric(tabCMF$NOM_pval)
topCMF <- topPath(fishMF ,"specificC", 0.5, 5)
tabCMF <- tabCMF[match(topCMF, rownames(tabCMF)),]
rownames(tabCMF) <- gsub("GOMF_", "",rownames(tabCMF)) 

#goLen<-lapply(goMF,length)
#goLen<-data.frame(t(data.frame(goLen)))
#tabCMF$Length <-goLen$t.data.frame.goLen..[match(rownames(tabCMF), rownames(goLen))]
#tabCMF <- subset(tabCMF, tabCMF$Length<1000)
#tabCMF <- head(tabCMF[order(tabCMF$NOM_pval,  decreasing= F),], n = 5)

tabFMF <- fishMF$specificF
tabFMF$fdr <- as.numeric(tabFMF$fdr)
tabFMF$OR.odds.ratio <- as.numeric(tabFMF$OR.odds.ratio) 
tabFMF$NOM_pval <- as.numeric(tabFMF$NOM_pval)
topFMF <- topPath(fishMF ,"specificF", 0.5, 5)
tabFMF <- tabFMF[match(topFMF, rownames(tabFMF)),]
rownames(tabFMF) <- gsub("GOMF_", "",rownames(tabFMF)) 

#tabFMF$Length <-goLen$t.data.frame.goLen..[match(rownames(tabFMF), rownames(goLen))]
#tabFMF <- subset(tabFMF, tabFMF$Length<1000)
#tabFMF <- head(tabFMF[order(tabFMF$NOM_pval,  decreasing= F),], n = 5)
##top5 pathway MF
tabCCC <- fishCC$specificC
tabCCC$fdr <- as.numeric(tabCCC$fdr)
tabCCC$OR.odds.ratio <- as.numeric(tabCCC$OR.odds.ratio) 
tabCCC$NOM_pval <- as.numeric(tabCCC$NOM_pval)
topCCC <- topPath(fishCC ,"specificC", 0.5, 5)
tabCCC <- tabCCC[match(topCCC, rownames(tabCCC)),]
rownames(tabCCC) <- gsub("GOCC_", "",rownames(tabCCC)) 


#goLen<-lapply(goCC,length)
#goLen<-data.frame(t(data.frame(goLen)))
#tabCCC$Length <-goLen$t.data.frame.goLen..[match(rownames(tabCCC), rownames(goLen))]
#tabCCC <- subset(tabCCC, tabCCC$Length<1000)
#tabCCC <- head(tabCCC[order(tabCCC$NOM_pval,  decreasing= F),], n = 5)

tabFCC <- fishCC$specificF
tabFCC$fdr <- as.numeric(tabFCC$fdr)
tabFCC$OR.odds.ratio <- as.numeric(tabFCC$OR.odds.ratio) 
tabFCC$NOM_pval <- as.numeric(tabFCC$NOM_pval)
topFCC <- topPath(fishCC ,"specificF", 0.5, 5)
tabFCC <- tabFCC[match(topFCC, rownames(tabFCC)),]
rownames(tabFCC) <- gsub("GOCC_", "",rownames(tabFCC)) 
#tabFCC$Length <-goLen$t.data.frame.goLen..[match(rownames(tabFCC), rownames(goLen))]
#tabFCC <- subset(tabFCC, tabFCC$Length<1000)
#tabFCC <- head(tabFCC[order(tabFCC$NOM_pval,  decreasing= F),], n = 5)

CBP <- ggplot(tabCBP, aes(x=reorder(rownames(tabCBP),-rank(NOM_pval)), y=-log10(NOM_pval))) + 
  geom_col(width=.65,fill ="olivedrab3") +theme(axis.text=element_text(size=8)) + scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) +
  labs(title = "Top 5 pathways associated with the BBB", subtitle = "GO Biological Process 2021") + xlab("") + coord_flip() + theme_bw()

FBP <- ggplot(tabFBP, aes(x=reorder(rownames(tabFBP),-rank(NOM_pval)), y=-log10(NOM_pval))) +
  geom_col(width=.65,fill='deepskyblue3') +theme(axis.text=element_text(size=8)) + scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) +
  xlab("")+ labs(title = "Top 5 pathways associated with the BMB", subtitle = "GO Biological Process 2021") + coord_flip()+theme_bw()

CMF <- ggplot(tabCMF, aes(x=reorder(rownames(tabCMF),-rank(NOM_pval)), y=-log10(NOM_pval))) + 
  geom_col(width=.65,fill ="olivedrab3") +theme(axis.text=element_text(size=8)) + scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) +
  xlab("")+ labs( subtitle = "GO Molecular Function 2021") + coord_flip()+theme_bw()

FMF <- ggplot(tabFMF, aes(x=reorder(rownames(tabFMF),-rank(NOM_pval)), y=-log10(NOM_pval))) + 
  geom_col(width=.65,fill='deepskyblue3') +theme(axis.text=element_text(size=8)) +scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) +
  xlab("")+ labs( subtitle = "GO Molecular Function 2021") + coord_flip()+theme_bw()

CCC <- ggplot(tabCCC, aes(x=reorder(rownames(tabCCC),-rank(NOM_pval)), y=-log10(NOM_pval))) + 
  geom_bar(stat = "identity",width=.65) +theme(axis.text=element_text(size=8)) + scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) +
  xlab("")+ labs( subtitle = "GO Cellular Component 2021") + coord_flip()+theme_bw()

FCC <- ggplot(tabFCC, aes(x=reorder(rownames(tabFCC),-rank(NOM_pval)), y=-log10(NOM_pval))) + 
  geom_bar(stat = "identity",width=.65) +theme(axis.text=element_text(size=8)) + scale_x_discrete(label = function(x) stringr::str_wrap(gsub("_", " ",x), width = 30)) +
  xlab("")+ labs( subtitle = "GO Cellular Component 2021") + coord_flip()+theme_bw()

ggarrange(CBP, FBP, CMF, FMF, #CCC, FCC, 
          ncol = 2, nrow = 2,align = c("hv"))


##GOBP_RESPONSE_TO_CYTOKINE in TabCBP & TabFBP
listSpecC = strsplit(tabCBP["GOBP_RESPONSE_TO_CYTOKINE",]$INTERSECT,",")[[1]]
listSpecF = strsplit(tabFBP["GOBP_RESPONSE_TO_CYTOKINE",]$INTERSECT,",")[[1]]
listSpecC1 = select(org.Hs.eg.db,keys=listSpecC,keytype='SYMBOL',columns='ENSEMBL')
listSpecF1 = select(org.Hs.eg.db,keys=listSpecF,keytype='SYMBOL',columns='ENSEMBL')
my_gene = c(listSpecC1[,2],listSpecF1[,2])
#pheatmap(E$E[rownames(E$E)%in%my_gene,])
#pheatmap(E$E[rownames(E$E)%in%my_gene,],scale='row')
sub = E$E[rownames(E$E)%in%my_gene,]
rownames(sub)=res.C[rownames(sub),'symbol']
#pheatmap(sub,annotation_col=meta[,c('group','condition')])
mg=subset(mg,my_gene%in%rownames(E$E))
rownames(mg)=res.C[mg[,1],'symbol']
#mg[which(mg$my_gene%in%listSpecF1[,2]),]
mg$symbol = rownames(mg)
#mg[which(mg$my_gene%in%listSpecF1[,2]),]
mg$specific = 'C'
mg[which(mg$my_gene%in%listSpecF1[,2]),'specific']='F'
#pheatmap(sub,annotation_col=meta[,c('group','condition')],annotation_row=mg[,3])
gene_annot=mg[,3]
names(gene_annot)=rownames(mg)
#pheatmap(sub,annotation_col=meta[,c('group','condition')],annotation_row=gene_annot)
#pheatmap(sub,annotation_col=meta[,c('group','condition')],annotation_row=mg[,2:3])
pheatmap(sub,annotation_col=meta[,c('group','condition')],annotation_row=mg[,2:3],scale='row')


### Comparaisons entre les conditions(groupM)
contr <- makeContrasts(
  A = groupM.conditionUntreated - groupM.conditionACM,
  B = groupM.conditionInflammation - groupM.conditionUntreated,
  C = groupM.conditionInflammation - groupM.conditionACM, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

res.A <- topTable(tmp, sort.by = "P", n = Inf, coef='A')
res.B <- topTable(tmp, sort.by = "P", n = Inf, coef='B')
res.C <- topTable(tmp, sort.by = "P", n = Inf, coef='C')

res.A$condition = 'UvA'
res.B$condition = 'UvI'
res.C$condition = 'AvI'

res.A$geneid = rownames(res.A)
res.B$geneid = rownames(res.B)
res.C$geneid = rownames(res.C)

res.conditionM = rbind(res.A,res.B,res.C)
res.conditionM$symbol = ref[res.conditionM$geneid,2]

venn.diagram(
  x=list(na.omit(subset(res.conditionM,condition=='UvA' & adj.P.Val<0.01 & abs(logFC)>1)$symbol), 
         na.omit(subset(res.conditionM,condition=='UvI' & adj.P.Val<0.01 & abs(logFC)>1)$symbol), 
         na.omit(subset(res.conditionM,condition=='AvI' & adj.P.Val<0.01 & abs(logFC)>1)$symbol)),
  filename='/Users/macbookpro/OneDrive/Documents/UDEM/Master/BulkData/Venn.conditionM.png',
  category.names=c('UvA','UvI','AvI'),
  fill=con_colors,
  col='black')

res.conditionM$sign = 'NS'
res.conditionM[which(res.conditionM$adj.P.Val<0.01 & res.conditionM$logFC>1),'sign']='UP'
res.conditionM[which(res.conditionM$adj.P.Val<0.01 & res.conditionM$logFC<(-1)),'sign']='DN'

volcano.plot.conditionM <- ggplot(res.conditionM,aes(x=logFC,y=-log10(P.Value),col=sign))+geom_point()+theme_bw()+scale_color_manual(values=c('deepskyblue','grey75','tomato'))+facet_wrap(~condition) + ggtitle("BMB")
volcano.plot.conditionM

## GSEA geneset approach: method #2
li.enr = list(common=common, common.inf=common.inf, specificA.inf=specificA.inf,specificU.inf=specificU.inf)
fish = fisher_enrichment(li.enr,go,unique(res.bbb_vs_bmb$symbol),'GO_BP')

## GSEA ranked approach: fastGSEA
res.B.noNA = subset(res.conditionH,!is.na(symbol) & condition == 'UvI')
tstat.B = res.B.noNA$t
names(tstat.B) = ref[as.character(res.B.noNA$geneid),2]

fgseaRes <- fgsea(pathways = go, 
                  stats    = tstat.B,
                  nperm = 2000)
head(fgseaRes[order(fgseaRes$padj),])
plotEnrichment(go[['GO_DNA_DEPENDENT_DNA_REPLICATION_MAINTENANCE_OF_FIDELITY']],tstat.B)



conv = select(org.Hs.eg.db,keys=rownames(res.C),keytype='ENSEMBL',columns='SYMBOL')
conv=conv[!duplicated(conv[,1]),]
rownames(conv)=conv[,1]
res.C$symbol = conv[rownames(res.C),'SYMBOL']
degs = subset(res.C,abs(logFC)>0.25&adj.P.Val <0.01)$symbol

## DEA avec DESeq2 --> http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
## Make R notebook https://bookdown.org/yihui/rmarkdown/notebook.html





#pca summarize all experiments
#BMB vs BBB
#volcano plot
#plot show enrichr associated with Down and upregulated genes
#boxplot of gene that of interresting to us 

##reverse regulation
test <- data.frame(E$E)
test1 <- rownames(subset(test,HI1-H1<(-1) & MI1-M1>1))
test2 <- ref[match(test1, ref$ENSEMBL),]

library(PCAtools)

groups = c(rep("H", 3), rep("HA", 3), rep("HI", 3),rep("M", 3), rep("MA", 3), rep("MI", 3))
cols = c('deepskyblue','olivedrab2','darkorange','deepskyblue','olivedrab2','darkorange')[factor(groups)]
cond = c("H","H","H","F","F","F")[factor(groups)]

#groups = c(rep("BBB_untreated", 3), rep("BBB_astrocyte", 3), rep("BBB_inflammation", 3),rep("BMB_untreated", 3), rep("BMB_astrocyte", 3), rep("BMB_inflammation", 3))

metadata <- data.frame(groups = groups,cond=cond)
rownames(metadata) <- colnames(E$E)
#rownames(metadata) <- paste0('BBB_', 1:nrow(metadata))
p <- pca(E$E, metadata = metadata)

plot1<-biplot(p, colby = 'groups', shape = 'cond'  ,hline = 0, vline = 0,legendPosition = 'right')
plot2<- plot1  +  scale_color_manual( name= "Condition",labels = c("Untreated", "ACM", "Inflammation","Untreated", "ACM", "Inflammation"), 
                                      values= c('deepskyblue','olivedrab2','darkorange','deepskyblue','olivedrab2','darkorange')) + 
  scale_shape_manual(name= "Groups", labels = c("BBB", "BMB"), 
                     values= c(16,18) )+
  theme(legend.title = element_text(size=10), legend.text = element_text(size=8,face="bold") )
#showLoadings = TRUE,labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
plot2


test <- data.frame(E$E)
test1 <- rownames(subset(test,HI1-H1<(-1) & MI1-M1>1))
test1 <- na.omit(ref[match(test1, ref$ENSEMBL),])

test2 <- rownames(subset(test,HI1-H1>1 & MI1-M1<(-1)))
test2 <- na.omit(ref[match(test2, ref$ENSEMBL),])

test3 <- rbind(test1,test2)


for (i in 1:nrow(test3)) {
  nam <- paste("a", i, sep = "")
  assign(nam, plot.gene(test3$SYMBOL[i],cols=con_colors))
}

liste="a1"
for(i in 2:nrow(test3)){
  nam <- paste("a", i, sep = "")
  liste=paste(liste, nam, sep = ",")
}
liste
plot_list = list(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48,a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63,a64,a65,a66,a67,a68,a69,a70,a71,a72,a73,a74,a75,a76,a77,a78,a79,a80,a81,a82,a83,a84,a85,a86) 
ggarrange(plotlist=plot_list, widths = c(5,18), align = c("hv"))

