#scr 0 vs 2
library(stringr)
counts_df=read.csv(file="./counts_igra_negative_scr_0hr_vs_igra_negative_b3_0hr.csv",header = TRUE,stringsAsFactors = FALSE)


#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

columns=colnames(counts_df)
columns[str_detect(columns,pattern = "SCO149")]="SCO149"
columns[str_detect(columns,pattern = "SCO68")]="SCO68"
columns[str_detect(columns,pattern = "SCO71")]="SCO71"
columns[1:3]=paste("scr_0hr_",columns[1:3],sep="")
columns[4:6]=paste("b3_0hr_",columns[4:6],sep="")
colnames(counts_df)=columns

head(counts_df)
counts_df = median_based_count_matrix_filtering(counts_df, expected_median = 10)

transposed_counts_df = t(counts_df)
head(transposed_counts_df)
colnames(transposed_counts_df)
rownames(transposed_counts_df)

library("WGCNA")
sampleTree = hclust(dist(transposed_counts_df), method = "average")
# Plot the sample tree: Open a graphic ouwindow of size 12 by 9 inches# The user should change the dimensions if the window is too large or too small.sizeGrWindow(12,9)#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=70, by=2))
sft = pickSoftThreshold(transposed_counts_df, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 50
adjacency_matrix = adjacency(transposed_counts_df, power = softPower)

TOM = TOMsimilarity(adjacency_matrix)
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM),method = "average")
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04)

minimum_module_size = 30 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minimum_module_size)
print(table(dynamicMods))

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)


sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(transposed_counts_df, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")


MEDissThres = 0.125
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(transposed_counts_df, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

# Define numbers of genes and samples
nGenes = ncol(transposed_counts_df)
nSamples = nrow(transposed_counts_df)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(transposed_counts_df, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
datTraits = read.csv("./traits_data_cross_sectional.csv", stringsAsFactors = FALSE)
rownames(datTraits) = datTraits[,1]
datTraits$sample_id = NULL
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)




#finding out module membership
sample_type = as.data.frame(datTraits %>% dplyr::select(type))

#name colours of the modules 
modNames = substring(names(MEs),3)
#calculating gene module membership 
geneModuleMembership = as.data.frame(cor(transposed_counts_df,MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples))
names(geneModuleMembership) = paste("MM",modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

geneTraitSignificance = as.data.frame(cor(transposed_counts_df,sample_type, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(sample_type), sep="");
names(GSPvalue) = paste("p.GS.", names(sample_type), sep="");

for(module in unique(moduleColors))
{
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
#sizeGrWindow(7, 7);
#par(mfrow = c(1,1));
  jpeg(paste("./generated_figures/module_",module,".jpeg", sep = ""))
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]),xlab = paste("Module Membership in", module, "module"),ylab = "Gene significance for sample type (SCR/B3)",main = paste("Module membership vs. gene significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
}


# par(mar = c(6, 8.5, 3, 3))
# # Display the correlation values within a heatmap plot
# labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(datTraits),yLabels = names(MEs),ySymbols = names(MEs),textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,zlim = c(-1,1),main = paste("Module-trait relationships"), colorMatrix = )
# 
# res = labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(datTraits),yLabels = names(MEs),ySymbols = names(MEs),colorLabels = FALSE,colors = greenWhiteRed(50),textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,zlim = c(-1,1),main = paste("Module-trait relationships"))
# 
# 
# 
# naColor = NULL
# if (length(naColor)==0) naColor = 0;
# signed = FALSE
# colors = greenWhiteRed(50)
# zlim = c(-1,1)
# if (signed) zlim = c(-max(abs(zlim)), max(abs(zlim)));
# colorMatrix = numbers2colors(moduleTraitCor, signed, colors = greenWhiteRed(50), lim = zlim, naColor = naColor)
