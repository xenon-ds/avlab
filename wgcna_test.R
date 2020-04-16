#scr 0 vs 2
library(stringr)
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/longitudinal_igra_negative/counts_igra_negative_scr_0hr_vs_igra_negative_scr_2hr.csv",header = TRUE,stringsAsFactors = FALSE)


#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

columns=colnames(counts_df)
columns[str_detect(columns,pattern = "SCO149")]="SCO149"
columns[str_detect(columns,pattern = "SCO68")]="SCO68"
columns[str_detect(columns,pattern = "SCO71")]="SCO71"
columns[1:3]=paste("scr_0hr_",columns[1:3],sep="")
columns[4:6]=paste("scr_2hr_",columns[4:6],sep="")
colnames(counts_df)=columns

head(counts_df)
counts_df = median_based_count_matrix_filtering(counts_df, expected_median = 10)

transposed_counts_df = t(counts_df)
head(transposed_counts_df)
colnames(transposed_counts_df)
rownames(transposed_counts_df)

library("WGCNA")
sampleTree = hclust(dist(transposed_counts_df), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches# The user should change the dimensions if the window is too large or too small.sizeGrWindow(12,9)#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(transposed_counts_df, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
