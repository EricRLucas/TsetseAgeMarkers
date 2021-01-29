library(stringr)
library(edgeR)
library(plotrix)
library(fdrtool)
library(WGCNA)
library(grid)

# As a reminder, here are what the different columns contain in the counts data:
#column 1: gene ID
#column 2: counts for unstranded RNA-seq
#column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
#column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
# The library used for these data ("Ultra Directional RNA library preparation kits") requires the "-s reverse"
# option

# First, load up all of the data tables that were produced as part of the STAR alignment process
for (i in 1:50){
	padi <- str_pad(i, 2, pad = '0')
	filepath <- paste('read_counts/ReadsPerGene_S', padi, '.out.tab', sep = '')
	these.counts <- read.table(filepath, row.names = 1, sep = '\t')
	if (i == 1){
		read.counts.all <- these.counts[,3, drop = F]
		colnames(read.counts.all) <- 'S01'
	}
	else{
		# Do a check on row names
		if (!all(rownames(these.counts) == rownames(read.counts.all)))
			stop(paste('Rownames for sample S', padi, ' do not match previous samples.\n', sep = ''))
		read.counts.all[[paste('S', padi, sep = '')]] <- these.counts[,3]
	}
}

# We remove the top four rows (which are global values) as well as genes for which expression is too low:
read.counts <- read.counts.all[-(1:4),]
read.counts <- read.counts[apply(read.counts, 1, sum) > 10, ]

# Now get the table that relates sequencing sample codes to the original sample codes
sample.conversion <- read.table('tables/sample_conversion.csv', sep = '\t', row.names = 1, stringsAsFactors = F)
colnames(sample.conversion) <- 'original.name'

# Now load up the metadata
meta <- read.table('tables/sample_metadata.csv', sep = '\t', header = T, row.names = 1, quote = '', comment.char = '#')
# Narrow it down to just the samples used for RNAseq. 
sample.conversion <- cbind(sample.conversion, meta[sample.conversion$original.name,])

# Check that the sample order is the same in the sample conversion table and the counts data
if (!all(rownames(sample.conversion) == colnames(read.counts)))
	stop('Order of samples in read counts and meta data should be the same.\n')

age <- sample.conversion$age
sex <- sample.conversion$sex
days.since.feed <- sample.conversion$days.since.feed

# Create the design matrix
design.matrix <- model.matrix(~age + sex + days.since.feed)
read.counts.dgelist <- DGEList(counts = read.counts)
read.counts.dgelist <- calcNormFactors(read.counts.dgelist)
read.counts.dgelist <- estimateGLMTrendedDisp(read.counts.dgelist, design.matrix)
read.counts.dgelist <- estimateGLMTagwiseDisp(read.counts.dgelist, design.matrix)
read.counts.glm <- glmFit(read.counts.dgelist, design.matrix)

# GLM test of age
age.test <- glmLRT(read.counts.glm, coef = 'age')
age.ps <- age.test$table$PValue
names(age.ps) <- rownames(age.test$table)
age.fdr <- fdrtool(age.ps, statistic='pvalue')$qval
ranked.diff.genes.age <- names(sort(age.fdr))[sort(age.fdr) < 0.01]

# GLM test of sex
sex.test <- glmLRT(read.counts.glm, coef = 'sexmale')
sex.ps <- sex.test$table$PValue
names(sex.ps) <- rownames(sex.test$table)
sex.fdr <- fdrtool(sex.ps, statistic='pvalue')$qval
ranked.diff.genes.sex <- names(sort(sex.fdr))[sort(sex.fdr) < 0.01]

# GLM test of days since blood meal
dsf.test <- glmLRT(read.counts.glm, coef = 'days.since.feed')
dsf.ps <- dsf.test$table$PValue
names(dsf.ps) <- rownames(dsf.test$table)
dsf.fdr <- fdrtool(dsf.ps, statistic='pvalue')$qval
ranked.diff.genes.dsf <- names(sort(dsf.fdr))[sort(dsf.fdr) < 0.01]

# Get the normalised read counts
normalised.totals <- read.counts.dgelist$counts
for (n in colnames(read.counts.dgelist$counts)){
	normalised.library.size <- read.counts.dgelist$samples[n, "lib.size"] * read.counts.dgelist$samples[n, "norm.factors"]
	normalised.totals[,n] <- 1000000 * read.counts.dgelist$counts[,n] / normalised.library.size
}

# Save the metadata table
output.table <- sample.conversion[, c('original.name', 'collection.date', 'age', 'sex', 'days.since.feed', 'box.uses')]
output.table <- output.table[order(as.numeric(sub('^T', '', output.table$original.name))), ]
colnames(output.table) <- c('sample.name', 'date', 'age', 'sex', 'days.since.feed', 'box.uses')
write.table(output.table, file = 'tables/RNAseq_sample_metadata.csv', sep = '\t', col.names = T, row.names = F)

# Set some plotting colours
age.colours <- c(young = 'yellow', old = 'brown')
agecolscale <- color.scale(1:61, extremes = age.colours)
sample.conversion$agecol <- agecolscale[sample.conversion$age-1]
sex.colours <- c(female = 'royalblue3', male = 'orangered')
sample.conversion$sexcol <- c(sex.colours)[as.numeric(sample.conversion$sex)]
sex.pch <- c(female = 24, male = 21)
sample.conversion$sexpch <- sex.pch[as.numeric(sample.conversion$sex)]
sexpch.noborder <- c(female = 17, male = 19)[as.numeric(sample.conversion$sex)]
dsf.colours <- c('1' = 'limegreen','2' = 'purple3', '3' = 'indianred3')
sample.conversion$dsfcol <- dsf.colours[as.numeric(sample.conversion$days.since.feed)]

# Specifically look for genes that change in old age
read.counts.old <- read.counts[, sample.conversion$age > 30]
sample.conversion.old <- subset(sample.conversion, age > 30)

age.old <- sample.conversion.old$age
sex.old <- sample.conversion.old$sex
days.since.feed.old <- sample.conversion.old$days.since.feed

# Create the design matrix
design.matrix.old <- model.matrix(~age.old + sex.old + days.since.feed.old)
read.counts.dgelist.old <- DGEList(counts = read.counts.old)
read.counts.dgelist.old <- calcNormFactors(read.counts.dgelist.old)
read.counts.dgelist.old <- estimateGLMTrendedDisp(read.counts.dgelist.old, design.matrix.old)
read.counts.dgelist.old <- estimateGLMTagwiseDisp(read.counts.dgelist.old, design.matrix.old)
read.counts.glm.old <- glmFit(read.counts.dgelist.old, design.matrix.old)

# GLM test of age
age.old.test <- glmLRT(read.counts.glm.old, coef = 'age.old')
age.old.ps <- age.old.test$table$PValue
names(age.old.ps) <- rownames(age.old.test$table)
age.old.fdr <- fdrtool(age.old.ps, statistic='pvalue')$qval
ranked.genes.age.old <- names(sort(age.old.fdr))

# GLM test of sex
sex.old.test <- glmLRT(read.counts.glm.old, coef = 'sex.oldmale')
sex.old.ps <- sex.old.test$table$PValue
names(sex.old.ps) <- rownames(sex.old.test$table)
sex.old.fdr <- fdrtool(sex.old.ps, statistic='pvalue')$qval
ranked.genes.sex.old <- names(sort(sex.old.fdr))

svg('diffexp_age_bysex.svg')
par(mfrow = c(3,3), mgp = c(1.5,0.5,0), mar = c(3,3,2,1))
for (i in 1:9){
	y <- normalised.totals[ranked.diff.genes.age[i], ]
	plot(age, y, main = i, col = 'grey70', bg = sample.conversion$sexcol, ylab = 'Normalised read counts', pch = sample.conversion$sexpch, lwd = 0.6, ylim = c(0,max(y)), cex = 1.2)
	# Add a legend
	if (i == 1){
		points(c(40,40), c(0.7, 0.8)*max(y), col = sex.colours, bg = sex.colours, pch = sex.pch)
		text(c(42,42), c(0.7, 0.8)*max(y), label = names(sex.colours), adj = 0)
	}
}
dev.off()

# There aren't really any genes strongly differentiated by age.
svg('diffexp_oldage_bysex.svg')
par(mfrow = c(3,3), mgp = c(1.5,0.5,0), mar = c(3,3,2,1))
for (i in 1:9){
	y <- normalised.totals[ranked.genes.age.old[i], ]
	plot(age, y, main = i, col = 'grey70', bg = sample.conversion$sexcol, ylab = 'Normalised read counts', ylim = c(0, max(y)), pch = sample.conversion$sexpch, lwd = 0.6, cex = 1.2)
	if (i == 9){
		points(c(40,40), c(0.8, 0.9)*max(y), col = sex.colours, pch = c(17, 19))
		text(c(42,42), c(0.8, 0.9)*max(y), label = names(sex.colours), adj = 0)
	}
}
dev.off()

# Let's look for the genes with the lowest BCV (these should be good housekeeping genes). 
ordered.dispersion <- order(read.counts.glm$dispersion)
ranked.dispersion <- rownames(read.counts)[ordered.dispersion]

svg('housekeeping.svg')
par(mfrow = c(3,3), mgp = c(1.5,0.5,0), mar = c(3,3,2,1))
for (i in 1:9){
	y <- normalised.totals[ranked.dispersion[i], ]
	plot(age, y, main = i, col = 'grey70', bg = sample.conversion$sexcol, ylab = 'Normalised read counts', ylim = c(0, max(y)), pch = sample.conversion$sexpch, lwd = 0.6, cex = 1.2)
	if (i == 1){
		points(c(40,40), c(0.1, 0.2)*max(y), col = sex.colours, pch = c(17, 19))
		text(c(42,42), c(0.1, 0.2)*max(y), label = names(sex.colours), adj = 0)
	}
}
dev.off()


##################
# WGCNA analysis #
##################

# Set the number of top genes that we want to keep
ntopgenes = 5000

# Select the top genes based on variance 
s <- rowMeans((normalised.totals - rowMeans(normalised.totals))^2)
genes.to.keep <- names(sort(s, decreasing = T))[1:ntopgenes]
normalised.totals.filtered <- normalised.totals[genes.to.keep, ]

# We take logs to normalised the variance, adding 1 before taking logs in order to avoid log(0)
log.normalised.totals.filtered <- log2(normalised.totals.filtered + 1)
# We turn the matrix around to be in the required format for the package
wgcna.exp.matrix <- t(log.normalised.totals.filtered)
nGenes <- ncol(wgcna.exp.matrix)
nSamples <- nrow(wgcna.exp.matrix)

# We do a quick clustering of samples out of interest. The distance is calculated as the euclidian distance of the
# samples to each other if they were located in an n-dimensional space where each dimension is a gene. The function
# "dist" outputs an object of class "dist", which is in essence a distance matrix, but instead of being a matrix, 
# it is a form of list, where each element is one of the distances.
sample.tree <- hclust(dist(wgcna.exp.matrix), method='average')
informative.sample.names <- paste(sample.conversion[sample.tree$labels, 'sex'], sample.conversion[sample.tree$labels, 'age'], sample.conversion[sample.tree$labels, 'days.since.feed'], sep = '.')
# We change the sample names so that the clustering tree has informative labels.
sample.tree$labels <- informative.sample.names

# Choose a value for the soft-thresholding power parameter
powers = 1:15
# we make some space by deleting garbage in the memory 
collectGarbage()
# Call the network topology analysis function. 
sft = pickSoftThreshold(wgcna.exp.matrix, powerVector = powers, verbose = 5)
# We clear space in the memory 
collectGarbage()

# Soft-thresholding as a function of soft-thresholding power
svg('Soft_thresholding.svg')
par(mfrow=c(1,2))
plot(sft$fitIndices$Power, -sign(sft$fitIndices$slope)*sft$fitIndices$SFT.R.sq, xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices$Power, -sign(sft$fitIndices$slope)*sft$fitIndices$SFT.R.sq, labels=powers, cex=0.9, col="red")
abline(h=0.90,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices$Power, sft$fitIndices$mean.k.,
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices$Power, sft$fitIndices$mean.k., labels=powers, cex=0.9, col="red")
dev.off()

# We choose a value of 8, which is where R^2 starts to plateau.
power.parameter <- 8
adjacency.matrix <- adjacency(wgcna.exp.matrix, power = power.parameter)
collectGarbage()

# We create the topological overlap matrix (TOM) which gives us our measures of similarity. Then we just 
# take 1- this to give us our dissimilarity
TOM = TOMsimilarity(adjacency.matrix)
dissTOM = 1-TOM
# Call the hierarchical clustering function
gene.tree = hclust(as.dist(dissTOM), method = "average")
# and plot the dendrogram
sizeGrWindow(12,9)
plot (gene.tree, xlab= "" , sub = "" , main = "Gene clustering on TOM-based dissimilarity" , labels = FALSE, hang = 0.04)

# Next, we build our clusters using this dendrogram
# We set a minimum module size 
minModuleSize = 30
# Module identification using dynamic tree cut, using the hybrid_4 method. 
dynamicMods_hybrid_4 = cutreeDynamic(dendro = gene.tree, distM = dissTOM, deepSplit = 4, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
# The following displays the number of genes assigned to each cluster. Colour "0" is for genes that were 
# not assigned a cluster
table(dynamicMods_hybrid_4)

# Assign a colour to each gene based on its module
dynamicColors_hybrid_4 = labels2colors(dynamicMods_hybrid_4) 
table(dynamicColors_hybrid_4)

# Clear some memory
collectGarbage()

# Calculate eigengenes (according to the "presentation.pdf" file by the authors of the package, an eigengene is 
# the first principle component of a module's expression matrix). These are the coordinates of each sample in 
# eigengene-space (where each dimension represents the eigengene for one of the modules). 
MEList_hybrid_4 = moduleEigengenes(wgcna.exp.matrix, colors = dynamicColors_hybrid_4)
MEs_hybrid_4 = MEList_hybrid_4$eigengenes
# Get the dissimilariy of the eigengenes. We take the absolute value of the correlation coefficient, rather than 
# its actual value. 
MEDiss_hybrid_4 = 1-abs(cor(MEs_hybrid_4))
# Build a dendrogram of module eigengenes. Here we are basically repeating the clustering analysis above, but
# applying it to the eigengenes instead of the whole dataset. 
METree_hybrid_4 = hclust(as.dist(MEDiss_hybrid_4), method = "average")
# Plot the result
X11()
sizeGrWindow(7, 6)
plot(METree_hybrid_4, main = "Clustering of hybrid (4) module eigengenes", xlab = "", sub = "")
abline(h = 0.2)

# Plot the eigengene values of each sample for each module
number.of.genes <- numeric()
for (module in colnames(MEs_hybrid_4))
	number.of.genes[module] <- sum(paste('ME', dynamicColors_hybrid_4, sep='') == module)
x11()
par(mfrow=c(6,7), mar=c(2,4,2,2))
for (module in colnames(MEs_hybrid_4))
	plot(age, MEs_hybrid_4[, module], ylab = module, main = number.of.genes[module], xlab = 'age', col = sample.conversion$sexcol, pch = 19)

# Even at a low cut point (0.2), some modules get merged that look quite similar, but others get merged despite
# looking quite different (eg: grey60 & tan or lightcyan & yellow). Since our main purpose here is to avoid choosing
# too many genes from the same module, it is better to do too much rather than too little merging. So let's merge 
# those.

MEDissThresHybrid_4 = 0.2
# We use an automatic merging function 
merge_hybrid_4 = mergeCloseModules(wgcna.exp.matrix, dynamicColors_hybrid_4, cutHeight = MEDissThresHybrid_4, useAbs = T, verbose = 3)
mergedColors_hybrid_4 = merge_hybrid_4$colors
names(mergedColors_hybrid_4) <- colnames(wgcna.exp.matrix)
mergedMEs_hybrid_4 = merge_hybrid_4$newMEs
rownames(mergedMEs_hybrid_4) <- rownames(wgcna.exp.matrix)

# Plot the new eigengene values 
svg('WGCNA_module_plots.svg', width = 15, height = 9)
number.of.genes.merged <- numeric()
for (module in colnames(mergedMEs_hybrid_4))
	number.of.genes.merged[module] <- sum(paste('ME', mergedColors_hybrid_4, sep='') == module)
par(mfrow=c(5,6), mar=c(2,4,2,2))
for (module in colnames(mergedMEs_hybrid_4))
	plot(age, mergedMEs_hybrid_4[, module], ylab = module, main = number.of.genes.merged[module], xlab = 'age', col = sample.conversion$sexcol, pch = 19)
dev.off()

# Nearly all of the top 30 age-related genes are in the turquoise cluster:
unique(mergedColors_hybrid_4[ranked.diff.genes.age[1:30]])
# GMOY008370 is from a different cluster and keeps declining in later ages (albeit only very slightly).

# The interesting candidates for genes that still change in older ages (GMOY007669, GMOY000919, GMOY011008, 
# GMOY001238, GMOY009985, GMOY001603, GMOY007469 are mostly different modules.

# Then we want two housekeeping genes. 
x11()
par(mfrow = c(5,6), mgp = c(1.5,0.5,0), mar = c(3,3,2,1))
for (i in 1:30){
	y <- normalised.totals[ranked.dispersion[i], ]
	plot(age, y, main = ranked.dispersion[i], col = sample.conversion$sexcol, ylab = 'Normalised read counts', pch = 19, ylim = c(0, max(y)))
}
# Unsurprisingly, most of housekeeping candidates don't have clusters assigned to them because by definition 
# they have very low variability so they got taken out of the WGCNA analysis.

write.table(normalised.totals, file = 'tables/normalised_totals.csv', sep = '\t', col.names = NA)

save.image('RNAseq_analysis.Rdata')
