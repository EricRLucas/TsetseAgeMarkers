library(stringr)
library(edgeR)
library(plotrix)
library(fdrtool)

# As a reminder, here are what the different columns contain in the counts data:
#column 1: gene ID
#column 2: counts for unstranded RNA-seq
#column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
#column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
# The library used for these data ("Ultra Directional RNA library preparation kits") requires the "-s reverse"
# option.

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
# Narrow it down to just the samples we have here
sample.conversion <- cbind(sample.conversion, meta[sample.conversion$original.name,])

# Check that the sample order is the same in the sample conversion table and the counts data
if (!all(rownames(sample.conversion) == colnames(read.counts)))
	stop('Order of samples in read counts and meta data should be the same.\n')

age <- sample.conversion$age
sex <- sample.conversion$sex
box.uses <- sample.conversion$box.uses
days.since.feed <- sample.conversion$days.since.feed

# Create the design matrix
design.matrix <- model.matrix(~age + sex + box.uses + days.since.feed)
read.counts.dgelist <- DGEList(counts = read.counts)
read.counts.dgelist <- calcNormFactors(read.counts.dgelist)
read.counts.dgelist <- estimateGLMTrendedDisp(read.counts.dgelist, design.matrix)
read.counts.dgelist <- estimateGLMTagwiseDisp(read.counts.dgelist, design.matrix)
read.counts.glm <- glmFit(read.counts.dgelist, design.matrix)

# GLM testing of age
age.test <- glmLRT(read.counts.glm, coef = 'age')
age.ps <- age.test$table$PValue
names(age.ps) <- rownames(age.test$table)
age.fdr <- fdrtool(age.ps, statistic='pvalue')$qval
ranked.diff.genes.age <- names(sort(age.fdr))[sort(age.fdr) < 0.01]
cat('AGE:\nOut of ', length(age.fdr), ', ', sum(age.fdr < 0.01), ' and ', sum(age.fdr < 0.05), ' genes were significantly associated with age after FDR control at 0.01 and 0.05 respectively.\n\n', sep = '')

# GLM testing of sex
sex.test <- glmLRT(read.counts.glm, coef = 'sexmale')
sex.ps <- sex.test$table$PValue
names(sex.ps) <- rownames(sex.test$table)
sex.fdr <- fdrtool(sex.ps, statistic='pvalue')$qval
ranked.diff.genes.sex <- names(sort(sex.fdr))[sort(sex.fdr) < 0.01]
cat('SEX:\nOut of ', length(sex.fdr), ', ', sum(sex.fdr < 0.01), ' and ', sum(sex.fdr < 0.05), ' genes were significantly associated with sex after FDR control at 0.01 and 0.05 respectively.\n\n', sep = '')
 
# GLM testing of says since blood meal
dsf.test <- glmLRT(read.counts.glm, coef = 'days.since.feed')
dsf.ps <- dsf.test$table$PValue
names(dsf.ps) <- rownames(dsf.test$table)
dsf.fdr <- fdrtool(dsf.ps, statistic='pvalue')$qval
ranked.diff.genes.dsf <- names(sort(dsf.fdr))[sort(dsf.fdr) < 0.01]
cat('DAYS SINCE FEEDING:\nOut of ', length(dsf.fdr), ', ', sum(dsf.fdr < 0.01), ' and ', sum(dsf.fdr < 0.05), ' genes were significantly associated with days since feeding after FDR control at 0.01 and 0.05 respectively.\n\n', sep = '')

# GLM testing of number of previous box uses
uses.test <- glmLRT(read.counts.glm, coef = 'box.uses')
uses.ps <- uses.test$table$PValue
names(uses.ps) <- rownames(uses.test$table)
uses.fdr <- fdrtool(uses.ps, statistic='pvalue')$qval
ranked.genes.uses <- names(sort(uses.fdr))
cat('PREVIOUS BOX USES:\nOut of ', length(uses.fdr), ', ', sum(uses.fdr < 0.01), ' and ', sum(uses.fdr < 0.05), ' genes were significantly associated with number of previous box uses after FDR control at 0.01 and 0.05 respectively.\n\n', sep = '')

# Plot some P-value histograms
pdf('P-hist_tests.pdf', height = 5)
par(mfrow = c(2,2), mgp = c(2,0.5,0), mar = c(2, 1, 2, 1), oma = c(1,2,0.5,0), xpd = NA)
hist(age.ps, ylim = c(0,4000), xaxt = 'n', yaxt = 'n', xlab = '', main = '', col = 'tan2')
mtext('Age', line = 0, font = 2)
axis(2, c(0, 1000, 2000, 3000, 4000), labels = c('0', '', '2000', '', '4000'))
#
hist(sex.ps, ylim = c(0,4000), xaxt = 'n', xlab = '', yaxt = 'n', ylab = '', main = '', col = 'tan2')
mtext('Sex', line = 0, font = 2)
#
hist(dsf.ps, ylim = c(0,4000), xlab = 'P-value', yaxt = 'n', main = '', col = 'tan2')
mtext('Days since blood meal', line = 0, font = 2)
axis(2, c(0, 1000, 2000, 3000, 4000), labels = c('0', '', '2000', '', '4000'))
#
hist(uses.ps, ylim = c(0,4000), xlab = 'P-value', yaxt = 'n', ylab = '', main = '', col = 'tan2')
mtext('Previous cold exposures', line = 0, font = 2)
dev.off()


