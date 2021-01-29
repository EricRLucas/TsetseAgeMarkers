library(glmnet)

load('RNAseq_analysis.Rdata')


# The basic edgeR analysis asks whether expression of each gene is affected by X (eg: age), so expression is
# the dependent variable. Here, we want to build a model that predicts age based on expression, so expression
# should be the explanatory variable, with age being the dependent variable.  

# Find the best lambda using cross-validation
set.seed(42)
cv.lasso <- cv.glmnet(t(normalised.totals), age, alpha = 1, family = "gaussian")
plot(cv.lasso)
# The smallest mean-square error is found at the smallest value of lambda, but we can see that MSE has 
# plateaued by this point, so carry on with this value. 

# With this model, we manage to eliminate 11226 variables, leaving 52 (plus the intercept):
all.variables <- coef(cv.lasso, cv.lasso$lambda.min)
kept.variables <- all.variables[all.variables[,1] != 0,1]

# Fit the final model on the training data.  
model <- glmnet(t(normalised.totals), age, alpha = 1, family = "gaussian", lambda = cv.lasso$lambda.min)
kept.genes <- model$beta[model$beta[,1] != 0, 1]
# This keeps 55 genes
# Plot the expression of those genes and see what they look like
plot.gene <- function(gene.name, ...){
	y <- normalised.totals[gene.name, ]
	plot(age, y, main = gene.name, col = sample.conversion$sexcol, ylab = 'Normalised read counts', ylim = c(0, max(y)), pch = 19, ...)
	if (i == 1){
		points(c(40,40), c(0.1, 0.2)*max(y), col = c('red', 'blue'), pch = 19)
		text(c(42,42), c(0.1, 0.2)*max(y), label = c('male', 'female'), adj = 0)
	}
}

plot.genes <- function(gene.list){
	for (i in 1:ceiling(length(gene.list)/12)){
		x11()
		par(mfrow = c(3,4), mgp = c(1.5,0.5,0), mar = c(3,3,2,1))
		if (i*12 > length(gene.list))
			these.genes <- names(gene.list[((i-1)*12+1):length(gene.list)])
		else
			these.genes <- names(gene.list[((i-1)*12+1):(i*12)])
		for (gene in these.genes){
			plot.gene(gene)
		}
	}
}

plot.genes(kept.genes)

# GMOY001603 looks promising for continuing to work in old age. 000155 and 000689 also looks good. 
# In fact, all of the chosen genes tend to show linear(ish) increases with age, perhaps not surprising given
# that we built a linear model. So a non-linear gene kay not work well. 

# What if we choose a smaller value of lambda to narrow it down to just 10 genes?
# A lambda of 7 works for this
model.10 <- glmnet(t(normalised.totals), age, alpha = 1, family = "gaussian", lambda = 7)
kept.genes.10 <- model.10$beta[model.10$beta[,1] != 0, 1]

plot.genes(kept.genes.10)

# 001603 and 000689 are still there. 008370 looks good for early to middle age. 


# Now what if we ask it to categorise individuals (<= 15 days or > 15 days) instead of actually estimate age? 
age.cat.15 <- as.numeric(age <= 20)
# With this function, if you run it several times, you get very different results, with a different minimum. 
# This wasn't so much the case when we ran it on the continuous data, but is very much the case here. I think
# it's to do with which samples end up getting used for validation, which is maybe set randomly and kept the 
# same for all values of lambda within a run, but not between runs. To make it more consistent, we run it
# with leave-one-out validation, which means that there is no longer any randomness in the validation step
cv.lasso.cat.15 <- cv.glmnet(t(normalised.totals), age.cat.15, alpha = 1, nfolds = ncol(normalised.totals), family = "binomial")
plot(cv.lasso.cat.15)
# We are left with 20 genes
model.cat.15.min <- glmnet(t(normalised.totals), age.cat.15, alpha = 1, family = "binomial", lambda = cv.lasso.cat.15$lambda.min)
kept.genes.cat.15.min <- model.cat.15.min$beta[model.cat.15.min$beta[,1] != 0, 1]
# Let's look at them
plot.genes(kept.genes.cat.15.min)

# 001603 still coming out strong!
# From a previous run that I can't recreate because I can't remember the seed, 009908 looked good too.

# What if we split the analysis by sex? 
normalised.males <- normalised.totals[, sex == 'male']
age.males <- age[sex == 'male']
cv.lasso.male <- cv.glmnet(t(normalised.males), age.males, alpha = 1, nfolds = ncol(normalised.males), family = "gaussian")
plot(cv.lasso.male)
# We have 21 genes
model.male <- glmnet(t(normalised.males), age.males, alpha = 1, family = "gaussian", lambda = cv.lasso.male$lambda.min)
kept.genes.male <- model.male$beta[model.male$beta[,1] != 0, 1]
plot.genes(kept.genes.male)
# 005053 looks like it could be useful for males but not females

normalised.females <- normalised.totals[, sex == 'female']
age.females <- age[sex == 'female']
cv.lasso.female <- cv.glmnet(t(normalised.females), age.females, alpha = 1, nfolds = ncol(normalised.females), family = "gaussian")
plot(cv.lasso.female)
# We have 30 genes
model.female <- glmnet(t(normalised.females), age.females, alpha = 1, family = "gaussian", lambda = cv.lasso.female$lambda.min)
kept.genes.female <- model.female$beta[model.female$beta[,1] != 0, 1]
plot.genes(kept.genes.female)
# 011979 and 003588 look good (for both sexes) but nothing is obviously crucial for females specifically.
# 000749 looks good and increases with age


# Things to bear in mind when choosing primers:
# 1. No evidence that the gene is differentially expressed by time since feeding. 
# 2. Ideally not from the same WGCNA module as another chosen gene
# 3. Is the coverage at the primer site representative of the gene?
# 4. Any evidence of variation in the primer site (within the colony, or relative to the genome)

# OK, let's make our shortlist
# From differential expression analysis:
# 005321/005590 looks perfect for around < 12 days in females and < 8 in males
# 010232/002920/008748 are good for slightly older (around < 20 days in females and < 12 in males)
# 008370 good for slightly older still (this one also came out from lasso analysis)
# 009985 maybe for old age? Has the advantage that it increases rather than decreases with age. 
# 007669 (differentiates middle age from other ages in females?)
# From lasso analysis:
# 001603 is consistently there so we want that one
# 000689
# 003090
# 009908
# 011979
# 003588
# 005053 for males (has an increase with age)
# 000749 for females (has an increase with age)
# 000155
# 003371 (increases with age)

shortlist <- c('GMOY005321', 'GMOY010232', 'GMOY008370', 'GMOY009985', 'GMOY007669', 'GMOY001603', 'GMOY000689', 'GMOY003090', 'GMOY009908', 'GMOY011979', 'GMOY003588', 'GMOY005053', 'GMOY000749', 'GMOY000155', 'GMOY003371')

# What about WGCNA modules? 
unique(mergedColors_hybrid_4[shortlist])
# We are representing four clusters overall (with two of the shortlist having not been assigned a cluster)
# Let's show the cluster membership through the plot title
par(mfrow = c(3,5), mgp = c(1.5,0.5,0), mar = c(3,3,2,1))
for (gene in shortlist){
	if (is.na(mergedColors_hybrid_4[gene]))
		title.col <- 'black'
	else if (mergedColors_hybrid_4[gene] == 'white')
		title.col <- 'yellow'
	else
		title.col <- mergedColors_hybrid_4[gene]
	plot.gene(gene, col.main = title.col)
}
# There are 4 modules represented (plus the grey one), but it will be hard to improve on this, as most modules 
# aren't affected by age. 

# After primer testing, GMOY010232 performed poorly, perhaps because it has low expression. Replacing it with 
# GMOY002920 in the shortlist
new.shortlist <- c('GMOY005321', 'GMOY002920', 'GMOY008370', 'GMOY009985', 'GMOY007669', 'GMOY001603', 'GMOY000689', 'GMOY003090', 'GMOY009908', 'GMOY011979', 'GMOY003588', 'GMOY005053', 'GMOY000749', 'GMOY000155', 'GMOY003371')

# Plot the expression
par(mfrow = c(3,5), mgp = c(1.5,0.5,0), mar = c(3,3,2,1))
for (gene in new.shortlist){
	if (is.na(mergedColors_hybrid_4[gene]))
		title.col <- 'black'
	else if (mergedColors_hybrid_4[gene] == 'white')
		title.col <- 'yellow'
	else
		title.col <- mergedColors_hybrid_4[gene]
	plot.gene(gene, col.main = title.col)
}
