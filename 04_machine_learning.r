library(caret)
library(plotrix)

  ###############################
  # Load and create data tables #
  ###############################

tsetse <- read.table('tables/tsetse_qpcr.csv', sep = '\t', row.names = 1, header = T)
# There is one sample where many genes failed, including the housekeeping genes, so remove this from
# the table
tsetse <- tsetse[rownames(tsetse) != 'T39', ]

# Scale the data 
preProcValues <- preProcess(tsetse, method = c('center', 'scale'))
tsetse.scaled <- predict(preProcValues, tsetse)

# Load the metadata
metadata <- read.table('tables/sample_metadata.csv', sep = '\t', header = T, row.names = 1, quote = '')[rownames(tsetse), ]

tsetse.scaled$age <- metadata[, 'age']
tsetse.scaled$sex <- metadata[, 'sex']

# Set up the plotting colours
tsetse$age <- metadata[, 'age']
age.colours <- c(young = 'yellow', old = 'brown')
agecolscale <- color.scale(1:61, extremes = age.colours)
tsetse$agecol <- agecolscale[tsetse$age-1]
tsetse$agecol.cutoff <- c('orange', 'purple')[(tsetse$age > 15) + 1]
tsetse$sex <- metadata[, 'sex']
sex.colours <- c(female = 'royalblue3', male = 'orangered')
tsetse$sexcol <- sex.colours[as.character(tsetse$sex)]
sex.pch <- c(female = 24, male = 21)
tsetse$sexpch <- sex.pch[as.character(tsetse$sex)]

hk.genes <- c('GMOY003952', 'GMOY010976')
age.genes <- c('GMOY000749', 'GMOY001603', 'GMOY002920', 'GMOY003090', 'GMOY003371', 'GMOY003588', 'GMOY005053', 'GMOY005321', 'GMOY009908', 'GMOY011979')
genes <- c(age.genes, hk.genes)
norm.age.genes <- paste(age.genes, 'n', sep = '')

# Let's split the data into males and females. 
tsetse.m <- subset(tsetse, sex == 'male')
tsetse.f <- subset(tsetse, sex == 'female')


  ##################################
  # Create test / train data split #
  ##################################

# Let's create a data partition at random. The function below, rather than just sampling randomly, will
# split the "age" variable into percentile bins and sample randomly within each of these bins. This makes
# it more likely that a representative sample of the age range ends up in the training set. The number of
# percentiles is set by the "groups" argument. To make sure the sexes are also equally represented, we 
# split the data by sex, then split each of those into training and testing, and then recombine. 
set.seed(42)
# Females
inTrain.f <- createDataPartition(y = tsetse.f$age, p = 0.75, groups = 10)$Resample1
tsetse.train.f <- tsetse.f[inTrain.f, ]
tsetse.train.f$age.categ <- as.factor(tsetse.train.f$age > 15)
tsetse.test.f <- tsetse.f[-inTrain.f, ]
tsetse.test.f$age.categ <- as.factor(tsetse.test.f$age > 15)
# Males
inTrain.m <- createDataPartition(y = tsetse.m$age, p = 0.75, groups = 10)$Resample1
tsetse.train.m <- tsetse.m[inTrain.m, ]
tsetse.train.m$age.categ <- as.factor(tsetse.train.m$age > 15)
tsetse.test.m <- tsetse.m[-inTrain.m, ]
tsetse.test.m$age.categ <- as.factor(tsetse.test.m$age > 15)

tsetse.train <- rbind(tsetse.train.m, tsetse.train.f)
tsetse.test <- rbind(tsetse.test.m, tsetse.test.f)


  ##########################
  # Functions and controls # 
  ##########################

# Write a function that will calculate the size of the residual band within which some percentile of 
# predictions fall.
estimate.percentile.accuracy <- function(pred, obs, perc = 0.95){
	if (perc < 0.5 | perc > 1)
		stop('perc should be between 0.5 and 1')
	if (length(pred) != length(obs))
		stop('Predicted and observed values should have the same length.')
	sorted.residuals <- sort(abs(pred - obs))
	# If the percentile is so small that there are not enough samples for at least one sample to fall
	# in the percentile, give a warning and report NA
	perc.cutoff <- floor((1-perc) * length(sorted.residuals))
	if (perc.cutoff == 0){
		cat('The sample size is not large enough for even one observation to fall within the ', perc, 
            ' percentile. Only reporting max values.\n', sep = '')
		upper.value <- NA
	}
	else{
		upper.value <- sorted.residuals[length(sorted.residuals) + 1 - perc.cutoff]
	}
	output <- c(upper.value, max(sorted.residuals))
	names(output) <- c('upper', 'max')
	output
}

# Write a function that will take a set of predictions and observations, plot the performance
# and output some metrics
assess.performance <- function(model.fit, test.set, perc.acc.quantile = 0.9, include.plot = T, title = 'Summary', sex.col = test.sex.colours, sex.pch = test.sex.pch, new.plot = T){
	# We allow the function to receive a list of model fits and test sets, which will be concantenated
	if (all(class(model.fit) == 'list')){
		if ((class(test.set) != 'list') | (length(model.fit) != length(test.set)))
			stop('If model.fit is passed as a list, then test.set needs to be a list of the same length.')
		else {
			predictions <- list()
			observations <- list()
			categ.obss <- list()
			var.imp <- list()
			for (i in 1:length(model.fit)){
				predictions[[i]] <- predict(model.fit[[i]], newdata = test.set[[i]])
				observations[[i]] <- test.set[[i]]$age
				categ.obss[[i]] <- test.set[[i]]$age.categ
				var.imp[[i]] <- varImp(model.fit[[i]])
			}
			prediction <- unlist(predictions)
			observation <- unlist(observations)
			categ.obs <- unlist(categ.obss)
			if (!is.null(names(model.fit)))
				names(var.imp) <- names(model.fit)
		}
	}
	else {
		prediction <- predict(model.fit, newdata = test.set)
		observation <- test.set$age
		categ.obs <- test.set$age.categ
		# Get the variable importance
		var.imp <- varImp(model.fit)
	}
	if (class(prediction) == 'numeric'){
		# Build logical vectors for age categories
		young.age <- 15
		mid.age <- 30
		young.prediction <- prediction <= young.age
		mid.prediction <- prediction <= mid.age
		young.observation <- observation <= young.age
		mid.observation <- observation <= mid.age
		# Get some measure of the quality of the fit, for the whole data, and for the data where the PREDICTED
		# age is below a certain value
		fit.qual <- postResample(pred = prediction, obs = observation)
		fit.qual.young <- postResample(pred = prediction[young.prediction], obs = observation[young.prediction])
		fit.qual.mid <- postResample(pred = prediction[mid.prediction], obs = observation[mid.prediction])
		# Get the percentile accuracy
		perc.acc <- estimate.percentile.accuracy(prediction, observation, perc.acc.quantile)
		perc.acc.young <- estimate.percentile.accuracy(prediction[young.prediction], observation[young.prediction], perc.acc.quantile)
		perc.acc.mid <- estimate.percentile.accuracy(prediction[mid.prediction], observation[mid.prediction], perc.acc.quantile)
		# See how well we could perform a categorical prediction based on these regression predictions
		confusion <- confusionMatrix(data = as.factor(!young.prediction), reference = categ.obs, positive = 'FALSE', mode = 'everything')
		confusion.summary <- c(confusion$overall['Accuracy'], confusion$byClass[c('Sensitivity', 'Specificity')])
		# Now do the same for individuals 30 days or over
		confusion.30 <- confusionMatrix(data = as.factor(!mid.prediction), reference = as.factor(!mid.observation), positive = 'FALSE', mode = 'everything')
		confusion.summary.30 <- c(confusion.30$overall['Accuracy'], confusion.30$byClass[c('Sensitivity', 'Specificity')])
		# Make the plot
		if (include.plot){
			if (new.plot)
				x11()
			par(mfrow = c(1,3))
			plot(observation, prediction, type = 'n', xaxt = 'n', yaxt = 'n', ylim = c(min(0, prediction), max(62, prediction)))
			axis(1, at = seq(0, 60, 15))
			axis(2, at = seq(0, 60, 15))
			# Get user coordinates
			lef <- par('usr')[1]; rig <- par('usr')[2]; bot <- par('usr')[3]; top <- par('usr')[4]
			# Draw rectangles to show age categories
			#rect(lef, bot, c(young.age, rig, young.age), c(top, young.age, young.age), col = c(rep('grey90', 2), 'grey80'), border = NA)
			rect(lef, bot, rep(c(rig, mid.age, young.age), 3)[2:9], rep(c(top, mid.age, young.age), each = 3)[2:9], col = c('grey90', 'grey80', 'grey90', 'grey80', 'grey70', 'grey80', 'grey70', 'grey60'), border = NA)
			points(observation, prediction, pch = sex.pch, col = 'grey70', bg = sex.col, lwd = 0.6)
			abline(0,1, col = 'orchid4', lwd = 2)
			label.x <- lef + (rig - lef) / 20
			label.ys <- top - (1:7) * ((top - bot) / 20)
			text(label.x, label.ys, label = paste(c(names(fit.qual[1:2]), 'RMSE <= 15do', 'RMSE <= 30do', paste(perc.acc.quantile*100, c('% Acc', '% Acc < 30do', '% Acc < 16do'), sep = '')), round(c(fit.qual[1:2], fit.qual.young[1], fit.qual.mid[1], perc.acc['upper'], perc.acc.mid['upper'], perc.acc.young['upper']), 2), sep = ': '), adj = 0, cex = 0.8)
			# We can also plot the residuals against age
			plot(observation, prediction - observation, pch = sex.pch, col = 'grey70', bg = sex.col, lwd = 0.6, xaxt = 'n')
			axis(1, at = seq(0, 60, 15))
			abline(h=0, col = 'orchid4', lwd = 2)
			plot(prediction, prediction - observation, pch = sex.pch, col = 'grey70', bg = sex.col, lwd = 0.6, xaxt = 'n', xlim = c(min(0, prediction), max(62, prediction)))
			axis(1, at = seq(0, 60, 15))
			abline(h=0, col = 'orchid4', lwd = 2)
			mtext(title, 3, -2, T, font = 2, cex = 1.3)
		}
		output.list <- list(fitted.age = prediction, fit.quality = fit.qual, fit.quality.mid = fit.qual.mid, fit.quality.young = fit.qual.young,
		                    perc.acc = perc.acc, perc.acc.mid = perc.acc.mid, perc.acc.young = perc.acc.young, 
		                    full.confusion = confusion, confusion.summary = confusion.summary, 
		                    full.confusion.30 = confusion.30, confusion.summary.30 = confusion.summary.30, var.imp = var.imp)
	}
	else {
		# We set "False" (ie: age <= 15) as "positive for the calculation of metrics)
		confusion <- confusionMatrix(data = prediction, reference = categ.obs, positive = 'FALSE', mode = 'everything')
		confusion.summary <- c(confusion$overall['Accuracy'], confusion$byClass[c('Sensitivity', 'Specificity')])
		output.list <- list(fitted.age = prediction, full.confusion = confusion, confusion.summary = confusion.summary, var.imp = var.imp)
	}
	cat('\n', title, ':\n\n', sep = '')
	print(confusion$table)
	cat('\n')
	output.list
}


# Set up a control that performs 10-fold cross-validation 3 times
kfold_10_3_ctrl <- trainControl(method = 'repeatedcv', 
                                # number of folds
                                number = 10,
								# number of repeats
                                repeats = 3
                                )

test.sex.colours <- sex.colours[as.character(tsetse.test$sex)]
test.f.colours <- sex.colours[as.character(tsetse.test.f$sex)]
test.m.colours <- sex.colours[as.character(tsetse.test.m$sex)]
test.sex.pch <- sex.pch[as.character(tsetse.test$sex)]
test.f.pch <- sex.pch[as.character(tsetse.test.f$sex)]
test.m.pch <- sex.pch[as.character(tsetse.test.m$sex)]


  ####################
  # Machine learning # 
  ####################

# Let's start by fitting a pls.
pls.fit <- train(age ~ ., 
                 data = tsetse.train[, c(norm.age.genes, 'sex', 'age')], 
                 # the ML method that will be used
                 method = 'pls', 
                 # the pre-processing that will be done
                 preProc = c('center', 'scale'),
                 # a value representing how many different values of each parameter to try
                 tuneLength = 15,
                 # the validation method
                 trControl = kfold_10_3_ctrl
                )

pls.fitted.age <- assess.performance(pls.fit, tsetse.test, title = 'Basic PLS')

# Random forest regression
rfr.fit <- train(age ~ .,
                 data = tsetse.train[, c(norm.age.genes, 'sex', 'age')], 
                 # the ML method that will be used
                 method = 'rf', 
                 # the pre-processing that will be done
                 preProc = c('center', 'scale'),
                 # a value representing how many different values of each parameter to try
                 tuneLength = 15,
                 # the validation method
                 trControl = kfold_10_3_ctrl, 
				 # Apparently with Random Forest you need to specify that you want to calculate the variable importance
                 importance = T
                 )

rfr.fitted.age <- assess.performance(rfr.fit, tsetse.test, title = 'RF regressor')

# Xgboost
xgbGrid <- expand.grid(nrounds = c(100, 200), 
                       max_depth = c(3, 10, 15),
                       colsample_bytree = seq(0.5, 0.9, length.out = 5),
                       eta = 0.1,
                       gamma = 0,
                       min_child_weight = 1,
                       subsample = 1)

xgb.fit <- train(age ~ .,
                 data = tsetse.train[, c(norm.age.genes, 'sex', 'age')], 
                 # the ML method that will be used
                 method = 'xgbTree', 
                 # the pre-processing that will be done
                 preProc = c('center', 'scale'),
                 # The grid of parameter values to explore
                 tuneGrid = xgbGrid,
                 # the validation method
                 trControl = kfold_10_3_ctrl
                )

xgb.fitted.age <- assess.performance(xgb.fit, tsetse.test, title = 'XGB regressor')

# Now the random forest regression split by sex
# Males
rfr.fit.m <- train(age ~ .,
                   data = tsetse.train.m[, c(norm.age.genes, 'age')], 
                   # the ML method that will be used
                   method = 'rf', 
                   # the pre-processing that will be done
                   preProc = c('center', 'scale'),
                   # a value representing how many different values of each parameter to try
                   tuneLength = 15,
                   # the validation method
                   trControl = kfold_10_3_ctrl, 
				   # Apparently with Random Forest you need to specify that you want to calculate the variable importance
                   importance = T
                   )

rfr.fitted.age.m <- assess.performance(rfr.fit.m, tsetse.test.m, sex.col = test.m.colours, sex.pch = test.m.pch, title = 'RF regressor for males')

# Females
rfr.fit.f <- train(age ~ .,
                   data = tsetse.train.f[, c(norm.age.genes, 'age')], 
                   # the ML method that will be used
                   method = 'rf', 
                   # the pre-processing that will be done
                   preProc = c('center', 'scale'),
                   # a value representing how many different values of each parameter to try
                   tuneLength = 15,
                   # the validation method
                   trControl = kfold_10_3_ctrl, 
				   # Apparently with Random Forest you need to specify that you want to calculate the variable importance
                   importance = T
                  )

rfr.fitted.age.f <- assess.performance(rfr.fit.f, tsetse.test.f, sex.col = test.f.colours, sex.pch = test.f.pch, title = 'RF regressor for females')

# We combine the male and female predictions and see how we did
rfr.fitted.age.mf <- assess.performance(list(male = rfr.fit.m, female = rfr.fit.f), list(tsetse.test.m, tsetse.test.f), title = 'RF regressor for combined sexes')

# Now the xgb split by sex
# Males
xgb.fit.m <- train(age ~ .,
                   data = tsetse.train.m[, c(norm.age.genes, 'age')], 
                   # the ML method that will be used
                   method = 'xgbTree', 
                   # the pre-processing that will be done
                   preProc = c('center', 'scale'),
                   # The grid of parameter values to explore
                   tuneGrid = xgbGrid,
                   # the validation method
                   trControl = kfold_10_3_ctrl
                   )

xgb.fitted.age.m <- assess.performance(xgb.fit.m, tsetse.test.m, sex.col = test.m.colours, sex.pch = test.m.pch, title = 'XGB regressor for males')

# Females
xgb.fit.f <- train(age ~ .,
                   data = tsetse.train.f[, c(norm.age.genes, 'age')], 
                   # the ML method that will be used
                   method = 'xgbTree', 
                   # the pre-processing that will be done
                   preProc = c('center', 'scale'),
                   # The grid of parameter values to explore
                   tuneGrid = xgbGrid,
                   # the validation method
                   trControl = kfold_10_3_ctrl
                  )

xgb.fitted.age.f <- assess.performance(xgb.fit.f, tsetse.test.f, sex.col = test.f.colours, sex.pch = test.f.pch, title = 'XGB regressor for females')

# We combine the male and female predictions and see how we did
xgb.fitted.age.mf <- assess.performance(list(male = xgb.fit.m, female = xgb.fit.f), list(tsetse.test.m, tsetse.test.f), title = 'XGB regressor for combined sexes')

# Now try at classifiers. See if we can group flies according to whether or not they are more than 15 days old. 
# First try a decision tree
tree.fit <- train(age.categ ~ .,
                  data = tsetse.train[, c(norm.age.genes, 'sex', 'age.categ')], 
                  # the ML method that will be used
                  method = 'rpart', 
                  # the pre-processing that will be done
                  preProc = c('center', 'scale'),
                  # a value representing how many different values of each parameter to try
                  tuneLength = 15,
                  # the validation method
                  trControl = kfold_10_3_ctrl
                 )

tree.fitted.age <- assess.performance(tree.fit, tsetse.test, title = 'Basic Tree')

# Random forest classifier
rfc.fit <- train(age.categ ~ .,
                 data = tsetse.train[, c(norm.age.genes, 'sex', 'age.categ')], 
                 # the ML method that will be used
                 method = 'rf', 
                 # the pre-processing that will be done
                 preProc = c('center', 'scale'),
                 # a value representing how many different values of each parameter to try
                 tuneLength = 15,
                 # the validation method
                 trControl = kfold_10_3_ctrl
                )

rfc.fitted.age <- assess.performance(rfc.fit, tsetse.test, title = 'RF classifier')

# XGB classifier
xgbc.fit <- train(age.categ ~ .,
                  data = tsetse.train[, c(norm.age.genes, 'sex', 'age.categ')], 
                  # the ML method that will be used
                  method = 'xgbTree', 
                  # the pre-processing that will be done
                  preProc = c('center', 'scale'),
                  # The grid of parameter values to explore
                  tuneGrid = xgbGrid,
                  # the validation method
                  trControl = kfold_10_3_ctrl
                 )
 
xgbc.fitted.age <- assess.performance(xgbc.fit, tsetse.test, title = 'XGB classifier')

# Save the output
save.image('machine_learning.Rdata')

