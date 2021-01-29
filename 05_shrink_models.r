library(caret)

load('machine_learning.Rdata')

set.seed(42)

test.train.split <- function(){
	inTrain.f <- createDataPartition(y = tsetse.f$age, p = 0.75, groups = 10)$Resample1
	train.f <- tsetse.f[inTrain.f, ]
	train.f$age.categ <- as.factor(train.f$age > 15)
	test.f <- tsetse.f[-inTrain.f, ]
	test.f$age.categ <- as.factor(test.f$age > 15)
	# Males
	inTrain.m <- createDataPartition(y = tsetse.m$age, p = 0.75, groups = 10)$Resample1
	train.m <- tsetse.m[inTrain.m, ]
	train.m$age.categ <- as.factor(train.m$age > 15)
	test.m <- tsetse.m[-inTrain.m, ]
	test.m$age.categ <- as.factor(test.m$age > 15)

	train.set <- rbind(train.m, train.f)
	test.set <- rbind(test.m, test.f)
	list('train' = train.set, 'test' = test.set)
}

# We run each model 20 times in order to get an idea of how variable the RMSE is. 
num.repeats <- 20

xgb.ranked.factors <- rownames(xgb.fitted.age$var.imp$importance)
rfr.ranked.factors <- rownames(rfr.fitted.age$var.imp$importance)
xgbc.ranked.factors <- rownames(xgbc.fitted.age$var.imp$importance)
rfc.ranked.factors <- rownames(rfc.fitted.age$var.imp$importance)
# Let's keep making ML models using fewer and fewer of these effects 
xgb.fitted.models <- list()
rfr.fitted.models <- list()
xgbc.fitted.models <- list()
rfc.fitted.models <- list()
for (i in 10:1){
	this.model.name <- paste('using', i, 'variables')
	xgb.fitted.models[[this.model.name]] <- list()
	rfr.fitted.models[[this.model.name]] <- list()
	xgbc.fitted.models[[this.model.name]] <- list()
	rfc.fitted.models[[this.model.name]] <- list()
	for (j in 1:num.repeats){
		this.test.train.split <- test.train.split()
		xgb.new.fit <- train(age ~ .,
		                     data = this.test.train.split$train[, c(xgb.ranked.factors[1:i], 'age')], 
		                     method = 'xgbTree', 
		                     preProc = c('center', 'scale'),
		                     tuneGrid = xgbGrid,
		                     trControl = kfold_10_3_ctrl
							)
		xgb.fitted.models[[this.model.name]][[j]] <- assess.performance(xgb.new.fit, this.test.train.split$test, include.plot = F)
		#
		rfr.new.fit <- train(age ~ .,
		                     data = this.test.train.split$train[, c(rfr.ranked.factors[1:i], 'age')], 
		                     method = 'rf', 
		                     preProc = c('center', 'scale'),
		                     tuneLength = 15,
		                     trControl = kfold_10_3_ctrl, 
		                     importance = T
							)
		rfr.fitted.models[[this.model.name]][[j]] <- assess.performance(rfr.new.fit, this.test.train.split$test, include.plot = F)
		# 
		xgbc.new.fit <- train(age.categ ~ .,
		                      data = this.test.train.split$train[, c(xgbc.ranked.factors[1:i], 'age.categ')], 
		                      method = 'xgbTree', 
		                      preProc = c('center', 'scale'),
		                      tuneGrid = xgbGrid,
		                      trControl = kfold_10_3_ctrl
							 )
		xgbc.fitted.models[[this.model.name]][[j]] <- assess.performance(xgbc.new.fit, this.test.train.split$test, include.plot = F)
		#
		rfc.new.fit <- train(age.categ ~ .,
		                     data = this.test.train.split$train[, c(rfc.ranked.factors[1:i], 'age.categ')], 
		                     method = 'rf', 
		                     preProc = c('center', 'scale'),
		                     tuneLength = 15,
		                     trControl = kfold_10_3_ctrl, 
		                     importance = T
							)
		rfc.fitted.models[[this.model.name]][[j]] <- assess.performance(rfc.new.fit, this.test.train.split$test, include.plot = F)
	}
}

all.xgb.RMSE <- sapply(xgb.fitted.models, function(L) sapply(L, function(x) x$fit.quality['RMSE']))
all.rfr.RMSE <- sapply(rfr.fitted.models, function(L) sapply(L, function(x) x$fit.quality['RMSE']))
all.xgbc.accuracy <- sapply(xgbc.fitted.models, function(L) sapply(L, function(x) x$full.confusion$overall['Accuracy']))
all.rfc.accuracy <- sapply(rfc.fitted.models, function(L) sapply(L, function(x) x$full.confusion$overall['Accuracy']))

save.image('shrunk_models.Rdata')



