# Just print out a summary of the important features of the two best prediction models

load('machine_learning.Rdata')

print(pls.fitted.age)

print(rfr.fitted.age)

print(xgb.fitted.age)

print(tree.fitted.age)

print(rfc.fitted.age)

print(xgbc.fitted.age)
