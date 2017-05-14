
ebay = read.table('eBayNumberOfBidderData.dat', header = TRUE)

glmModel = glm(nBids ~ 0 + ., data=ebay, family=poisson)
glmModel

# Significant:
# verify, sealed, majblem, minbidsshare