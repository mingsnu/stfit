### pidx0.1
res = readRDS("./pidx0.1/res.rds")
RMSEmat1 = matrix(unlist(lapply(res, function(x)x[1])), 20)
RMSEmat2 = readRDS("../simulation2/pidx0.1/RMSEmat2.rds")
RMSEmat1 < RMSEmat2
table(c(RMSEmat1 < RMSEmat2))
# FALSE  TRUE 
# 13    82 
table(c(RMSEmat1.0 < RMSEmat2))
# FALSE  TRUE 
# 26    69 
RMSEmat1.0 = readRDS("../simulation2/pidx0.1/RMSEmat1.rds")
RMSEmat1 < RMSEmat1.0
table(c(RMSEmat1 < RMSEmat1.0))

res = readRDS("./pidx0.4_0.6/res.rds")
RMSEmat1 = matrix(unlist(lapply(res, function(x)x[1])), 20)
RMSEmat2 = readRDS("../simulation2/pidx0.4_0.6/RMSEmat2.rds")
RMSEmat1 < RMSEmat2
table(c(RMSEmat1 < RMSEmat2))
# FALSE  TRUE 
# 5    90 
table(c(RMSEmat1.0 < RMSEmat2))
# FALSE  TRUE 
# 24    71 
RMSEmat1.0 = readRDS("../simulation2/pidx0.4_0.6/RMSEmat1.rds")
RMSEmat1 < RMSEmat1.0
table(c(RMSEmat1 < RMSEmat1.0))

res = readRDS("./pidx0.8_0.95/res.rds")
RMSEmat1 = matrix(unlist(lapply(res, function(x)x[1])), 20)
RMSEmat2 = readRDS("../simulation2/pidx0.8_0.95/RMSEmat2.rds")
RMSEmat1 < RMSEmat2
table(c(RMSEmat1 < RMSEmat2))
# FALSE  TRUE 
# 5    90 
table(c(RMSEmat1.0 < RMSEmat2))
# FALSE  TRUE 
# 19    76
RMSEmat1.0 = readRDS("../simulation2/pidx0.8_0.95/RMSEmat1.rds")
RMSEmat1 < RMSEmat1.0
table(c(RMSEmat1 < RMSEmat1.0))
