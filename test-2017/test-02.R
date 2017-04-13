

varCoef.list <- llply(1:12, function(i) ldply(VSI, function(x) x[[i]][["varCoef"]], .id = NULL), .progress = "text") %>% set_names(1:12)

sd.list   <- llply(1:12, function(i) ldply(VSI,   function(x) x[[i]]$sdDetVec, .id = "Id"), .progress = "text") %>% set_names(1:12)
mean.list <- llply(1:12, function(i) ldply(VSI, function(x) x[[i]]$meanVec, .id = "Id"), .progress = "text") %>% set_names(1:12)
## need to test monthly

i = 5

# X$Id %<>% as.numeric()



