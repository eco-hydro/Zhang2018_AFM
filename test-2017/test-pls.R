valid_pls <- function(X, Y, minVar = 10, method = c("SSX", "VIP"), plsModel = FALSE){
  method <- toupper(method[1])
  
  nvar <- rev(minVar:ncol(X)); vars <- colnames(X)
  errorId <- list()
  
  fits <- foreach(i = nvar, k = icount()) %do%{
    #acording to PRESS and Q2 select comps
    plsr_fit <- tryCatch(
      plsreg1(X, Y, comps = NULL),
      error = function(e) {
        errorId[[k]] <<- k
        # message(sprintf("[%d]error: ", i))
        message(sprintf("[%d] %s", i, e))
        plsreg1(X, Y, comps = 2)
      }
    )
    #delete unimportant variables based on "SSX" or "pls"
    if (method == "SSX"){
      SSXi <- SSX_plsreg1(plsr_fit, X)
      Id_del <- which.min(SSXi)
    }else if(method == "VIP"){
      Id_del <- VIP_plsreg(plsr_fit) %>% {which.min(.[nrow(.), ])}
    }else{
      stop("method should be 'SSX' or 'VIP'!\n")
    }
    
    X <- X[, -Id_del, drop = F]
    plsr_fit
  }
  fits %<>% set_names(nvar)
  Id_del <- unlist(errorId)
  # if (length(Id_del)) fits <- fits[-Id_del]
  
  Q2 <- llply(fits, `[[`, "Q2") %>% ldply(function(x) cbind.data.frame(predI = 1:nrow(x), x), .id = "nvar")
  if (nrow(Q2) > 0) Q2$nvar %<>% {as.numeric(as.character(.))}
  
  #according to Q2 select significant variables
  get_sign <- function(x){
    limq2 <- 0.05
    Id <- which(x$Q2 > limq2)
    x[Id, ]
  }
  info <- ddply(Q2, .(nvar), get_sign)
  if (nrow(info) > 0) info %<>% mutate(PRESS_adj = PRESS/(nrow(X) - predI - 1))
  #delete nvars which not all components Q2 is significant
  if (nrow(info) > 0){
    Id_nvar <- split(info$predI, info$nvar) %>% {which(laply(., function(x) 1 %in% x))}
    info <- subset(info, nvar %in% unique(info$nvar)[Id_nvar]) %>% set_rownames(NULL)
    # info <- info[Id_var, ]
  }
  
  # if (plsModel) {
  Id.best <- which(nvar == info$nvar[which.min(info$PRESS_adj)])
  plsr_fit <- fits[[Id.best]]
  vip <- VIP_plsreg(plsr_fit) %>% {.[nrow(.), ]}
  coef <- plsr_fit$std.coefs
  COEF <- VIP <- {numeric(length(vars))*NA} %>% set_names(vars)
  Id.var <- match(names(vip), vars)
  VIP[Id.var] <- vip
  COEF[Id.var] <- coef
  coef = data.frame(vip = VIP, coef = COEF)#RETURN
  # }else{
  # info#RETURN
  # }
  return(list(coef = coef, validInfo = info))
}
