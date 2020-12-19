## get only the top x percent of loadings components from a pcaObject if needed

getComponents <- function(
  pcaobj,
  components = NULL)
{
  if (!is.null(components)) {
    return(pcaobj$components[components])
  } else if (is.null(components)) {
    return(pcaobj$components)
  }
}

components = getComponents(p, c(1:3))
{
  # filter in the variables in the top percentRetain of the loadings range
  x <- pcaobj$loadings[,components]
  retain <- c()
  for (i in seq_along(components)) {
    # for each PC, based on the loadings range, calculate the rangeRetain
    # fraction of the range
    offset <- (max(x[,i]) - min(x[,i])) * rangeRetain
    
    # to obtain upper and lower cut-offs, offset max and min by the offset
    uppercutoff <- max(x[,i]) - offset
    lowercutoff <- min(x[,i]) + offset
    
    # build a consensus list of variables that will be included
    retain <- unique(c(retain,
                       which(x[,i] >= uppercutoff),
                       which(x[,i] <= lowercutoff)))
  }
  message('-- variables retained:')
  message(paste0(rownames(x)[retain], collapse = ', '))
  x <- x[retain,]
  
rangeRetain = 0.01
    # filter in the variables in the top percentRetain of the loadings range
    x <- p$loadings[,components]
    retain <- c()
    for (i in seq_along(components)) {
      # for each PC, based on the loadings range, calculate the rangeRetain
      # fraction of the range
      offset <- (max(x[,i]) - min(x[,i])) * rangeRetain
      
      # to obtain upper and lower cut-offs, offset max and min by the offset
      uppercutoff <- max(x[,i]) - offset
      lowercutoff <- min(x[,i]) + offset
      
      # build a consensus list of variables that will be included
      retain <- unique(c(retain,
                         which(x[,i] >= uppercutoff),
                         which(x[,i] <= lowercutoff)))
    }
    message('-- variables retained:')
    message(paste0(rownames(x)[retain], collapse = ', '))
    x <- x[retain,]
dim(x)






