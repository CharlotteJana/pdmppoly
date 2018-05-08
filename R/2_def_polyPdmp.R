library(prodlim) # needed for function "row.match"

##### Todo #####

# validity schreiben (was kommt da genau rein?, siehe validObject(.Object))
# arity(spray) == length(obj@init) für dynpolys und ratepolys
# ist arity = number cont variables oder arity = length(obj@init)? (entsprechend polyGenerator umändern)
# eventuell redefineDynpolys so ändern, dass der θ-Anteil in 'overall' wegfällt
# generator <-> jumptypes: which(jumpfunc(t,z,parms,jtype) == discVar)???
# matchRow: multipleMatches notwendig? brauche ich die funktion überhaupt noch?



##### methods ####

# needed for "ratepolys←" and "dynpolys←" :
redefineRatepolys <- function(x, obj, where = parent.frame()){
  # check if 'x' is defined correctly, numbers are coerced to 'spray'-objects.
  # (This is still an object of type 'language')
  
  eva <- with(as.list(obj@parms), eval(x)) #  only the names of the parameters are important
  dim <- length(obj@init)
  if(is.numeric(eva)) {bquote(.(x)*one(.(dim)))} # turn numbers into 'spray'-objects
  else if(is.spray(eva)) {x}
  else if(is.list(eva)) {as.call(lapply(x, redefineRatepolys, obj = obj, where = where))}
  else if(is.name(x) && identical(x, as.name("list"))) {x}
  else { stop("Input '", x, "' in 'ratepolys' is not correct.")}
}
redefineDynpolys <- function(x, obj, where = parent.frame(), overall = 0){
  # check if 'x' is defined correctly, numbers are coerced to 'spray'-objects,
  # all entrys of 'x' are coerced to the form 'variant 1' (description: see "dynpolys←")
  # (This is still an object of type 'language')
  
  eva <- with(as.list(obj@parms), eval(x)) #  only the names of the parameters are important
  eoa <- with(as.list(obj@parms), eval(overall))
  dim <- length(obj@init)
  
  # if 'overall' == number: turn it into a 'spray'-object
  if(is.numeric(eoa)) { 
    overall <- bquote(.(overall)*one(.(dim)));
    eoa <- with(as.list(obj@parms), eval(overall))
  }
  
  # if x == value (number or 'spray'): add 'overall' to x and turn x into a 'spray'-object (when necessary)
  if(is.numeric(eva)) {
    if(is.zero(eoa)) bquote(.(x)*one(.(dim))) # if 'overall' == 0: don't add anything, turn x into 'spray'
    else  bquote(.(overall)+.(x)*one(.(dim)))
  }
  else if(is.spray(eva)){
    if(is.zero(eoa)) x
    else  bquote(.(overall)+.(x))
  }
  
  # if x == 'list':  rerun 'redefineDyn'
  else if(is.name(x) && identical(x, as.name("list"))) {x}
  else if(is.list(eva)) {
    if(!is.null(x$overall) & is.null(x$specific)) {
      overall <- x$overall
      zeroList <- as.list(rep(0, length(obj@discDomain)+1));
      zeroList[[1]] <- "list"
      x <- do.call(call, zeroList)
    } 
    if(!is.null(x$overall)) {overall <- x$overall; x$overall <- NULL}
    if(!is.null(x$specific)) {x <- x$specific}
    as.call(lapply(x, redefineDynpolys, obj = obj, where = where, overall = overall))
  }
  
  # if x == something else: give error
  else {stop("Input ", x, " for 'dynpolys' is not correct.")}
}  

getIndex <- function(var, vect){
  index <- which(vect == var) # index of vect that corresponds to var
  if(length(index) == 0){ # if former calculation didnt work (because of numerical issues): check for "nearly" equality (this is less efficient)
    #print("all.equal was used")
    comp <- mapply(function(x) {isTRUE(all.equal(x, var, check.names = FALSE))}, vect) 
    index <- which(comp == TRUE)
  }
  return(index)
}      # get index of variable „var“ in vector „vect“
# matchRow <- function(row, matrix, multipleMatches = TRUE){
#   matchings <- apply(matrix, 1, identical, row)
#   print(matchings)
#   if(sum(matchings) == 0) return(NA)
#   if(multipleMatches) which(matchings)
#   else which(matchings)[1]
# } # get rowindex of a specific row in a matrix

# matchRow <- function (x, table, nomatch = NA){
#   if (class(table) == "matrix") 
#     table <- as.data.frame(table)
#   if (is.null(dim(x))) 
#     x <- as.data.frame(matrix(x, nrow = 1))
#   cx <- do.call("paste", c(x[, , drop = FALSE], sep = "\r"))
#   ct <- do.call("paste", c(table[, , drop = FALSE], sep = "\r"))
#   print(sapply(ct, function(y) identical(cx, y)))
#   match(cx, ct, nomatch = nomatch)
# }
ratespraysToMatrix <- function(obj){
  # compute a 'ratematrix' out of ratesprays with ratematrix[[i]][[j]] = rate from discrete state i to discrete state j
  # (ratematrix is a nested list of sprays)
  
  n <- length(obj@init) - 1 # number of continuous variables
  nj <- length(obj@ratesprays) # number of jumptypes
  k <- length(obj@discDomain) # number of different discrete states
  z <- obj@init
  
  ratematrix <- rep(list(rep(list(NULL), length(obj@discDomain))),length(obj@discDomain))
  for(jtype in 1:nj){ 
    for(i in 1:k){
      oldDiscVar <- obj@discDomain[[i]]
      z[length(z)] <- oldDiscVar
      for(j in 1:k){
        newDiscVar <- obj@discDomain[[j]]
        if(obj@jumpfunc(1, z, obj@parms, jtype)[length(z)] == newDiscVar){
          ratematrix[[i]][[j]] <- obj@ratesprays[[jtype]][[i]]
        }
      } 
    }
  }
  return(ratematrix)
}  # compute nested list with entry [[i]][[j]] = rate i → j





