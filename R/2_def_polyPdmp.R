library(prodlim) # needed for function "row.match"

##### Todo #####

# validity schreiben (was kommt da genau rein?, siehe validObject(.Object))
# arity(spray) == length(obj@init) für dynpolys und ratepolys
# ist arity = number cont variables oder arity = length(obj@init)? (entsprechend polyGenerator umändern)
# eventuell redefineDynpolys so ändern, dass der θ-Anteil in 'overall' wegfällt
# generator <-> jumptypes: which(jumpfunc(t,z,parms,jtype) == discVar)???
# matchRow: multipleMatches notwendig? brauche ich die funktion überhaupt noch?



##### methods ####



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

