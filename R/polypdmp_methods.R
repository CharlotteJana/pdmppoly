
blowupSpray <- function(obj, spray){
  # This function blows the last variable of the spray object 
  # (which stands for the discrete variable θ) up to several different 
  # indicator varialbes θ₁,…,θₖ, where k = length(obj@discDomain).
  
  # example: x² + 2*θ + θy becomes x² + 0*θ₁ + 2*1*θ₂ + 0*θ₁*y + 1*θ₂*y, if discDomain = c(0,1)
  #          blowupPoly(polyModel9, linear(c(1,0,0),2) + linear(c(0,0,2),1) + product(c(0,1,1)))
  
  if(is.zero(spray)) return(spray)
  
  n <- length(obj@init) # number of continous variables
  k <- length(obj@discDomain) # number of different discrete states
  spray <- increase_arity(spray, n - 1 + k)
  
  # part of spray that is independent of θ:
  core <- subs(spray, n, 0, keepArity = TRUE) 
  blowedSpray <- core
  
  # add the specific parts (depend on θ):
  for(i in 1:k){
    discVar <- obj@discDomain[i]
    if(!is.zero(spray-core)){
      specific <- subs(spray-core, n, discVar, keepArity = TRUE)*lone(n-1+i, n-1+k)
      blowedSpray <- blowedSpray + specific
    }
  }
  return(blowedSpray)
}  # change discrete variable θ to indicator variables θ₁,…,θₖ


##### output methods ####

setMethod(f="print",
          signature="polyPdmpModel",
          definition=function(x, ...)
          {
            .local <- function (x, all = FALSE, ...) 
            {
              if (all) {
                print.default(x, all = TRUE, ...)
              }
              else {
                cat("An S4-object of class", class(x)[1], "\n\n")
                slotnames <- slotNames(x)
                for (slotname in slotnames) {
                  slotcontent <- slot(x, slotname)
                  if (!is.null(slotcontent)) {
                    if (slotname == "main" || 
                        slotname == "dynfunc" ||
                        slotname == "dynsprays" ||
                        slotname == "ratesprays" ||
                        slotname == "ratefunc") next
                    cat("Slot ", dQuote(slotname), ":\n", sep = "")
                    if (slotname == "out")    cat("  outputs exist ...\n")
                    else if (slotname == "parms")  str(slotcontent)
                    else {
                      print(slotcontent)
                    }
                    cat("\n")
                  }
                }
                cat("Hint: use print(x, all=TRUE) to see all polyPdmpModel slots.\n")
              }
            }
            .local(x, ...)
          }
)
