##' .. content for \description{} (no empty lines) ..
##'  For matched pairs, the singleton is just the treatment group member, ie pairs are grouped with one-many matched sets.
##' .. content for \details{} ..
##' @title Identify the one among one-many and many-one matched sets
##' @param omatch An optmatch object
##' @return list of named vectors theOne.name, a character vector;
##' theOne.position, a vector of integers; and isOneMany, a logical vector.
##' Names of these vectors are just the levels of omatch.  Entries of theOne.name
##' are names of the omatch vector, where theOne.position are corresponding positions in that vector.
##' @author Ben B Hansen
findMatchedSetSingletons <- function(omatch)
  { # NOT YET TESTED!
stopifnot(inherits(omatch, "optmatch"),
          !is.null(theTx <- attr(omatch, "contrast.group"))
          )
nlev <- nlevels(omatch)
levs <- levels(omatch)
isOneMany <- (table(omatch[theTx, drop=FALSE]) <2) #ATTN PLEASE
### probably need to coerce
### this to vector/remove weird table attribute

theOne.pos <- integer(nlev)
theOne.name <- character(nlev)
names(theOne.pos) <- levs

for (lev in levs) theOne.pos[lev] <- which(omatch==lev && theTx==isOneMany[lev])
theOne.name <- names(omatch)[theOne.pos]
names(theOne.name) <- levs

return(list(theOne.name=theOne.name, theOne.position=theOne.pos, isOneMany=isOneMany))
### Can immediately coerce the above to data frame if desired.
### This offers user more flexibility than returning as a data frame,
### however, since pulling a column out of a d.f. tends to miss the row names.
  }
