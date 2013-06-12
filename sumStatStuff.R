"sumStatVar" <-
  function(qvec, tvar, mid)
{
### qbars by matched set
  qbar <- unsplit( lapply( split(qvec, mid), 
                          function(x) {rep(mean(x), length(x))} ) ,
                  mid )
  
  sum(
### stratum scaling for variance
      tapply(tvar, mid,
             function(x) {
               ifelse(length(x)>1 ,sum(x)*sum(!x)/(length(x)*(length(x) -1)), 0) }
             ) *
### stratum ss qs
      tapply((qvec - qbar), mid, function(x){sum(x^2)})
### close call to sum()
      )
}
"sumStatOMinusE" <-
  function(qvec, tvar, mid)
{
### qbars by matched set
  qbar <- unsplit( lapply( split(qvec, mid), 
                          function(x) {rep(mean(x), length(x))} ) ,
                  mid )
  sum(tvar*(qvec-qbar))
}
