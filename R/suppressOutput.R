suppressOutput <- function(...,verbose=0){
  if(verbose > 0) return(eval(...)) else return(capture.output(eval(...)))
}