

#' threshold the number of species, ranked by abundance
#' @export
thrspp <- function(datl,thr){

  sdi <- data.frame(t =  seq(20000,2000000,20000), di = 0)
  for(rr in 1:length(datl)){
    ag <- datl[[rr]]
    ag <- ag[ag$cum<thr,]

    if(nrow(ag)>0){
      un <- ag$Freq*(ag$Freq-1)
      sdi$di[rr] <- 1 - (sum(un)/(sum(ag$Freq)*(sum(ag$Freq)-1)))
    }
  }
  return(sdi)
}


#' Plot Simpson's diversity index
#' @export
diversityplot <- function(cd){

  lin <- thrspp(cd,1)
  plot(x=lin$t,y=lin$di,type="l",xlab="Timesteps",ylim=c(0,1),xlim=c(0,2000000),ylab="Simpson's diversity")

  tvals <- c(0.75,0.5,0.25,0.1)
  cols <- rainbow(length(tvals))
  ci <- 1
  for(th in tvals){
    #message(sprintf("th is %0.2f",th))
    lin <- thrspp(cd,th)
    lines(x=lin$t,y=(lin$di-(th*0)),col=cols[ci])
    ci <- ci+1
  }

  legend("bottomleft",legend=c("100%","75%","50%","25%","10%"),lty=1,col=c("black",cols))
}
