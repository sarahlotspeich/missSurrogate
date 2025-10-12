kf=function(x, h){return(dnorm(x/h) / h)}

VTM<-function(vc, dm){
  #takes vc and makes it the repeated row of a matrix, repeats it dm times
  matrix(vc, ncol=length(vc), nrow = dm, byrow = T)
}

Kern.FUN <- function(zz,zi,bw) ## returns an (n x nz) matrix ##
{
  out = (VTM(zz,length(zi))- zi)/bw
  return(dnorm(out)/bw)

}

pred.smooth <-function(zz,zi.one, bw=NULL,y1, weight=NULL) {
  if(is.null(bw)) { bw = bw.nrd(zz)/((length(zz))^(.10))}
  if(is.null(weight)) {weight = rep(1, length(y1))}
  return(sum(Kern.FUN(zz,zi.one,bw=bw)*y1*(weight))/sum(Kern.FUN(zz,zi.one,bw=bw)*weight))
}

delta.s.single = function(sone,szero,yone,yzero, h.select = NULL, weight.1 = NULL, weight.0=NULL, n0.all=NULL) {
  #we can talk about the bandwidth later, this default should work ok
  if(is.null(h.select)) {h.select = bw.nrd(sone)*(length(sone)^(-0.25))
  }
  if(is.null(weight.1) & is.null(weight.0)){
    mu.1.s0 = sapply(szero,pred.smooth,zz=sone, bw=h.select, y1=yone)
    if(sum(is.na(mu.1.s0))>0){
      print(paste("Note: ", sum(is.na(mu.1.s0)), " values extrapolated."))
      c.mat = cbind(szero, mu.1.s0)
      for(o in 1:length(mu.1.s0)) {
        if(is.na(mu.1.s0[o])){
          distance = abs(szero - szero[o])
          c.temp = cbind(c.mat, distance)
          c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where mean is not na
          new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
          mu.1.s0[o] = new.est[1]   #in case there are multiple matches
        }
      }}
    delta.s = mean(mu.1.s0) - mean(yzero)
  }
  if(!is.null(weight.1) & !is.null(weight.0)){
    mu.1.s0 = sapply(szero,pred.smooth,zz=sone, bw=h.select, y1=yone, weight=weight.1)
    if(sum(is.na(mu.1.s0))>0){
      print(paste("Note: ", sum(is.na(mu.1.s0)), " values extrapolated."))
      c.mat = cbind(szero, mu.1.s0)
      for(o in 1:length(mu.1.s0)) {
        if(is.na(mu.1.s0[o])){
          distance = abs(szero - szero[o])
          c.temp = cbind(c.mat, distance)
          c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where mean is not na
          new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
          mu.1.s0[o] = new.est[1]   #in case there are multiple matches
        }
      }}
    delta.s = sum((weight.0)*mu.1.s0)/n0.all-sum((weight.0)*yzero)/n0.all
  }
  return(delta.s)
}

expit = function(x){
  return((exp(x))/(1+exp(x)))
}
