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

delta.s.single.ipw = function(sone,szero,yone,yzero, h.select = NULL, weight.1, weight.0) {
	s1.observed = sone[!is.na(sone)]
  	y1.observed = yone[!is.na(sone)] 
  	s0.observed = szero[!is.na(szero)]
 	w1.observed = weight.1[!is.na(sone)]
 	w0.observed = weight.0[!is.na(szero)]
 	
  	if(is.null(h.select)) {h.select = bw.nrd(s1.observed)*(length(s1.observed)^(-0.25))
  }
  
  mu.1.s0 = sapply(s0.observed,pred.smooth,zz=s1.observed, bw=h.select, y1=y1.observed, weight=w1.observed)
    if(sum(is.na(mu.1.s0))>0){
      print(paste("Note: ", sum(is.na(mu.1.s0)), " values extrapolated."))
      c.mat = cbind(s0.observed, mu.1.s0)
      for(o in 1:length(mu.1.s0)) {
        if(is.na(mu.1.s0[o])){
          distance = abs(s0.observed - s0.observed[o])
          c.temp = cbind(c.mat, distance)
          c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where mean is not na
          new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
          mu.1.s0[o] = new.est[1]   #in case there are multiple matches
        }
      }}
    delta.s = sum(w0.observed*mu.1.s0)/sum(w0.observed)- mean(yzero)
  return(delta.s)
}

expit = function(x){
  return((exp(x))/(1+exp(x)))
}
