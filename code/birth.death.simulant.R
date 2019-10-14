#
# Forward simulation of process

#note, condCounts 
birth.death.simulant = function(t,X0=1,lambda=1,mu=2,nu=1,
                                condCounts=NULL) {
  hx = list()
  hx.times  = c()
  hx.states = c()
  cTime = 0
  cState = X0
  if (is.null(condCounts)){
    while( cTime <= t ) { 
      hx.times  = c(hx.times,cTime)
      hx.states = c(hx.states,cState)
      if( nu==0 && cState == 0 ) {
        break()
      }
      birth = lambda * cState
      death = mu * cState
      rate = nu + birth + death
      cTime = cTime + rexp(n=1,rate=rate)
      event = sample(3,1,prob=c(birth,death,nu)/rate)
      if( event == 1 || event == 3 ) {
        cState = cState + 1
      } else {
        cState = cState - 1
      }
    }
  }
  else {
    while( cTime <= t ) { 
      hx.times  = c(hx.times,cTime)
      hx.states = c(hx.states,cState)
      {#quit if you violate conditional assumptions
        nis <- NijBD(hx.states);
        cNminus <- sum(nis[1,]); #current Nminus
        cNplus <- sum(nis[2,]);
        if (cNminus > condCounts["Nminus"] || cNplus > condCounts["Nplus"]) {
          break();
        }
      }
      if( nu==0 && cState == 0 ) {
        break()
      }
      birth = lambda * cState
      death = mu * cState
      rate = nu + birth + death
      cTime = cTime + rexp(n=1,rate=rate)    
      event = sample(3,1,prob=c(birth,death,nu)/rate)
      if( event == 1 || event == 3 ) {
        cState = cState + 1
      } else {
        cState = cState - 1
      }
    }
  }
  new("BDMC", times=hx.times,states=hx.states,T=t)
}
