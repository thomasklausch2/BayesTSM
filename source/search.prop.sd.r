# Implements heuristic search of Metropolis proposal sd
search.prop.sd = function(m, mc = 3000, succ.min = 3, acc.bounds.X =c(0.2,0.3), acc.bounds.S =c(0.3,0.4)){
  found.X = found.S = F
  it = 1; succ = 0;
  while(succ!=succ.min){
    print(paste('Iteration',it))
    if(it == 1) { ac.X.cur = mean(m$ac.X) 
    prop.sd.X = m$prop.sd.X}
    if(it >1) {
      m = bts_survreg( prev.run = m, mc = mc, prop.sd.X = prop.sd.X, prop.sd.S = NULL)
      ac.X.cur = mean(m$ac.X.cur)
    }
    acc.bounds.mean = (acc.bounds.X[2]-acc.bounds.X[1])/2 + acc.bounds.X[1]
    ac.X = mean(m$ac.X)
    print(c(ac.X, ac.X.cur))
    ac.X = ac.X.cur
    if((ac.X > acc.bounds.X[1] & ac.X < acc.bounds.X[2]) ){ found.X = T ; print(paste('found.X:',prop.sd.X))} else{
      if(ac.X < acc.bounds.X[1]) {
        dif =  1-(acc.bounds.X[1] - ac.X)/acc.bounds.X[1]
        prop.sd.X=prop.sd.X * dif
      }
      if(ac.X > acc.bounds.X[2]) {
        dif =  1+(ac.X - acc.bounds.X[2])/acc.bounds.X[2]
        prop.sd.X=prop.sd.X * dif
      }
      found.X = F
      print(paste('prop.sd.X now', prop.sd.X))
    }
    it=it+1
    if(found.X ){ mc = mc*2; succ = succ +1; found.X = found.S = F
    print(paste('Doubling number of MC iterations:',mc))}
  }
  ret= list()
  ret$prop.sd.X = prop.sd.X
  ret$ac.X      = ac.X
  ret
}

# For use in the simulation only
search.prop.sd_seq = function(m, mc = 3000, succ.min = 3, acc.bounds.X =c(0.2,0.3), acc.bounds.S =c(0.3,0.4)){
  found.X = found.S = F
  it = 1; succ = 0;
  while(succ!=succ.min){
    print(paste('Iteration',it))
    if(it == 1) { ac.X.cur = mean(m$ac.X) 
    prop.sd.X = m$prop.sd.X}
    if(it >1) {
      m = bts_survreg_seq( prev.run = m, mc = mc, prop.sd.X = prop.sd.X, prop.sd.S = NULL)
      ac.X.cur = mean(m$ac.X.cur)
    }
    acc.bounds.mean = (acc.bounds.X[2]-acc.bounds.X[1])/2 + acc.bounds.X[1]
    ac.X = mean(m$ac.X)
    print(c(ac.X, ac.X.cur))
    ac.X = ac.X.cur
    if((ac.X > acc.bounds.X[1] & ac.X < acc.bounds.X[2]) ){ found.X = T ; print(paste('found.X:',prop.sd.X))} else{
      if(ac.X < acc.bounds.X[1]) {
        dif =  1-(acc.bounds.X[1] - ac.X)/acc.bounds.X[1]
        prop.sd.X=prop.sd.X * dif
      }
      if(ac.X > acc.bounds.X[2]) {
        dif =  1+(ac.X - acc.bounds.X[2])/acc.bounds.X[2]
        prop.sd.X=prop.sd.X * dif
      }
      found.X = F
      print(paste('prop.sd.X now', prop.sd.X))
    }
    it=it+1
    if(found.X ){ mc = mc*2; succ = succ +1; found.X = found.S = F
    print(paste('Doubling number of MC iterations:',mc))}
  }
  ret= list()
  ret$prop.sd.X = prop.sd.X
  ret$ac.X      = ac.X
  ret
}
