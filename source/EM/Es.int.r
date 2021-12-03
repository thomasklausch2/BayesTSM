# Helper functions for the E step
Es.int = function(s, left, right, d, cur.par.S, cur.par.X, ps.int){
  p = ps(s, left = left, right= right, d = d, 
     cur.par.S = cur.par.S, cur.par.X = cur.par.X)
  p *s/ps.int}
Es2.int = function(s, left, right, d, cur.par.S, cur.par.X, ps.int){
  p = ps(s, left = left, right= right, d = d, 
     cur.par.S = cur.par.S, cur.par.X = cur.par.X)
  p*(s^2)/ps.int}

