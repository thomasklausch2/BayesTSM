# Helper functions for the E step
Ex.int = function(x, left, right, d, cur.par.S, cur.par.X, px.int){
  p = px(x, left = left, right= right, d = d, 
     cur.par.S = cur.par.S, cur.par.X = cur.par.X)
  p*x/px.int}
Ex2.int = function(x, left, right, d, cur.par.S, cur.par.X, px.int){
  p = px(x, left = left, right= right, d = d, 
     cur.par.S = cur.par.S, cur.par.X = cur.par.X)
  p*(x^2)/px.int}
