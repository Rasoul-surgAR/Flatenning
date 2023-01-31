function c = lsconformal(x,A,b)
  c =  (A * x - b)'* (A * x - b );
end
