function yinit = guess(x)
global infinity
  yinit = [ (x / infinity) ^ 0.5 
             0 
             0 %(x / infinity) ^ 0.5 
             1
             0 % 0.4%(x / infinity) ^ 0.5 
             0
             0
           ];  
          