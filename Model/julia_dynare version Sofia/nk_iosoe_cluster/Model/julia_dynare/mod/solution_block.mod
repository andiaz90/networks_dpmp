shocks;
var eps_om=0.1;
var eps_i=0.1;
var epschi=0.0;
@#for z in 1:nsec
   var epsA_@{z}=.01;
@#endfor
end;
stoch_simul(order=1,periods=150,nograph);
