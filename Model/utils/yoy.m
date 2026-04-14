% yoy

function y = yoy(x)



T = length(x) ;



y=x*0;



for t = 1 : T



if t==1; y(t)=x(t); end
if t==2; y(t)=x(t)+x(t-1); end
if t==3; y(t)=x(t)+x(t-1)+x(t-2); end
if t>=4; y(t)=x(t)+x(t-1)+x(t-2)+x(t-3); end



end



end