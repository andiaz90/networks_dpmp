function [ y, maxyrel, tmaxy ] = sim_ar2(x,nper,iplot)

eps = x(3);
rho1 = x(1);
rho2 = x(2);

y(1) = eps  ;
y(2) = rho1*y(1) ;

for t=3:nper
y(t) = rho1*y(t-1) - rho2*y(t-2) ;
end

if iplot==1
    plot(y)
end

[ maxyrel, tmaxy ]= max(y/eps);


% disp('Maximum y')
% disp(max(y))
% disp(' ')
% disp('Maximum y relative to shock')
% disp(max(y)/eps)