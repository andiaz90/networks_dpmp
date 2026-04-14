function [U,U_g,U_s,U_gseb,leontief] = upstream(modbeta,modalpha,goods,M_,oo_)

nsec_val = size(modbeta,1);

for i=1:M_.endo_nbr
   eval([ deblank(char(M_.endo_names(i,:))) '_ss_lev = exp(oo_.steady_state(i)); ']);
end

for i=1:M_.param_nbr
    eval([ deblank(char(M_.param_names(i,:))) ' = M_.params(i);']);
end


m_mat=zeros(nsec_val,nsec_val);
for i=1:nsec_val
   for j=1:nsec_val
   m_mat(j,i)=modbeta(j,i)*eval(sprintf('(PM_%d_ss_lev/P_%d_ss_lev)^(epsM_%d) *M_%d_ss_lev',j,i,i,j));   
   end
end


cg_vec=zeros(nsec_val,1);
cs_vec=zeros(nsec_val,1);
y_vec=zeros(nsec_val,1);

for i=1:nsec_val 
y_vec(i)=eval(sprintf('Y_%d_ss_lev ',i));
gammag(i)=eval(sprintf('gammag_%d',i));
gammas(i)=eval(sprintf('gammas_%d',i));
end


for i=1:nsec_val 
cg_vec(i)=gammag(i)*C_g_ss_lev*p_g_ss_lev / eval(sprintf('P_%d_ss_lev ',i));
cs_vec(i)=gammas(i)*C_s_ss_lev*p_s_ss_lev / eval(sprintf('P_%d_ss_lev ',i));
end


c_vec=cs_vec+cg_vec;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPSTREAMNESS AND DOWNSTREAMNESS Antras et al. (2012) method
m_mat_t=m_mat';  %adjust matrix to be consistent with Antras forumla
d_mat=zeros(nsec_val,nsec_val);

for i=1:nsec_val
d_mat(i,:)=m_mat_t(i,:)./y_vec';        
end

%Compute upstreamness with the two different formulas in Antras (they give
%same result)
U=( ((eye(nsec_val)-d_mat)^-2)*c_vec)./y_vec;
U_g=( ((eye(nsec_val)-d_mat)^-2)*cg_vec)./y_vec;
U_s=( ((eye(nsec_val)-d_mat)^-2)*cs_vec)./y_vec;



leontief = inv(eye(nsec_val)-modbeta'.*modalpha);
% % U_g2 = leontief'*(goods);
% % U_g2 = (leontief'*(y_vec.*goods))./(leontief'*y_vec);
% % U_g2 = (leontief'*( goods))./(leontief'*ones(Nsec,1));
% U_gseb = (leontief'*(c_vec.*goods))./(leontief'*c_vec);
U_gseb = (modbeta*(c_vec.*goods))./(modbeta*c_vec);


end

