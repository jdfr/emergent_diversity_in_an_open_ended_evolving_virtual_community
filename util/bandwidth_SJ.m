function [hout]=bandwidth_SJ(xi,K)
%
% Bandwidth selection with Sheater and Jones(1991) 
% and Wand and Jones (1995) pp.74-75

syms uh
if K=='norm'
    fk=1/sqrt(2*pi)*exp(-0.5*uh^2);
    interK1=-inf;interK2=inf;
elseif K=='biwt'    
    fk=15/16*(1-uh^2)^2;
    interK1=-1;interK2=1;
elseif K=='expn'
    fk=exp(-uh);
    interK1=0;interK2=inf;
end

fk_der4_0=subs(diff(fk,4),0);
fk_der6_0=subs(diff(fk,6),0);
  
RK=int(fk^2,interK1,interK2);
RKr=double(RK);
mu2=int(uh^2*fk,interK1,interK2);
mu2r=double(mu2);
nL=length(xi);
%Psi6_NS=-15/(16*pi^0.5*std(xi)^7);
Psi8_NS=105/(32*pi^0.5*std(xi)^9);
g1=(-2*fk_der6_0/(mu2r*Psi8_NS*nL))^(1/9);
Psi6=f_Psi_Lrg(fk,6,g1,xi,interK1,interK2);

g2=(-2*fk_der4_0/(mu2r*Psi6*nL))^(1/7);
Psi4=f_Psi_Lrg(fk,4,g2,xi,interK1,interK2);

hout=(RKr/(mu2r^2*Psi4*nL))^(1/5);


function [outPsi]=f_Psi_Lrg(L,r,g,xx,interK1,interK2)
fL_r_der=diff(L,r);
nL=length(xx);
sum1=0;
for iL1=1:nL
        uu1=(xx-xx(iL1))./g;
        uu2=uu1(uu1~=0 & (uu1>interK1 & uu1<interK2));
        sum1=sum1+sum(subs(fL_r_der,uu2));
end

outPsi=(nL*(nL-1))^(-1)*g^(-r-1)*sum1;
%outPsi=(nL*(nL))^(-1)*sum1;
