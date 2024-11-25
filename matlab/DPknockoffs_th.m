% Power/FDR Theorical value eastimation
function [fdr_est_th,fdr_th,power_th]= DPknockoffs_th(n,p,mu,r,eps,del,b,T,lambda,gamma,c0)

B = sqrt(b^2*p); % bound on row norm of X
w_sq = B^2*(1+(2/eps+2/log(4/del))*(sqrt(2*r*log(4/del))+log(4/del)));
eta = @(x) 1/(1+w_sq/(gamma*n))*sign(x).*max((abs(x)-lambda/gamma),0);
rho = sqrt((n+r+1)/(n*r*gamma)+w_sq/(gamma^2*n*r));

% monte carlo 
P1 = 0;
P2 = 0;
P3 = 0;
mc_iter = 10000;
for i = 1: mc_iter
    Z = randn(1,6);
    if abs(eta(mu+rho*Z(1)))-abs(eta(rho*Z(2)))<=-T
        P1 = P1+1;
    end
    if abs(eta(mu+rho*Z(3)))-abs(eta(rho*Z(4)))>=T
        P2 = P2+1;
    end
    if abs(eta(rho*Z(5)))-abs(eta(rho*Z(6)))<=-T
        P3 = P3+1;
    end
end
P1 = P1/mc_iter;
P2 = P2/mc_iter;
P3 = P3/mc_iter;
if (c0*P2+(1-c0)*P3)>0
    fdr_est_th = (c0*P1+(1-c0)*P3)/(c0*P2+(1-c0)*P3);
    fdr_th = (1-c0)*P3/(c0*P2+(1-c0)*P3);
else
    fdr_est_th = 0;
    fdr_th = 0;
end

power_th = P2;
end