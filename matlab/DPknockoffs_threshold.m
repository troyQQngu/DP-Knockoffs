% choose threshold t based on target fdr q
function T = DPknockoffs_threshold(n,p,mu,r,eps,del,b,lambda,gamma,s0,q)

% initialization
c0 = s0/p;
B = sqrt(b^2*p); % bound on row norm of X
w_sq = B^2*(1+(2/eps+2/log(4/del))*(sqrt(2*r*log(4/del))+log(4/del)));

eta = @(x) 1/(1+w_sq/(gamma*n))*sign(x).*max((abs(x)-lambda/gamma),0);
rho = sqrt((n+r+1)/(n*r*gamma)+w_sq/(gamma^2*n*r));

mc_size = 100000; % monte carlo size

% simulate unbiased theta
s0 = floor(c0*mc_size);

signals = rho*randn(mc_size,1)+[mu*ones(s0,1);zeros(mc_size-s0,1)];
knockoffs = rho*randn(mc_size,1);

W = abs(eta(signals))-abs(eta(knockoffs));

t = sort( [0;abs(W(W~=0))]);
ratio = zeros(1,length(t));
for i = 1:length(t)
    ratio(i) = sum(W<=-t(i))/max(1,sum(W>=t(i)));
end

index = find(ratio <=q,1,'first');

if isempty(index)
    T = inf;
else
    T = t(index);
end

end