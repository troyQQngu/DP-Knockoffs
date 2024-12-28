function [S_hat] = selection(W,q)
    t = sort([0;abs(W(W~=0))]);
    ratio = zeros(1,length(t));
    for j = 1:length(ratio)
        ratio(j) = (1+sum(W<=-t(j)))/sum(W>=t(j));
    end
    index = find(ratio<=q,1,'first');
    if isempty(index)
        T = Inf;
    else
        T = t(index);
    end

    % selection

    S_hat = find(W>=T); % selection
end