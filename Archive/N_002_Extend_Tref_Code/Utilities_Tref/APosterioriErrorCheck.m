function [err2, maxerr, Z2, G2, wt2] = APosterioriErrorCheck(nw,tt,pt,g,w,scl,u,rel,err)
Z2 = []; G2 = [];
for k = 1:nw
    newtt = mean([tt{k}(1:end-1) tt{k}(2:end)],2);
    newpts = pt{k}(newtt);
    Z2 = [Z2; newpts];
    G2 = [G2; g{k}(newpts)];
end
M2 = length(Z2);
K2j = zeros(M2,1);
for j = 1:M2
    dd2 = abs(Z2(j)-w); K2j(j) = find(dd2==min(dd2),1); % nearest corner to Z2j
end
if rel
    wt2 = abs(Z2-w(K2j))/scl;                           % weights to measure error
else
    wt2 = ones(M2,1);
end
err2 = norm(wt2.*(G2-u(Z2)),inf);
maxerr = max(err,err2);                                % estimated max error

end