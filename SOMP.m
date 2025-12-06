function X_hat = SOMP(Y,A,At,K)
L=size(Y,2); N=size(A,2); r=Y; X_hat=zeros(N,L);sup=[];
for k = 1:K
cr=At*r; [~,idx]=max(sum(cr.^ 2,2));  
sup = [sup; idx];A_sup=A(:,sup); X_est=A_sup\Y; 
r=Y-A_sup*X_est;
if norm(r,'fro') < 1e-6
    break;
end
end
X_hat(sup, :) = X_est;
end
