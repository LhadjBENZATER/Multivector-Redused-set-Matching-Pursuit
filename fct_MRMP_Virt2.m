function X_hat = fct_MRMP_Virt2(A, At, Y, K, bk, a_cor)
N=size(A,2); [M,L]=size(Y); Supp=[]; Y_v1=Y*~eye(L)/(L-1);
C2=nchoosek(1:L,L-2); L2=size(C2,1); Y_v2=zeros(M, L2);
for i = 1:L2
    Y_v2(:,i) = mean(Y(:,C2(i,:)),2);
end
Y_virt=[Y_v1 (L-2)/(L-1)*Y_v2];r=Y_virt;er=norm(r);X_hat=zeros(N,L+L2);
for iter = 1:K
cr=abs(At*r);  cr=cr./max(cr,[],1); 
[~, top_k] = maxk(cr(:), round(bk * K)); [ind,~]=ind2sub(size(cr),top_k);
Supp=unique([Supp;ind(cr(ind)>a_cor)]);
X_hat(Supp,:)=lsqminnorm(A(:,Supp),Y_virt);
[~,idx]=maxk(abs(X_hat),K,1);
msk=false(size(X_hat));
msk(sub2ind(size(X_hat),idx(:),repelem(1:size(X_hat, 2),K)')) = true;
X_hat(~msk)=0;

r = Y_virt - A * X_hat;
Supp = find(sum(X_hat ~= 0, 2)); err=norm(r,'fro');
if abs(er-err)< 1e-6
    break;
end
er = err;
end
end
