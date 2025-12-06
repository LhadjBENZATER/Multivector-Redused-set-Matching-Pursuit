function X_hat = fct_MRMP_Virt(A, At, Y, K, bk, a)
N=size(A,2); L=(size(Y,2)-1)/2; Sp_r=[];Sp_v=[];X_hat=zeros(N,2*L+1); 
Y_r = Y(:,1:L); Y_v = Y(:, L+1:2*L); Y_avg = Y(:,end);  r_r = Y_r;  
r_v = Y_v; er_r = norm(r_r);  er_v = norm(r_v);
    for iter = 1:K
cr_r=abs(At*r_r);  cr_r=cr_r./max(cr_r,[],1);  
cr_v=abs(At*r_v);  cr_v=cr_v./max(cr_v,[],1);  
[~,sor_idr]=sort(cr_r(:),'descend');  [~,sor_idv]=sort(cr_v(:),'descend');
top_r=sor_idr(1:round(bk*K)); top_v=sor_idv(1:round(bk*K));  
[ind_r,~]=ind2sub(size(cr_r),top_r);[ind_v,~]=ind2sub(size(cr_v),top_v);
ind_r=ind_r(cr_r(ind_r)>a); ind_v=ind_v(cr_v(ind_v)>a);  
Sp_v=unique([Sp_v;ind_v]); Sp_r=unique([Sp_r; ind_r]);  
Supp = unique([Sp_r;Sp_v]);


X_hat(Supp,1:L)=A(:,Supp)\Y_r;X_hat(Supp, L+1:2*L)= A(:,Supp)\Y_v;
X_hat(Supp, end)=A(:,Supp)\Y_avg; 
temp_X_hat=zeros(N,2*L+1);
    for i = 1:2*L+1
    [~, X_ind] = sort(abs(X_hat(:, i)), 'descend');
    temp_X_hat(X_ind(1:K), i) = X_hat(X_ind(1:K), i);
    end
X_hat = temp_X_hat;
r_r=Y_r-A*X_hat(:, 1:L);
r_v=Y_v-A*X_hat(:,L+1:2*L);
if (er_r - norm(r_r) < 1e-6) && (er_v - norm(r_v) < 1e-6)
    break;
end
er_r = norm(r_r); er_v = norm(r_v);
    end
X_hat=X_hat(:,1:L);
end
