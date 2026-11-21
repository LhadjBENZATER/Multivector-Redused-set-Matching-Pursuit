function X = somp_baseline(A,Y,K)
% SOMP baseline: l2 row-correlation, common support, restricted LS.
N=size(A,2); X=zeros(N,size(Y,2)); R=Y; support=[];
for it=1:K
 score=sum(abs(A'*R).^2,2); score(support)=-Inf;
 [~,idx]=max(score); support=[support;idx];
 coefficients=A(:,support)\Y;
 R=Y-A(:,support)*coefficients;
 if norm(R,'fro')<1e-10, break; end
end
X(support,:)=coefficients;
end
