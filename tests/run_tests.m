function run_tests
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'examples'));
m=support_metrics([1 1 0 0],[1 0 1 0]);assert(max(abs(m-[0.5 0.5 0.5]))<1e-12);
m=support_metrics([1 1 0 0],[1 1 0 0]);assert(isequal(m,[1 0 1]));
m=support_metrics([1 1 0 0],[0 0 0 0]);assert(isequal(m,[0 0 0]));
A=eye(16);Y=zeros(16,3);Y([2 6 11],:)=repmat([1+1i;2-1i;3+2i],1,3);
[X,support]=fct_MRMP(A,A',Y,3,0.65,0.7);
assert(norm(X-Y,'fro')<1e-10);assert(isequal(support,[2;6;11]));
[Z,support]=fct_MRMP(A,A',zeros(16,3),3,0.65,0.7);assert(~any(Z(:))&&isempty(support));
X=somp_baseline(A,Y,3);assert(norm(X-Y,'fro')<1e-10);
% Independent-noise VDS variance and coherent-signal preservation.
rng(7);W=randn(10000,4);V=(ones(4)-eye(4))/3;
ratio=var(W*V,0,1)./var(W,0,1);assert(max(abs(ratio-1/3))<0.03);
assert(norm(ones(10,4)*V-ones(10,4),'fro')<1e-12);
fprintf('All numerical checks passed.\n');
end
