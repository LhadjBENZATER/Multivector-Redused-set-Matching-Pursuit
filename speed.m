t=zeros(3,1);
for d=1:100
tic;fftshift(abs(fct_MRMP(A,At,y(:,1:L(j)),K,bk,a)));t(1)=t(1)+toc;
tic;fftshift(abs(SOMP(y(:,1:L(j)),A,At,K)));t(2)=t(2)+toc;
tic;
for i=1:L(j)
fftshift(abs(fct_MRMP(A,At,y(:,i),K,bk,a)));
end
t(3)=t(3)+toc;
end
t(2)/t(1)
t(3)/t(1)