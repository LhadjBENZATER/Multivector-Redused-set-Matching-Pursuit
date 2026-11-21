function result = demo_mrmp(makePlot)
% Self-contained synthetic MMV example using the supplied fct_MRMP solver.
if nargin<1, makePlot=true; end
root=fileparts(mfilename('fullpath')); addpath(root,fullfile(root,'examples'));
rng(20261121,'twister'); N=128; M=64; K=6; L=4; MC=20;
snrDb=[-15 -10 -5 0 5 10]; labels={'MRMP','SOMP','MRMP + VDS','SOMP + VDS'};
A=randn(M,N); A=A./sqrt(sum(abs(A).^2,1));
truth=false(N,1); truth([14:16 77:79])=true;
% Coherent equal-gain channels permit leave-one-out averaging without
% signal cancellation. Independent receiver AWGN is added after A*X.
x=zeros(N,1); x(truth)=exp(1i*linspace(0,pi/3,K)).';
X=repmat(x,1,L); clean=A*X; power=mean(abs(clean(:)).^2);
scores=zeros(numel(snrDb),4,3); elapsed=zeros(numel(snrDb),4);
for s=1:numel(snrDb)
 for trial=1:MC
  noise=sqrt(power/10^(snrDb(s)/10)/2)*(randn(M,L)+1i*randn(M,L));
  Y=clean+noise; Yv=Y*(ones(L)-eye(L))/(L-1);
  for method=1:4
   data=Y; if method>2, data=Yv; end
   timer=tic;
   if mod(method,2)==1, estimate=fct_MRMP(A,A',data,K,0.65,0.7);
   else, estimate=somp_baseline(A,data,K); end
   elapsed(s,method)=elapsed(s,method)+toc(timer)/MC;
   % Identical reporting policy for both solvers: top K row-l2 scores.
   % This reporting step does not alter the returned solver coefficients.
   energy=sum(abs(estimate).^2,2); [~,order]=sort(energy,'descend');
   decision=false(N,1); decision(order(1:K))=true;
   scores(s,method,:)=scores(s,method,:)+reshape(support_metrics(truth,decision),1,1,3)/MC;
  end
 end
end
result=struct('snrDb',snrDb,'metrics',scores,'seconds',elapsed,'seed',20261121,'trials',MC);
out=fullfile(root,'results'); if ~exist(out,'dir'), mkdir(out); end
rows=[]; for m=1:4, rows=[rows;snrDb(:),repmat(m,numel(snrDb),1),squeeze(scores(:,m,:)),elapsed(:,m)]; end
write_demo_csv(fullfile(out,'demo_metrics.csv'),'snr_db,method,pd,pf,f1,seconds_per_recovery',rows);
save(fullfile(out,'demo_results.mat'),'result');
fprintf('MRMP demo: N=%d M=%d K=%d L=%d, %d paired trials/SNR.\n',N,M,K,L,MC);
disp(array2table(squeeze(scores(end,:,:)),'VariableNames',{'Pd','Pf','F1'},'RowNames',labels));
if makePlot
 fig=figure('Visible','off','Color','w','Position',[100 100 1100 650]);
 for metric=1:4
  subplot(2,2,metric); hold on; grid on;
  for method=1:4
   values=elapsed(:,method); if metric<4, values=scores(:,method,metric); end
   plot(snrDb,values,'-o','LineWidth',1.5);
  end
  xlabel('Measurement SNR (dB)'); names={'P_d','P_f','F1','Seconds per recovery'}; ylabel(names{metric});
  if metric<4, ylim([0 1]); end
 end
 subplot(2,2,1); legend(labels,'Location','best');
 sgtitle('Synthetic MMV example: identical inputs and bin scoring; not paper-figure reproduction');
 print(fig,fullfile(root,'docs','preview.svg'),'-dsvg'); close(fig);
end
end
