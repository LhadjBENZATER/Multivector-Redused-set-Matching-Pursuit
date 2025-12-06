%% Joint Sparse Recovery
clear; clc; close all; set(0,'DefaultAxesFontSize',16); N=1024; 
M=round(0.25*N); a=0.7; bk=0.65; Psi=fft(eye(N))/sqrt(N); tic;
%W=[0 1 0 1 1 0 0];[x,xf,index,K,f,t,ind_shft]=Sig(2.9e9,4e7,N,1,W);
[x,xf,index,K,f,t,ind_shft]=LFM(2.9e9,4e7,N,1);
Phi=Phi_selection(N,M,2,2e7); A=Phi/Psi; At=A'; Y = Phi*x;
th=0.1; L=[5,10,20]; SNR=-20:5:20; MC=1;l=length(L);tech=4;
snr=length(SNR); Pd=zeros(tech,l,snr); Pf=Pd; RMSE=Pd; ncor=zeros(l,snr);  
for s = 1:snr
fprintf('SNR(dB) = %d\n', SNR(s));  
Pd_loc=zeros(l,tech); Pf_loc=Pd_loc; RMSE_loc=Pd_loc;
for j = 1:l
for frame = 1:MC
y = Virtual(diff_chn(Y, L(j), SNR(s), 0)); Xf_hat = zeros(N,tech); 
Xf_hat(:,1)=prod(fftshift(abs(fct_MRMP(A,At,y(:,1:L(j)),...
      K,bk,a))),2).^(1/L(j));%Physical Ch
Xf_hat(:,2)=prod(fftshift(abs(fct_MRMP(A,At,y(:,L(j)+1:...
     2*L(j)),K,bk,a))),2).^(1/L(j));%MRMP L-1 Virt 
Xf_hat(:,3)=prod(fftshift(abs(SOMP(y(:,1:L(j)),A,At,K))),2).^(1/L(j));
Xf_hat(:,4)=prod(fftshift(abs(SOMP(y(:,1+L(j):2*L(j)),A,At,K))),2).^(1/L(j));
Xf_hat=Xf_hat./max(Xf_hat,[],1);
Pd_loc(j,:)=Pd_loc(j,:)+sum(Xf_hat(index,:)>th)/K;
Pf_loc(j,:)=Pf_loc(j,:)+(sum(Xf_hat>th)-sum(Xf_hat(index,:)>th))/(N-K);
RMSE_loc(j,:)=RMSE_loc(j,:)+sqrt(mean((xf-Xf_hat).^2, 1)); 
end
end
Pd(:,:,s)=Pd_loc'/MC; Pf(:,:,s)=Pf_loc'/MC; RMSE(:,:,s)=RMSE_loc'/MC;
end
disp(['Time=', datestr(now, 'HH:MM:SS'), ', Duration(min) = ', num2str(toc / 60)]);
[L_grid,SNR_grid] = meshgrid(L,SNR);
%Qd_R=zeros()
Qd_R(:,:)=Pd(1,:,:);Qf_R(:,:)=Pf(1,:,:);RMSE_R(:,:)=RMSE(1,:,:);
Qd_V(:,:)=Pd(2,:,:);Qf_V(:,:)=Pf(2,:,:);RMSE_V(:,:)=RMSE(2,:,:);
Qd_SOMP(:,:)=Pd(3,:,:);Qf_SOMP(:,:)=Pf(3,:,:);RMSE_SOMP(:,:)=RMSE(3,:,:);
Qd_VSOMP(:,:)=Pd(4,:,:);Qf_VSOMP(:,:)=Pf(4,:,:);RMSE_VSOMP(:,:)=RMSE(4,:,:);

subplot(1,3, 1);    
surf(L_grid, SNR_grid, Qd_R','FaceAlpha',0.2,'EdgeColor','blue','LineWidth',2);hold on; grid on;    
surf(L_grid, SNR_grid, Qd_V','LineStyle','--','FaceColor','none','FaceAlpha',0.3,'EdgeColor', 'blue','LineWidth',2);
surf(L_grid, SNR_grid, Qd_SOMP','FaceAlpha', 0.2,'EdgeColor','green','LineWidth',2);
surf(L_grid, SNR_grid, Qd_VSOMP','LineStyle','--','FaceColor','none','FaceAlpha', 0.3,'EdgeColor','green','LineWidth',2);
zlabel('Recovered Support'); ylabel('SNR (dB)'); xlabel('L'); 
set(gca,'XScale','log');
subplot(1,3,2);
surf(L_grid, SNR_grid, Qf_R','FaceAlpha',0.2,'EdgeColor','blue','LineWidth',2);hold on; grid on;    
surf(L_grid, SNR_grid, Qf_V','LineStyle','--','FaceColor','none','FaceAlpha',0.3,'EdgeColor', 'blue','LineWidth',2);
surf(L_grid, SNR_grid, Qf_SOMP','FaceAlpha', 0.2,'EdgeColor','green','LineWidth',2);
surf(L_grid, SNR_grid, Qf_VSOMP','LineStyle','--','FaceColor','none','FaceAlpha', 0.3,'EdgeColor','green','LineWidth',2);
legend('MRMP_{Re}','MRMP_{Virt}','SOMP_{Re}','SOMP_{Virt}', 'Location', 'best'); 
set(gca,'ZScale','log'); zlabel('Erroneous Support'); ylabel('SNR (dB)'); xlabel('L');
set(gca,'XScale','log');
subplot(1,3,3);
surf(L_grid, SNR_grid, RMSE_R','FaceAlpha',0.2,'EdgeColor','blue','LineWidth',2);hold on; grid on;    
surf(L_grid, SNR_grid, RMSE_V','LineStyle','--','FaceColor','none','FaceAlpha',0.3,'EdgeColor', 'blue','LineWidth',2);
surf(L_grid, SNR_grid, RMSE_SOMP','FaceAlpha', 0.2,'EdgeColor','green','LineWidth',2);
surf(L_grid, SNR_grid, RMSE_VSOMP','LineStyle','--','FaceColor','none','FaceAlpha', 0.3,'EdgeColor','green','LineWidth',2);
zlabel('RMSE'); ylabel('SNR (dB)'); xlabel('L'); set(gca, 'ZScale', 'log');set(gca, 'XScale', 'log');
set(gca,'XScale','log');
