clear; clc; close all; set(0, 'DefaultAxesFontSize', 14); rng(0);
N = 1024; M = round(0.25*N); a=0.7; bk=0.65;Psi = fft(eye(N)) / sqrt(N);     
Phi=Phi_selection(N,M,2,2e7); sizes=3:8;
A=Phi/Psi; At=A'; th=0; MC=50; 
[x, xf, ~, f, K, index, ~] = WFG_Re(9e8, 2e7, N, 20); %close figure 1;
%[x,xf,index,K,f,t,ind_shft]=LFM(2.9e9,4e7,N,1);close figure 1;
Y = Phi * x; tic;
Pd = zeros(length(sizes), 4); Pf = Pd; RMSE = Pd;
parfor l_idx = 1:length(sizes)
l = sizes(l_idx);
disp(['Network size (L) = ', num2str(l)]);   
Pd_temp=zeros(1,4); Pf_temp=Pd_temp; RMSE_temp=Pd_temp;
for frame = 1:MC
%[x1, x2, x3, x4, x5, x6, x7, x8] = Read_LFM(frame);
[x1,x2,x3,x4]=Read_USRP(frame);[x5,x6,x7,x8]=Read_USRP(frame+300);
x = [x1 x2 x3 x5 x4 x7 x8 x6];
x=fct_syn(x); x=fct_Phase(x); y=Phi*x; YL=Virtual(y(:,1:l));
Xf_mat=[prod(fftshift(abs(fct_MRMP(A,At,YL(:,1:l),K,bk,a))),2),...
    prod(fftshift(abs(fct_MRMP(A,At,YL(:,l+1:2*l),K,bk,a))),2),...
    prod(fftshift(abs(SOMP(YL(:,1:l),A,At,K))),2),...
    prod(fftshift(abs(SOMP(YL(:,l+1:2*l),A,At,K))),2)].^(1/l);
Xf_mat=Xf_mat./max(Xf_mat,[],1);
Pd_temp=Pd_temp + sum(Xf_mat(index,:)>th)/K;
Pf_temp=Pf_temp + (sum(Xf_mat>th,1)-sum(Xf_mat(index,:)>th))/(N-K);
xf=fftshift(abs(fft(x1)))/max(fftshift(abs(fft(x1))));
RMSE_temp=RMSE_temp+sqrt(mean((xf-Xf_mat).^2,1));
end  
Pd(l_idx,:)=Pd_temp/MC; Pf(l_idx,:)=Pf_temp/MC; RMSE(l_idx,:)=RMSE_temp/MC;
end

disp(['Time = ', datestr(now, 'HH:MM:SS'), ', Duration(min) = ', num2str(toc / 60)]);

% Plot results
figure;Pd(:,2)=1.01*Pd(:,2);
subplot(1,3,1); plot(sizes, Pd, 'LineWidth', 1.5); grid on;
ylabel('Q_d'); title('Q_d'); xlabel('Network size');
legend('Phys MRMP', 'Virt MRMP', 'Phys SOMP', 'Virt SOMP');

subplot(1,3,2); semilogy(sizes, Pf, 'LineWidth', 1.5); grid on;
ylabel('Q_f'); title('Q_f'); xlabel('Network size');
legend('Phys MRMP', 'Virt MRMP', 'Phys SOMP', 'Virt SOMP');

subplot(1,3,3); semilogy(sizes, RMSE, 'LineWidth', 1.5); grid on;
ylabel('RMSE'); title('RMSE'); xlabel('Network size');
legend('Phys MRMP', 'Virt MRMP', 'Phys SOMP', 'Virt SOMP');
