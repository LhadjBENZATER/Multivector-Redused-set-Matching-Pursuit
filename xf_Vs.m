set(0,'DefaultAxesFontSize',14);
plot(f,xf,'LineWidth',1.5);grid on; hold on
plot(f,Xf_hat(:,2),'LineWidth',2,'LineStyle','--');
plot(f,Xf_hat(:,4),'LineWidth',2,'LineStyle','-.')
legend('Original spectrum','MRMP','SOMP','Location','best');
ylabel('Amplitude (V/Hz)'); xlabel('Frequency (GHz)'); 