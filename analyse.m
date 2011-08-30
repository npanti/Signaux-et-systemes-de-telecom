format long;


%% EbNo
% EbNo=0:0.5:16;
% SNRc = zeros(1,length(EbNo));
% SNRo = zeros(1,length(EbNo));
% FOM = zeros(1,length(EbNo));
% BER_FSK = zeros(1,length(EbNo));
% BER_FM = zeros(1,length(EbNo));
% 
% for k=1:length(EbNo)
%     [SNRc_tmp,SNRo_tmp,FOM_tmp,BER_FSK_tmp,BER_FM_tmp] =  modul_fm_analyse(EbNo(k));
%     SNRc(k) = SNRc_tmp;
%     SNRo(k) = SNRo_tmp; 
%     FOM(k) = FOM_tmp; 
%     BER_FSK(k) = BER_FSK_tmp;
%     BER_FM(k) = BER_FM_tmp;
% end
% figure(1);
% %plot(EbNo,SNRc,'*-','Color','Green');
% hold on; grid on;
% %plot(EbNo,SNRo,'*-','Color','Red');
% plot(EbNo,FOM,'*-');
% xlabel('Eb/No');
% ylabel('[dB]');
% %legend('SNRc','SNRo','FOM');
% hold on; grid on;
% 
% figure(2);
% plot(EbNo,BER_FM,'*-');
% hold on; grid on;
% plot(EbNo,BER_FSK,'*-','Color','Green');
% xlabel('Eb/No');
% ylabel('BER');
% legend('FM','FSK');
% 
% clear; clc;

%% f

parametre=100:100:2000;
SNRc = zeros(1,length(parametre));
SNRo = zeros(1,length(parametre));
FOM = zeros(1,length(parametre));
BER_FSK = zeros(1,length(parametre));
BER_FM = zeros(1,length(parametre));

for k=1:length(parametre)
    [SNRc_tmp,SNRo_tmp,FOM_tmp,BER_FSK_tmp,BER_FM_tmp] =  modul_fm_analyse(parametre(k));
    SNRc(k) = SNRc_tmp;
    SNRo(k) = SNRo_tmp; 
    FOM(k) = FOM_tmp; 
    BER_FSK(k) = BER_FSK_tmp;
    BER_FM(k) = BER_FM_tmp;
end
figure(1);
%plot(EbNo,SNRc,'*-','Color','Green');
hold on; grid on;
%plot(EbNo,SNRo,'*-','Color','Red');
plot(parametre,FOM,'*-');
xlabel('f');
ylabel('FOM [dB]');
%legend('SNRc','SNRo','FOM');
hold on; grid on;

figure(2);
plot(parametre,BER_FM,'*-');
hold on; grid on;
plot(parametre,BER_FSK,'*-','Color','Green');
xlabel('f');
ylabel('BER');
legend('Démodulateur FM','Démodulateur FSK');

clear; clc;

%% fc

% parametre=3000:100:9000;
% SNRc = zeros(1,length(parametre));
% SNRo = zeros(1,length(parametre));
% FOM = zeros(1,length(parametre));
% BER_FSK = zeros(1,length(parametre));
% BER_FM = zeros(1,length(parametre));
% 
% for k=1:length(parametre)
%     [SNRc_tmp,SNRo_tmp,FOM_tmp,BER_FSK_tmp,BER_FM_tmp] =  modul_fm_analyse(parametre(k));
%     SNRc(k) = SNRc_tmp;
%     SNRo(k) = SNRo_tmp; 
%     FOM(k) = FOM_tmp; 
%     BER_FSK(k) = BER_FSK_tmp;
%     BER_FM(k) = BER_FM_tmp;
% end
% figure(1);
% %plot(EbNo,SNRc,'*-','Color','Green');
% hold on; grid on;
% %plot(EbNo,SNRo,'*-','Color','Red');
% plot(parametre,FOM,'*-');
% xlabel('f');
% ylabel('FOM [dB]');
% %legend('SNRc','SNRo','FOM');
% hold on; grid on;
% 
% figure(2);
% plot(parametre,BER_FM,'*-');
% hold on; grid on;
% plot(parametre,BER_FSK,'*-','Color','Green');
% xlabel('f');
% ylabel('BER');
% legend('FM','FSK');
% 
% clear; clc;

% %% phi
% 
% parametre=0:2*pi/50:2*pi;
% SNRc = zeros(1,length(parametre));
% SNRo = zeros(1,length(parametre));
% FOM = zeros(1,length(parametre));
% BER_FSK = zeros(1,length(parametre));
% BER_FM = zeros(1,length(parametre));
% 
% for k=1:length(parametre)
%     [SNRc_tmp,SNRo_tmp,FOM_tmp,BER_FSK_tmp,BER_FM_tmp] =  modul_fm_analyse(parametre(k));
%     SNRc(k) = SNRc_tmp;
%     SNRo(k) = SNRo_tmp; 
%     FOM(k) = FOM_tmp; 
%     BER_FSK(k) = BER_FSK_tmp;
%     BER_FM(k) = BER_FM_tmp;
% end
% figure(1);
% %plot(EbNo,SNRc,'*-','Color','Green');
% hold on; grid on;
% %plot(EbNo,SNRo,'*-','Color','Red');
% plot(parametre,FOM,'*-');
% xlabel('\phi');
% ylabel('FOM [dB]');
% %legend('SNRc','SNRo','FOM');
% hold on; grid on;
% 
% figure(2);
% plot(parametre,BER_FM,'*-');
% hold on; grid on;
% plot(parametre,BER_FSK,'*-','Color','Green');
% xlabel('\phi');
% ylabel('BER');
% legend('FM','FSK');
% 
% clear; clc;