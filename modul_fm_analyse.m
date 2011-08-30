function [SNR_c,SNR_o,FOM,BER_FSK,BER_FM] = modul_fm_analyse(parametre)


%% Initialisation
 
bitsSend = 10000;   %Nombre de bits émis [1000]
%OSF = f_s/f;        %Oversampling factor [40]
OSF = 40;
f = parametre;      %Fréquence symbole [500 Hz]
f_s = OSF*f;        %Fréquence d'échantillonage [20000 Hz]
T_s = 1/f_s;        %Période d'échantillonage
k_f = f;          %Sélectivité fréquencielle [500 Hz]
f_c = 6000;         %Fréquence porteuse [6000 Hz]

%Condition sur OSF
% l = f_c*T_s/0.5;
% if mod(l,2) == 0 || mod(1,2) == 1
%     while f_c*T_s > 0.46*l
%         fprintf('bla');
%         OSF = OSF + 1;
%         f_s = OSF*f;
%         T_s = 1/f_s;
%     end
% end
if f == 100
    OSF = 131;
elseif f == 200;
    OSF = 66;
elseif f == 300
    OSF = 44;
end
    f_s = OSF*f;
    T_s = 1/f_s;

t = 0:T_s:OSF*T_s*bitsSend-T_s; %Temps d'échantillonage
phi = 0;            %Déphasage
EbN0_ratio = 10;    %defaut 10

rate = 10;
n_T = 40;
bits = randi(2,1, bitsSend) - 1; %Séquence aléatoire

%% Création du message
m = ones(1,bitsSend*OSF);

k=1;
for i=1:bitsSend
    if bits(i) == 0
        for j=1:OSF
            m(k+j-1) = -1;
        end   
    end
    k = k+OSF;
end

%% Modulation FM

e_s = exp(1i*2*pi*k_f.*cumsum(m)*T_s);

%% Singal RF

s = real(e_s).*cos(2*pi*f_c*t) - imag(e_s).*sin(2*pi*f_c*t);

%% Bruit blanc

%Paramètres du bruit
Ps = sum(s.^2)/(bitsSend*OSF);
Eb = Ps/f;
N0 = Eb/(10^(EbN0_ratio/10));
sigma_n = sqrt(N0*f_s/2);

%Génération du bruit
n = sigma_n*randn(1,bitsSend*OSF);

%Ajout du bruit
r = s + n;

%% Filtre passe-bas

%Réponse impulsionnel du filtre passe bas
h = rcosfir(0,n_T,rate,1/(3*f_c));
enlever = (length(h)-1)/2;

%Filtre appliqué juste sur le message sans bruit
r_cos = s.*cos(2*pi*f_c*t+phi);
r_sin = s.*sin(2*pi*f_c*t+phi);

e_r = conv(r_cos,h) + 1i*conv(r_sin,h);
e_r = e_r(1+enlever:end-enlever);

%Filtre appliqué sur le message avec bruit
r_cos = r.*cos(2*pi*f_c*t+phi);
r_sin = r.*sin(2*pi*f_c*t+phi);

e_r_bruit = conv(r_cos,h) + 1i*conv(r_sin,h);
e_r_bruit = e_r_bruit(1+enlever:end-enlever);

%% Démodulateur FSK

e_s0 = exp(-1i*2*pi*k_f*t);
e_s1 = exp(1i*2*pi*k_f*t);

m_fsk = zeros(1,bitsSend);
for k=1:bitsSend
    
    s0_sum = 0;
    s1_sum = 0;
     
    for i=1:OSF
        s0_sum = s0_sum + e_r_bruit((k-1)*OSF+i)*conj(e_s0((k-1)*OSF+i));
        s1_sum = s1_sum + e_r_bruit((k-1)*OSF+i)*conj(e_s1((k-1)*OSF+i));
    end
    
    if abs(s0_sum) > abs(s1_sum)
        m_fsk(k) = 1;
    end
end

%% Démodulateur FM
g = 1/(2*pi*k_f*rate);
B_T = 2*k_f + 2*f;

%Démodulation sans bruit
d_er = gradient(e_r,T_s);  %derivee par rapport au temps (1x400)

s1 = g.*(d_er+1i*pi*B_T.*e_r);    % apres 2 filtres a pente
s2 = g.*(-d_er+1i*pi*B_T.*e_r);    % apres 2 filtres a pente 

m_r = abs(s2)-abs(s1);

m_fm = zeros(1,bitsSend*OSF);
for k=1:bitsSend*OSF
    if m_r(k) >= 0
        m_fm(k) = 1;
    end
end

%Démodulation avec bruit
d_er_bruit = gradient(e_r_bruit,T_s);  %derivee par rapport au temps (1x400)

s1 = g.*(d_er_bruit+1i*pi*B_T.*e_r_bruit);    % apres 2 filtres a pente
s2 = g.*(-d_er_bruit+1i*pi*B_T.*e_r_bruit);    % apres 2 filtres a pente

m_r_bruit = abs(s2)-abs(s1);

m_fm_bruit = zeros(1,bitsSend*OSF);
for k=1:bitsSend*OSF
    if m_r_bruit(k) >= 0
        m_fm_bruit(k) = 1;
    end
end

% Interprétation du message avec bruit
m_fm_final = zeros(1,bitsSend);
for k=1:bitsSend

    for l=1:OSF
        m_fm_final(k) = m_fm_final(k) + m_fm_bruit((k-1)*OSF+l);
    end
    m_fm_final(k) = round(m_fm_final(k)/OSF);
end

% j=1;
% for k=round(OSF/2):OSF:OSF*bitsSend-round(OSF/2)
%     if(m_fm_bruit(k) > 0)
%         m_fm_final(j) = 1;
%     end
%     j=j+1;
% end

%% BER

incorrectBits = 0;
for k=1:bitsSend
   if(m_fm_final(k) ~= bits(k))
       incorrectBits = incorrectBits + 1;
   end
end

BER_FM = incorrectBits/bitsSend;

incorrectBits = 0;
for k=1:bitsSend
   if(m_fsk(k) ~= bits(k))
       incorrectBits = incorrectBits + 1;
   end
end

BER_FSK = incorrectBits/bitsSend;

%% PSD
% L = 512;
% D = 250;
% [PSD_m, PSD_freq] = welch(m,t,L,D);            %PSD du message
% PSD_es = welch(e_s,t,L,D);                     %PSD du signal modulé en FM
% PSD_s = welch(s,t,L,D);                        %PSD du signal RF
% %PSD_er = welch(e_r,t,L,D);                     %PSD du signal après le FPB
% %PSD_er = welch(m_fm,t,L,D);                    %PSD du singal démodulé FM
% 
% plot(PSD_freq, PSD_m,'Color','Green');
% grid on;
% hold on;
% plot(PSD_freq, PSD_es);
% plot(PSD_freq, PSD_s, 'Color', 'Red');
% legend('Message','Signal FM','Signal RF');
% xlabel('Fréquence [Hz]');
% ylabel('PSD');

%% SNR & FOM

Ps_c = sum(s.^2)/(OSF*bitsSend);
Ps_o = sum(m_r.^2)/(OSF*bitsSend);
Pn_o = sum((m_r_bruit - m_r).^2)/(OSF*bitsSend);
Pn_c = N0*f*OSF/2;

SNR_c = 10*log10(Ps_c/Pn_c);
SNR_o = 10*log10(Ps_o/Pn_o);
FOM = SNR_o - SNR_c; 
