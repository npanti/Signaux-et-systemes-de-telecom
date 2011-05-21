%% Initialisation
clc;
clear;

F = 500;            %Fréquence symbole
T = 1/F;
Fs = 20000;         %Fréquence d'échantillonnage Matlab
Ts = 1/Fs;
OSF = Fs/F;
bits_size = 5;   %Nombre de symbole émis
Kf = 500;           %Sélectivité fréquentielle
Fc = 6000;          %Fréquence porteuse
Tc = 1/Fc;
vector_size = OSF*bits_size;

t =0:Ts:OSF*Ts*bits_size-Ts;

%% Création des bits
bits = randi(2,1, bits_size) - 1;   %-1 pour avoir 0 ou 1 en aléatoire

%% Création du message
msg = [];   %Définition de notre matrisse message = m(t)

for i=1:bits_size
   
    bit = bits(1,i);
    if bit == 0         %Si le bit vaut 0 on écrit -1
        bit = -1;
    end
    
    bit_msg = bit;      %On duplique la valeur du bit OSF fois
    bit_msg = repmat(bit_msg,1,OSF);
    
    msg = horzcat(msg,bit_msg);
    
end

figure(1);
plot(msg,'-');

%% Modulation FM
%Es = fmmod(msg,Fc,Fs,50); ou modulate
% somme = spalloc(1,40000,40000);
% for n=1:40000
%     somme(1,n) = sum(msg(1:n))/Fs;
% end

%Cumsum est la somme cumulée, c'est l'intégrale discrète de notre message
e_s = exp(1i*2*pi*Kf*cumsum(msg)*Ts);

% figure(2);
% plot(t,real(e_s),'-'); hold on;
% plot(t,imag(e_s),'-r');
% plot(t,msg,'xg','MarkerSize',10);


%% Signal RF

s = real(e_s).*cos(2*pi*Fc*(t-1)*1/Fs) - imag(e_s).*sin(2*pi*Fc*(t-1)*1/Fs);

% figure(3);
% plot(t,s); hold on;
% plot(t,msg,'xg','MarkerSize',10);

%s2 = fmmod(msg,Fc,Fs,Kf);

%plot(t,s,'-or');
% hold on;
% plot(n(1:plt),s2(1:plt),'-xb');

%% Bruit blanc


%% Filtre passe-bas
%Optention des deux signal qui vont rentrer dans le filtre passe-bas

r_cos = s.*cos(2*pi*Fc*(t-1)/Fs);
r_sin = s.*sin(2*pi*Fc*(t-1)/Fs);


% subplot(211); hold on;
% subplot(212); hold on;
% rcosfir(0,20,10,1/(3*Fc),'r');

h = rcosfir(0,20,10,1/(3*Fc));

e_r = conv(r_cos,h) + 1i*conv(r_sin,h);

e_r = e_r(20*10+1:vector_size+20*10);


%% Démodulation FSK non cohérente

e_s0 = exp(-1i*2*pi*Kf*t);
e_s1 = exp(1i*2*pi*Kf*t);

for k=1:bits_size
    
    sum_s0 = 0;
    sum_s1 = 0;
    
    for i=1:OSF
        sums_s0 = sum_s0 + e_r((k-1)*OSF+i)*conj(e_s0((k-1)*OSF+i));
        sums_s1 = sum_s1 + e_r((k-1)*OSF+i)*conj(e_s1((k-1)*OSF+i));
    end
    
    if abs(sums_s0) < abs(sums_s1)
        m(1,k) = 1;
    else
        m(1,k) = 0;
    end
    
%     if  abs(sum(e_r((k-1)*OSF+1:k*OSF).*conj(e_s0((k-1)*OSF+1:k*OSF)))) > abs(sum(e_r((k-1)*OSF+1:k*OSF).*conj(e_s1((k-1)*OSF+1:k*OSF))))
%         m(k) = 0;
%     else
%         m(k) = 1;
%     end
end

%m = abs(cumsum(e_r.*conj(e_s1))) - abs(cumsum(e_r.*conj(e_s0)));
% for t2 = 1:size(m,2)
%     if m(1,t2) > 0
%         m(1,t2) = 1;
%     else
%         m(1,t2) = 0;
%     end
% end



%% Démodulation FM

