%% Initialisation
clc;
clear;
plt = 50;

F = 500;            %Fréquence symbole
Fs = 20000;         %Fréquence d'échantillonnage Matlab
OSF = Fs/F;
bits_size = 10;   %Nombre de symbole émis
Kf = 500;           %Sélectivité fréquentielle
Fc = 6000;          %Fréquence porteuse

t=0:1/Fs:OSF*bits_size*1/Fs-1/Fs;

longueur_chaine = bits_size*OSF;

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
plot(msg,'-','LineWidth',3);

%% Modulation FM

%Cumsum est la somme cumulée, c'est l'intégrale discrète de notre message
e_s = exp(1i*2*pi*Kf*cumsum(msg)/Fs);


%% Signal RF

s = real(e_s).*cos(2*pi*Fc*t) - imag(e_s).*sin(2*pi*Fc*t);

% s2 = fmmod(msg,Fc,Fs,Kf);
% 
% plot(t,s,'-or');
% hold on;
% plot(t,s2,'-xb');

%% Bruit blanc

n = zeros(1, longueur_chaine);

r = s + n;

%% Filtre passe bas

r_cos = r.*cos(2*pi*Fc*t);
r_sin = r.*sin(2*pi*Fc*t);

h = rcosfir(0,10,20,1/(3*Fc));

e_r = conv(r_cos,h) + 1i*conv(r_sin,h);

e_r = e_r(1+10*20:longueur_chaine+20*10);

%% Démodulateur FSK

e_s0 = exp(-1i*2*pi*Kf*t);
e_s1 = exp(1i*2*pi*Kf*t);

for k=1:bits_size
    
    s0_sum = 0;
    s1_sum = 0;
     
    for i=1:OSF
        s0_sum = s0_sum + e_r((k-1)*OSF+i)*conj(e_s0((k-1)*OSF+i));
        s1_sum = s1_sum + e_r((k-1)*OSF+i)*conj(e_s1((k-1)*OSF+i));
    end
    
    if abs(s0_sum) > abs(s1_sum)
        m(1,k) = 1;
    else
        m(1,k) = 0;
    end
end

%% Comparaison des résultats

msg = [];
for i=1:bits_size
   
    bit = m(1,i);
    if bit == 0         %Si le bit vaut 0 on écrit -1
        bit = -1;
    end
    
    bit_msg = bit;      %On duplique la valeur du bit OSF fois
    bit_msg = repmat(bit_msg,1,OSF);
    
    msg = horzcat(msg,bit_msg);
    
end

hold on;
plot(msg,'-r');