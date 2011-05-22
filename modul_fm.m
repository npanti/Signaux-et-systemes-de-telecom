%% Initialisation
clc;
clear;

F = 500;            %Fréquence symbole
Fs = 20000;         %Fréquence d'échantillonnage Matlab
Ts = 1/Fs;
OSF = Fs/F;
bits_size = 1000;   %Nombre de symbole émis
Kf = 500;           %Sélectivité fréquentielle
Fc = 6000;          %Fréquence porteuse

t=0:1/Fs:OSF*bits_size*1/Fs-1/Fs;

longueur_chaine = bits_size*OSF;
tic();
%% Création des bits
bits = randi(2,1, bits_size) - 1;   %-1 pour avoir 0 ou 1 en aléatoire

%% Création du message
msg = zeros(1,longueur_chaine);   %Définition de notre matrisse message = m(t)

k=0;
for i=1:bits_size
    
    if bits(1,i) == 0
        for j=1:OSF
            msg(1,k+j) = -1;
        end
    else
        for j=1:OSF
            msg(1,k+j) = 1;
        end
    end
    k = k+OSF;
end

%% Modulation FM

%Cumsum est la somme cumulée, c'est l'intégrale discrète de notre message
e_s = exp(1i*2*pi*Kf*cumsum(msg)/Fs);


%% Signal RF

s = real(e_s).*cos(2*pi*Fc*t) - imag(e_s).*sin(2*pi*Fc*t);

% s2 = fmmod(msg,Fc,Fs,Kf);

%% Bruit blanc

%Paramètre du bruit
N0 = sum(s.^2)/(10*F*OSF*longueur_chaine);
sigma_n = sqrt(N0*Fs/2);

%Génération du bruit
n = sigma_n*randn(1,longueur_chaine);

%n = zeros(1, longueur_chaine*OSF);

r = s + n;

% figure(1);
% hold on;
% %plot(t,r,'-k');
% plot(t,s);
% plot(t,n,'-r');

SNRc = (sum(r.^2)/longueur_chaine)/(N0*F*OSF);

fprintf('SNRc = %f \n',SNRc);
%% Filtre passe bas

r_cos = r.*cos(2*pi*Fc*t);
r_sin = r.*sin(2*pi*Fc*t);

h = rcosfir(0,10,20,1/(3*Fc));

e_r = conv(r_cos,h) + 1i*conv(r_sin,h);

e_r = e_r(1+10*20:longueur_chaine+20*10);

%% Démodulateur FSK

e_s0 = exp(-1i*2*pi*Kf*t);
e_s1 = exp(1i*2*pi*Kf*t);

m = zeros(1,bits_size);
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

% figure(2);
% plot(msg,'-','LineWidth',3);
% 
% msg2 = zeros(1,longueur_chaine);
% k=0;
% for i=1:bits_size
%     
%     if m(1,i) == 0
%         for j=1:OSF
%             msg2(1,k+j) = -1;
%         end
%     else
%         for j=1:OSF
%             msg2(1,k+j) = 1;
%         end
%     end
%     k = k+OSF;
% end
% 
% hold on;
% plot(msg2,'-r');

%% Calcul de pourcentage d'erreur
erreur = 0;
for i=1:bits_size
    if bits(1,i) ~= m(1,i)
        erreur = erreur + 1;
    end
end

pourcentage_erreur = erreur/bits_size*100;

fprintf('Pourcentage d''erreur: %f %% \n Il y a %i bits erronés sur %i \n', pourcentage_erreur, erreur, bits_size);
toc();