function [FOM, SNRc, SNR0] = modul_fm(F,Fs,Kf,Fc,Eb_N0)
%% Initialisation

clc;
%clear;
% 
% F = 500;            %Fr�quence symbole
% Fs = 20000;         %Fr�quence d'�chantillonnage Matlab
Ts = 1/Fs;
OSF = Fs/F;
bits_size = 10000;   %Nombre de symbole �mis
% Kf = 250;           %S�lectivit� fr�quentielle
% Fc = 6000;          %Fr�quence porteuse
% Eb_N0 = 10;         %en db

t=0:1/Fs:OSF*bits_size*1/Fs-1/Fs;

longueur_chaine = bits_size*OSF;
tic();

%% Cr�ation des bits

bits = randi(2,1, bits_size) - 1;   %-1 pour avoir 0 ou 1 en al�atoire

%% Cr�ation du message

fprintf('Cr�ation du message\n');

msg = zeros(1,longueur_chaine);   %D�finition de notre matrisse message = m(t)

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

fprintf('Modulation FM\n');

%Cumsum est la somme cumul�e, c'est l'int�grale discr�te de notre message
e_s = exp(1i*2*pi*Kf*cumsum(msg)/Fs);


%% Signal RF

fprintf('Signal RF\n');

s = real(e_s).*cos(2*pi*Fc*t) - imag(e_s).*sin(2*pi*Fc*t);

%s = fmmod(msg,Fc,Fs,Kf);

%% Bruit blanc

fprintf('Bruit Blanc\n');

%Param�tre du bruit
N0 = sum(s.^2)/(10^(Eb_N0/10)*F*OSF*longueur_chaine);
sigma_n = sqrt(N0*Fs/2);

%G�n�ration du bruit
n = sigma_n*randn(1,longueur_chaine);

r = s + n;

% figure(1);
% plot(t,r); hold on;
% plot(t,s,'r');

%% Filtre passe bas

fprintf('FPB\n');

%Filtre appliqu� juste sur le message sans bruit
r_cos = s.*cos(2*pi*Fc*t);
r_sin = s.*sin(2*pi*Fc*t);

h = rcosfir(0,10,20,1/(3*Fc));

e_r = conv(r_cos,h) + 1i*conv(r_sin,h);

e_r = e_r(1+10*20:longueur_chaine+20*10);

%Filtre appliqu� sur le message avec bruit
r_cos = r.*cos(2*pi*Fc*t);
r_sin = r.*sin(2*pi*Fc*t);

h = rcosfir(0,10,20,1/(3*Fc));

e_r_bruit = conv(r_cos,h) + 1i*conv(r_sin,h);

e_r_bruit = e_r_bruit(1+10*20:longueur_chaine+20*10);

%% D�modulateur FSK

fprintf('D�modulateur FSK\n');

e_s0 = exp(-1i*2*pi*Kf*t);
e_s1 = exp(1i*2*pi*Kf*t);

m_fsk = zeros(1,bits_size);
for k=1:bits_size
    
    s0_sum = 0;
    s1_sum = 0;
     
    for i=1:OSF
        s0_sum = s0_sum + e_r_bruit((k-1)*OSF+i)*conj(e_s0((k-1)*OSF+i));
        s1_sum = s1_sum + e_r_bruit((k-1)*OSF+i)*conj(e_s1((k-1)*OSF+i));
    end
    
    if abs(s0_sum) > abs(s1_sum)
        m_fsk(1,k) = 1;
    else
        m_fsk(1,k) = 0;
    end
end



%% D�modulateur FM

fprintf('D�modulateur FM\n');

Bt = 2*(Kf+Fc); %largeur de bande
a = 1/(2*pi*Kf); % valeur de la pente pour retomber sur le message d'origine

% d_er = diff(e_r)./diff(t);   %derivee par rapport au temps (1x399)
% %probleme de dimension e_r=1x400

%D�modulation sur le message sans bruit

d_er = gradient(e_r,1/Fs);  %derivee par rapport au temps (1x400)

s1 = a.*d_er+1i*pi*a*Bt.*e_r;    % apres 2 filtres a pente
s2 = -a.*d_er+1i*pi*a*Bt.*e_r;    % apres 2 filtres a pente 

m_r = zeros(1,longueur_chaine);

for k=1:longueur_chaine
    if abs (s2(k))-abs(s1(k))<0
        m_r(k)=0;
    else
        m_r(k)=1;
    end
end

%D�modulation sur le message avec bruit
d_er_bruit = gradient(e_r_bruit,1/Fs);  %derivee par rapport au temps (1x400)

s1 = a.*d_er_bruit+1i*pi*a*Bt.*e_r_bruit;    % apres 2 filtres a pente
s2 = -a.*d_er_bruit+1i*pi*a*Bt.*e_r_bruit;    % apres 2 filtres a pente 

m_r_bruit = zeros(1,longueur_chaine);

for k=1:longueur_chaine
    if abs (s2(k))-abs(s1(k))<0
        m_r_bruit(k)=0;
    else
        m_r_bruit(k)=1;
    end
end


%% Interpr�tation du message
fprintf('Interpr�tation du message apr�s d�modulation FM\n');

m_fm = zeros(1,bits_size);
for k=1:bits_size

    for l=1:OSF
        m_fm(k)=(m_fm(k)+m_r_bruit((k-1)*OSF+l));
    end
    m_fm(k)=round(m_fm(k)/OSF);
end

%% Comparaison des r�sultats

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
fprintf('Calcul d''erreur\n');

erreur_FSK = 0;
erreur_FM = 0;
for i=1:bits_size
    if bits(1,i) ~= m_fsk(1,i)
        erreur_FSK = erreur_FSK + 1;
    end
    
    if bits(1,i) ~= m_fm(1,i)
        erreur_FM = erreur_FM + 1;
    end
end

pourcentage_erreur_fsk = erreur_FSK/bits_size*100;
pourcentage_erreur_fm = erreur_FM/bits_size*100;

clc;
fprintf('Pourcentage d''erreur FSK: %f %% \n Il y a %i bits erron�s sur %i \nPourcentage d''erreur FM: %f %% \n Il y a %i bits erron�s sur %i \n', pourcentage_erreur_fsk, erreur_FSK, bits_size, pourcentage_erreur_fm, erreur_FM, bits_size);

%% Calcul des diff�rents SNRc SNRo, FOM

fprintf('Calcul de SNRc, SNRo et FOM\n');

%Puissance d'un bit avant d�modulation
Ps_c = sum(r.^2)/longueur_chaine;
Pn_c = (N0*F*OSF)/2;

%Puissance d'un bit apr�s d�modulation
Ps_0 = sum(m_r_bruit.^2)/longueur_chaine;
Pn_0 = sum((m_r_bruit-m_r).^2)/longueur_chaine;

SNRc = 10*log10(Ps_c/Pn_c);
SNR0 = 10*log10(Ps_0/Pn_0);
FOM = SNR0-SNRc;

fprintf('SNRc = %.2f \nSNR0 = %.2f \nFOM = %.2f \nEb/N0 = %.0f \n',SNRc,SNR0,FOM, Eb_N0);

%% Densit� de puissance par unit� de fr�quence

[psd_m,f_m] = welch(msg,t);
[psd_e_s, f_e_s] = welch(e_s,t);
[psd_s, f_s] = welch(s,t);
[psd_r, f_r] = welch(r,t);
[psd_e_r, f_e_r] = welch(e_r,t);
[psd_m_rb, f_m_rb] = welch(m_r_bruit,t);

figure(2);

subplot(3,2,1);
plot(f_m,psd_m);
title('Densit� spectrale de m');

subplot(3,2,2);
plot(f_e_s,psd_e_s);
title('Densit� spectrale de Es');

subplot(3,2,3);
plot(f_s,psd_s);
title('Densit� spectrale de s');

subplot(3,2,4);
plot(f_r,psd_r);
title('Densit� spectrale de r');

subplot(3,2,5);
plot(f_e_r,psd_e_r);
title('Densit� spectrale de Er');

subplot(3,2,6);
plot(f_m_rb,psd_m_rb);
title('Densit� spectrale de Mr');
clear;


toc();

end