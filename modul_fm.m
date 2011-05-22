%% Initialisation
clc;
clear;
plt = 50;

F = 500;            %Fr�quence symbole
Fs = 20000;         %Fr�quence d'�chantillonnage Matlab
OSF = Fs/F;
bits_size = 10;   %Nombre de symbole �mis
Kf = 500;           %S�lectivit� fr�quentielle
Fc = 6000;          %Fr�quence porteuse

t=0:1/Fs:OSF*bits_size*1/Fs-1/Fs;

longueur_chaine = bits_size*OSF;

%% Cr�ation des bits
bits = randi(2,1, bits_size) - 1;   %-1 pour avoir 0 ou 1 en al�atoire

%% Cr�ation du message
msg = [];   %D�finition de notre matrisse message = m(t)

for i=1:bits_size
   
    bit = bits(1,i);
    if bit == 0         %Si le bit vaut 0 on �crit -1
        bit = -1;
    end
    
    bit_msg = bit;      %On duplique la valeur du bit OSF fois
    bit_msg = repmat(bit_msg,1,OSF);
    
    msg = horzcat(msg,bit_msg);
    
end

figure(1);
plot(msg,'-','LineWidth',3);

%% Modulation FM

%Cumsum est la somme cumul�e, c'est l'int�grale discr�te de notre message
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

%% D�modulateur FSK

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

%% Comparaison des r�sultats

msg = [];
for i=1:bits_size
   
    bit = m(1,i);
    if bit == 0         %Si le bit vaut 0 on �crit -1
        bit = -1;
    end
    
    bit_msg = bit;      %On duplique la valeur du bit OSF fois
    bit_msg = repmat(bit_msg,1,OSF);
    
    msg = horzcat(msg,bit_msg);
    
end

hold on;
plot(msg,'-r');