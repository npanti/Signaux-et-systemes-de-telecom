%% Initialisation
clc;
clear;

F = 500;            %Fr�quence symbole
Fs = 20000;         %Fr�quence d'�chantillonnage Matlab
OSF = Fs/F;
bits_size = 1000;   %Nombre de symbole �mis
Kf = 500;           %S�lectivit� fr�quentielle
Fc = 6000;          %Fr�quence porteuse

t = 1/Fs:1/Fs:OSF*1/Fs*bits_size;

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

%% Modulation FM
%Es = fmmod(msg,Fc,Fs,50); ou modulate
% somme = spalloc(1,40000,40000);
% for n=1:40000
%     somme(1,n) = sum(msg(1:n))/Fs;
% end

%Cumsum est la somme cumul�e, c'est l'int�grale discr�te de notre message
Es = exp(1i*2*pi*Kf*cumsum(msg)/Fs);


%% Signal RF
%t=1:40000;

s = real(Es).*cos(2*pi*Fc*(t-1)*1/Fs) - imag(Es).*sin(2*pi*Fc*(t-1)*1/Fs); 

%s2 = fmmod(msg,Fc,Fs,Kf);

%plot(t,s,'-or');
% hold on;
% plot(n(1:plt),s2(1:plt),'-xb');

%Optention des deux signal qui vont rentrer dans le filtre passe-bas

r_cos = s.*cos(2*pi*Fc*(t-1)/Fs);
r_sin = s.*sin(2*pi*Fc*(t-1)/Fs);

%% Filtre passe-bas

b = rcosfir(0,20,10,1/(3*Fc));
conv_cos = conv(r_cos,b);
conv_sin = conv(r_sin,b);

e_r = conv_cos + 1i*conv_sin;

e_r(:,1:200) = [];
e_r(:,size(e_r,2)-199:size(e_r,2)) = []; 


%% D�modulation FSK non coh�rente

Ac = 1;

e_s0 = exp(-1i*2*pi*Kf*t);
e_s1 = exp(1i*2*pi*Kf*t);

m = abs(cumsum(e_r.*conj(e_s1))) - abs(cumsum(e_r.*conj(e_s0)));

for t2 = 1:size(m,2)
    if m(1,t2) > 0
        m(1,t2) = 1;
    else
        m(1,t2) = 0;
    end
end

%% D�modulation FM

