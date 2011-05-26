% 
F = 500;            %Fréquence symbole
Fs = 20000;         %Fréquence d'échantillonnage Matlab
Kf = 500;           %Sélectivité fréquentielle
Fc = 6000;          %Fréquence porteuse
Eb_N0 = 10;         %en db

result = zeros(3,nombre_tests);

 %nombre_tests = 30;
%Eb_N0 = 0:1:29;

%nombre_tests = 10;
%Kf = 100:100:1000;

% nombre_tests = 5;
% F = [250 500 1000 2000 2500];

%nombre_tests = 40;
% Fs = 20000; %Fs = 40000;
% Fc = 1000:1000:40000;

% nombre_tests = 31;
% Fs = 10000:1000:40000;

%Paramètres optimisés en fonction des maximums de chaque essais éffectués
%plus haut
nombre_tests = 1;
F = 250;            %Fréquence symbole
Fs = 40000;         %Fréquence d'échantillonnage Matlab
Kf = 400;           %Sélectivité fréquentielle
Fc = 9000;          %Fréquence porteuse
Eb_N0 = 1;         %en db



for i=1:nombre_tests

    [result(1,i), result(2,i), result(3,i)] = modul_fm(F,Fs,Kf,Fc(i),Eb_N0);
    
    
end

figure(1)
plot(Fc,result);
legend('FOM','SNRc','SNR0');

% figure(2)
% 
% plot(F,result);
% Title('SNRc');
% 
% figure(3)
% 
% plot(F,result);
% Title('SNR0');
