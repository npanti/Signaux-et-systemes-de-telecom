% 
F = 500;            %Fréquence symbole
Fs = 20000;         %Fréquence d'échantillonnage Matlab
Kf = 400;           %Sélectivité fréquentielle
Fc = 6000;          %Fréquence porteuse
Eb_N0 = 10;         %en db
nombre_tests =1;

% nombre_tests = 59;
% Eb_N0 = 0:0.5:29;

% nombre_tests = 19;
% Kf = 100:50:1000;

% nombre_tests = 5;
% F = [250 500 1000 2000 2500];

% nombre_tests = 40;
% %Fs = 20000; 
% Fs = 40000;
% Fc = 1000:1000:40000;

% nombre_tests = 31;
% Fs = 10000:1000:40000;

%Paramètres optimisés en fonction des maximums de chaque essais éffectués
%plus haut
% nombre_tests = 1;
% F = 250;            %Fréquence symbole
% Fs = 40000;         %Fréquence d'échantillonnage Matlab
% Kf = 400;           %Sélectivité fréquentielle
% Fc = 9000;          %Fréquence porteuse
% Eb_N0 = 10;         %en db

result = zeros(5,nombre_tests);

for i=1:nombre_tests

    [result(1,i), result(2,i), result(3,i), result(4,i), result(5,i)] = modul_fm(F,Fs(i),Kf,Fc,Eb_N0);
    
    
end

% figure(1);
% plot(Eb_N0,result(1:3,:));
% legend('FOM','SNRc','SNR0');

figure(2);
plot(Fs, result(4:5,:));
legend('% erreur démodulateur FSK', '% erreur démodulateur FM');

% figure(2)
% 
% plot(F,result);
% Title('SNRc');
% 
% figure(3)
% 
% plot(F,result);
% Title('SNR0');
