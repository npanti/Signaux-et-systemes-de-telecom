% 
F = 500;            %Fr�quence symbole
Fs = 20000;         %Fr�quence d'�chantillonnage Matlab
Kf = 500;           %S�lectivit� fr�quentielle
Fc = 6000;          %Fr�quence porteuse
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

%Param�tres optimis�s en fonction des maximums de chaque essais �ffectu�s
%plus haut
nombre_tests = 1;
F = 250;            %Fr�quence symbole
Fs = 40000;         %Fr�quence d'�chantillonnage Matlab
Kf = 400;           %S�lectivit� fr�quentielle
Fc = 9000;          %Fr�quence porteuse
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
