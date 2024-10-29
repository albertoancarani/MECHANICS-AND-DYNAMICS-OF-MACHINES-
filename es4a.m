clc
clear

% Una ruota cilindrica a denti dritti con Z = 23 denti, viene tagliata con una dentiera avente
% modulo m = 1.5 mm. Determinare il diametro della circonferenza di testa nel caso in cui
% il taglio avvenga con una correzione negativa con coefficiente di spostamento x = −0.25.
% Soluzione: Diametro della circonferenza di testa = 36.75 mm

%% dati
z=23; % denti
mo=1.5; % [mm] modulo dentiera
x= -0.25; % coefficiente di spostamento

%% calcoli

R= mo*z/2; % [mm] 
v=x*mo;;

e1=mo-v; %addendum pignone
e2=mo+v; %addendum ruota
Re1=R+e1; 
Re2=R+e2; % raggio di testa della ruota
De2=Re2*2