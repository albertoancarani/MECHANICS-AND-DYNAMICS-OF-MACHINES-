clc 
clear


% Una ruota cilindrica a denti dritti, i cui fianchi hanno profilo ad evolvente, viene corretta
% durante il taglio con dentiera. Con i dati a disposizione, determinare il raggio di curvatura
% del profilo dei denti in corrispondenza della circonferenza primitiva e della circonferenza
% di testa (si ipotizzi l’assenza di smussi o raggi di raccordo).

% soluzione:
% raggio di curvatura sulla circonferenza primitiva = 10.603 mm
% raggio di curvatura sulla circonferenza di testa = 14.411 mm

%% dati

m0 = 2;% [mm] Modulo della dentiera generatrice 
alpha0 = 20; % [deg] Angolo di pressione della dentiera generatrice
z = 31; % Numero di denti della ruota
x = -0.25; %Coefficiente di spostamento 

%% calcoli

R= m0*z/2; % [mm] 
v=x*m0; % vano [mm]
rho= R*sind(alpha0); % raggio di curvatura sulla circonferenza primitiva



e1=m0-v; %addendum pignone
e2=m0+v; %addendum ruota
Re1=R+e1; % raggio di testa del pignone
Re2=R+e2; % raggio di testa della ruota




