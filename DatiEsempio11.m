% Modello di rettifica per sbarbatori

% Parametri del motore elettrico
% R 		[ohm]
% L 		[Henry=s*Volt/A]
% Kc 		[Nm/A]
% Kb 		[Volt/(rad/s)]

% Parametri del controllo di corrente
% Kpc [V/A]
% Tic [s]

% Parametri del controllo di velocita'
% Kpv [Nm/(rad/s)]
% Tiv [s]

% Parametri del controllo di posizione
% Kpp [1/s]
% Tip [s]

% Motore elettrico
Kc=5;
Kb=Kc;
R=0.4;
L=0.003;

% Controllo corrente PI
Kpc=8;
Tic=0.002;

% Controllo di velocita' PI
Tiv = 0.1;
Kpv=95;

% Controllo di posizione PI o P a seconda del valore di Tip
Kpp=72;
Tip=1000;		% equivale ad un controllo P
%Tip=0.1;

% Legge di moto
load LeggeEsempio11

% Velocita' di rotazione [giri/min]
rpm=20;

% Periodo di rotazione [s]
per=60/rpm;

% Parametri del sistema meccanico (J in kg*m^2, k in Nm/rad):
Jm=0.6;
J2=0.085;
J3=0.085;

k1=1.15e6;
k2=1.15e5;

qq=1e-05;
c1=qq*k1;
c2=qq*k2;

% Gioco
%g=1.e-6;
g=0;


