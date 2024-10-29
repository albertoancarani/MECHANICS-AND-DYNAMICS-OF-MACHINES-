clc
clear

%% TESTO

% Un punto P di biella di un quadrilatero articolato deve assumere tre 
% posizioni assegnate da occuparsi in tre istanti prefissati.
% In particolare, a partire dall’istante t1 = 0 in cui il punto di biella è 
% nella posizione P1, i punti P2 e P3 devono essere raggiunti 
% rispettivamente negli istanti t2 e t3. Il movente del quadrilatero 
% (la prima asta A0A della diade di sinistra) ruota con velocità angolare 
% pari a n.
% Mediante sintesi cinematica, individuare il quadrilatero.
% In particolare, con riferimento alla Fig. 2, occorre determinare:
% 
% 1. le coordinate degli assi delle coppie rotoidali fisse (A0 e B0)
% 2. le lunghezze delle aste accoppiate a telaio (A0A e B0B)
% 3. la distanza tra gli assi delle coppie rotoidali di biella (AB)
% 4. l’angolo che l’asta A0A forma con la congiungente gli assi delle
% coppie fisse quando il punto P è nella prima posizione

% Si chiede inoltre di disegnare il quadrilatero ottenuto in opportuna scala.

% Affinché con i dati forniti il problema abbia una soluzione unica 
% (in modo da consentire il controllo dei risultati da parte del docente), 
% si assuma che quando il punto P occupa rispettivamente le posizioni P2 e P3:

% 1. la prima asta (B0B) della diade di destra abbia compiuto,
% a partire dalla posizione iniziale, rotazioni assegnate ψ2 e ψ3;
% 2. un segmento di biella abbia compiuto, a partire dalla posizione iniziale,
% rotazioni assegnate γ2 e γ3.
% Oltre alla soluzione così ottenuta, si chiede di riportare anche altre 
% soluzioni (almeno due) scegliendo arbitrariamente gli angoli ψj e γj . 
% I quadrilateri corrispondenti a tali ulteriori soluzioni vanno disegnati
% in scala opportuna e comparati con quello ottenuto dalla prima soluzione.


u=0; % è il resto della divisione per 4 del numero di matricola

%% DATI

P1 = 200-75i % [mm]
P2 = (25 * u + 60)+(25*u - 151)*1i; % [mm]
P3 = (25 * u + 100)+ (25 * u - 305)*1i; % [mm]
t1=0; %[s]
t2 = 1.05; %[s]
t3 = 2.1; %[s]
n = 20; %[rpm]
psi2 = 20*pi/180; %[rad]
psi3 = 50*pi/180; % [rad]
gamma2 = 15*pi/180; % [rad]
gamma3 = -30*pi/180; % [rad]

%% ANGOLI DI INCLINAZIONE DELL'ASTA AoA 

%L'asta AoA si porta nelle posizione A2 e A3

delta2= P2- P1; %[mm]
delta3= P3- P1; %[mm]

w_asta_AoA= n*2*pi/60; %[rad/s]
phi1= w_asta_AoA*t1; %[rad]
phi2= w_asta_AoA*t2; %[rad] A-->A2
phi3= w_asta_AoA*t3; %[rad] A-->A3

%% DIADE DI SINISTRA

% Ho 2 equazioni vettoriali in 6 incognite. Devo perciò scegliere
% arbitrariamente 2 valori e avere così 4 equazioni in 4 incognite
% scelgo phi2 e phi3 --> ho come incognite W e Z

% 1) W*(e^i*phi2 - 1) + Z*( e^i*gamma2 - 1) = delta2
% 2) W*(e^i*phi3 - 1) + Z*( e^i*gamma3 - 1) = delta3 

% ricorda e^i*x= cosx+isinx)
% CRAMER :

numW = [ delta2 , (cos(gamma2) +1i*sin(gamma2) - 1) ; delta3 , (cos(gamma3) +1i*sin(gamma3) - 1) ];
denW = [(cos(phi2) +1i*sin(phi2) - 1) , (cos(gamma2) +1i*sin(gamma2) - 1) ; (cos(phi3) +1i*sin(phi3) - 1), (cos(gamma3) +1i*sin(gamma3) - 1) ];

W= det(numW) / det(denW); %[mm]

numZ= [(cos(phi2) +1i*sin(phi2) - 1) , delta2 ; (cos(phi3) +1i*sin(phi3) - 1), delta3 ];
denZ= [(cos(phi2) +1i*sin(phi2) - 1) , (cos(gamma2) +1i*sin(gamma2) - 1); (cos(phi3) +1i*sin(phi3) - 1), (cos(gamma3) +1i*sin(gamma3) - 1) ];

Z= det(numZ) / det(denZ); %[mm]

%% CATENA DI DESTRA

% 1) W_asterisco*(e^i*psi2 - 1) + Z_asterisco*( e^i*gamma2 - 1) = delta2
% 2) W_asterisco*(e^i*psi3 - 1) + Z_asterisco*( e^i*gamma3 - 1) = delta3 

% CRAMER :

numW_asterisco = [ delta2 , (cos(gamma2) +1i*sin(gamma2) - 1) ; delta3 , (cos(gamma3) +1i*sin(gamma3) - 1) ];
denW_asterisco = [(cos(psi2) +1i*sin(psi2) - 1) , (cos(gamma2) +1i*sin(gamma2) - 1) ; (cos(psi3) +1i*sin(psi3) - 1), (cos(gamma3) +1i*sin(gamma3) - 1) ];

W_asterisco= det(numW_asterisco) / det(denW_asterisco); %[mm]

numZ_asterisco= [(cos(psi2) +1i*sin(psi2) - 1) , delta2 ; (cos(psi3) +1i*sin(psi3) - 1), delta3 ];
denZ_asterisco= [(cos(psi2) +1i*sin(psi2) - 1) , (cos(gamma2) +1i*sin(gamma2) - 1); (cos(psi3) +1i*sin(psi3) - 1), (cos(gamma3) +1i*sin(gamma3) - 1) ];

Z_asterisco= det(numZ_asterisco) / det(denZ_asterisco); %[mm]

%% COPPIE ROTOIDALI FISSE A0 B0

A0 = P1-Z-W
B0 = P1-Z_asterisco-W_asterisco

%% LUNGHEZZA ASTE ACCOPPIATE A TELAIO

A0A=sqrt( real(W)^2 + imag(W)^2 );
B0B=sqrt( real(W_asterisco)^2 + imag(W_asterisco)^2 );

%% DISTANZA TRA LE COPPIE ROTOIDALI DI BIELLA 
OA = P1-Z ;
OB = P1-Z_asterisco;
AB = sqrt((real(OA)-real(OB))^2+(imag(OA)-imag(OB))^2);
A=OA
B=OB
%% angolo tra asta A0A E A0B0
alfa = atand((imag(A0)-imag(B0))/(real(A0)-real(B0)));
beta = atand((imag(OA)-imag(A0))/(real(OA)-real(A0)));
gamma = abs(alfa)+abs(beta);

















