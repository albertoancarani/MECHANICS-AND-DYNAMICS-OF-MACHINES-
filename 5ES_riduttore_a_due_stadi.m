clc 
clear

%% TESTO

% In Fig. 8 è rappresentato un riduttore coassiale a due stadi le cui 
% ruote dentate sono cilindriche a denti dritti.
% Volendo garantire la coassialità tra albero di ingresso 
% (su cui è calettata la ruota 1) e quello di uscita (su cui è montata 
% la ruota 4), i dati impongono di operare una correzione.
% Si sceglie di intervenire sul secondo stadio, per il quale si chiede di determinare:

% 1. l’angolo di pressione di lavoro α3,4 (espresso in gradi);
% 2. la somma dei coefficienti di spostamento x3 + x4;
% 3. la somma degli spostamenti di profilo v3 + v4 (espressa in mm).

u=0; % u è il resto della divisione per 4 del numero di matricola

%% DATI

M1_2 = 2.5; % [mm] Modulo del primo stadio 
M3_4 = 3; % [mm] Modulo del secondo stadio 
% Numero di denti delle ruote:
Z1 = 50;
Z2 = 160 +u ;
Z3 = 59 +u ;
Z4 = 118 - u;
alpha = 20*pi/180; % [rad] Angolo di pressione di taglio 

%% RAGGI PRIMITIVI DELLE RUOTE DENTATE

R1=M1_2*Z1/2; %[mm]
R2=M1_2*Z2/2; %[mm]
R3=M3_4*Z3/2; %[mm]
R4=M3_4*Z4/2; %[mm]

%% INTERASSI DI TAGLIO

% a= interassee di taglio

a1_2= R1+R2; %[mm] 1 stadio
a3_4= R3+R4; %[mm] 2 stadio

%% COASSIALITA' TRA I DUE ALBERI

% correggo il 2° stadio 

% condizione di coassialità : a3_4_primo = a1_2  
% (a3_4_primo va calcolato dopo la correzione)

a3_4_primo = a1_2; % a_primo=interasse di lavoro

%% ANGOLI DI PRESSIONE

% a3_4_primo = a3_4 * cos(alpha) / cos (alpha3_4)
alpha3_4 = acos( a3_4 * cos(alpha) /  a3_4_primo ) % [rad]

%% COEFFICIENTI DI SPOSTAMENTO

%cerco la somma dei coefficienti di spostamento X3+X4

% teoria: involute(alpha_primo)= involute(alpha0) + 2*tan(alpha0) * (x1+x2) / (z1+z2)

% 2*tan(alpha) * (X3+X4) + (Z3+Z4)*(involute(alpha) - involute(alpha3_4) ) = 0

Somma_X3_X4= (Z3+Z4)*( - involute(alpha) + involute(alpha3_4) ) / (2*tan(alpha) )

%% SPOSTAMENTI DI PROFILO 

%cerco la somma degli spostamenti di profilo v3+v4 [mm]
% v= m*x 

Somma_v3_v4 = Somma_X3_X4 * M3_4

