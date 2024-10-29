clc,clear

%% TESTO

% Determinare la prima frequenza propria di vibrazione flessionale della 
% colonna con serbatoio elevato mostrata in Fig. 19,
% supponendo che la sezione tubolare della colonna sia costante.
% Si esprima il risultato in Hz impiegando almeno cinque cifre significative.
% uv ultime due cifre del numero di matricola
% (matricola = 0000####uv).

u=8;
v=4;
%% DATI

D = 3 + u/10; % [m]  diametro esterno della colonna
rho = 2400 + v^2 + u; %[kg/m3] massa volumica del materiale della colonna
d = 2.45 + v/30; % [m]  diametro interno della colonna
E = 2.8 * 10^10 ; % [N/m2]  modulo di elasticità del materiale della colonna
l = 90 + (u^2) /5 - v; % [m]  lunghezza della colonna
Q = (2.7 + (u^2) /100 + u * v/50) * 10^6; %[N]  peso del serbatoio
g=9.81; %[m/s^2]
%% SEZIONE DELLA TRAVE

%la prima frequenza di vibrazione si studia ipotizzando la prima forma
%modale alla linea elastica. Il mio sistema è approssimabile a una trave
%incastrata con carico all'estremità

x=[0:1:l]; % posizione nella trave
I=pi*(D^4-d^4)/64; %momento di inerzia trave incastrata
% P= carico all'estremità
% y(x)=P*x^2 / (6*E*I) * (3*l-x); %equazione linea elastica
%ymax =y(x=l) =P*(l^3) / 3*E*I;


% k= rigidezza flessionale .
% K si trova imponendo P=1--> k= 1/ymax
k = 1 / ( (l^3) / 3*E*I );

% S= sezione della trave
S=pi * (D^2-d^2)/4;
%% Metodo energetico di Rayleigh Tmax= Vmax

% funzione vibratoria
% v(t,x)= f(t)*y(x);
% f(t)= funzione del tempo
% y(x) = funzione della cordinaata spaziale

% A= ampiezza
% w1 = approssimazione 1 pulsazione
%f(t)= A*cos(w1*t)
%df(t)/dt= -w1*A*sin(w1*t) 

% T= energia cinetica. T= Tcolonna + Tserbatoio

% vedi passaggi integrali Tc e Ts su word
% Tc=1/8*rho*S*( (df(t)/dt)^2 ) * ymax^2 *33/35 *l
% Ts=1/2 * Q/g * ( ( (df(t)/dt)^2 ) * ymax^2 )

% T = ( ( (df(t)/dt)^2 ) * yma x^2 ) * (1/8 * rho*S*l *33/35 + 1/2 Q/g )

% V= 1/2 * k * (v(t,x=l)^2 =1/2 * 3*E*I / (l^3 ) *f(t)^2 *ymax^2

% Tmax = w1^2 * A^2* max( (sin(w1*t))^2 ) * ymax^2 *(1/8 *rho*S*l*33/35 + 1/2 Q/g )
% max( (sin(w1*t))^2 ) =1

% Vmax= 1/2 * k  *A^2 *max( (cos(w1*t))^2 ) *ymax^2
% max( (cos(w1*t))^2 )

% Tmax= Vmax
% w1^2 * A^2 * ymax^2 *(1/8 *rho*S*l*33/35 + 1/2 Q/g ) = 1/2 * k  *A^2 *ymax^2
% w1^2  = 1/2 * k / (33/280 *rho*S*l + 1/2 Q/g )
w1= sqrt( ( (3/2 *E*I ) / (l^3 ) )  / (33/280 *rho*S*l + Q /(2*g) ) ); % rad/s

% f1 = prima frequenza del sistema
f1= w1/( 2*pi) % 1/s