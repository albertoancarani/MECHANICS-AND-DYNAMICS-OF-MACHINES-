clc
clear

% Determinare l’interasse con il quale montare un ingranaggio costituito
% da due ruote cilindriche a denti diritti non corrette, 
% in modo che sulle primitive di lavoro le dentature
% presentino un assegnato gioco laterale necessario 
% per una adeguata lubrificazione.

% Soluzione: Interasse di montaggio = 150.407 mm

%% Dati

m0 = 2; % [mm] Modulo della dentiera generatrice 
alpha0 = 20; % [deg] Angolo di pressione della dentiera generatrice 
z1 = 40; % Numero di denti della prima ruota 
z2 = 110; % Numero di denti della seconda ruota 
delta = 0.3; % [mm] Gioco laterale sulla primitiva di lavoro 

alpha0=alpha0*pi/180; %[rad]

%% traccia 

% % R= raggio della primitiva di taglio
% S = pi*m0/2; % spessore sulla primitiva di taglio
% 
% %primitiva di lavoro = primo
% fdelta= ( cos(alpha0)/ cos(alphaprimo) );
% Rprimo=R*fdelta;
% aprimo=a*fdelta;
% mprimo=m0*fdelta;
% pprimo=p*fdelta;
% 
% % Sprimo = Rprimo * ( S/R - 2( involute(alphaprimo)- involute(alpha0) )); 
% Sprimo = S*fdelta- 2*Rprimo* ( involute(alphaprimo)- involute(alpha0) );
% % spessore sulla primitiva di lavoro
% 
% %somma spessori denti sulla primitiva di lavoro
% SommaS1primoS2primo=pi*mprimo*-2*aprimo* ( involute(alphaprimo)- involute(alpha0) )


%% cerco lo zero della funzione  f(α′, δ, α0, a)

R1= m0*z1/2;
R2= m0*z2/2;
a= R1+R2;

% delta - 2*a*( cos(alpha0)/ cos(alphaprimo) )* (involute(alphaprimo) - involute(alpha0) ) = 0
alphaprimovettore= [ 1:1:360 ].*(pi/180);
alpha0vettore= ones(1,length(alphaprimovettore) )*alpha0;
deltavettore= ones(1,length(alphaprimovettore) )*delta;


fun4b= delta - ( cos(alpha0vettore)/ cos(alphaprimovettore) )* (involute(alphaprimovettore) - involute(alpha0vettore) ).*(2*a);

plot(alphaprimovettore.*(180/pi),fun4b);
grid on
%x=fzero(fun4b,alphaprimo);
%alphaprimo=x;

%fdelta= ( cos(alpha0)/ cos(alphaprimo) );
%aprimo=a*fdelta;

