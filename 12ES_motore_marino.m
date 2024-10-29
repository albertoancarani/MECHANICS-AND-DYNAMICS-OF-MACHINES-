clc
clear

%% TESTO

% motore marino connesso all’elica mediante
% un riduttore ad ingranaggi ad uno stadio. 
% Noti i momenti di inerzia del volano, del motore, delle ruote dentate, 
% dell’elica e le dimensioni degli alberi, trovare le frequenze naturali
% e i modi di vibrare torsionali del sistema.

% Si trascuri l’inerzia degli alberi e si esprima il risultato 
% utilizzando almeno cinque cifre significative. Inoltre:
% 1. si esprimano le frequenze naturali in Hz;
% 2. indicata con θ la rotazione dell’asse motore e con φ la rotazione
% dell’asse dell’elica,esprimere i modi di vibrare nella seguente forma:
% r1 = { PHI2/ THETA1 }|1     r2 = { PHI2 / THETA1 }|2

%% DATI

u = 4; % ultima cifra del numero di matricola

Jv= 35000; % [kg · m^2] momento di inerzia del volano 
Jm = 1000 - 5 * u^2; %momento di inerzia del motore [kg · m^2]
Jz_1 = 250 - u;% [kg · m^2] momento di inerzia della ruota 1 
Jz_2 = 150 + 2 * u;% [kg · m^2] momento di inerzia della ruota  2
Je = 2000 + 20 * u^2;%[kg · m2] momento di inerzia dell’elica 
G = 8 * 10^10 ;%[N/m^2] modulo di elasticità tangenziale acciaio 
lm = 0.8; %[m] lunghezza albero motore
le = 1; %[m] lunghezza albero elica
dm =0.1; %[m] diametro albero motore
de =0.15; %[m] diametro albero elica
z1= 40; % denti_ruota1
z2= 20; % denti_ruota2

%% MODELLO A PARAMETRI CONCENTRATI

% 1° asse uso cordinata theta

% thetav= gdl rotazione volano;
% thetam= gdl rotazione motore;
% k1= molla torsionale con rigidezza uguale a quella dell'albero motore

% 2° asse uso cordinata phi

% k2= molla torsionale con rigidezza uguale a quella dell'albero elica
% phi_2= posizione angolare dell'elica
% phiZ_2=posizione angolare dell'inerzia Jz_2

% GDL

% thetav= gdl rotazione volano;
% thetam= gdl rotazione motore;
% phi_2= posizione angolare dell'elica

% phiZ_2 non è un gdl perchè dipende da thetam
% tau = z1/z2 = omega2/omega1 
% omega2 = phi_punto vel angolare a valle
% omega1 = theta_punto vel angolare a monte

Coppia_torcente=1; % coppia torcente unitaria per trovare K torsionale
Ip_1= pi*dm^4/32; %momento di inerzia polare di sezione dell'albero motore
Ip_2= pi*de^4/32; %momento di inerzia polare di sezione dell'albero elica

deltaTHETA_1= ( Coppia_torcente*lm ) / (G*Ip_1); % deformazione angolare causata da M
K1= 1/deltaTHETA_1; %K torsionale relativa all'albero motore

deltaTHETA_2= ( Coppia_torcente*le ) / (G*Ip_2); % deformazione angolare causata da M
K2= 1/deltaTHETA_2; %K torsionale relativa all'albero elica

tau= z1/z2;
%% IPOTESI 

% 1)
% ipotizzo l'albero tra motore e ruota 1 infinitamente rigido
% torsionalmente dato che è di lunghezza << lm e le .
% grazie a questa ipotesi posso affermare che Jz_1 e Jm siano rigidamente collegati 

% 2)
% Stesso ragionamento per Jz_2 e Je anche se ho un albero diverso su un
% asse diverso --> Jz_2 e Je sono rigidamente collegati 

% Ho 3 gdl --> mi aspetto quindi 3 pulsazioni naturali
% Osservando che il sistema non è collegato a telaio mi aspetto 1 Moto
% rigido e 2 moti vibratori ( w1=0 ; w2=\= 0 ; w3 =\=0 )

% 3)
% Volano perfetto --> Inerzia infinita --> theta 1° sezione = costante
% Studio un modello equivalente con un incastro in corrispondenza di Jv
% perchè Jv  molto maggiore rispetto alle altre inerzie ( Jv >> Je > Jm )
% avrò quindi solo 2 pulsazioni naturali ma molto simili a quelle del
% modello senza volano perfetto

%% MODELLO A 2 GDL

% GDL
% thetam= gdl rotazione motore
% phi_2= posizione angolare dell'elica

%% RIDUZIONE AL 1° ASSE

% _primo = riportato all'asse motore

% nuove variabili
% Jz_2_primo
% Je_primo
% K2_primo

%faccio delle equivalenze dinamiche per cercare le variabili dinamiche

% _punto = derivata prima
% _2punti = derivata seconda

% Energia cinetica:

% T (Jz_2) : 1/2 *Jz_2 * phi_punto^2 = 1/2* Jz_2_primo *theta_punto^2
% --> Jz_2_primo= Jz_2 * ( phi_punto^2 / theta_punto^2 ) = Jz_2 * tau^2

Jz_2_primo= Jz_2 * tau^2; % [kg * m^2]

% T (Je) : 
Je_primo = Je * tau^2; % [kg * m^2]

% energia potenziale : 
% 1/2 K2* (delta_phi)^2 = 1/2 K2_primo* (delta_theta)^2 

K2_primo= K2*tau^2; % [N*m]

%% EQUAZIONI DI EQUILIBRIO SISTEMA RIDOTTO

%theta2=theta elica

% 1) ( Jm + Jz_1 + Jz_2_primo) * thetam_2punti + K1*thetam +K2_primo*thetam - K2_primo*theta2 =0 

% 2) Je_primo * theta2_2punti + K2_primo*theta2 - K2_primo*thetam = 0

M = [ ( Jm + Jz_1 + Jz_2_primo) , 0 ; 0 Je_primo ];
K = [ (K1 + K2_primo) , -K2_primo ; -K2_primo , K2_primo];

%% RISULTATI

[V,D]=eig(K,M); 
% D = matrice dinamica diagonale
% V= matrice avente come colonne gli autovettori
theta11=V(1,1);
theta12=V(1,2);
theta21=V(2,1);
theta22=V(2,2);

w1= sqrt(D(1,1) ); % 1° pulsazione naturale
w2= sqrt(D(2,2) ); % 2° pulsazione naturale

% frequenze naturali
f1= w1/(2*pi) %[s^-1]
f2= w2/(2*pi) %[s^-1]

%modi di vibrare

% r1 = { PHI2 / THETA1 }|1 = (tau*phi2 /theta1 )|1 
r1 = tau * theta21/theta11

% r2 = { PHI2 / THETA1 }|2 = (tau*phi2 /theta1 )|2
r2 = tau * theta22/theta12