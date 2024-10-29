clc,clear

%% TESTO

% In Fig. 25 è rappresentato un sistema a 3 gdl. 
% Noti i valori delle masse e delle rigidezze,calcolare:
% 1. le 3 pulsazioni naturali del sistema (esprimerle in rad/s);
% 2. le 3 forme modali (eseguire la normalizzazione in modo che la prima
% componente sia unitaria).
% 3. introdotte nel sistema le modifiche strutturali indicate nei dati,
% calcolare il nuovo valore della seconda pulsazione propria del sistema 
% impiegando il quoziente di Rayleigh

%% DATI

u =1; %ultima cifra del numero di matricola

m = 1 + u/10; % [kg]
k = 1 - u/10; % [N/m]
m1 = 2 * m;
m2 = 3 * m;
m3 = 2 * m;
k1 = 4 * k;
k2 = 3 * k;
k3 = 5 * k;

% Modifiche strutturali:

% variazioni dovute alle modifiche strutturali del modello
deltam3 = 0.4 * m;
deltak2 = 0.7 * k;

%% EQUAZIONI DI EQUILIBRIO

% 1) m1*x1_2punti + k1*x1 +k2*x1 - k2*x2 = 0
% 2) m2*x2_2punti + k2*x2 + k3*x2 -k2*x1 -k3*x3 = 0
% 3) m3*x3_2punti + k3*x3 - k3*x2 = 0

% forma matriciale [M]{X_2punti} + [K]{X}={0}
[M]=[m1 0 0; 0 m2 0; 0 0 m3];
[K]=[(k1+k2) -k2 0; -k2 (k2+k3) -k3; 0 -k3 k3];

%% PULSAZIONI NATURALI 

[Mmod,D]=eig(K,M);
% D matrice diagonale con pulsazioni al quadrato 
% Mmod matrice con forme modali
% Mmod= [{X}1 {X}2 {X}3 ]

w1=sqrt(D(1,1)) %rad/s
w2=sqrt(D(2,2)) %rad/s
w3=sqrt(D(3,3)) %rad/s

%% NORMALIZZAZIONE FORME MODALI

p = Mmod(1,:);
P = ones(1,length(Mmod)).*p;
Mmod_normalizzata = Mmod./P;

prima_forma_modale_normalizzata_trasposta=(Mmod_normalizzata(:,1)).';
seconda_forma_modale_normalizzata_trasposta=(Mmod_normalizzata(:,2)).';
terza_forma_modale_normalizzata_trasposta=(Mmod_normalizzata(:,3)).';

%% INTRODUZIONE MODIFICHE STRUTTURALI

% modificando masse e rigidezze variano anche i modi e le frequenze
% {X}j* = {X}j + delta{X}j
% wj* = wj +deltawj

[deltaK] = [ deltak2 -deltak2 0; -deltak2 deltak2 0; 0 0 0];
[deltaM] = [ 0 0 0; 0 0 0; 0 0 deltam3 ];

%% QUOZIENTE DI RAYLEIGH
% cerco la seconda frequenza propria dopo le modifiche 

% wj*^2= ( {X}j_trasposta*[K+deltaK]*{X}j ) / ( {X}j_trasposta*[M+deltaM]*{X}j )
% --> wj*^2= ( {X}j_trasposta*[K]*{X}j + {X}j_trasposta*[deltaK]*{X}j ) /
%          / ( {X}j_trasposta*[M]*{X}j + {X}j_trasposta*[deltaM]*{X}j )

% ricorda : wj^2 = ( {X}j_trasposta*[K]*{X}j ) / ({X}j_trasposta*[M]*{X}j)
% raccolgo wj^2 per semplificare 
% wj^2 = wj^2 * (Num/Den)
% Num = 1 + ( {X}j_trasposta*[deltaK]*{X}j ) / ( {X}j_trasposta*[K]*{X}j )
% Den = 1 + ( {X}j_trasposta*[deltaM]*{X}j ) / ( {X}j_trasposta*[M]*{X}j )

seconda_forma_modale_trasposta=(Mmod(:,2)).';
seconda_forma_modale=Mmod(:,2);

Num = 1 + ( seconda_forma_modale_trasposta * deltaK * seconda_forma_modale ) / ( seconda_forma_modale_trasposta * K * seconda_forma_modale );
Den = 1 + ( seconda_forma_modale_trasposta * deltaM * seconda_forma_modale ) / ( seconda_forma_modale_trasposta * M * seconda_forma_modale );

w2_asterisco =sqrt( w2^2 * (Num/Den) )





