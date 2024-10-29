clc
clear
%% DATI

u=4;

%coppie rotoidali [mm]
O6= [ 0   , 0  ];
O1= [ 500 , 0  ];
O2= [ 385 , 15 ];
O3= [ 120 , -78 ];
O4= [ 80  , 115];
O5= [ -15 , 50 ];

%lunghezze aste [mm]
AB=100;
O2B=75-2*u;
BC=100;
CD=125;
O3D=150;
O3E=150;
EF=50;
O1A=42;
O4F=75;
O4G=85;
O4I=115;
GH=75;
O5H=72;
IL=50;
O6L=42;

%angoli [rad]
ang_DO3E=60*pi/180;
ang_FO4G=195*pi/180;
ang_FO4I=265*pi/180;
ang_ABC=180*pi/180;
ang_AO1O6=30*pi/180;

%% GRUPPO DRIVER cerco A

vett_O1A=O1A.*[-cos(ang_AO1O6) , sin(ang_AO1O6)];
A=O1+vett_O1A;

%% DIADE RRR A O2 AB O2B +1

O2A = distanza(O2,A);
%sdr locale
lambda1 = 1/2 * (1 + (AB^2 - O2B^2)/(O2A^2) );
mu1 = sqrt ( ( AB^2 / O2A^2 ) - lambda1^2 );
vettloc_AB = [lambda1*O2A , mu1*O2A ]; 

%sdr con origine in A

beta = acos(vettloc_AB(1)/AB ); %[rad]
alpha = atan( (A(2) - O2(2) )/ O2A ); %[rad]
gamma = alpha + beta ;
vett_AB= AB * [cos(gamma) , sin(gamma)];

B=A-vett_AB;

%% ASTA ABC

AC= AB+BC;
vett_AC= AC * [cos(gamma) , sin(gamma)];
C = A - vett_AC;

%% DIADE RRR O3 C O3D CD +1

teta=  - atan( ( C(2) - O3(2) )/ ( C(1) - O3(1) ) ); % [rad]
O3C = ( C(1) - O3(1) ) / cos(teta) ;

%sdr locale

lambda2 = 1/2 * (1 + (O3D^2 - CD^2)/(O3C^2) );
mu2 = sqrt ( ( O3D^2 / O3C^2 ) - lambda2^2 );
vettloc_O3D = [lambda2*O3C , mu2*O3C ]; 

%sdr originale

phi=atan( vettloc_O3D(2) / vettloc_O3D(1) ); % [rad]
delta= pi/2 - phi + teta; % [rad]
D = O3 + O3D * [sin(delta) , cos(delta)];
%nota bene sin e cos invertiti vedi disegno

%% TRIANGOLO EQUILATERO O3DE

psi = phi - teta;
vett_O3E=O3E*[cos(pi/3+psi) , sin(pi/3+psi)];
E=O3+vett_O3E;

%% DIADE RRR O4 E O4F EF +1

EO4=distanza(O4,E);

lambda3 = 1/2 * (1 + (O4F^2 - EF^2)/(EO4^2) );
mu3 = sqrt ( ( O4F^2 / EO4^2 ) - lambda3^2 );

Fx= lambda3*( E(1) - O4(1) ) + mu3* ( O4(2) - E(2) ) + O4(1);
Fy= lambda3*( E(2) - O4(2) ) + mu3* ( E(1) - O4(1) ) + O4(2);
F = [Fx , Fy];

%% DIADE RRR F O4 FI O4I +1

chi=2*pi - ang_FO4I;
IF= sqrt( O4F^2 + O4I^2 -2*O4F*O4I*cos(chi) );

lambda4 = 1/2 * (1 + (IF^2 - O4I^2)/(O4F^2) );
mu4 = sqrt ( ( IF^2 / O4F^2 ) - lambda4^2 );

Ix= lambda4*( O4(1) - F(1) ) + mu4* ( F(2) - O4(2) ) + F(1);
Iy= lambda4*( O4(2) - F(2) ) + mu4* ( O4(1) - F(1) ) + F(2);
I = [Ix , Iy];

%% DIADE RRR O6 I O6L IL +1

O6I=distanza(O6,I);

lambda5 = 1/2 * (1 + (O6L^2 - IL^2)/(O6I^2) );
mu5 = sqrt ( ( O6L^2 / O6I^2 ) - lambda5^2 );

Lx= lambda5*( I(1) - O6(1) ) + mu5* ( O6(2) - I(2) ) + O6(1);
Ly= lambda5*( I(2) - O6(2) ) + mu5* ( I(1) - O6(1) ) + O6(2);
L = [Lx , Ly]


