clc, clear

%% dati

h=8; % alzata al netto delle rampe di ripresa gioco [mm]
alpha=240; %[deg]
h1= 0.3; %gioco [mm]
u=270; % max vel d'urto [mm/s]
n=6000; % massimo regime di rotazione motore [rpm]
theta1=38; %[deg]
theta2=42; %[deg]
deltamin=5; %minimo raggio di curvatura del profilo camma [mm]
Fmin=140; %[N]
k=8.8; % [N/mm]
massa= 0.19; %[kg] massa equipaggio mobile della valvola


%% calcoli

%theta
tq= alpha/4-theta2; % [deg]
thetaR= alpha/4-theta1; % [deg]
thetaS= alpha/4; % [deg]

omega_camma= ( n/2 )* (2*pi/60 ) ; %rad/s
% omegamax = 2pi*n/(60*2) = pi*n/60 [rpm] --> 
% omegamax = pi*n/60 * 180/pi=3n [deg/s]
omegamax = 3*n; %[deg/s]



b0=h; %coefficiente legge di alzata [mm]


% yprimo= u/omegamax %[mm/s / deg/s] --> [mm/deg]
yprimo=u/omegamax; %[mm/deg]

thetaP = -h1/yprimo; %[deg]

format shortE
A=[20*tq^3 12*tq^2 6*tq  0 0; thetaR^5 thetaR^4 thetaR^3 -(thetaR-thetaS)^4 -(thetaR-thetaS)^2; 5*thetaR^4  4*thetaR^3 3*thetaR^2 -4*(thetaR-thetaS)^3 -2*(thetaR-thetaS); 20*thetaR^3 12*thetaR^2 6*thetaR -12*(thetaR-thetaS)^2 -2; 60*thetaR^2 24*thetaR 6 -24*(thetaR-thetaS) 0];
b=[0; b0-(yprimo*thetaR); -yprimo; 0; 0];
x=A\b;

a5=x(1);
a4=x(2);
a3=x(3);
b4=x(4);
b2=x(5); 
a0=0;
a1=yprimo;
a2=0;


%% andamento di y in funzione di theta 

% lineare da P a O 
% pol 5 grado da O a R
% Pol 4 grado da R a S 
% speculare da S in poi

% P - O
% m = y' = u/ omegamax
% y = mx+q 
% q = 0
% ciclo for che fa variare theta da thetaP a thetaO(0,0)
m= yprimo;
thetaPO = thetaP : 0.1 : 0;
yPO= thetaPO.*m;

% O - R
% y = a5*theta^5 + a4*theta^4 + a3*theta^3 + a2*theta^2 + a1*theta + a0
thetaOR = 0 : 0.1 : thetaR ;
a0vettore= ones(1,length(thetaOR) ) *a0;
yOR =  a5.*thetaOR.^5 + a4.*thetaOR.^4 + a3.*thetaOR.^3 + a2.*thetaOR.^2 + a1.*thetaOR + a0vettore ;

% R - S
% y= b4*(thetaRS-thetaS)^4 + b2*(thetaRS-thetaS)^2+ b0
thetaRS = thetaR : 0.1 : thetaS;
b0vettore= ones(1,length(thetaRS) ) *b0;
thetaSvettore= ones(1,length(thetaRS) ) *thetaS;
yRS= b4.*(thetaRS-thetaSvettore).^4 + b2.*(thetaRS-thetaSvettore).^2+ b0vettore;

thetaPS = [thetaPO thetaOR thetaRS];
yPS= [yPO yOR yRS];
% simmetria
% chiamo U il valore speculare a P

thetaSasse= ones(1,length(thetaPS) ) *thetaS;
distanza = thetaSasse - thetaPS;
thetaUS = thetaPS + distanza.*2;
yUS=yPS;

figure(1)

plot(thetaPS,yPS)
hold on 
plot(thetaUS,yUS)
grid on
title('andamento di y in funzione di theta ')
xlabel('theta [deg]')
ylabel('y(theta) [mm]')
% figure(1)
% plot(thetaPO,yPO);
% hold on
% plot(thetaOR,yOR);
% hold on 
% plot(thetaRS,yRS);

%% andamento di y' in funzione di theta 

% P - O
%y'=dy/dtheta
%yPO= thetaPO.*m
%y'=m
coefvelPO=ones(1,length(thetaPO) )*m;


% O - R
% y = a5*theta^5 + a4*theta^4 + a3*theta^3 + a2*theta^2 + a1*theta + a0
% y'= 5*a5*theta^4 + 4*a4*theta^3 + 3*a3*theta^2 + 2*a2*theta + a1

a1vettore= ones(1,length(thetaOR) ) *a1;
coefvelOR= (5*a5).*thetaOR.^4 + (4*a4).*thetaOR.^3 + (3*a3).*thetaOR.^2 + (2*a2).*thetaOR + a1vettore;

% R - S
% y= b4*(thetaRS-thetaS)^4 + b2*(thetaRS-thetaS)^2+ b0
% y'= 4*b4*(thetaRS-thetaS)^3 + 2*b2*(thetaRS-thetaS)
coefvelRS = (4*b4).*((thetaRS-thetaSvettore).^3) + (2*b2).*(thetaRS-thetaSvettore);
% unisco le 3 parti
coefvelPS= [coefvelPO coefvelOR coefvelRS];

% antisimmetria 
% chiamo U il valore speculare a P

coefvelUS=(-1).*coefvelPS;


figure(2)
plot(thetaPS,coefvelPS)
hold on 
plot(thetaUS,coefvelUS)
grid on
title('andamento di yprimo in funzione di theta ')
xlabel('theta [deg]')
ylabel('yprimo(theta) [mm/deg]')
% la derivata è antisimmetrica

% plot(thetaPO,coefvelPO)
% hold on 
% plot(thetaOR,coefvelOR)
% hold on
% plot(thetaRS,coefvelRS)

%% andamento di y'' in funzione di theta 

% derivata di una antisimmetrica è simmetrica


% P - O
%yPO= thetaPO.*m
%y'=m
%y''=0

coefaccPO=zeros(1,length(thetaPO) );

% O - R
% y = a5*theta^5 + a4*theta^4 + a3*theta^3 + a2*theta^2 + a1*theta + a0
% y'= 5*a5*theta^4 + 4*a4*theta^3 + 3*a3*theta^2 + 2*a2*theta + a1
% y''= 20*a5*theta^3 + 12*a4*theta^2 + 6*a3*theta + 2*a2
a2vettore= ones(1,length(thetaOR) ) *a2;
coefaccOR= (20*a5).*thetaOR.^3 + (12*a4).*thetaOR.^2 + (6*a3).*thetaOR+ 2.*a2vettore;

% R - S
% y= b4*(thetaRS-thetaS)^4 + b2*(thetaRS-thetaS)^2+ b0
% y'= 4*b4*(thetaRS-thetaS)^3 + 2*b2*(thetaRS-thetaS)
% y''=  12*b4*(thetaRS-thetaS)^2 + 2*b2
b2vettore= ones(1,length(thetaRS) ) *b2;
coefaccRS = (12*b4).*((thetaRS-thetaSvettore).^2) + 2.*b2vettore;



% unisco le 3 parti
coefaccPS= [coefaccPO coefaccOR coefaccRS];

% simmetria 
% chiamo U il valore speculare a P

coefaccUS=coefaccPS;


figure(3)
plot(thetaPS,coefaccPS)
hold on 
plot(thetaUS,coefaccUS)
grid on
title('andamento di ysecondo in funzione di theta ')
xlabel('theta [deg]')
ylabel('ysecondo(theta) [mm/deg^2]')
%% determinare il raggio base della camma 

% il minimo raggio di curvatura è pari deltamin

% rho = Rb + h1 + y + y'' >= deltamin
% z(theta)= y + y''
% cerco il minimo di rho --> d( z(theta) ) / dtheta = 0

% ho rischio maggiore di avere un rho basso nella zona che è più appuntita
% perchè la circonferenza osculatrice è più piccola. La zona più appuntita
% è quella rappresentata dal polinomio di quarto grado R - S
% chiamo yb yRS per non sovrascrivere le variabili

% yb= b4*(thetaRS-thetaS)^4 + b2*(thetaRS-thetaS)^2+ b0 
% yb'= 4*b4*(thetaRS-thetaS)^3 + 2*b2*(thetaRS-thetaS)
% yb''=  12*b4*(thetaRS-thetaS)^2 + 2*b2
% yb'''= 24*b4*(thetaRS-thetaS)

% d(yb + yb'')/dtheta = yb'+yb'''= 
% =4*b4*(thetaRS-thetaS)^3 + 2*b2*(thetaRS-thetaS) + 24*b4*(thetaRS-thetaS) 
% = (thetaRS-thetaS)*[ 4*b4*(thetaRS-thetaS)^2 + 2*b2 + 24*b4 ] 
% devo imporre la derivata nulla ho due soluzioni:

% 1) thetaRS = thetaS

% 2) 4*b4*(thetaRS-thetaS)^2 + 2*b2 + 24*b4 = 0
% 4*b4*thetaRS^2 + 4*b4*thetaS^2 -8*b4*thetaRS*thetaS + 2*b2 + 24*b4 = 0
% [ 4*b4*thetaRS^2 ] + [ - 8*b4*thetaS*thetaRS ] + 
% +[ b4(4*thetaS^2 + 24)+ 2b2 ]= 0

% thetaRS1,2 = { 8*b4*thetaS +- sqrt{64*(b4^2)*(thetaS^2) - 
% - 4*[(b4*(4*thetaS^2 + 24) + 2*b2 )]*4b4} } / (8*b4)

% delta = 64*(b4^2)*(thetaS^2) - 4*[(b4*(4*thetaS^2 + 24) + 2*b2 )]*4b4 =
% = 64*(b4^2)*(thetaS^2) - 64*(b4^2)*(thetaS^2) - 384*b4^2 - 32*b2*b4 =
% = - 384*b4^2 - 32*b2*b4 
% il delta è negativo quindi avrei valori non reali perciò escludo la
% seconda soluzione

% rhomin = Rb + h1 + z(thetaRS = thetaS) = deltamin
% Metto a sistema:
% 1) Rb =  deltamin - h1 - z(thetaRS = thetaS)  
% 2) z = y + y''
% 3) yb=b4*(thetaRS-thetaS)^4 + b2*(thetaRS-thetaS)^2+ b0 
% 4) yb''=  12*b4*(thetaRS-thetaS)^2 + 2*b2

% Rb =  deltamin - h1 - yb(thetaS) - yb''(thetaS)  =
% Rb = deltamin - h1 - b0 - 2*b2
% [mm] - [mm] - [mm] - [mm/deg^2]
% devo moltiplicare b2 per *(360/2pi)^2 per avere la stessa unità di misura

Rb = deltamin - h1 - b0 - 2*b2*( (360/(2*pi))^2 ); % [mm]

%%  disegnare il profilo della camma in scala 2:1

% conviene utilizzare le cordinate polari OE , theta + delta
% OE^2 = OD^2 + OC^2
% OC = y'
% OD = Rb + h1
% tan(delta) = OC / OD

% devo tenere conto di 3 zone:

% 1) U - P (198°)
% RAGGIO BASE : OD = Rb ; OC = 0;
figure(4)
OD1=Rb;
OC1=0;
OE1 = sqrt(OD1^2+OC1^2); %OE^2 = OD^2 + OC^2
delta1=atan(OC1/OD1); % arctan(0/Rb)
%thetaU=thetaS+(thetaS-thetaR)+(thetaR-thetaP);
%thetaN=2*thetaS;%thetaS+(thetaS-thetaR)+(thetaR-thetaP)+thetaP=2*thetaS;
thetaU=thetaS+(thetaS-thetaR)+(thetaR-thetaP);
thetaPnormalizzato=thetaP+360; 
thetaUP=thetaU : 0.1 : thetaPnormalizzato;
theta1=thetaUP*pi/180; %[rad]
OE1vettore= ones(1,length(theta1) )*OE1;
delta1vettore= ones(1,length(theta1) )*delta1;
angolo1=theta1+delta1;
% 2) P - O + N - U (18°+18°)
% TRATTO DI RECUPERO GIOCO : 
% OD = Rb + h1 + y0' 
% y0'= OC = costante

% P - O
Rbvettore= ones(1,length(thetaPO) ) *Rb;
h1vettore= ones(1,length(thetaPO) ) *h1;
OD2a= Rbvettore + h1vettore + coefvelPO.*thetaPO;
OC2a= coefvelPO*180/pi;
delta2a=atan(OC2a./OD2a);
OE2a = sqrt(OD2a.^2+OC2a.^2);
theta2a = thetaPO*pi/180; %[rad]
angolo2a= theta2a + delta2a;

% N - U

%thetaT=2*(thetaS-thetaR) + thetaR;
%thetaTvettore=ones(1,length(thetaOR) ).*thetaT;
%thetaTN = thetaOR  + thetaTvettore;

thetaO=0;
thetaN= thetaS +(thetaS - thetaO);
thetaNU = thetaPO + (thetaN - thetaP);% non va normalizzato
%coefvelTN=(-1).*coefvelOR;
coefvelNU=(-1).*fliplr(coefvelPO); % non va normalizzato
OD2b= Rbvettore + h1vettore + coefvelNU.*thetaNU;  
OC2b= coefvelNU*180/pi;
OE2b = fliplr(OE2a);
delta2b=atan(OC2b./OD2b);
theta2b = thetaNU*pi/180; %[rad]
angolo2b= theta2b + delta2b;

%  continuità
OE12= [OE2b OE1vettore OE2a];
angolo12 = [angolo2b angolo1 angolo2a];

% 3) O - N (126°)
% TRATTO ATTIVO :
% OD = Rb + h1 + y
% OC = y'

% yST=fliplr(yRS);
% yRT=[yRS yST];
yOS=[yOR yRS(2:end)];
ySN=fliplr(yOS);
yON=[yOS ySN(2:end)];
% coefvelST=fliplr(coefvelRS).*(-1);
% coefvelRT= [coefvelRS coefvelST];
coefvelOS=[coefvelOR coefvelRS(2:end)];
coefvelSN=fliplr(coefvelOS).*(-1);
coefvelON= [ coefvelOS coefvelSN(2:end)];
Rbvettore3= ones(1,length(yON) ).*Rb;
h1vettore3= ones(1,length(yON) ).*h1;
% thetaST= thetaRS + (thetaS - thetaR);
% thetaRT=[thetaRS thetaST];
thetaON=thetaO:0.1:thetaN;
OD3= Rbvettore3 + h1vettore3 + yON; 
OC3= coefvelON*180/pi;
delta3=atan(OC3./OD3);
OE3 = sqrt(OD3.^2+OC3.^2);
theta3 = thetaON*pi/180; %[rad]
angolo3= theta3 + delta3;

% continuità
OE= [OE3 OE12];
angolo = [angolo3 angolo12];

% figure(4)
% polarplot(angolo1,2*OE1vettore,'LineWidth',2);
% hold on
% polarplot(angolo2a,2*OE2a,'LineWidth',1);
% polarplot(angolo12,2*OE12,'LineWidth',2);
% hold on
% polarplot(angolo2b,2*OE2b,'LineWidth',1);
% hold on
% polarplot(angolo3,2*OE3,'LineWidth',3);

figure(4)
polarplot(angolo,2*OE,'LineWidth',3);
title('profilo camma in scala 2:1')

%% bicchierino

% determino la massima distanza dall'asse valvola del contatto
% camma-bicchierino 

%devo avere sempre contatto fra camma e bicchierino quindi il bicchierino
%deve avere un ampiezza almeno uguale a y'max

OCmax=max([OC1 OC2a OC2b OC3]); %[mm] 

% ricorda : OC = y'=dy/dtheta quindi devo correggere l'unità di misura se
% voglio una distanza in mm


% determino il diametro minimo d del bicchierino essendo pari a d/2 lo
% spessore della camma

%d/2= y'max / cos(30) 
%y'max= OCmax

d= 2 * OCmax / cos(pi/6);

%% precarico molla

% F = Fel + Fin  + T;
% T = Precarico
% Fel = Felastica
% Fel = - ky  
% Fin = Finerziale
% Fin = - m*d^2y/dt^2 = -m*coefacc*omega^2 

% Cerco la Fmin Fmin =  (Fel + Fin)min  + T;
% il minimo si raggiunge sul polinomio yb (zona R - S )
% (Fel + Fin)min = - (ky + m*coefacc*omega^2 )
% [N/mm] * [mm] + [kg * mm / s^2 * 1/1000]

% devo trovare il minimo di - (k*yRS + m*coefaccRS*omegamax^2 )
% d( (Fel + Fin) ) / dtheta = 0
% d ( (k*yRS + m*coefaccRS*omegamax^2 ) ) = 0
% k*yb' + m* yb''' * omegamax^2 = 0

% ricavate precedentemente (zona R - S) : 
% yb= b4*(thetaRS-thetaS)^4 + b2*(thetaRS-thetaS)^2+ b0 
% yb'= 4*b4*(thetaRS-thetaS)^3 + 2*b2*(thetaRS-thetaS)
% yb''=  12*b4*(thetaRS-thetaS)^2 + 2*b2
% yb'''= 24*b4*(thetaRS-thetaS)
 
% k*(4*b4*(thetaRS-thetaS)^3 + 2*b2*(thetaRS-thetaS) ) +
% + m*omegamax^2*24*b4*(thetaRS-thetaS)= 0 

%(thetaRS-thetaS)*[4k*b4*(thetaRS-thetaS)^2 +2k*b2 + 24m*omegamax^2*b4] = 0
% ho due soluzioni.( thetaRS = theta sono letteralmente la stessa cosa )

% 1) theta = thetaS

% 2) 4k*b4*(theta-thetaS)^2 +2k*b2 + 24m*omegamax^2*b4 = 0
% 4k*b4*theta^2 + 4k*b4*thetaS^2 - 8k*b4*theta*thetaS + 2k*b2 + 24m*omegamax^2*b4 = 0
% 4k*b4*(theta^2) - 8k*b4*thetaS*(theta) + (4k*b4*thetaS^2 + 2k*b2 + 24m*omegamax^2*b4 ) = 0

% theta1,2=  8k*b4*thetaS +-
% +- sqrt (64*k^2*b4^2*thetaS^2 -4*(4k*b4*thetaS^2+2k*b2+24m*omegamax^2*b4 )*(4k*b4))

% delta = 64*k^2*b4^2*thetaS^2 -4*(4k*b4*thetaS^2+2k*b2+24m*omegamax^2*b4)*(4k*b4) =
% = 64*k^2*b4^2*thetaS^2 - 4^3*b4^2*k^2*thetaS^2 - 4^2*2*k^2*b2*b4 -
% - 4^2*24*m*omegamax^2*b4^2 =
% =  - 32*k^2*b2*b4 - 384*m*omegamax^2*b4^2 < 0 
% delta < 0 --> le soluzioni non sono reali
% utilizzo perciò la 1 soluzione

% 1) theta = thetaS
% T = Fmin - k*yb(thetaS) - m*coefaccRS(thetaS)*omegamax^2

% T = Fmin - k*bo - m*2b2*omegamax^2

% Fmin = [N]
% k = [N /mm]
% bo = [mm]
% omegamax=3*n;  [mm/deg^2]
% b2 = [deg^2/s^2]
% m= massa equipaggio mobile della valvola [kg]

% [N] - [N /mm] * [mm] - [kg]*[deg^2/s^2]*[mm/deg^2]
% [N] - [N ] - [kg*mm/s^2]
% [N] - [N ] - [N*1/1000]

T = Fmin - k*b0 - (massa*2*b2*(omegamax^2))/1000;

