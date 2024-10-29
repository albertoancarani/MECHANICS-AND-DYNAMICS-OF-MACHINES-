function [yp] = fun04a(t,y,flag,A1,B1,A2,B2)
% dwevo scrivere una function che è una matrice che ha 2 righe 
% una per ogni funzione da integrare

yp=[ A1+B1*y(1) ; A2+B2*y(2) ];


end

