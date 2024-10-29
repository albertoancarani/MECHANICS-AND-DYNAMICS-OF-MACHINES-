clc
clear

m=2;
k=2*100*(pi^2);
w1= sqrt(k/m);
w2= sqrt(3*k/m);

M=[ m 0 ; 0 m];
K=[2*k -k ;-k 2*k];

[V,D]= eig(K,M);
om2=diag(D);
om2Hz=om2/(2*pi);
[om2s,ind]=sort(om2);
Vs=V(:,ind);
