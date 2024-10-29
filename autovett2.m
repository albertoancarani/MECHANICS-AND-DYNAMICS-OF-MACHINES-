clc
clear

m=2;
k=2*100*(pi^2);
A=100;
%% sistema a 3 gdl


M=[ A*m 0 0; 0 m 0; 0 0 m];
K=[k -k 0; -k 2*k -k; 0 -k k];

[V,D]= eig(M\K);
om2=diag(D);
% om2Hz=om2/(2*pi);
[om2,ind]=sort(om2);
Vs=V(:,ind);
%normalizzazione
X1=V(:,1);
X2=V(:,2);
X3=V(:,3);
m1=(max(abs(X1)));
m2=max(abs(X2));
m3=max(abs(X3));
X1norm=X1./m1;
X2norm=X2./m2;
X3norm=X3./m3;

%% sistema a 2 gdl


M2=[ m 0; 0 m];
K2=[2*k -k; -k k];

[V2,D2]= eig(M2\K2);

om22=diag(D2);
% om2Hz2=om22/(2*pi);
[om2s,ind2]=sort(om22);
V2s=V2(:,ind2);

