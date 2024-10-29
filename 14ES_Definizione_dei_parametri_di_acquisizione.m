clc, clear

%% TESTO

% Si vogliono effettuare rilievi sperimentali di vibrazione sulla 
% struttura rappresentata schematicamente in Fig. 26. 
% La misura va condotta all’interno del campo di frequenze 0 ÷ f∗
% e, ai fini dell’analisi, occorre ottenere una risoluzione spettrale 
% massima pari a ∆f.

% Determinare:
% 1. La frequenza di taglio del filtro passa basso anti-aliasing
% 2. La frequenza di campionamento minima
% 3. Il numero di punti da elaborare tenendo conto che si desidera
% utilizzare l’algoritmo FFT (Fast Fourier Transform)

%% DATI
u=4; %ultimo numero matricola
f_asterisco = 4000 - 200 * u; %[Hz] estremo superiore del campo delle frequenze
deltaf = 10 + 0.5 * u^2; %[Hz] massima risoluzione spettrale

%% FREQUENZA DI TAGLIO FILTRO PASSA BASSO

% La misura va condotta all'interno del campo del campo di frequenze 
% [0;f*] Hz. Selezione un filtro passa basso con frequenza di taglio Ft= f*

Ft = f_asterisco % frequenza di taglio filtro passa basso [Hz]

%% FREQUENZA DI CAMPIONAMENTO MINIMA

% la frequenza di campionamento minima deve essere almeno 2 volte la
% frequenza massima (teorema di Shannon) per stare sul sicuro la prendo 2.5
% volte superiore

Fs_minima = 2.5 * Ft % [Hz]

%% NUMERO DI PUNTI DA ELABORARE PER FFT

T_asterisco=1/deltaf; %tempo di osservazione 
deltat=1/Fs_minima; %intervallo di campionamento
N_perfetto= T_asterisco/deltat; %numero di campioni necessari

% Ricordo che l'algoritmo FFT elabora solo potenze di 2

vettore_potenze_2 = [1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 ];
indice=find(vettore_potenze_2 >= N_perfetto);
N= vettore_potenze_2(min(indice)) %numero di campioni minimo per usare FFT

