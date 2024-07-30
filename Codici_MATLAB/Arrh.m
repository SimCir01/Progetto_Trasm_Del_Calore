% La funzione "Arr" per il calcolo della percentuale di danno all'interno
% di un tessuto generico

function perc_danno= Arrh(T,j, dt, nodi, tessuto)

% Calcolo passo temporale
L = 36.01e-3;                   % Lunghezza del sistema [m] 
x = linspace(0,L,nodi);         % [m]
% Numero minimo di nodi 12 con dx= 0.0327 m
dx = x(2)- x(1);        % [m]

% Contributi spaziali-temporali necessari per la leggi di Arrhenius
Omega_x_t=zeros(size(T,1), j);
Omega_t=zeros(1,j);
% Costanti dei gas, vettore attivazione dell'energia cinetica e vettore
% energia di attivazione
R_g= 8.31; %[J/mK]
vet_A=[5.60e63 1.18e46 1.18e36 1.18e46 1.80e50 1.80e50 1.80e50];
vet_H=[4.30e5 3.02e5 2.38e5 3.02e5 3.27e5 3.27e5 3.27e5];

% Contributi nel tempo
for index_t=1:j
    for i=1:size(T,1)
         [nodo_piu, nodo, nodo_meno] = cond(i, dx);
         if nodo==tessuto
             Omega_x_t(i,index_t)= vet_A(nodo)*exp(-vet_H(nodo)/(R_g*(T(i,index_t)+273)));
         else
             Omega_x_t(i,index_t)=0;
         end
    end
    Omega_t(1,index_t)= sum(Omega_x_t(:,index_t))*dt;
end

Omega= sum(Omega_t);
perc_danno=(1-exp(-Omega))*100;

end
