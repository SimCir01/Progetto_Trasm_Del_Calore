% La funzione "matrix" serve per il calcolo delle matrice A 
% (matrice coefficienti) e b (termini % noti) e 
% il campo di temperatura T(x,t) 

function [A, b, T]= matrix(nodi, tTot, dt,csi, corr, probe, obbiettivo)

%% Dati Tessuti %%

% ARIA
T_inf_c = 20 ; % [K]
h_inf_c = 13.5;     % ma dipendenza dalla T [w/ m^2K]
u_inf_c = 0.2;   % [m/s]

% EPIDERMIDE - DERMA - SOTTOCUTE
rho_ep = 1000;            % [Kg /m^3]
k_ep = 0.37;              % [W/mK]
c_ep = 3391;              % [J/Kg k]
Arr_ep = 1.80e51;         % [1/s]
Delta_H_ep = 3.27e5;      % [J/mol]
u_met_ep = 400;           % [W/m^3]
s_ep = 3.5e-3;               % [m]
s_d = 4e-3;                  % [m]
s_s = 3.5e-3;                % [m]

% STRATO GHIANDOLARE
rho_g = 1041;             % [Kg /m^3]
k_g = 0.33;               % [W/mK]
c_g = 2960;               % [J/Kg k]
Arr_g = 1.18e44;          % [1/s]
Delta_H_g = 3.02e5;       % [J/mol]
u_met_g = 700;            % [W/m^3]
s_g_sup = 8*1e-3;         % [m]
s_g_inf = 13e-3;          % [m]

% TUMORE
rho_t = 1060;           % [Kg /m^3]
k_t = 0.415;            % [W/mK]
c_t = 2726;             % [J/Kg k]
Arr_t = 1.18e36;        % [1/s]
Delta_H_t = 2.38e5;     % [J/mol]
u_met_t = 1400;         % [W/m^3]
s_t = 4e-3;                  % [m]

% CAPILLARE
rho_c = 1102;                   % [Kg /m^3]
k_c = 0.460;                    % [W/mK]
c_c = 3306;                     % [J/Kg k]
Arr_c = 5.60e63;                % [1/s]
Delta_H_c = 4.30e5;             % [J/mol]
u_met_c = 400;                  % [W/m^3]
s_c = 10e-6;                    % [m]

% SANGUE
u_inf_b = 2e-4;          % [W/m^3]
rho_b = 1050;            % [Kg /m^3]
c_b = 3617;              % [J/Kg k]
Arr_b = 5.60e63;         % [1/s]
Delta_H_b = 4.30e5;      % [J/mol]
% Parametri sangue c.c.
T_inf_bA = 37 ;         %[C]       
A_b = 1.2 ;             %[C]
periodo_b = 0.8;


% DEFINIZIONI PASSO SPAZIALE E TEMPORALE
% Settaggio parametri per la soluzione numerica nello spazio
L = 36.01e-3;                   % Lunghezza del sistema [m] 
x = linspace(0,L,nodi);         % [m]
% Numero minimo di nodi 12 con dx= 0.0327 m
dx = x(2)- x(1);        % [m]
dx1 = L/(nodi-1);       % [m]

if abs(dx-dx1)>0
    display('error')
    return
end 

% Numero step (arrotondato per eccesso)
steps = ceil(tTot/dt);
% Vettore dei tempi
t= linspace(0,tTot, steps+1);
% Settaggio parametri per la soluzione numerica nel tempo
% Condizione iniziale
Tini =37 ; % [C]
T = ones(nodi,steps) * Tini;
T_err = ones(nodi,steps) * Tini;

% Vettori utili per la costruzione delle matrici
% Vettore delle densità
vet_rho=[1102; 1041; 1060; 1041; 1000; 1000; 1000]';
% Vettore dei calori specifici
vet_c=[3306; 2960; 2726; 2960; 3391; 3391; 3391]';
% Vettore dei termini metabolici
u_metabolici=[400; 700; 1400; 700; 400; 400; 400]';
% Vettore delle conducibilità
vet_k= [0.460; 0.33 ; 0.415; 0.33; 0.37; 0.37 ;0.37 ]';

% COSTRUZIONE MATRICE TERMINI NOTI E DEI COEFFICIENTI
% Matrice dei coefficienti A
A = zeros(nodi,nodi, steps);
% Vettore dei termini noti b
b = zeros(nodi,1,steps);

for j =1:steps
    
    it=1;
    err_1(it)=1;
    err_end(it)=1;
    tol = 1e-8;

    while err_1(it)>tol && err_end(it)>tol

    % Perfusione sanguigna
    omega_b = zeros(nodi,1);
     for i=1:nodi
        omega_b(i) = 0.000021*T(i,j) + 0.00035;
     end   

    % Temperatura del sangue
    T_b_inf = T_inf_bA + A_b*sin((2*pi/periodo_b)*((j-1)*dt));
    

    % Coefficiente Convettivo sangue
    h_s(1,j,it+1)= h_sangue(T_err(1,j,it), T_b_inf, u_inf_b);
    % Coefficiente Convettivo sangue
    h_a(1,j,it+1) = h_aria(T_err(end,j,it), T_inf_c);
    

    % Nodo interfaccia superiore (EPIDERMIDE-ARIA) 
    Fo_ep = (k_ep*dt)/(rho_ep*c_ep*(dx^2));
    Bi_ep = h_a(j)*dx/k_ep;

    A(end, end-1,j) = 2*Fo_ep;
    A(end, end,j) = -2*Fo_ep*(1 + Bi_ep + (omega_b(end)*rho_b*c_b*dx^2)/(2*k_ep));
    b(end,1,j) = 2*Fo_ep *(Bi_ep*T_inf_c + (dx^2/(2*k_ep))*(u_met_ep + corr*probe*SAR(0)*rho_ep  + omega_b(end)*rho_b*c_b*T_b_inf));

    % Nodi interni (STRATO-STRATO)
    for ir = 2:nodi-1    
        [nodo_piu, nodo, nodo_meno] = cond(ir, dx);
        Sigma = sigma(ir,dx,dt);

        for ic = 1:nodi
            if ir==ic
                A(ir,ic,j) = - dt*Sigma*(vet_k(nodo_meno)/dx + vet_k(nodo_piu)/dx + omega_b(ir)*rho_b*c_b*dx);
            elseif ir-ic==1
                A(ir,ic,j) =dt* Sigma*vet_k(nodo_meno)/dx;
            elseif ir-ic==-1
                A(ir,ic,j) = dt*Sigma*vet_k(nodo_piu)/dx;
            end
        end
 
        b(ir,1,j) = dt*Sigma*dx*(u_metabolici(nodo) + corr*probe*SAR(L-(ir-1)*dx)*vet_rho(nodo) + omega_b(ir)*rho_b*c_b*T_b_inf);
    end

    % Nodo interfaccia inferiore (SANGUE-CAPILLARE)
   [nodo_piu, nodo, nodo_meno] = cond(1, dx);
   Sigma = sigma(1,dx,dt);

    A(1, 1,j) = - dt*Sigma*(vet_k(nodo_piu)/dx + h_s(j) + (omega_b(1)*rho_b*c_b*dx)/2);
    A(1, 2,j) = dt*Sigma*(vet_k(nodo_piu)/dx);
    b(1,1,j) = dt*Sigma*( h_s(j)*T_b_inf + (dx/2)*( u_met_c + corr*probe*SAR(L)*rho_c + omega_b(1)*rho_b*c_b*T_b_inf));

    if j==1 && csi==0 && (obbiettivo == 0 || obbiettivo == 1 )
            dt_critico= zeros(nodi,1);
            for i=1:nodi
                dt_critico(i,1)= -dt/A(i,i,j);
            end
            dt_limite= min(dt_critico); 
            while dt>dt_limite
                fprintf('\rERRORE INSTABILITA''\rPasso temporale troppo grande \r')
                dt = input('Inserire un passo temporale più piccolo:'); 
            end 
    end 

% IMPOSTAZIONE SISTEMA MATRICIALE  
% Matrice dei coefficienti e dei termini noti
Am = speye(size(A(:,:,j))) - A(:,:,j)*csi;
bm = speye(size(A(:,:,j)))*T(:,j) + (1-csi)*A(:,:,j)*T(:,j) + b(:,1,j); 
% Risoluzione campo di temeperatura
T(:,j+1) = pinv(Am(:,:))*bm(:,1);
T_err(:,j+1,it+1) = pinv(Am(:,:))*bm(:,1);
% Condizione da rispettare per spegnere il probe raggiunti i 55° lungo la
% superficie del tumore
if obbiettivo ~=1
    probe = T_sup_tumore(T(:,j+1), x);
end
err_1(it+1) = abs(T_err(1,j,it+1)-T_err(1,j,it));
err_end(it+1) =  abs(T_err(end,j,it+1)-T_err(end,j,it));
it = it+1;
end
    
end

end
