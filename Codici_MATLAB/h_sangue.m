% La funzione "h_sangue" consente il calcolo del coefficiente convettivo
% nell'interazione termica con il sangue (convezione forzata)

function h_s = h_sangue(T, T_b_inf, u_inf_b)

T= T +273.15;
T_b_inf= T_b_inf+273.15;
T_film= (T+T_b_inf)/2;
% Dimensione caratteristica
D = 0.1; %[m], lunghezza 


% Conducibilità terminica [W/mK]
x = [37+273.14, 100+273.15 ]; % Punti sull'ascissa (Temperatura K)
y =[0.54, 0.35]; 
% Interpolazione polinomiale di grado 1
degree = 1;
p = polyfit(x, y, degree);
% Sostituzione del valore di x nell'espressione simbolica
k_b = polyval(p, T_film);
k_b = fillmissing(k_b, 'previous');

% Viscosità [m^2/s]
x = [35+273.14, 36+273.14, 37+273.14, 38+273.14, 39+273.14, 40+273.14, 41+273.14, 42+273.14, 80+273.14 ];
y =[2.74, 2.70, 2.65, 2.57, 2.51, 2.43, 2.33, 2.24, 1.20].*10^-6; 
% Interpolazione polinomiale di grado 4
degree = 1;
p = polyfit(x, y, degree);
% Sostituzione del valore di x nell'espressione simbolica
ni_b = polyval(p, T_film);
ni_b = fillmissing(ni_b, 'previous');
% Viscosità superficie
ni_s = polyval(p, T);
ni_s = fillmissing(ni_s, 'previous');
% Viscosità nel sangue
ni_inf = polyval(p, T_b_inf);
ni_inf = fillmissing(ni_inf, 'previous');

% Numero di Prandtl
% Dati Pr [-]
rho_b = 1050;            % [Kg /m^3]
c_b = 3617;              % [J/Kg k]

Pr= (rho_b*ni_b*c_b)/k_b;

Re= u_inf_b*D/ ni_b;


% Flusso laminare
if Re<2e5
    if Pr < 0.6 

        Nu= 1.13*(Pr^0.5)*(Re^0.5);
        h_s = k_b*Nu/D;

    elseif Pr >= 0.6 && Pr<= 10

        Nu= 0.664*(Pr^(1/3))*(Re^0.5);
        h_s = k_b*Nu/D;
    elseif Pr > 10

        Nu= 0.678*(Pr^(1/3))*(Re^0.5);
        h_s = k_b*Nu/D;
    end
    
% Flusso turbolento
else
    if Pr >= 0.7 && Pr<= 380

        Nu= 0.036*(((Re^0.80)*(Pr^0.43)-17400) + 289*(Pr^(1/3)))*(((rho_b*ni_inf)/ (rho_b*ni_s))^(1/4));
        h_s = k_b*Nu/D;
    else
        h_s = 8;
    end

end


end