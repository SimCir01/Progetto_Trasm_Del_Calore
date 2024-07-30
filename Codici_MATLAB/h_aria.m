% La funzione "h_aria" consente il calcolo del coefficiente convettivo
% nell'interazione termica con l'aria (convezione naturale)

function h_a = h_aria(T, T_inf_c)

D = 1; %[m], raggio esterno

% Conversione della temperatura da Celsius a Kelvin
T = T +273.15;
T_inf_c= T_inf_c+273.15;
% Calcolo temperatura di film
T_film= (T+ T_inf_c)/2;
% Numero di Prandtl
% Dati Pr [-]
x = linspace(290,400, 12); % Punti sull'ascissa (Temperatura K)
y =[0.715, 0.713, 0.712, 0.710, 0.710, 0.708, 0.707, 0.706, 0.705, 0.705, 0.704, 0.704]; % Punti sull'ordinata (SAR in W/Kg)
% Interpolazione polinomiale di grado 6
degree = 4;
p = polyfit(x, y, degree);
% Sostituzione del valore di x nell'espressione simbolica
Pr = polyval(p, T_film);
Pr = fillmissing(Pr, 'previous');

% Conducibilità terminica [W/mK]
y =[2.538, 2.614, 2.687, 2.759, 2.830, 2.900, 2.970, 3.039, 3.107, 3.173, 3.239, 3.305]*1e-2; % Punti sull'ordinata (SAR in W/Kg)
% Interpolazione polinomiale di grado 4
degree = 4;
p = polyfit(x, y, degree);
% Sostituzione del valore di x nell'espressione simbolica
k = polyval(p, T_film);
k = fillmissing(k, 'previous');

% Fattore moltipicativo [1/Km^3]
y =[154.032, 131.092, 113.723, 98.488, 85.772, 74.921, 65.983, 58.015, 51.328, 45.557, 40.622, 36.271]*1e6; % Punti sull'ordinata (SAR in W/Kg)
% Interpolazione polinomiale di grado 4
degree = 4;
p = polyfit(x, y, degree);
% Sostituzione del valore di x nell'espressione simbolica
gbv = polyval(p, T_film);
gbv = fillmissing(gbv, 'previous');

Gr= gbv*(T- T_inf_c)*D;
Ra= Gr*Pr;

% Numero di Nusselt ----> Trovare la relazione giusta 
if Ra>2e8

    Nu= 0.14*(Ra^(1/3));
    h_a= Nu*k/D;

elseif Ra >1e6 && Ra >1e11

    Nu= 0.58*(Ra^(1/5));
    h_a= Nu*k/D;

else
    h_a= 13.5;

end




end