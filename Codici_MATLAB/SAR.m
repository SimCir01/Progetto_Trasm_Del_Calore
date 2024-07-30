% La funzione "SAR" consente di ottenere, scelta in input la posizione 
% spaziale, il corrispettivo valore SAR in base all'andamento di
% riferimento.

function SAR = SAR(depth)

% Dati 
x = [0, 3, 6, 9, 12, 14, 18, 21, 24, 27, 30, 33, 36, 40]*1e-3; % Punti sull'ascissa (Insertion depth in m)
y =[0, 0, 2, 6, 13, 20, 22, 19, 17, 15, 10, 8, 3, 1]; % Punti sull'ordinata (SAR in W/Kg)

% Interpolazione polinomiale di grado 6
degree = 7;
p = polyfit(x, y, degree);

% Sostituzione del valore di x nell'espressione simbolica
SAR = polyval(p, depth);

if SAR<0
    SAR=0;
end

%disp(['Per x = ', num2str(depth), ', y = ', num2str(SAR)]);

end