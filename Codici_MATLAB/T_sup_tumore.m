% La funzione "T_sup_tumore" serve per il calcolo delle temperatura lungo
% la superficie esterna del tumore e per l'impostazione della 
% condizione temporale sui suoi bordi

function [probe]=T_sup_tumore(T,x)

% Interpolazione del vettore delle temperature
degree = 6;
p = polyfit(x, T, degree);

% Sostituzione del valore di x nell'espressione simbolica
T_piu_x1=polyval(p,0.01701);
T_meno_x1=polyval(p,0.01301);

if T_piu_x1>=55 || T_meno_x1>=55
    probe=0;
else
    probe=1;
end

end