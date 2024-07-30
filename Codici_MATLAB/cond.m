% La funzione "cond" consente, conoscendo la posizione nel dominio
% spaziale (quindi i e dx), di definire in quali tessuti vengono a trovarsi
% i nodi i-esimo, (i-1)-esimo e (i+1)-esimo


function [nodo_piu, nodo, nodo_meno] = cond(i, dx)

% Posizione del nodo iesimo
dx= dx*1e3;
pos = (i-1)*(dx);


% CAPILLARE
if pos>=0 && pos<= 0.01
    nodo = 1;
    % Nodo superiore
    if pos+dx<=0.01 
        nodo_piu = 1;
    else
        nodo_piu = 2;
    end
    % Spessori inferiori
    if pos == 0
        nodo_meno = 0;
    else
        nodo_meno = 1;
    end    

% STRATO GHIANDOLARE INF
elseif pos> 0.01 && pos<= 13.01
    nodo = 2;
    % Nodo inferiore
    if pos-dx>= 0.01
        nodo_meno = 2; 
    else
        nodo_meno = 1;
    end
    % Nodo superiore
    if pos+dx<= 13.01
        nodo_piu = 2;
    else
        nodo_piu = 3;
    end

% TUMORE 
elseif pos>13.01 && pos<= 17.01
    nodo = 3;
    % Nodo inferiore
    if pos-dx>= 13.01
        nodo_meno = 3; 
    else
        nodo_meno = 2;
    end
    % Nodo superiore
    if pos+dx<= 17.01
        nodo_piu = 3;
    else
        nodo_piu = 4;
    end

% STRATO GHIANDOLARE SUP
elseif pos>17.01 && pos<= 25.01
    nodo = 4;
    % Nodo inferiore
    if pos-dx>= 17.01
        nodo_meno = 4; 
    else
        nodo_meno = 3;
    end

    % Nodo superiore
    if pos+dx<= 25.01
        nodo_piu = 4;
    else
        nodo_piu = 5;
    end


% SOTTOCUTE
elseif pos>25.01 && pos<= 28.51
    nodo = 5;
    % Nodo inferiore
    if pos-dx>= 25.01
        nodo_meno = 5; 
    else
        nodo_meno = 4;
    end
    % Nodo superiore
    if pos+dx<=28.51
        nodo_piu = 5;
    else
        nodo_piu = 6;
    end

% DERMA
elseif pos>28.51 && pos<= 32.51
    nodo = 6;
    % Nodo inferiore
    if pos-dx>= 28.51
        nodo_meno = 6; 
    else
        nodo_meno = 5;
    end
    % Nodo superiore
    if pos+dx<=32.51
        nodo_piu = 6;
    else
        nodo_piu = 7;
    end

% EPIDERMIDE
elseif pos>32.51 && pos<= 36.015
    nodo = 7;
    % Nodo inferiore
    if pos-dx>=32.51
        nodo_meno = 7; 
    else
        nodo_meno = 6;
    end
    % Spessore superiore
    if pos+dx >= 36.015
        nodo_piu=0;
    else
        nodo_piu = 7;
    end

end

dx= dx*1e-3;





end