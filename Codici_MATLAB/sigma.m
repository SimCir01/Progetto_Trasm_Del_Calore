% La funzione "Sigma" consente, conoscendo la posizione nel dominio
% spaziale (quindi i e dx), di definire gli eventuali spessori
% dei nodi (i-1)-esimo e (i+1)-esimo, in modo alla fine
% da poter ricavare il parametro Sigma necessario nei bilanci

function Sigma = sigma(i,dx,dt)

% Vettori degli pessori
s_piu = zeros(7,1);
s_meno = zeros(7,1);

% Vettore delle densità
vet_rho=[1102; 1041; 1060; 1041; 1000; 1000; 1000];
% Vettore dei calori specifici
vet_c=[3306; 2960; 2726; 2960; 3391; 3391; 3391];


% Posizione del nodo iesimo
dx = dx/2;
pos = (i-1)*(dx*1e3);

% CAPILLARE
if pos>=0 && pos<= 0.01
    % Spessori superiori
    if pos+dx<0.01 
        s_piu(1) = (pos + dx)- pos;
    elseif pos==0.01
        s_piu(1)=0;
        s_piu(2)= (pos+dx)-pos;
    elseif pos+dx>0.01 
        s_piu(1) = 0.01;
        if pos+dx <= 13.01
            s_piu(2) = (pos+dx) - 0.01;
        else
            s_piu(2) = 13;
            if pos+dx<=17.01
                s_piu(3) = (pos+dx) - 13.01;
            else
                s_piu(3) = 4;
                s_piu(4) = (pos+dx)-17.01;
            end    
        end
    end
    % Spessori inferiori
    if pos-dx>=0  && pos-dx< 0.01
        s_meno(1) = pos - (pos- dx); 
    end

% STRATO GHIANDOLARE INF
elseif pos> 0.01 && pos<= 13.01
    % Spessori inferiori
    if pos-dx>= 0.01
        s_meno(2) = pos - (pos - dx); 
    else
        s_meno(2) = pos - 0.01;
        s_meno(1) = 0.01 - (pos - dx);
    end
    % Spessori superiori
    if pos+dx<= 13.01
        s_piu(2) = (pos + dx) - pos;
    else
        s_piu(2) = 13.01 - pos;
        if pos+dx<= 17.01
            s_piu(3) = (pos + dx) -13.01;
        else
            s_piu(3) = 4;
            if pos+dx <=25.01
                s_piu(4) = (pos + dx) - 17.01;
            else
                s_piu(4)= 8;
                s_piu(5)= (pos+dx)-25.01;
            end
        end
    end

% TUMORE 
elseif pos>13.01 && pos<= 17.01
    % Spessori inferiori
    if pos-dx>= 13.01
        s_meno(3) = pos - (pos - dx); 
    else
        s_meno(3) = pos - 13.01;
        if pos-dx>=0.01
            s_meno(2)= 13.01 -(pos - dx);
        else
            s_meno(2)= 13;
            s_meno(1) = 0.01 - (pos - dx);
        end 
    end
    % Spessori superiori
    if pos+dx<= 17.01
        s_piu(3) = (pos + dx) - pos;
    else
        s_piu(3) = 17.01 - pos;
        if pos+dx<= 25.01
            s_piu(4) = (pos + dx) -17.01;
        else
            s_piu(4) = 8;
            if pos+dx <=28.51
                s_piu(5) = (pos + dx) - 25.01;
            else
                s_piu(5)= 3.5;
                if pos+dx<=32.51
                    s_piu(6)= (pos+dx)-28.51;
                else
                    s_piu(6)=4;
                    s_piu(7)= (pos+dx)-32.51;
                end
            end
        end
    end

% STRATO GHIANDOLARE SUP
elseif pos>17.01 && pos<= 25.01
    % Spessori inferiori
    if pos-dx>= 17.01
        s_meno(4) = pos - (pos - dx); 
    else
        s_meno(4) = pos - 17.01;
        if pos-dx>=13.01
            s_meno(3)= 17.01 -(pos - dx);
        else
            s_meno(3)= 4;
            if pos-dx>= 0.01
                s_meno(2)= 13.01-(pos-dx);
            else 
                s_meno(2)=13;
                s_meno(1) = 0.01 - (pos - dx);
            end
        end 
    end

    % Spessori superiori
    if pos+dx<= 25.01
        s_piu(4) = (pos + dx) -pos;
    else
        s_piu(4) = 25.01-pos;
        if pos+dx <=28.51
            s_piu(5) = (pos + dx) - 25.01;
        else
            s_piu(5)=3.5;
            if pos+dx<=32.51
                s_piu(6)=(pos+dx)-28.51;
            else
                s_piu(6)=4;
                s_piu(7)= (pos + dx) - 32.51;
            end
        end
            
    end


% SOTTOCUTE
elseif pos>25.01 && pos<= 28.51
    % Spessori inferiori
    if pos-dx>= 25.01
        s_meno(5) = pos - (pos - dx); 
    else
        s_meno(5) = pos - 25.01;
        if pos-dx>=17.01
            s_meno(4) = 25.01 -(pos - dx);
        else
            s_meno(4) = 8;
            s_meno(3) = 17.01 - (pos - dx);   
        end 
    end
    % Spessori superiore
    if pos+dx<=28.51
        s_piu(5)= (pos+dx)-pos;
    else
        s_piu(5)= 28.51-pos;
        if pos+dx<=32.51
            s_piu(6) = (pos+dx)-28.51;
        else
            s_piu(6) = 4;
            s_piu(7) = (pos+dx)-32.51;
        end
    end

% DERMA
elseif pos>28.51 && pos<= 32.51
    % Spessori inferiori
    if pos-dx>= 28.51
        s_meno(6) = pos - (pos - dx) ;
    else
        s_meno(6) = pos - 28.51;
        if pos-dx>=25.01
            s_meno(5) = 28.51 -(pos - dx) ;
        else
            s_meno(5) = 3.5;
            s_meno(4) = 25.01 - (pos - dx) ;
        end 
    end
    % Spessori superiori
    if pos+dx<=32.51
        s_piu(6) = (pos+dx)-pos ;
    else
        s_piu(6) = 32.51 -pos ;
        s_piu(7) = (pos+dx)-32.51 ;
    end

% EPIDERMIDE
elseif pos>32.51 && pos<= 36.015
    % Spessori inferiori
    if pos-dx>=32.51
        s_meno(7) = pos - (pos - dx) ;
    else
        s_meno(7) = pos- 32.51;
        if pos-dx >= 28.51
            s_meno(6) = 32.51 - (pos-dx) ;
        else
            s_meno(6) = 4;
            if pos-dx >= 25.01
                s_meno(5) = 28.51 - (pos-dx) ;
            else
                s_meno(5) = 3.5;
                s_meno(4) = 25.01 - (pos-dx) ;
            end
        end
    end
    % Spessore superiore
    if pos+dx>32.51 && pos+dx<=36.015
        s_piu(7)= (pos+dx)-pos ;
    end

end

Sigma= 1/(sum(vet_rho.*vet_c.*((s_piu+s_meno)*1e-3))) ;     %m

dx = 2*dx;

end 