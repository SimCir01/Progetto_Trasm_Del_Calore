
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                                                                           %%%%%%%%%%%%%
%%%%%%%%%%%%%                 Analisi termica del transitorio 1D                        %%%%%%%%%%%%%
%%%%%%%%%%%%%           relativo ad una massa tumorale nel tessuto mammario             %%%%%%%%%%%%%
%%%%%%%%%%%%%            a seguito di una termoablazione a microonde                    %%%%%%%%%%%%%
%%%%%%%%%%%%%                                                                           %%%%%%%%%%%%%
%%%%%%%%%%%%%                  CORSO DI TRASMISSIONE DEL CALORE                         %%%%%%%%%%%%%
%%%%%%%%%%%%%              IN APPLICAZIONI BIOMEDICALI (A.A. 2023-24)                   %%%%%%%%%%%%%                   Numerical Anylysis cource project                       %%%%%%%%%%%%%
%%%%%%%%%%%%%                    Masone Benedetta (matr. 177470)                        %%%%%%%%%%%%%
%%%%%%%%%%%%%                    Cirnelli Simone (matr. 177084)                         %%%%%%%%%%%%%
%%%%%%%%%%%%%                                                                           %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all
clc

%% Definizione dello studio
% Scelta tra le seguenti opzioni in base a quello che si intende fare
fprintf('OBBIETTIVO DELLO STUDIO\rInserire uno dei seguenti valori numeri in base all''opzione che si intende selezionare:\r');
fprintf('[0] CAMPO DI TEMPERATURA \r');
fprintf('[1] STIMA DELLA MASSIMA POTENZA FISSATA UNA PERCENTUALE DI DANNO\r');
fprintf('[2] VERIFICA CONDIZIONI DI STABILITA'' PER IL METODO ESPLICITO \r');

obbiettivo=input('');
while obbiettivo ~= 0 && obbiettivo ~= 1 && obbiettivo ~= 2
    fprintf('\rERRORE\rValore inserito non valido\r')
    obbiettivo = input('Inserire nuovamente il valore prescelto:'); 
end 
clc



% CASO [0]-[1]-[2]
% Parametri per la risoluzione del campo di temperatura

% CASO 0
if obbiettivo == 0
    % RICHIESTA ALL'UTENTE DI INSERIRE DEI PARAMETRI NUMERICI
    fprintf('SETTAGGIO PARAMETRI SPAZIALI E TEMPORALI \r')
    % Scelta del metodo
    fprintf('Metodo da utilizzare per la risoluzione\rInserire uno dei seguenti valori numeri in base al metodo scelto\r');
    fprintf('[0] METODO ESPLICITO  \r');
    fprintf('[1] METODO IMPLICITO\r');
    fprintf('[2] METODO DI CRANK NICOLSON\r');
    fprintf('[3] TUTTI I METODI\r');
    csi=input('');

    while csi ~= 0 && csi ~= 1 && csi ~= 2 && csi ~= 3 
        fprintf('\rERRORE\rMetodo inserito non valido\r')
        csi = input('Inserire nuovamente il metodo prescelto:'); 
    end

    if csi == 2
         csi =0.5;
    end 

    % Numero di nodi
    clc
    nodi = input('Inserire il numero di nodi:');
    while nodi<12
        fprintf('\rERRORE\rNumero di nodi non sufficiente a discretizzare il problema\r')
        nodi = input('Inserire nuovamente il numero di nodi:'); 
    end 
    % Tempo totale
    tTot = input('Inserire la durata della simulazione (espressa in secondi):'); 
    % Passo temporale
    dt =input('Inserire il passo temporale (espresso in secondi):');

% CASO 1
elseif  obbiettivo == 1
    % Parametri per la ricerca del coefficiente correttivo nel CASO [1]
    % Numero di nodi
    nodi = input('Inserire il numero di nodi:');
    while nodi<12
        fprintf('\rERRORE\rNumero di nodi non sufficiente a discretizzare il problema\r')
        nodi = input('Inserire nuovamente il numero di nodi:'); 
    end 
    % Tempo totale
    tTot = input('Inserire la durata della simulazione (espressa in secondi):'); 
    % Passo temporale
    dt =input('Inserire il passo temporale (espresso in secondi):');

    fprintf('\rSpecificare la percentuale massima di danno:\r')
    perc_max_danno=input('');

% CASO 2
elseif obbiettivo == 2
    csi=0;
    fprintf('PARAMETRI PER L''IMPOSTAZIONE METODO ESPLICITO\r')

    fprintf('Inserire un vettore rappresentante il numero di nodi\r')
    nodi= ones(5,1);
    for i= 1:size(nodi,1)
         fprintf('Inserire %d-iesimo valore:', i)
         nodi(i) = input('');
    end

    % Tempo totale
    tTot = input('\rInserire la durata della simulazione (espressa in secondi):'); 
    % Passo temporale
    fprintf('Inserire un vettore rappresentante il passo temporale (espresso in secondi)\r')
    dt= ones(4,1);
    for i= 1:size(dt,1)
         fprintf('Inserire %d-iesimo valore:', i)
         dt(i) = input('');
    end
    dt=sort(dt);

end

clc
% Coefficiente correttivo SAR
corr=1;
probe=1;


%% RISOLUZIONE NUMERICA 
% Campo di temperatura
if obbiettivo == 0

    if csi ==0 || csi==1 || csi==0.5

        [A, b, T]= matrix(nodi, tTot, dt, csi, corr, probe, obbiettivo);

        % Vettore percentuale danno
        perc_danno_capillare=zeros(1,size(T, 2));
        for j= 1: size(T,2)

            perc_danno_capillare(1,j)= Arrh(T,j, dt, nodi,1);

        end

        % Paremetri per il plot
        % DEFINIZIONI PASSO SPAZIALE E TEMPORALE
        % Settaggio parametri per la soluzione numerica nello spazio
        L = 36.01e-3;                   % Lunghezza del sistema [m] 
        x = linspace(0,L,nodi);         % [m]
    
        % Numero step (arrotondato per eccesso)
        steps = ceil(tTot/dt);
        % Vettore dei tempi
        t= linspace(0,tTot, steps+1);

return
%% ANDAMENTI TEMPORALI E SPAZIALI DELLA TEMPERATURA
figure
%% 1. Progressione nel tempo della distribuzione spaziale di temperatura
subplot(2, 2, 1);
% 
% % Initialize video
% myVideo = VideoWriter('video_T_x_varT_corr'); %open video file
% myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
% open(myVideo)

hold on;
for j = 1:steps
    plot(x-0.0150, T(:, j), 'b', 'LineWidth', 1.5);
    title('Progressione nel tempo della distribuzione spaziale di temperatura');
    xlabel('x [m]', 'FontSize', 12);
    ylabel('T [°C]', 'FontSize', 12);
    grid on;
    note = sprintf('Istante di tempo %d', j);  
    if exist('note_handle', 'var') && isvalid(note_handle)
        delete(note_handle); % Rimuove nota ad ogni istante di tempo
    end
    note_handle = text('Units', 'normalized', 'Position', [0.5, 0.95], ...
        'String', note, 'FontSize', 10, 'HorizontalAlignment', 'center', ...
        'Color', 'black');
    drawnow; % Aggiorna la figura
    pause(0.2); % Ritardo di 0.1 secondi 

%     frame = getframe(gcf); %get frame
%     writeVideo(myVideo, frame);
end
% close(myVideo)

hold off;


%% 2. Progressione nello spazio della distribuzione temporale di temperatura
subplot(2, 2, 2);
%figure
% Initialize video
% myVideo = VideoWriter('video_T_t_varx.avi'); %open video file
% myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
% open(myVideo)
% 
tessuti = ["Capillare", "Strato Ghiandolare Inf", "Tumore", "Strato Ghiandolare Sup", "Sottocute", "Derma", "Epidermide"];
colori = [0.8500 0.3250 0.0980; 0.827 0.827 0.827; 1 0 0; 0.827 0.827 0.827;0.9290 0.6940 0.1250;  0.933  0.914  0.839; 1 0.753 0.796 ];

hold on;
grid on;


for j = 1:steps+1
for i = 1:nodi    
    
    [~, nodo,~ ] = cond(i, (x(2)-x(1)));

    if nodo == 1 
       colore = [0.8500 0.3250 0.0980]; % capillare (arancio scuro)
     
    elseif nodo ==2 || nodo == 4
        colore = [0.827, 0.827, 0.827]; % strato ghiand. (grigio)
      
    elseif nodo ==3
        colore = [1 0 0]; % tumore (rosso)
        
    elseif nodo ==5
        colore = [0.9290 0.6940 0.1250]; % sottocute (giallo)
 
    elseif nodo ==6
        colore = [0.933  0.914  0.839]; % derma (beige chiaro)

    elseif nodo ==7
        colore = [1 0.753 0.796]; % epidermide (rosa)

    end
    

    if j~=steps+1
        plot(t(1:j), T(i,1:j), 'Color', colore, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        xlim([0,10])
        title('Progressione nello spazio della distribuzione temporale di temperatura');
        xlabel('tempo [s]', 'FontSize', 12);
        ylabel('T [°C]', 'FontSize', 12);
    else
         p(i,j) = plot(t(1:j), T(i,1:j), 'Color', colore, 'LineWidth', 1.5, 'DisplayName', [convertStringsToChars(tessuti(nodo))]);
         xlim([0,10]);
         title('Progressione nello spazio della distribuzione temporale di temperatura');
        xlabel('tempo [s]', 'FontSize', 12);
        ylabel('T [°C]', 'FontSize', 12);
         %legend
    end
end
    drawnow; % Aggiorna la figura
    pause(0.1); % Ritardo di 0.1 secondi 

%     frame = getframe(gcf); %get frame
%     writeVideo(myVideo, frame);

end

% close(myVideo)
hold off;



%% 3. Plot 3D della temperatura nello spazio e nel tempo
subplot(2, 2, 3);
%figure
[X, Y] = meshgrid(t, x-0.0150);
surf(X, Y, T);
title('Plot 3D della temperatura nello spazio e nel tempo');
xlabel('tempo [s]', 'FontSize', 12);
ylabel('spazio [m]', 'FontSize', 12);
zlabel('T [°C]', 'FontSize', 12);
grid on;

%% 4. Rappresentazione grafica strati di tessuti
subplot(2, 2, 4);
title('Rappresentazione grafica strati di tessuti');
plot_tissues = tissue(nodi); % Supponiamo che tissue sia una funzione definita che restituisce il plot dei tessuti



% %% Variabilità di dx per T 
% Nnodi = 12:10:1000;
% 
% figure
% count = 1
% %Tumore
% for n=1:length(Nnodi)
%     [~, ~, T]= matrix(Nnodi(n), tTot, dt, csi, corr, probe, obbiettivo);
%     
%     x = linspace(0,L,Nnodi(n));
%     dx = x(2)-x(1);
%     DX(n)=dx;
%     
%     % Posizione del tumore
% %     pos=0.0150;
%     pos = 15.01e-3;
%     i = ceil((pos/dx)+1);
%     T_tum(n)=T(i,end);
%     
% end
% count=2
% 
% %Epidermide
% for n=1:length(Nnodi)
%     [~, ~, T]= matrix(Nnodi(n), tTot, dt, csi, corr, probe, obbiettivo);
%     
%     x = linspace(0,L,Nnodi(n));
%     dx = x(2)-x(1);
%     DX(n)=dx;
%     
%     pos = 35.0e-3;
%     i = ceil((pos/dx)+1);
%     T_ep(n)=T(i,end);
%     
% end
% 
% count = 3
% % Strato Ghiandolare inferiore
% for n=1:length(Nnodi)
%     [~, ~, T]= matrix(Nnodi(n), tTot, dt, csi, corr, probe, obbiettivo);
%     
%     x = linspace(0,L,Nnodi(n));
%     dx = x(2)-x(1);
%     DX(n)=dx;
%     
%     pos = 12.0e-3;
%     i = ceil((pos/dx)+1);
%     T_gh(n)=T(i,end);
%     
% end
% 
%     hold on;
%     grid on
%     plot(Nnodi, T_tum, 'r', 'LineWidth', 2, 'DisplayName', 'Tumore');
%     plot(Nnodi, T_ep, 'Color', [1 0.753 0.796], 'LineWidth', 2, 'DisplayName', 'Epidermide');
%     plot(Nnodi, T_gh, 'Color', [0.827 0.827 0.827],  'LineWidth', 2, 'DisplayName', 'Strato Ghiandolare Inf');
%     % Nodi = 12, tTot = 10, dt = 1
% %   title('Variabilità di dx per T nel tumore all ''istante finale');
%     legend
%     title('Temperatura in funzione del passo spaziale all''istante finale');
%     xlabel('Nodi [a.u.]', 'FontSize', 12);
%     ylabel('T [°C]', 'FontSize', 12);
%     set(gca, 'FontSize', 12, 'FontWeight', 'bold');
%     hold off;
% 
% %% Variabilità di dt per T nel tumore
% 
% dt=0.0001:0.01:10.0001;
% time_tot = 10.0001;
% Nnodi=12;
% 
% 
% figure
% 
% count = 1
% % Tumore
% for it=1:length(dt)
%     [~, ~, T]= matrix(Nnodi, time_tot, dt(it), csi, corr, probe, obbiettivo);
% 
%     step = ceil(time_tot/dt(it));
%     time= linspace(0,time_tot, step+1);
%     DT=dt(it);
% 
%     x = linspace(0,L,Nnodi);
%     dx = x(2)-x(1);
%    
% 
%     % Posizione del tumore
%     pos=0.0150;
%     i = ceil((pos/dx)+1);
%     T_tum(it)=T(i,end);
%    
% end
% 
% count = 2
% % Epidermide
% for it=1:length(dt)
%     [~, ~, T]= matrix(Nnodi, time_tot, dt(it), csi, corr, probe, obbiettivo);
% 
%     step = ceil(time_tot/dt(it));
%     time= linspace(0,time_tot, step+1);
%     DT=dt(it);
% 
%     x = linspace(0,L,Nnodi);
%     dx = x(2)-x(1);
%    
% 
%     % Posizione 
%     pos=36e-3;
%     i = ceil((pos/dx)+1);
%     T_ep(it)=T(i,end);
%    
% end
% 
% count = 3
% % Strato Ghiandolare inf
% for it=1:length(dt)
%     [~, ~, T]= matrix(Nnodi, time_tot, dt(it), csi, corr, probe, obbiettivo);
% 
%     step = ceil(time_tot/dt(it));
%     time= linspace(0,time_tot, step+1);
%     DT=dt(it);
% 
%     x = linspace(0,L,Nnodi);
%     dx = x(2)-x(1);
%    
% 
%     % Posizione del tumore
%     pos=12e-3;
%     i = ceil((pos/dx)+1);
%     T_gh(it)=T(i,end);
%    
% end
% 
%     hold on;
%     grid on;
%     plot(dt, T_tum, 'r', 'LineWidth', 2, 'DisplayName', 'Tumore');
%     plot(dt, T_ep, 'Color', [1 0.753 0.796], 'LineWidth', 2, 'DisplayName', 'Epidermide');
%     plot(dt, T_gh, 'Color', [0.827 0.827 0.827],  'LineWidth', 2, 'DisplayName', 'Strato Ghiandolare Inf');
%     % Nodi = 12, tTot = 10, dt = 1
%     title('Temperatura in funzione del passo temporale all ''istante finale');
%     xlabel('dt [s]', 'FontSize', 12);
%     ylabel('T [°C]', 'FontSize', 12);
%     legend
%     set(gca, 'FontSize', 12, 'FontWeight', 'bold');
%     
% 
%     hold off;

% Plot dell'andamento nel tempo delle percentuali di danno per i vari tessuti
    
% Valutazione del danno effettuata a fattore di correzione unitario

% num_tessuti=7;
% tessuti = ["Capillare", "Strato Ghiandolare Inf", "Tumore", "Strato Ghiandolare Sup", "Sottocute", "Derma", "Epidermide"];
% colori = [0.8500 0.3250 0.0980; 0.827 0.827 0.827; 1 0 0; 0.827 0.827 0.827;0.9290 0.6940 0.1250;  0.933  0.914  0.839; 1 0.753 0.796 ];
% 
% [~, ~, T]= matrix(nodi, tTot, dt, csi, corr, probe, obbiettivo);
% 
% 
% % Vettore percentuale danno
% perc_danno_tessuti=zeros(num_tessuti,size(T, 2));
% 
% figure
% hold on;
% grid on;
% 
% % Calcolo percentuale di danno all'interno del capillare
% for tess = 1:size(perc_danno_tessuti,1)
% 
%     for j= 1:size(perc_danno_tessuti,2)
%         perc_danno_tessuti(tess,j)= Arrh(T,j, dt, nodi,tess);       
%     end
% 
%     plot(t, perc_danno_tessuti(tess,:), 'Color', colori(tess,:), 'LineWidth', 2, 'DisplayName', [convertStringsToChars(tessuti(tess))]);
%     
% end
% 
% legend('Location', 'best');
% 
% title('Percentuale di danno nel tempo per i vari tessuti a fattore correttivo SAR unitario');
% xlabel('Tempo [s]', 'FontSize', 15);
% ylabel('Percentuale di danno [%]', 'FontSize', 15);
% set(gca, 'FontSize', 12, 'FontWeight', 'bold');
% hold off;


% %% Plot delle temperatura nel tessuto ghiandolare superiore e inferiore 
% hold on;
% grid on;
% 
% k=1;
% h=1;
% for i = 1:nodi
%     [nodo_piu, nodo, nodo_meno] = cond(i, dx);
% 
%     %Tessuto ghiandolare inferiore
%     if nodo ==2 
% 
%          T_ghiand_inf(k,:)= T(i,:);
%          k=k+1;
%          
%     % Tessuto ghiandolare superiore
%     elseif nodo ==4
% 
%         T_ghiand_sup(h,:)= T(i,:);
%         h=h+1;
% 
%     end
% end
% 
% for j=1:steps+1
%     T_ghiand_inf_media(j)= mean(T_ghiand_inf(:,j));
%     T_ghiand_sup_media(j)= mean(T_ghiand_sup(:,j));
% end 
% plot(t, T_ghiand_inf_media, 'b', 'LineWidth', 2, 'DisplayName', ['Strato ghiandolare inf.'])
% plot(t, T_ghiand_sup_media, 'g',  'LineWidth', 2, 'DisplayName', ['Strato ghiandolare sup.'])
% title('Progressione nel tempo della temperatura nel tessuto ghiandolare');
% xlabel('tempo [s]', 'FontSize', 12);
% ylabel('Temperatura [°C]', 'FontSize', 12);
% legend('Location', 'best');
% set(gca, 'FontSize', 12, 'FontWeight', 'bold');
% hold off;
% 

% % %% Plot delle temperature nel tumore e lungo i suoi bordi
% figure
% hold on;
% grid on;
% 
% T_centro_tumore=[1, steps+1];
% T_sup_tumore=[1, steps+1];
% T_inf_tumore=[1, steps+1];
% 
% for j = 1:steps+1
% 
%     p = polyfit(x, T(:,j), 6);
% 
%     T_centro_tumore(1,j)=polyval(p,0.0150);
%     T_sup_tumore(1,j)=polyval(p, 0.01701);
%     T_inf_tumore(1,j)=polyval(p, 0.01301);
%    
% end
% 
% plot(t, T_centro_tumore, 'r', 'LineWidth', 2, 'DisplayName', ['Centro del tumore'])
% plot(t, T_sup_tumore, 'g',  'LineWidth', 2, 'DisplayName', ['Bordo superiore del tumore'])
% plot(t, T_inf_tumore, 'y', 'LineWidth', 2, 'DisplayName', ['Bordo inferiore del tumore'])
% title('Progressione nel tempo della temperatura nel tessuto tumorale');
% xlabel('t [s]', 'FontSize', 12);
% ylabel('T [°C]', 'FontSize', 12);
% legend('Location', 'best');
% set(gca, 'FontSize', 12, 'FontWeight', 'bold');
% hold off;



    else 
         csi=[0, 1, 0.5];


        % Paremetri per il plot
        % DEFINIZIONI PASSO SPAZIALE E TEMPORALE
        % Settaggio parametri per la soluzione numerica nello spazio
        L = 36.01e-3;                   % Lunghezza del sistema [m] 
        x = linspace(0,L,nodi)-0.0150;         % [m]
    
        % Numero step (arrotondato per eccesso)
        steps = ceil(tTot/dt);
        % Vettore dei tempi
        t= linspace(0,tTot, steps+1);

        soluzione= ones(nodi, steps+1, 3);


        % Per vedere anche il metodo esplicito in condizioni di instabilità
        % occorre commentare la parte di codice in "matrix"
        % che regola il reinserimento del passo temporale
        % fino a che non se ne sceglie uno che porti alla stabilità
        
        % Confronto tra i metodi a istante temporale fissato
        figure
        hold on;
        grid on;

        for i=1:size(csi,2)

            [~, ~, T]= matrix(nodi, tTot, dt, csi(i), corr, probe, obbiettivo);
            csi(i)
            % Plot
            plot(x,T(:,end), 'LineWidth',2, 'DisplayName', ['csi = ' num2str(csi(i))]);
            % Nodi = 20,   Ttot = 1,  dt = 2e-2
        end
        legend('Location', 'best');

        title('Confronto tra i metodi (esplicito, implicito, Crank-Nicolson) fissato un valore temporale');
        xlabel('Spazio [m]', 'FontSize', 15);
        ylabel('Temperatura [°C]', 'FontSize', 15);
        hold off
        
        % Confronto tra i metodi a posizione spaziale fissata (tumore)
        figure
        hold on;
        grid on;

        for i=1:size(csi,2)
            csi(i)
            [A, b, T]= matrix(nodi, tTot, dt, csi(i), corr, probe, obbiettivo);
            % Plot
            plot(t(1:end),T(8,1:end), 'LineWidth',2, 'DisplayName', ['csi = ' num2str(csi(i))]);
            % Nodi = 20,   Ttot = 10,  dt = 1e-1
        end
        legend('Location', 'best');

        title('Confronto tra i metodi (esplicito, implicito, Crank-Nicolson) fissato la posizione spaziale (tumore)');
        xlabel('Tempo [s]', 'FontSize', 15);
        ylabel('Temperatura [°C]', 'FontSize', 15);
        hold off



    end

%% RICERCA DEL COEFFICIENTE CORRETTIVO
elseif obbiettivo == 1
    csi=0.5;
    max_perc_danno_capillare=[0 ];
    i=2;
    
    while max_perc_danno_capillare(i-1)< perc_max_danno+80


        
        [A, b, T]= matrix(nodi, tTot, dt, csi, corr, probe, obbiettivo);
        % Vettore percentuale danno
        perc_danno_capillare=zeros(1,size(T, 2));

   
        % Calcolo percentuale di danno all'interno del capillare
        for j= 1: size(T,2)
            
            perc_danno_capillare(1,j)= Arrh(T,j, dt, nodi,1);

        end

       % Incremento coefficiente corretivo per il calcolo della massima
       % potenza
       max_perc_danno_capillare(i)= max(perc_danno_capillare);
       corr = corr + 1
       i=i+1;
    
    end

    figure
    grid on
    hold on

    x_corr=linspace(1, corr, corr);
    plot(x_corr, max_perc_danno_capillare, 'r', 'LineWidth', 2 )
    yline(20, '--', 'LineWidth', 2 )
    
    title('Varianzione della percentuale di danno nel capillare in funzione del coefficiente correttivo del SAR')
    xlabel('Coefficiente correttivo [a.u.]')
    ylabel('Percentuale danno [ % ]')
    ylim([0 100])
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');
    hold off


%% STUDIO DELLA STABILITÀ
elseif obbiettivo == 2

    % Matrice contenente l'andamento del dt_limite al variare del passo
    % spaziale e temporale
    dt_limite= zeros(size(nodi,1), size(dt,1));

    Tempo_simulazioni = zeros(size(nodi,1),size(dt,1));

    for stab_dx=1:size(nodi,1)
        for stab_dt=1:size(dt,1)

            % Inizio conteggio
            tic;
            
            [A, b, T]= matrix(nodi(stab_dx), tTot, dt(stab_dt), csi, corr, probe, obbiettivo);

             Tempo_simulazioni(stab_dx, stab_dt) = toc;

            dt_critico= zeros(nodi(stab_dx),1);
            for i=1:nodi(stab_dx)
                dt_critico(i,1)= -dt(stab_dt)/A(i,i,1);
            end
            dt_limite(stab_dx, stab_dt) = min(dt_critico);
            
            fprintf('nodi= %f, passo temporale = %f --->  passo critico %f\r\r', nodi(stab_dx), dt(stab_dt), dt_limite(stab_dx, stab_dt))

        end
    end

% figure
% hold on;
% grid on
% 
% 
% dt_2 = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]';
% dt_limite = repmat(dt_limite,1,5);
% 
% % Definisci una palette di colori
% colori = {'b', 'g', 'r', 'c', 'm', 'y', 'k'};
% 
% 
% for stab_dx = 1:size(nodi, 1)
%     % Seleziona un colore dalla palette
%     colore = colori{mod(stab_dx, length(colori)) + 1};
% 
%     % Plot delle linee con etichetta
%     plot(dt_2, dt_limite(stab_dx,:), colore, 'LineWidth', 2, 'DisplayName', ['nodi = ' num2str(nodi(stab_dx))]);
% 
%     % Aggiunta di cerchi sui punti senza etichetta
%     plot(dt_2, dt_limite(stab_dx,:), 'o', 'Color', colore, 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
% 
% end
% 
%     plot(dt_limite(:,1), dt_limite(:,1), '--',  'LineWidth', 2, 'HandleVisibility', 'off');
% 
% % % Aggiunta della legenda solo per il primo plot
% legend('Location', 'best');
% 
% title('Zona di stabilità per il metodo esplicito');
% xlabel('Passo temporale [s]', 'FontSize', 12);
% ylabel('Passo temporale critico [°C]', 'FontSize', 12);
% hold off;

% Grafico Tempo simulazione 
figure
hold on;
grid on
% Definisci una palette di colori
colori = {'b', 'g', 'r', 'c', 'm', 'y', 'k'};


for stab_dx = 1:size(dt, 1)
    % Seleziona un colore dalla palette
    colore = colori{mod(stab_dx, length(colori)) + 1};

    % Plot delle linee con etichetta
    plot(nodi, Tempo_simulazioni(:,stab_dx), colore, 'LineWidth', 2, 'DisplayName', ['dt = ' num2str(dt(stab_dx))]);

    % Aggiunta di cerchi sui punti senza etichetta
    plot(nodi, Tempo_simulazioni(:,stab_dx, :), 'o', 'Color', colore, 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');

end


% % Aggiunta della legenda solo per il primo plot
legend('Location', 'best');

title('Stima del tempo impiegato in fase si simulazione variando la discretizzazione del sistema');
xlabel('Nodi [a.u.]', 'FontSize', 12);
ylabel('Tempo simulazione [s]', 'FontSize', 12);
hold off;





end















