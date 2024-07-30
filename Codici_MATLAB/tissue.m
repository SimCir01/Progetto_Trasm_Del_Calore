% La funzione "tissue" consente di plottare l'immamgine dei vari tessuti
% che compongono il dominio spaziale con la discretizzazione impostata
% dei nodi

function plot_tissues = tissue(nodi)

% Definizione del raggio del tumore
raggio_tumore = 2e-3; % 2 mm

% Creazione della figura
%figure;
hold on;
axis equal;

% Calcolo delle estremità del rettangolo
y_rec = [-15.1e-3, -15.01e-3, -15.0e-3, -15.1e-3]; % Coordinate x dei vertici
x_rec = [-72e-3 , -72e-3 , 72e-3 , 72e-3 ]; % Coordinate y dei vertici

% Disegno del rettangolo sotto gli ovali
fill(x_rec, y_rec, hex2rgb('#D95319')); % Colore arancione


% Definizione dei parametri per l'ovale rosa
centro_or_x = 0;
centro_or_y = -15.01e-3; % Centrato a -15.01 mm lungo y
asse_maggiore_or = 72e-3;
semiasse_minore_or = 36.01e-3;

% Creazione di 100 punti per l'ovale rosa
theta_or = linspace(0, pi, 100);
x_or = centro_or_x + asse_maggiore_or * cos(theta_or);
y_or = centro_or_y + semiasse_minore_or * sin(theta_or);

% Colorazione della parte superiore dell'ovale rosa
fill(x_or, y_or, [1, 0.753, 0.796]); % Rosa

% Definizione dei parametri per l'ovale beige
centro_ob_x = 0;
centro_ob_y = -15.01e-3; % Centrato a -15.01 mm lungo y
asse_maggiore_ob = 65e-3;
semiasse_minore_ob = 32.51e-3;

% Creazione di 100 punti per l'ovale beige
theta_ob = linspace(0, pi, 100);
x_ob = centro_ob_x + asse_maggiore_ob * cos(theta_ob);
y_ob = centro_ob_y + semiasse_minore_ob * sin(theta_ob);

% Colorazione della parte superiore dell'ovale beige
fill(x_ob, y_ob, [0.933, 0.914, 0.839]); % Beige




% Definizione dei parametri per l'ovale giallo
centro_oy_x = 0;
centro_oy_y = -15.01e-3; % Centrato a -15.01 mm lungo y
asse_maggiore_oy = 57e-3;
semiasse_minore_oy = 28.51e-3;

% Creazione di 100 punti per l'ovale giallo
theta_oy = linspace(0, pi, 100);
x_oy = centro_oy_x + asse_maggiore_oy * cos(theta_oy);
y_oy = centro_oy_y + semiasse_minore_oy * sin(theta_oy);

% Colorazione della parte superiore dell'ovale giallo
fill(x_oy, y_oy, [0.929, 0.694, 0.125]); % Giallo

% Definizione dei parametri per l'ovale grigio
centro_og_x = 0;
centro_og_y = -15.01e-3; % Centrato a -15.01 mm lungo y
asse_maggiore_og = 50e-3;
semiasse_minore_og = 25.01e-3;

% Creazione di 100 punti per l'ovale grigio
theta_og = linspace(0, pi, 100);
x_og = centro_og_x + asse_maggiore_og * cos(theta_og);
y_og = centro_og_y + semiasse_minore_og * sin(theta_og);

% Colorazione della parte superiore dell'ovale grigio
fill(x_og, y_og, [0.827, 0.827, 0.827]); % Grigio chiaro

% Disegno del cerchio del tumore (cerchio rosso)
rectangle('Position', [-raggio_tumore, -raggio_tumore, 2*raggio_tumore, 2*raggio_tumore], 'Curvature', [1, 1], 'FaceColor', [1, 0, 0], 'EdgeColor', 'none');


% Disegno dell'asse verticale x
line([0 0], [-15.01e-3 21e-3], 'Color', 'k', 'LineWidth', 1);

% Impostazione del titolo e degli assi
xlabel('x (m)');


% Impostazione dei limiti degli assi
axis([-15.01e-3*1.2, 15.01e-3*1.2, -15.01e-3*1.2, 21e-3*1.2]);

% Attiva la modalità cursore dati
datacursormode on;

% Imposta le opzioni del cursore dati
dcm_obj = datacursormode(gcf);
set(dcm_obj, 'UpdateFcn', @myupdatefcn);

% Imposta i limiti dell'asse y
limite_inf_y = -15.01e-3;
limite_sup_y = 21e-3;

% Genera i punti lungo l'asse y
y_points = linspace(limite_inf_y, limite_sup_y, nodi); % Genera 10 punti equidistanti

% Disegna i pallini rossi ai punti generati
hold on; % Mantieni il grafico esistente
plot_tissues = plot(zeros(size(y_points)), y_points, 'ro', 'MarkerSize', 8); % 'ro' indica rosso e cerchio



function rgb = hex2rgb(hex)
    % Converti il colore esadecimale in RGB
    hex = hex(2:end); % Rimuovi il carattere '#'
    r = hex2dec(hex(1:2)) / 255;
    g = hex2dec(hex(3:4)) / 255;
    b = hex2dec(hex(5:6)) / 255;
    rgb = [r, g, b];
end

function txt = myupdatefcn(~, event)
    pos = get(event,'Position');
    txt = {['X: ', num2str(pos(1))], ['Y: ', num2str(pos(2))]};
end



end