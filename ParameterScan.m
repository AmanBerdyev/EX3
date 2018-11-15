% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.
%

%% Parametres %%
%%%%%%%%%%%%%%%%

sizeX = 10;
sizeY = 10;


%{
    colors = [  45, 149, 191; %Blue
            240, 196, 25; %Yellow
            241, 90, 90;  %Red
            149, 91, 165; %Purple
            78, 186, 111; %Green
             ]/255;
             
%}

colors = [166,206,227; %light blue
31,120,180; % dark blue
178,223,138; % light green
51,160,44; % dark green
251,154,153; % light red
227,26,28; %dark red
253,191,111; %light orange
255,127,0; % dark orange
202,178,214; %light purple
106,61,154; %dark purple
]/255;

set(groot, 'DefaultFigureResize',               'on'               );
set(groot, 'DefaultFigurePaperUnits',           'centimeters'       );
set(groot, 'DefaultFigureUnits',                'centimeters'       );
set(groot, 'DefaultFigurePaperSize',            [sizeX, sizeY]      );
set(groot, 'DefaultFigureInvertHardcopy',       'on'                );
set(groot, 'DefaultFigurePaperPosition',        [0, 0, sizeX, sizeY]);
set(groot, 'DefaultFigurePosition',             [10,10,sizeX,sizeY] );

set(groot, 'DefaultAxesColorOrder',             colors          );
set(groot, 'DefaultLineLineWidth',              0.25            );
set(groot, 'DefaultLineMarker',                     'o'         );
set(groot, 'DefaultLineMarkerEdgeColor',            colors(1,:) );
set(groot, 'DefaultLineMarkerFaceColor',            colors(1,:) );
set(groot, 'DefaultLineMarkerSize',                 3           );

set(groot, 'DefaultTextInterpreter',            'LaTeX' );
set(groot, 'DefaultAxesTickLabelInterpreter',   'LaTeX' );
set(groot, 'DefaultAxesFontName',               'LaTeX' );
set(groot, 'DefaultAxesFontSize',               11      );
set(groot, 'DefaultAxesBox',                    'off'   );
set(groot, 'DefaultAxesXGrid',                  'on'    );
set(groot, 'DefaultAxesYGrid',                  'on'    );
set(groot, 'DefaultAxesGridLineStyle',          ':'     );
set(groot, 'DefaultAxesLayer',                  'top'   );
set(groot, 'DefaultLegendInterpreter',          'LaTeX' );


repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base



%Exercice:

Ex = '3a';

switch Ex

case {'3b','3a'}

    w0        = sqrt(g/L);
    % Parametres physiques :
    tFin      = 250;
    d         = 0.03;
    Omega     = w0;
    kappa     = 0;
    m         = 0.1;
    g         = 9.81;
    L         = 0.1;
    theta0    = 0;
    thetadot0 = 0.01;

    % Parametres numeriques :
    Dt       = 0.001;

    config = sprintf(  ['%s=%.15g %s=%.15g %s=%.15g %s=%.15g %s=%.15g'...
                    ' %s=%.15g %s=%.15g %s=%.15g %s=%.15g %s=%.15g'] , ...
                    'tFin'      ,tFin       ,...
                    'd'         ,d          ,...
                    'Omega'     ,Omega      ,...
                    'kappa'     ,kappa      ,...
                    'm'         ,m          ,...
                    'g'         ,g          ,...
                    'L'         ,L          ,...
                    'theta0'    ,theta0     ,...
                    'thetadot0' ,thetadot0  ,...
                    'dt'        ,Dt    );


    %% Simulations %%
    %%%%%%%%%%%%%%%%%

    switch Ex
    case '3a'

        name = [Ex,'.out'];
        if (exist(name,'file') ~= 2)
        cmd = sprintf('%s%s %s %s output=%s', repertoire, executable, input,config ,name);
        system(strcat("wsl ",cmd));
        end

        data = load(name); % Chargement du fichier de sortie pour Omega = w0

        t        = data(:,1);
        theta    = data(:,2);
        thetaDot = data(:,3);
        emec     = data(:,4);
        pnc      = data(:,5);

        demec =diff(emec)/Dt;
        figure
        c = linspace(1,100,length(demec));
        scatter(demec,pnc(1:end-1),[],c)

        %plot(demec,pnc(1:end-1),'o');
        hold on;
        id = plot(demec,demec,'-');
        id.LineWidth = 1.5;

        axis tight;


        figure
        plot(t(2:end),1*demec,'.');
        hold on;
        plot(t,pnc,'.');

        figure
        plot(t(2:end),pnc(2:end)-demec);

        figure
        plot(t,emec);

    case '3b'

        nsimul = 21; % Nombre de simulations a faire

        voisinage = linspace(0.98,1.02,nsimul);
        Omega     = w0.*voisinage;

        paramstr = 'Omega'; % Nom du parametre a scanner (changer ici 'dt' ou 'Omega' ou autre)
        param    = Omega;

        output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie

        for i = 1:nsimul
            output{i} = [Ex,paramstr, '=', num2str(param(i)), '.out'];
            % Execution du programme en lui envoyant la valeur a scanner en argument
            cmd = sprintf('%s%s %s %s %s=%.15g output=%s', repertoire, executable, input,config ,paramstr, param(i), output{i});
            disp(cmd)
            system(strcat("wsl ",cmd));
        end


        f=figure;
        for i = 1:nsimul % Parcours des resultats de toutes les simulations
            data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation

            t        = data(:,1);
            theta    = data(:,2);
            thetaDot = data(:,3);
            emec     = data(:,4);
            pnc      = data(:,5);

            figure(f)
            plot3(ones(size(t))*Omega(i),t,emec-emec(1));
            hold on;

            Emax(i) = max(emec-emec(1)); 
        end

        %% Figures %%
        %%%%%%%%%%%%%

        figure
        plot(Omega, Emax, 'k-+')
        xlabel('\Omega [rad/s]')
        ylabel('max(E_{mec}(t)) [J]')
        grid on

        for i = 1:nsimul
        cmd = sprintf('rm %s', output{i});
        disp(cmd)
        system(strcat("wsl ",cmd));
        end
    end


otherwise




nsimul = 20; % Nombre de simulations a faire
nsteps = logspace(3,5,nsimul);

dt     = tFin./nsteps; % TODO: Choisir des valeurs de dt pour faire une etude de convergence
Omega  = linspace(w0*0.98,w0*1.02,nsimul); % TODO: Choisir des valeurs de Omega pour trouver la resonance
theta0 = linspace(0.01,pi-0.01,nsimul);

traverse = true;
name = 'C';

paramstr = 'Omega'; % Nom du parametre a scanner (changer ici 'dt' ou 'Omega' ou autre)
param = Omega; % Valeurs du parametre a scanner (changer ici dt ou Omega ou autre)


Step   = 500;
    Start  = 1e+3;
    nsteps = Start:step:(Start+nsimul*Step);

    dt     = tFin./nsteps;


config = sprintf(  ['%s=%.15g %s=%.15g %s=%.15g %s=%.15g %s=%.15g' ...
                   ' %s=%.15g %s=%.15g %s=%.15g %s=%.15g %s=%.15g'] , ...
                    'tFin'      ,tFin       ,...
                    'd'         ,d          ,...
                    'Omega'     ,w0         ,...
                    'kappa'     ,kappa      ,...
                    'm'         ,m          ,...
                    'g'         ,g          ,...
                    'L'         ,L          ,...
                    'theta0'    ,theta0     ,...
                    'thetadot0' ,thetadot0  ,...
                    'dt'        ,Dt          ...
);

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie

if traverse
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i)), '.out'];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s %s=%.15g output=%s', repertoire, executable, input,config ,paramstr, param(i), output{i});
    disp(cmd)
    system(strcat("wsl ",cmd));
end
else
    cmd = sprintf('%s%s %s %s output=%s', repertoire, executable, input,config ,name);
    system(strcat("wsl ",cmd));
end
%% Analyse %%
%%%%%%%%%%%%%


if traverse
if strcmp(paramstr, 'dt')
    error = zeros(1,nsimul);
elseif strcmp(paramstr, 'Omega')
    Emax = zeros(1,nsimul);
end
f=figure;
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation

    t        = data(:,1);
    theta    = data(:,2);
    thetaDot = data(:,3);
    emec     = data(:,4);
    pnc      = data(:,5);

    %this is a tess!!!!

    theta_th = thetadot0./w0 *sin(w0*t) + theta0 * cos(w0*t);
    
        figure(f)

    plot3(ones(size(t))*Omega(i),t,emec-emec(1));
    hold on;
    
    if strcmp(paramstr, 'dt')
        error(i) = theta(end);%max(abs(theta_th-theta));%sqrt(sum((theta_th-theta).^2)/(length(theta)-1)); % TODO: Calculer l'erreur a partir de l'output
    elseif strcmp(paramstr, 'Omega')
        Emax(i) = max(emec-emec(1)); % TODO: Calculer le maximum de l'energie
    end
end

%% Figures %%
%%%%%%%%%%%%%

if strcmp(paramstr, 'dt')
    figure
    plot(dt, error, 'k+')
    xlabel('\Delta t')
    ylabel('Erreur sur \theta(t_{fin}) [rad]')
    grid on
elseif strcmp(paramstr, 'Omega')
    figure
    plot(Omega, Emax, 'k-+')
    xlabel('\Omega [rad/s]')
    ylabel('max(E_{mec}(t)) [J]')
    grid on
end

else

    data = load(name); % Chargement du fichier de sortie de la i-ieme simulation

    t        = data(:,1);
    theta    = data(:,2);
    thetaDot = data(:,3);
    emec     = data(:,4);
    pnc      = data(:,5);

    demec =diff(emec)/Dt;
    figure

    plot(demec,pnc(1:end-1),'--');
    axis equal;

    figure
    plot(t(2:end),1*demec,'.');
    hold on;
    plot(t,pnc,'.');

    figure
    plot(t(2:end),pnc(2:end)-demec);

    figure
    plot(t,emec);

end



end


