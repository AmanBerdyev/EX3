% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'inputName
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
inputName = 'configuration.in'; % Nom du fichier d'entree de base



%Exercice:
promptEx = 'Exercise : ';
Ex = input(promptEx,'s');
deleteAfter = input('Do you want to delete the generated files? [Y/N] ','s');

%%% --- Global parameters   --- %%%
g         = 9.81;
L         = 0.1;
m         = 0.1;
w0        = sqrt(g/L);


switch Ex

case {'3b','3a'}


    % Parametres physiques :
    tFin      = 250;
    d         = 0.03;
    kappa     = 0;
    m         = 0.1;
    theta0    = 0;
    thetadot0 = 0.01;

    Omega     = w0;

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

            % Name of output file to generate
            name = [Ex,'.out'];


            %%%%%  --- SIMULATION ---   %%%%%
            if (exist(name,'file') ~= 2) %test if the file exists
                %if not:
                cmd = sprintf('%s%s %s %s output=%s', repertoire, executable, inputName,config ,name);
                system(strcat("wsl ",cmd)); % Wsl to compile using gcc on the wsl (windows subsystem for linux)
            end

            data = load(name); % Load generated file

            %%%%%   --- Load data   --- %%%%%

            t        = data(:,1);
            theta    = data(:,2);
            thetaDot = data(:,3);
            emec     = data(:,4);
            pnc      = data(:,5);

            demec =diff(emec)/Dt; %Derivative of mechanical energy

            %%%%%   --- Plot        --- %%%%%

            %%  - Mechanical energy theorem verification    - %%
            mecEnThmFig = figure;
            mecEnThmax  = axes(mecEnThmFig);
            mecEnThmScatter = scatter(mecEnThmax,demec,pnc(1:end-1),2       ,...
                'Marker'         , 'o'          ,...
                'MarkerFaceColor', colors(5,:)  ,...
                'MarkerEdgeColor', colors(5,:)  ,...
                'MarkerFaceAlpha', .01          ,...
                'MarkerEdgeAlpha', .01);
            hold on;
            %id = plot(mecEnThmax,demec,demec,'-');
            %id.LineWidth = 1.5;
            axis tight;

            %%  - Other verification method of the theorem  - %%
            dmecEnPncFig = figure;
            dmecEnPncax  = axes(dmecEnPncFig);

            dmecEnScatter = scatter(dmecEnPncax,t(1:end-1),1*demec,2 ,...
                'LineWidth'      , 3            ,...
                'Marker'         , 'o'          ,...
                'MarkerFaceColor', colors(2,:)  ,...
                'MarkerEdgeColor', colors(2,:)  ,...
                'MarkerFaceAlpha', .1          ,...
                'MarkerEdgeAlpha', .1);

            hold on;

            PncScatter = scatter(dmecEnPncax,t,pnc,2,...
                'LineWidth'      , 3            ,...
                'Marker'         , 'o'          ,...
                'MarkerFaceColor', colors(6,:)  ,...
                'MarkerEdgeColor', colors(6,:)  ,...
                'MarkerFaceAlpha', .1          ,...
                'MarkerEdgeAlpha', .1);

            %%%%%   --- Delete after simulation?    --- %%%%%

            if deleteAfter == 'Y'
            cmd = sprintf('rm %s', name);
            disp(cmd)
            system(strcat("wsl ",cmd));
            end

        case '3b'

            %%%%%  --- SIMULATION ---   %%%%%
            nsimul = 21; % Nombre de simulations a faire

            voisinage = linspace(0.98,1.02,nsimul);
            Omega     = w0.*voisinage;

            paramstr = 'Omega'; % Nom du parametre a scanner (changer ici 'dt' ou 'Omega' ou autre)
            param    = Omega;

            output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie

            for i = 1:nsimul
                output{i} = [Ex,paramstr, '=', num2str(param(i)), '.out'];
                if (exist(output{i},'file') ~= 2) %test if the file exists
                %if not:
                % Execution du programme en lui envoyant la valeur a scanner en argument
                cmd = sprintf('%s%s %s %s %s=%.15g output=%s', repertoire, executable, inputName,config ,paramstr, param(i), output{i});
                disp(cmd)
                system(strcat("wsl ",cmd));
                end
            end

            %%%%%   --- Load data   --- %%%%%
            % and dynamically plot the energy

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

            %%%%%   --- Figures      --- %%%%%

            resFig = figure;
            resAx  = axes(resFig);
            plot(resAx,Omega, Emax, 'k-+')
            xlabel('$\Omega$ [rad/s]')
            ylabel('$max(E_{mec}(t))$ [J]')
            grid off

            %%%%%   --- Delete after simulation?    --- %%%%%

            if deleteAfter == 'Y'
                for i = 1:nsimul
                    cmd = sprintf('rm %s', output{i});
                    disp(cmd)
                    system(strcat("wsl ",cmd));
                end
            end

    end
case {'4a','4b'}

    % Parametres physiques :
    tFin      = 100;
    d         = 0.005;
    kappa     = 0.05;
    m         = 0.1;
    theta0    = 0;
    thetadot0 = 0.01;

    Omega     = 2*w0;

    % Parametres numeriques :
    Dt       = 0.0001;

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
        case '4a'

            % Name of output file to generate
            name = [Ex,'.out'];

            if (exist(name,'file') ~= 2) %test if the file exists
                %if not:
                cmd = sprintf('%s%s %s %s output=%s', repertoire, executable, inputName,config ,name);
                system(strcat("wsl ",cmd)); % Wsl to compile using gcc on the wsl (windows subsystem for linux)
            end

            data = load(name); % Load generated file

            t        = data(:,1);
            theta    = data(:,2);
            thetaDot = data(:,3);
            emec     = data(:,4);
            pnc      = data(:,5);

            demec =diff(emec)/Dt; %Derivative of mechanical energy


            mecEnThmFig = figure;
            mecEnThmax  = axes(mecEnThmFig);
            mecEnThmScatter = scatter(mecEnThmax,demec,pnc(1:end-1),10       ,...
                'Marker'         , 'o'          ,...
                'MarkerFaceColor', colors(5,:)  ,...
                'MarkerEdgeColor', colors(5,:)  ,...
                'MarkerFaceAlpha', .1          ,...
                'MarkerEdgeAlpha', .1);
            hold on;
            id = plot(mecEnThmax,demec,demec,'-');
            id.LineWidth = 1.5;
            axis tight;


            dmecEnPncFig = figure;
            dmecEnPncax  = axes(dmecEnPncFig);
            dmecEnScatter = scatter(dmecEnPncax,t(1:end-1),1*demec,10 ,...
                'LineWidth'      , 1            ,...
                'Marker'         , 'o'          ,...
                'MarkerFaceColor', colors(1,:)  ,...
                'MarkerEdgeColor', colors(1,:)  ,...
                'MarkerFaceAlpha', .1          ,...
                'MarkerEdgeAlpha', .1);
            hold on;
            PncScatter = scatter(dmecEnPncax,t,pnc,10 ,...
                'LineWidth'      , 1            ,...
                'Marker'         , 'o'          ,...
                'MarkerFaceColor', colors(5,:)  ,...
                'MarkerEdgeColor', colors(5,:)  ,...
                'MarkerFaceAlpha', .1          ,...
                'MarkerEdgeAlpha', .1);

            if deleteAfter == 'Y'
            cmd = sprintf('rm %s', name);
            disp(cmd)
            system(strcat("wsl ",cmd));
            end
        case '4b'

            nsimul = 21; % Nombre de simulations a faire

            voisinage = linspace(0.98,1.02,nsimul);
            Omega     = 2*w0.*voisinage;

            paramstr = 'Omega'; % Nom du parametre a scanner (changer ici 'dt' ou 'Omega' ou autre)
            param    = Omega;

            output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie

            for i = 1:nsimul
                output{i} = [Ex,paramstr, '=', num2str(param(i)), '.out'];
                if (exist(output{i},'file') ~= 2) %test if the file exists
                %if not:
                % Execution du programme en lui envoyant la valeur a scanner en argument
                cmd = sprintf('%s%s %s %s %s=%.15g output=%s', repertoire, executable, inputName,config ,paramstr, param(i), output{i});
                disp(cmd)
                system(strcat("wsl ",cmd));
                end
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

            resFig = figure;
            resAx  = axes(resFig);
            plot(resAx,Omega, Emax, 'k-+')
            xlabel('$\Omega$ [rad/s]')
            ylabel('$max(E_{mec}(t))$ [J]')
            grid off

            if deleteAfter == 'Y'
                for i = 1:nsimul
                    cmd = sprintf('rm %s', output{i});
                    disp(cmd)
                    system(strcat("wsl ",cmd));
                end
            end
    end

case {'5a','5b','5c','5d','5e'}

    % Parametres physiques :
    d         = 0.005;
    kappa     = 0.05;
    theta0    = 1e-6;
    thetadot0 = 0.01;
    omega0    = sqrt(g/L);
    n         = 100;

    Omega     = 2*w0;

    % Parametres numeriques :
    dt = tfin ./ (1e3:450:1e4-450);

    switch Ex
        case '5a'

            %%%%%  --- SIMULATION ---   %%%%%
            nsimul = 20; % Nombre de simulations a faire

            % TODO: Choisir des valeurs de dt pour faire une etude de convergence
            dt = logspace(-2, -4, nsimul);
            % Convergence (e) : dt = tfin ./ (1e3:450:1e4-450);
            thetavar0 = linspace(pi/(nsimul+1), pi - pi/(nsimul+1), nsimul);

            paramstr = 'dt'; % Nom du parametre a scanner (changer ici 'dt' ou 'Omega' ou autre)
            param = dt; % Valeurs du parametre a scanner (changer ici dt ou Omega ou autre)

            tfin   = 20;

            config = sprintf(  ['%s=%.15g %s=%.15g %s=%.15g %s=%.15g %s=%.15g'...
                        ' %s=%.15g %s=%.15g %s=%.15g %s=%.15g'] , ...
                        'tFin'      ,tFin       ,...
                        'd'         ,d          ,...
                        'Omega'     ,Omega      ,...
                        'kappa'     ,kappa      ,...
                        'm'         ,m          ,...
                        'g'         ,g          ,...
                        'L'         ,L          ,...
                        'theta0'    ,theta0     ,...
                        'thetadot0' ,thetadot0  );


            output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie

            for i = 1:nsimul
                output{i} = [Ex,paramstr, '=', num2str(param(i)), '.out'];
                if (exist(output{i},'file') ~= 2) %test if the file exists
                %if not:
                % Execution du programme en lui envoyant la valeur a scanner en argument
                cmd = sprintf('%s%s %s %s %s=%.15g output=%s', repertoire, executable, inputName,config ,paramstr, param(i), output{i});
                disp(cmd)
                system(strcat("wsl ",cmd));
                end
            end

            %% Etude de convergence (a) %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            error = zeros(1,nsimul);

            for i = 1:nsimul % Parcours des resultats de toutes les simulations
                data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
                time  = data(:,1);
                theta = data(:,2);
                thetadot = data(:,3);

                error(i) = max(abs(theta - theta0*cos(omega0*time))); % TODO: Calculer l'erreur a partir de l'output

            end

            %%  - Other verification method of the theorem  - %%
            convFig = figure;
            convax  = axes(convFig);

            convScatter = scatter(convax,dt.^2, error,5 ,...
                'LineWidth'      , 3            ,...
                'Marker'         , 'o'          ,...
                'MarkerFaceColor', colors(2,:)  ,...
                'MarkerEdgeColor', colors(2,:) );

            %%%%%   --- Delete after simulation?    --- %%%%%

            if deleteAfter == 'Y'
                for i = 1:nsimul
                    cmd = sprintf('rm %s', output{i});
                    disp(cmd)
                    system(strcat("wsl ",cmd));
                end
            end

        case '5b'

            T      = ones(1,nsimul);
            Tth    = 4/omega0 * ellipticK(sin(thetavar0/2).*sin(thetavar0/2));
            erreur = ones(1,nsimul);

            last    = 250;
            data    = load(output{20});
            t       = data(:,1);
            theta   = data(:,2);

            figure
            plot(t(1:last), theta(1:last), 'b+--', [t(1) t(last)], [0, 0], 'k')
            ylabel('$\theta$ [rad]')
            xlabel('$t$ [s]')
            grid on
            hold on

            sgn   = sign(theta(1));
            for j = 2:last
                if not(sgn == sign(theta(j)))
                    plot(t(j), theta(j), 'r+')
                    sgn = sign(theta(j));
                end
            end


            for i = 1:nsimul % Parcours des resultats de toutes les simulations
                data    = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
                t       = data(:,1);
                theta   = data(:,2);
                sgn     = sign(theta(1));
                nchange = 0;
                for j = 2:length(theta)
                    if not(sgn == sign(theta(j)))
                        
                        nchange = nchange + 1;
                        tchange = (t(j-1) + t(j))/2;
                        sgn = sign(theta(j));
                    end
                end
                T(i)      = 4*tchange/(2*nchange-1);
                erreur(i) = T(i) - Tth(i);
            end
            
            hold off
            
            figure
            plot(thetavar0,T,'b+')
            xlabel('$\theta_0$ [rad]')
            ylabel('Period $T(\theta_0)$ [s]')
            grid on
            
            figure
            plot(thetavar0,Tth,'r+')
            xlabel('$\theta_0$ [rad]')
            ylabel('Period $T(\theta_0)$ [s]')
            grid on
            
            figure
            plot(thetavar0,erreur,'k+')
            xlabel('$\theta_0$ [rad]')
            ylabel('Error [s]')
            grid on
    
    end





otherwise


%{

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
    cmd = sprintf('%s%s %s %s %s=%.15g output=%s', repertoire, executable, inputName,config ,paramstr, param(i), output{i});
    disp(cmd)
    system(strcat("wsl ",cmd));
end
else
    cmd = sprintf('%s%s %s %s output=%s', repertoire, executable, inputName,config ,name);
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

%}
end