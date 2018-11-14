% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
% 
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.
%

set(groot, 'DefaultTextInterpreter',            'LaTeX' );

%% Donnees %%
%%%%%%%%%%%%%

m      = 0.1;
g      = 9.81;
L      = 0.1;
theta0 = 1e-6;
omega0 = sqrt(g/L);
tfin   = 20; % 20 dans (a), (b) ; 40 * pi / omega0 dans (e).
% 20000 * pi / omega0
n      = 100;

%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configurationAB.in'; % Nom du fichier d'entree de base

nsimul = 20; % Nombre de simulations a faire

% SUPPRIMER LES SIMULATIONS ANTERIEURES
mode = 0; % Simulation, si mode = 1; suppression, si mode = -1.
tache = 'anperB'; % Tache a executer
% Taches : model ; convA ; compsolA ; anperB ; convE ; poincarE

dt = tfin ./ (1e3:450:1e4-450); % TODO: Choisir des valeurs de dt pour faire une etude de convergence
% Convergence (a) : dt = logspace(-2, -4, nsimul);
% Convergence (e) : dt = tfin ./ (1e3:450:1e4-450);
Omega = ones(1,nsimul); % TODO: Choisir des valeurs de Omega pour trouver la resonance
thetavar0 = linspace(pi/(nsimul+1), pi - pi/(nsimul+1), nsimul);

paramstr = 'theta0'; % Nom du parametre a scanner (changer ici 'dt' ou 'Omega' ou autre)
param = thetavar0; % Valeurs du parametre a scanner (changer ici dt ou Omega ou autre)

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i)), '.out'];
end

%% Simulations %%
%%%%%%%%%%%%%%%%%

if mode == 1
    for i = 1:nsimul
        % Execution du programme en lui envoyant la valeur a scanner en argument
        cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i});
        disp(cmd)
        system(cmd);
    end
end

%% Suppression %%
%%%%%%%%%%%%%%%%%

if mode == -1
    for i = 1:nsimul
        cmd = sprintf('rm %s', output{i});
        disp(cmd)
        system(cmd);
    end
    return
end

%% Modelisation %%
%%%%%%%%%%%%%%%%%%

if strcmp(tache, 'model')
    data = load('a.out');
    time  = data(:,1);
    theta = data(:,2);
    
    figure
    plot(L*sin(theta),-L*cos(theta),'b.')
    axis equal
    grid on
    figure
    plot3(L*sin(theta),time,-L*cos(theta),'b.')
    axis equal
    grid on
end

%% Etude de convergence (a) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(tache, 'convA')
    if strcmp(paramstr, 'dt')
        error = zeros(1,nsimul);
    elseif strcmp(paramstr, 'Omega')
        Emax = zeros(1,nsimul);
    end

    for i = 1:nsimul % Parcours des resultats de toutes les simulations
        data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
        time  = data(:,1);
        theta = data(:,2);
        thetadot = data(:,3);

        if strcmp(paramstr, 'dt')
            error(i) = max(abs(theta - theta0*cos(omega0*time))); % TODO: Calculer l'erreur a partir de l'output
        elseif strcmp(paramstr, 'Omega')
            Emax(i) = max(m*L*L*thetadot*thetadot/2 - m*g*l*cos(theta)); % TODO: Calculer le maximum de l'energie
        end
    end

    if strcmp(paramstr, 'dt')
        figure
        loglog(dt, error, 'k+')
        xlabel('\Deltat')
        ylabel('Erreur sur \theta(t_{fin}) [rad]')
        grid on
    elseif strcmp(paramstr, 'Omega')
        figure
        plot(Omega, Emax, 'k-+')
        xlabel('\Omega [rad/s]')
        ylabel('max(E_{mec}(t)) [J]')
        grid on
    end
end


%% Theta, thetadot et emec (a) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(tache, 'compsolA')
    figure(1)
    xlabel('t')
    ylabel('\theta [rad]')
    grid on
    hold on
    
    figure(2)
    xlabel('t')
    ylabel('thetadot [rad/s]')
    grid on
    hold on
    
    figure(3)
    xlabel('t')
    ylabel('E_{mec} [J]')
    grid on
    hold on
    
    for j = 1:nsimul
        data     = load(output{j});
        t        = data(:,1);
        theta    = data(:,2);
        thetadot = data(:,3);
        emec     = data(:,4);
        pnc      = data(:,5);

        figure(1)
        plot(t,theta,'.-')
        figure(2)
        plot(t,thetadot,'.-')
        figure(3)
        plot(t,emec)
    end
    thetath    = theta0*cos(omega0*t);
    thetadotth = -theta0*omega0*sin(omega0*t);
    emecth     = m*L*L*thetadotth.*thetadotth/2 - m*g*L*cos(thetath);
    
    figure(1)
    plot(t,thetath,'r--')
    hold off
    figure(2)
    plot(t,thetadotth,'r--')
    hold off
    figure(3)
    plot(t,emecth,'r--')
    hold off
end

%% Analyse de la periode (b) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(tache, 'anperB')
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

%% Etude de convergence (e) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(tache, 'convE')
    thetafin    = ones(1,nsimul);
    thetadotfin = ones(1,nsimul);
    
    for i = 1:nsimul
        data           = load(output{i});
        thetafin(i)    = data(end,2);
        thetadotfin(i) = data(end,3);
    end
    
    figure
    plot(dt.*dt,thetafin,'k+')
    xlabel('\Deltat 2')
    ylabel('\thetạa_{fin} [rad]')
    grid on
    figure
    plot(dt.*dt,thetadotfin,'k+')
    xlabel('\Deltat 2')
    ylabel('dot\theta_{fin} [rad]')
    grid on
end

%% Sections de Poincare (e) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(tache, 'poincarE')
    
    thetavar0    = -3*pi/4:pi/4:3*pi/4;
    thetadotvar0 = 0;
    k = 0.01;
    for p = 1:3
        thetavardot0 = horzcat(thetadotvar0, [1*k 3*k 5*k]);
        k = 10*k;
    end
    erreur       = 1e-8;
    
    name = cell(length(thetavar0), length(thetadotvar0));
    for i = 1:length(thetavar0)
        for j = 1:length(thetadotvar0)
            name{i,j} = ['theta0=', num2str(thetavar0(i)), '_thetadot0=', num2str(thetadotvar0(j)), '.out'];
            cmd = sprintf('%s%s %s theta0=%s thetadot0=%s dt=%s output=%s', ...
            repertoire, executable, input, num2str(thetavar0(i)), num2str(thetadotvar0(j)), num2str(2 * pi / (omega0 * n)), name{i,j});
            disp(cmd)
            system(cmd);
        end
    end
    
    data     = load(name{1,1});
    theta    = data(101:end,2);
    thetadot = data(101:end,3);
    
    figure
    plot(theta, thetadot, 'b.')
    title(sprintf('theta_0 = %g ; thetadot_0 = %d', theta0, thetadot0))
    xlabel('\theta [rad]')
    ylabel('dot \theta [rad]')
    grid on
    
    figure
    plot(L*sin(theta),-L*cos(theta),'b.')
    axis equal
    grid on
end

