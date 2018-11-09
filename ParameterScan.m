% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
% 
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.
%

%% Donnees %%
%%%%%%%%%%%%%

m      = 0.1;
g      = 9.81;
L      = 0.1;
theta0 = 1e-6;
omega0 = sqrt(g/L);
tfin   = 20000 * pi / omega0; % 20 dans (a), (b) ; 40 * pi / omega0 dans (e).
n      = 5;

%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configurationE.in'; % Nom du fichier d'entree de base

nsimul = 20; % Nombre de simulations a faire

% SUPPRIMER LES SIMULATIONS ANTERIEURES
mode = 0; % Simulation, si mode = 1; suppression, si mode = -1.
tache = 'poincarE'; % Tache a executer
% Taches : model ; convA ; compsolA ; anperB ; convE ; poincarE

dt = tfin ./ (1e3:450:1e4-450); % TODO: Choisir des valeurs de dt pour faire une etude de convergence
% Convergence (a) : dt = logspace(-2, -4, nsimul);
% Convergence (e) : dt = tfin ./ (1e3:450:1e4-450);
Omega = ones(1,nsimul); % TODO: Choisir des valeurs de Omega pour trouver la resonance
thetavar0 = linspace(pi/(nsimul+1), pi - pi/(nsimul+1), nsimul);

paramstr = 'dt'; % Nom du parametre a scanner (changer ici 'dt' ou 'Omega' ou autre)
param = dt; % Valeurs du parametre a scanner (changer ici dt ou Omega ou autre)

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
    data = load(output{1});
    time  = data(:,1);
    theta = data(:,2);
    
    figure
    plot(L*sin(theta),-L*cos(theta),'b.')
    grid on
    figure
    plot3(L*sin(theta),time,-L*cos(theta),'b.')
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
    figure
    plot(thetavar0,T,'b+',thetavar0,Tth,'r+')
    xlabel('\theta_0 [rad]')
    ylabel('Periode T(\theta_0)')
    grid on
    
    figure
    plot(thetavar0,erreur,'k+')
    xlabel('\theta_0 [rad]')
    ylabel('Erreur sur la periode')
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

%% Sections de Poincare %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(tache, 'poincarE')
    cmd = sprintf('%s%s %s dt=%s output=%s', repertoire, executable, input, num2str(2 * pi / (omega0 * n)), 'a.out');
    disp(cmd)
    system(cmd);
    
    data     = load('a.out');
    theta    = data(101:end,2);
    thetadot = data(101:end,3);
    
    figure
    plot(theta, thetadot, 'b.')
    grid on
end

