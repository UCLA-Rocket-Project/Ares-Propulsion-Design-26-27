% Sizing Script for Coaxial Swirl Injector Elements
% Written for internal mixing types (stage 1 nozzle recessed into stage 2)
% Stage 1, OX; Stage 2, Fuel

clc; clear; close all;

%% Setup

% Read Charts

muACorrelation = readmatrix('muACorrelation.csv');
A = muACorrelation(:, 1); 
mu = muACorrelation(:, 2);
[Araw, idx] = sort(A);
Chart.mu = mu(idx);
[Chart.A, ia] = unique(Araw);          
Mu_u = Chart.mu(ia);                 
Chart.AQ = linspace(min(Chart.A), max(Chart.A), 500);
Chart.mu = interp1(Chart.A, Mu_u, Chart.AQ, 'pchip');

% Design Parameters

Param.numElements = 9;

    % Imperial
Param.PcUS = 370; % Psia
Param.mdotTotUS = 6.639; % lbm/s
Param.OF = 1.4;
Param.massFrac = 0.75;
Param.mdotFUS = Param.mdotTotUS * Param.massFrac; % lbm/s
Param.mdotOXUS = Param.mdotTotUS * (1 - Param.massFrac); % lbm/s

    % SI
Param.Pc = convpres(Param.PcUS, 'psi', 'Pa'); % Pa
Fuel.mdot = convmass(Param.mdotFUS, 'lbm', 'kg') / Param.numElements; % kg/s
Ox.mdot = convmass(Param.mdotOXUS, 'lbm', 'kg') / Param.numElements; % kg/s
Fuel.stiffness = 0.2;
Ox.stiffness = 0.2;
Fuel.dP = Fuel.stiffness*Param.Pc;
Ox.dP = Ox.stiffness*Param.Pc;
Fuel.manifoldPressure = Param.Pc + Fuel.dP;
Ox.manifoldPressure = Param.Pc + Ox.dP;

% Chosen Parameters

Fuel.sprayConeAngle = 90; % 90-120 deg
Ox.sprayConeAngle = 60; %60-80 deg
Fuel.numInPass = 2;
Ox.numInPass = 2;
Fuel.inLengthFactor = 3; % 3-6, factor to mult inlet radius for inlet length
Ox.inLengthFactor = 3; 
Fuel.nozLengthFactor = 0.5; % 0.5-2 factor to mult nozzle radius for nozzle length
Ox.nozLengthFactor = 0.5;
Fuel.chamLengthFactor = 3; % >2 factor to mult inlet radius vortex-chamber length
Ox.chamLengthFactor = 3;
Fuel.RIn = 3; % Radial distance from centerline to inlet center ( = 3 for closed, 0.7 <= x <= 0.8 for open)
Ox.RInFactor = 0.7;

% Fluid Properties

%Fuel.rho
%Ox.rho

% Array Initialization
Array.AOx = [];
Array.muOx = [];
Array.RnOx = [];

%% Main

stageOneParam(Ox, Chart);

%% Functions
function stageOneParam(Ox, Chart)
    
    Atol = 1;
    AError = realmax;
    AGuess = 5;% Geometric Parameter A
    while (abs(AError) > Atol)
        Ox.mu = interp1(Chart.AQ, Chart.mu, AGuess, 'pchip'); % Mass flow coefficient
        Ox.Rn = 0.475 * sqrt(Ox.mdot / (Ox.mu...
            * sqrt(Ox.rho * Ox.dP))); % Nozzle radius (103)
        Ox.rIn = sqrt((Ox.RIn * Ox.Rn) / (Ox.numInPass * Ox.A)); % Inlet passage radius (104)
        Ox.inletLength = Ox.rIn * Ox.inLengthFactor; % Inlet length
        Ox.nozzleLength = Ox.Rn * Ox.nozLengthFactor; % Nozzle length
        Ox.chamLength = Ox.rIn * Ox.chamLengthFactor; % Vortex chamber length
    
        Ox.Rin = Ox.RinFactor * Ox.Rn; % Radial distance from centerline to inlet center
        Ox.Rs = Ox.Rin + Ox.rIn; % Vortex chamber radius
    
        Ox.vIn = (Ox.mdot/Ox.numInPass) / (Ox.rho * pi * Ox.rIn^2); % Inlet fluid velocity
        Ox.ReIn = 0.637 * ((Ox.mdot/Ox.numInPass) / (sqrt(Ox.numInPass) * Ox.rIn...
            * Ox.rhoOx * Ox.vIn)); % Reynolds in inlet passages
        Ox.lambda = 0.3164 / (Ox.ReIn)^0.25; % Friction factor

        Ox.AEq = (Ox.Rin * Ox.Rn) / (Ox.numInPass * Ox.rIn^2 +...
            Ox.lambda * Ox.Rin * (Ox.Rin - Ox.Rn)); % A corrected for viscous losses (100)
        Ox.tilt = 90 - atand(Ox.Rs / Ox.inletLength); % Inlet tilting angle
        Ox.muEq = interp1(Chart.AQ, Chart.mu, Ox.AEq, 'pchip'); % Mass flow coefficient corrected for viscous losses
        Ox.alphaEq = NaN;% Spray cone angle coefficient corrected for viscous losses
        Ox.xiIn = -(0.4 * Ox.tilt)/60 + 0.9; % Hydrolic-loss coefficient (Fig. 25)
        Ox.xi = Ox.xiIn + Ox.lambda * (Ox.inletLength / (2 * Ox.rIn));
    
        Ox.mu = Ox.muEq / sqrt(1 +...
            Ox.xi * Ox.muEq^2 * AGuess^2 / Ox.RInFactor^2);
        Ox.Rn = 0.475 * sqrt(Ox.mdot / (Ox.mu...
            * sqrt(Ox.rho * Ox.dP)));
        Ox.A = (Ox.Rin * Ox.Rn) / (Ox.numInPass * Ox.rIn^2);

        AError = Ox.A - AGuess;
        AGuess = Ox.A + 0.1*AError;
    end
    
end

