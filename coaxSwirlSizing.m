% Sizing Script for Coaxial Swirl Injector Elements
% Written for internal mixing types (stage 1 nozzle recessed into stage 2)
% Stage 1, OX; Stage 2, Fuel

clc; clear; close all;

%% Setup

% Read Charts

muACorrelation = readmatrix('muACorrelation.csv');
A1 = muACorrelation(:, 1); 
mu = muACorrelation(:, 2);
[Araw, idx] = sort(A1);
Chart.mu = mu(idx);
[Chart.Amu, ia] = unique(Araw);          
MuU = Chart.mu(ia);                 
Chart.AQmu = linspace(min(Chart.Amu), max(Chart.Amu), 500);
Chart.mu = interp1(Chart.Amu, MuU, Chart.AQmu, 'pchip');

AalphaCorrelation = readmatrix('AalphaCorrelation.csv');
A2 = AalphaCorrelation(:, 1); 
alpha = AalphaCorrelation(:, 2);
[Araw, idx] = sort(A2);
Chart.alpha = alpha(idx);
[Chart.Aalpha, ia] = unique(Araw);          
alphaU = Chart.alpha(ia);                 
Chart.AQalpha = linspace(min(Chart.Aalpha), max(Chart.Aalpha), 500);
Chart.alpha = interp1(Chart.Aalpha, alphaU, Chart.AQalpha, 'pchip');

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

Ox.sprayConeHalfAngle = 40; %30-40 deg
Fuel.numInPass = 3;
Ox.numInPass = 3;
Fuel.inLengthFactor = 3; % 3-6, factor to mult inlet radius for inlet length
Ox.inLengthFactor = 3; 
Fuel.nozLengthFactor = 0.5; % 0.5-2 factor to mult nozzle radius for nozzle length
Ox.nozLengthFactor = 0.5;
Fuel.chamLengthFactor = 3; % >2 factor to mult inlet radius vortex-chamber length
Ox.chamLengthFactor = 3;
Fuel.RIn = 3; % Radial distance from centerline to inlet center
Ox.RInFactor = 3; % May want to iterate to get desired spray angle

% Fluid Properties

Fuel.temp = 293; % Change depending on regen script
Ox.temp = 90; % Change depending on predicted LOX temp
Fuel.viscosity = py.CoolProp.CoolProp.PropsSI('V', 'T', Fuel.temp, 'P', Fuel.manifoldPressure, 'Ethanol');
Ox.viscosity = py.CoolProp.CoolProp.PropsSI('V', 'T', Ox.temp, 'P', Ox.manifoldPressure, 'Oxygen');
Ox.rho = py.CoolProp.CoolProp.PropsSI('D', 'T', Ox.temp, 'P', Ox.manifoldPressure, 'Oxygen'); 

% Array Initialization
Array.AOx = [];
Array.muOx = [];
Array.RnOx = [];

%% Main

Ox = stageOneParam(Ox, Chart);
disp(Ox);

%% Functions
function Ox = stageOneParam(Ox, Chart)
    
    Atol = 10^-10;
    AError = realmax;
    AGuess = interp1(Chart.alpha, Chart.AQalpha, Ox.sprayConeHalfAngle, 'pchip');
    while (abs(AError) > Atol)
        Ox.mu = interp1(Chart.AQmu, Chart.mu, AGuess, 'pchip'); % Mass flow coefficient
        Ox.Rn = 0.475 * sqrt(Ox.mdot / (Ox.mu...
            * sqrt(Ox.rho * Ox.dP))); % Nozzle radius (103)
        Ox.Rin = Ox.RInFactor * Ox.Rn; % Radial distance from centerline to inlet center
        Ox.rIn = sqrt((Ox.Rin * Ox.Rn) / (Ox.numInPass * AGuess)); % Inlet passage radius (104)
        Ox.inletLength = Ox.rIn * Ox.inLengthFactor; % Inlet length
        Ox.nozzleLength = Ox.Rn * Ox.nozLengthFactor; % Nozzle length
        Ox.chamLength = Ox.Rin * Ox.chamLengthFactor; % Vortex chamber length
    
        Ox.Rs = Ox.Rin + Ox.rIn; % Vortex chamber radius

        Ox.ReIn = 0.637 * (Ox.mdot / (sqrt(Ox.numInPass) * Ox.rIn...
            * Ox.viscosity)); % Reynolds in inlet passages

        Ox.lambda = 0.3164 / (Ox.ReIn)^0.25; % Friction factor

        Ox.AEq = (Ox.Rin * Ox.Rn) / (Ox.numInPass * Ox.rIn^2 +...
            Ox.lambda * Ox.Rin * (Ox.Rin - Ox.Rn)*0.5); % A corrected for viscous losses (100)
        Ox.tilt = 90 - atand(Ox.Rs / Ox.inletLength); % Inlet tilting angle
        Ox.muEq = interp1(Chart.AQmu, Chart.mu, Ox.AEq, 'pchip'); % Mass flow coefficient corrected for viscous losses
        Ox.alphaEq = interp1(Chart.AQalpha, Chart.alpha, Ox.AEq, 'pchip');% Spray cone angle coefficient corrected for viscous losses
        Ox.xiIn = -(0.4 * Ox.tilt)/60 + 0.9; % Hydrolic-loss coefficient (Fig. 25)
        Ox.xi = Ox.xiIn + Ox.lambda * (Ox.inletLength / (2 * Ox.rIn));
    
        Ox.mu = Ox.muEq / sqrt(1 +...
            Ox.xi * Ox.muEq^2 * AGuess^2 / Ox.RInFactor^2);
        Ox.Rn = 0.475 * sqrt(Ox.mdot / (Ox.mu...
            * sqrt(Ox.rho * Ox.dP)));
        Ox.A = (Ox.Rin * Ox.Rn) / (Ox.numInPass * Ox.rIn^2);

        AError = AGuess - Ox.A;
        AGuess = Ox.A + 0.1*AError;
    end
end