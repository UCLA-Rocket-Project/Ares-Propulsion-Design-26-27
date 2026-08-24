% Sizing Script for Coaxial Swirl Injector Elements
% Written for externally mixing types (spray cones do not mix before
% nozzle)

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

Param.numElements = 4;
Param.minFlangeThickness = 0;
Param.wallThickness = 0;

    % Imperial
Param.PcUS = 370; % Psia
Param.mdotTotUS = 6.639; % lbm/s
Param.OF = 1.4;
Param.mdotFUS = Param.mdotTotUS/(1+Param.OF); % lbm/s
Param.mdotOXUS = (Param.mdotTotUS*Param.OF)/(1+Param.OF); % lbm/s

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

% Fluid Properties

Fuel.temp = 293; % Change depending on regen script
Ox.temp = 90; % Change depending on predicted LOX temp
Fuel.viscosity = py.CoolProp.CoolProp.PropsSI('V', 'T', Fuel.temp, 'P', Fuel.manifoldPressure, 'Ethanol');
Ox.viscosity = py.CoolProp.CoolProp.PropsSI('V', 'T', Ox.temp, 'P', Ox.manifoldPressure, 'Oxygen');
Fuel.rho = py.CoolProp.CoolProp.PropsSI('D', 'T', Fuel.temp, 'P', Fuel.manifoldPressure, 'Ethanol'); 
Ox.rho = py.CoolProp.CoolProp.PropsSI('D', 'T', Ox.temp, 'P', Ox.manifoldPressure, 'Oxygen'); 

% Chosen Parameters

Ox.sprayConeHalfAngle = 50;
Fuel.sprayConeHalfAngle = Ox.sprayConeHalfAngle - 7.5;
Fuel.numInPass = 2;
Ox.numInPass = 8;
Fuel.inLengthFactor = 3; % 3-6, factor to mult inlet radius for inlet length
Ox.inLengthFactor = 3; 
Fuel.nozLengthFactor = 0.5; % 0.5-2 factor to mult nozzle radius for nozzle length
Ox.nozLengthFactor = 2;
Fuel.chamLengthFactor = 2; % >2 factor to mult inlet radius vortex-chamber length
Ox.chamLengthFactor = 2;
Fuel.RInFactor = 3; % Radial distance from centerline to inlet center; affects nozzle/chamber length ratio; 3 for closed, 0.7-0.8 for open
Ox.RInFactor = 0.7;

%% Main

report = 'swirlOutputGeo.txt';
Ox = stageOneParam(Ox, Chart);
Fuel = stageOneParam(Fuel, Chart);
checkReportGeometry(Ox, Fuel, report, Param);

%% Functions
function species = stageOneParam(species, Chart)
    
    Atol = 0.1;
    AError = realmax;
    AGuess = interp1(Chart.alpha, Chart.AQalpha, species.sprayConeHalfAngle, 'pchip');
    while (abs(AError) > Atol)
        species.mu = interp1(Chart.AQmu, Chart.mu, AGuess, 'pchip'); % Mass flow coefficient
        species.Rn = 0.475 * sqrt(species.mdot / (species.mu...
            * sqrt(species.rho * species.dP))); % Nozzle radius (103)
        species.RIn = species.RInFactor * species.Rn; % Radial distance from centerline to inlet center
        species.rIn = sqrt((species.RIn * species.Rn) / (species.numInPass * AGuess)); % Inlet passage radius (104)
        species.inletLength = species.rIn * species.inLengthFactor; % Inlet length
        species.nozzleLength = species.Rn * species.nozLengthFactor; % Nozzle length
        species.chamLength = species.RIn * species.chamLengthFactor; % Vortex chamber length
    
        species.Rs = species.RIn + species.rIn; % Vortex chamber radius

        species.ReIn = 0.637 * (species.mdot / (sqrt(species.numInPass) * species.rIn...
            * species.viscosity)); % Reynolds in inlet passages

        species.lambda = 0.3164 / (species.ReIn)^0.25; % Friction factor

        species.AEq = (species.RIn * species.Rn) / (species.numInPass * species.rIn^2 +...
            species.lambda * species.RIn * (species.RIn - species.Rn)*0.5); % A corrected for viscous losses (100)
        species.tilt = 90 - atand(species.Rs / species.inletLength); % Inlet tilting angle
        species.muEq = interp1(Chart.AQmu, Chart.mu, species.AEq, 'pchip'); % Mass flow coefficient corrected for viscous losses
        species.alphaEq = interp1(Chart.AQalpha, Chart.alpha, species.AEq, 'pchip');% Spray cone angle coefficient corrected for viscous losses
        species.xiIn = -(0.4 * species.tilt)/60 + 0.9; % Hydrolic-loss coefficient (Fig. 25)
        species.xi = species.xiIn + species.lambda * (species.inletLength / (2 * species.rIn));
    
        species.mu = species.muEq / sqrt(1 +...
            species.xi * species.muEq^2 * AGuess^2 / species.RInFactor^2);
        species.Rn = 0.475 * sqrt(species.mdot / (species.mu...
            * sqrt(species.rho * species.dP)));
        ANextGuess = (species.RIn * species.Rn) / (species.numInPass * species.rIn^2);

        AError = AGuess - ANextGuess;
        AGuess = ANextGuess;
    end

    species.length = species.nozzleLength + species.chamLength;
end

function checkReportGeometry(Ox, Fuel, report, Param)
    reason = {};

    if (Ox.Rn + Param.wallThickness < Fuel.Rn)
        speciesStageOne = 'Ox';
        speciesStageTwo = 'Fuel';
        if (Ox.length - Fuel.length - 2*Ox.rIn < Param.minFlangeThickness)
            pass = false;
            reason = [reason, {'Flange Thickness Violated'}];
        else
            pass = true;
        end
    elseif (Fuel.Rn + Param.wallThickness < Ox.Rn)
        speciesStageOne = 'Fuel';
        speciesStageTwo = 'Ox';
        if (Fuel.length - Ox.length - 2*Fuel.rIn < Param.minFlangeThickness)
            pass = false;
            reason = [reason, {'Flange Thickness Violated'}];
        else
            pass = true;
        end
    else
        pass = false;
        reason = [reason, {'Nozzle Collision'}];
        speciesStageOne = 'N/A';
        speciesStageTwo = 'N/A';
    end

    reportWrite = fopen(report, 'w');

    if (pass == false)
        for i = 1:numel(reason)
            fprintf(reportWrite, 'NO PASS:\n%s\n' ,reason{i});
        end
        fprintf(reportWrite, '\n');
    end
    fprintf(reportWrite, 'Stage One: %s\nStage Two: %s\n', speciesStageOne,...
        speciesStageTwo);
    fprintf(reportWrite, '\nOX DIMENSIONS\nNozzle Radius: %f m\nChamber Radius: %f m\nInlet Radius: %f m\nInlet Tilt: %f deg\nNumber of Inlet Passes: %f\nInlet Position: %f m\nInlet Length: %f m\nNozzle Length: %f m\nChamber Length: %f m\nA Parameter: %f\n', ...
        Ox.Rn, Ox.Rs, Ox.rIn, Ox.tilt, Ox.numInPass, Ox.RIn, ...
        Ox.inletLength, Ox.nozzleLength, Ox.chamLength, Ox.AEq);
    fprintf(reportWrite, '\nFUEL DIMENSIONS\nNozzle Radius: %f m\nChamber Radius: %f m\nInlet Radius: %f m\nInlet Tilt: %f deg\nNumber of Inlet Passes: %f\nInlet Position: %f m\nInlet Length: %f m\nNozzle Length: %f m\nChamber Length: %f m\nA Parameter: %f\n', ...
        Fuel.Rn, Fuel.Rs, Fuel.rIn, Fuel.tilt, Fuel.numInPass, Fuel.RIn, ...
        Fuel.inletLength, Fuel.nozzleLength, Fuel.chamLength, Fuel.AEq);
end