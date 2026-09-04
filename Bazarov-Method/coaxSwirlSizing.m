% Sizing Script for Coaxial Swirl Injector Elements
% Written for externally mixing types (spray cones do not mix before
% nozzle)

clc; clear; close all;

%% Setup

% Read for 2D Interpolations

angleRbarA1 = readmatrix('2alphaRbarA1.csv');
angleRbarA2 = readmatrix('2alphaRbarA2.csv');
angleRbarA3 = readmatrix('2alphaRbarA3.csv');
angleRbarA4 = readmatrix('2alphaRbarA4.csv');
Chart.angleRbarA = {angleRbarA1, angleRbarA2, angleRbarA3, angleRbarA4};

muRbarA1 = readmatrix('muRbarA1.csv');
muRbarA2 = readmatrix('muRbarA2.csv');
muRbarA3 = readmatrix('muRbarA3.csv');
muRbarA4 = readmatrix('muRbarA4.csv');
muRbarA6 = readmatrix('muRbarA6.csv');
Chart.muRbarA = {muRbarA1, muRbarA2, muRbarA3, muRbarA4, muRbarA6};

% Read for 1D Interpolations

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
Param.minFlangeThickness = 0.01;
Param.wallThickness = 0.0005;

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

Fuel.temp = 400; % Change depending on regen script
Ox.temp = 90; % Change depending on predicted LOX temp
Fuel.viscosity = py.CoolProp.CoolProp.PropsSI('V', 'T', Fuel.temp, 'P', Fuel.manifoldPressure, 'Ethanol');
Ox.viscosity = py.CoolProp.CoolProp.PropsSI('V', 'T', Ox.temp, 'P', Ox.manifoldPressure, 'Oxygen');
Fuel.rho = py.CoolProp.CoolProp.PropsSI('D', 'T', Fuel.temp, 'P', Fuel.manifoldPressure, 'Ethanol'); 
Ox.rho = py.CoolProp.CoolProp.PropsSI('D', 'T', Ox.temp, 'P', Ox.manifoldPressure, 'Oxygen'); 

% Chosen Parameters

Fuel.inLengthFactor = 3; % 3-6, factor to mult inlet radius for inlet length
Ox.inLengthFactor = 3; 
Fuel.nozLengthFactor = 0.5; % 0.5-2 factor to mult nozzle radius for nozzle length
Ox.nozLengthFactor = 2;
Fuel.chamLengthFactor = 2; % >2 factor to mult inlet radius vortex-chamber length
Ox.chamLengthFactor = 6;
Ox.sprayConeHalfAngle = 50;
Fuel.sprayConeHalfAngle = Ox.sprayConeHalfAngle - 7.5;

%% Main

report = 'swirlOutputGeo.txt';
findValidGeometry(Ox, Fuel, Chart, report, Param);

%% Functions
function species = stageOneParam(species, Chart)
    
    Atol = 0.001;
    AError = realmax;
    AGuess = interp1(Chart.alpha, Chart.AQalpha, species.sprayConeHalfAngle, 'pchip');
    species.isConverged = false;
    while (abs(AError) > Atol)
        if (~isreal(AGuess))
            return;
        end
        % Prepare guessed geometry
        species.mu = interp1(Chart.AQmu, Chart.mu, AGuess, 'pchip'); % Mass flow coefficient
        species.Rn = 0.475 * sqrt(species.mdot / (species.mu...
            * sqrt(species.rho * species.dP))); % Nozzle radius (103)
        species.RIn = species.RInFactor * species.Rn; % Radial distance from centerline to inlet center
        species.rIn = sqrt((species.RIn * species.Rn) / (species.numInPass * AGuess)); % Inlet passage radius (104)
        species.inletLength = species.rIn * species.inLengthFactor; % Inlet length
        species.nozzleLength = species.Rn * species.nozLengthFactor; % Nozzle length
        species.chamLength = species.RIn * species.chamLengthFactor; % Vortex chamber length
        species.Rs = species.RIn + species.rIn; % Vortex chamber radius

        % Calculate relevant fluid properties
        species.ReIn = 0.637 * (species.mdot / (sqrt(species.numInPass) * species.rIn...
            * species.viscosity)); % Reynolds in inlet passages
        species.lambda = 0.3164 / (species.ReIn)^0.25; % Friction factor

        % Correct geometry for viscous losses
        species.AEq = (species.RIn * species.Rn) / (species.numInPass * species.rIn^2 +...
            species.lambda * species.RIn * (species.RIn - species.Rn)*0.5); % A corrected for viscous losses (100)
        species.tilt = 90 - atand(species.Rs / species.inletLength); % Inlet tilting angle
        species.muEq = interp1(Chart.AQmu, Chart.mu, species.AEq, 'pchip'); % Mass flow coefficient corrected for viscous losses
        species.alphaEq = interp1(Chart.AQalpha, Chart.alpha, species.AEq, 'pchip');% Spray cone angle coefficient corrected for viscous losses
        species.xiIn = -(0.4 * species.tilt)/60 + 0.9; % Hydrolic-loss coefficient (Fig. 25)
        species.xi = species.xiIn + species.lambda * (species.inletLength / (2 * species.rIn));

        % Calculate next guess using corrected geometry
        species.mu = species.muEq / sqrt(1 +...
            species.xi * species.muEq^2 * AGuess^2 / species.RInFactor^2);
        species.Rn = 0.475 * sqrt(species.mdot / (species.mu...
            * sqrt(species.rho * species.dP)));
        ANextGuess = (species.RIn * species.Rn) / (species.numInPass * species.rIn^2);

        AError = AGuess - ANextGuess;
        AGuess = ANextGuess;
    end
    if ~isnan(species.Rn)
        species.isConverged = true;
    end
end

function species = stageTwoParam(species, speciesOne, Chart, Param)
    species.isConverged = false;

    if (~isreal(speciesOne.Rn) || isnan(speciesOne.Rn))
        return;
    end
    oneExR = speciesOne.Rn + Param.wallThickness;

    nozzleSpacingGuess = 0.0003;
    stepSize = 0.000005;
    tol = 6000; % Pa
    dPError = realmax;

    while (abs(dPError) > tol || isnan(species.A) || isnan(species.mu)) && nozzleSpacingGuess < 0.01
        nozzleSpacingGuess = nozzleSpacingGuess + stepSize;
        species.Rn = oneExR + nozzleSpacingGuess;
        species.RBar = species.Rn/oneExR;
        species.A = thirtySevenInterpA(Chart, species.RBar, species.sprayConeHalfAngle); % Fig. 37, alpha must be < 50, > 10
        species.mu = thirtyEightInterp(Chart, species.RBar, species.A);
        dP = 0.05*((species.mdot^2)/(species.mu^2 * species.rho * species.Rn^4));
        dPError = dP - species.dP;
    end
    
    if (nozzleSpacingGuess ~= 0.01)
        species.isConverged = true;
    end

    species.RIn = species.RInFactor * species.Rn;
    species.rIn = sqrt((species.RIn * species.Rn) / (species.numInPass * species.A));
    species.alphaCalc = thirtySevenInterpAlpha(Chart, species.RBar, species.A);
    species.inletLength = species.rIn * species.inLengthFactor; % Inlet length
    species.nozzleLength = species.Rn * species.nozLengthFactor;
    species.chamLength = species.RIn * species.chamLengthFactor;
    species.Rs = species.RIn + species.rIn;
    species.tilt = 90 - atand(species.Rs / species.inletLength);
end

function reportGeometry(Ox, Fuel, reportWrite, Param, numValid)
    reason = {};
    Ox.length = Ox.chamLength + Ox.nozzleLength;
    Fuel.length = Fuel.chamLength + Fuel.nozzleLength;

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

    if (Fuel.Rs < Ox.Rs)
        pass = false;
        reason = [reason, {'Chamber Violation'}];
    end

    if (pass == false)
        for i = 1:numel(reason)
            fprintf(reportWrite, 'NO PASS:\n%s\n' ,reason{i});
        end
        fprintf(reportWrite, '\n');
    else
        fprintf(reportWrite, 'YES PASS\n\n');
    end
    
    fprintf(reportWrite, 'GEOMETRY #%d\nStage One: %s\nStage Two: %s\n', numValid, speciesStageOne,...
        speciesStageTwo);
    fprintf(reportWrite, '\nOX DIMENSIONS\nNozzle Radius: %f m\nChamber Radius: %f m\nInlet Radius: %f m\nInlet Tilt: %f deg\nNumber of Inlet Passes: %f\nInlet Position: %f m\nInlet Length: %f m\nNozzle Length: %f m\nChamber Length: %f m\nA Parameter: %f\nSpray Cone Angle: %f deg\n', ...
        Ox.Rn, Ox.Rs, Ox.rIn, Ox.tilt, Ox.numInPass, Ox.RIn, ...
        Ox.inletLength, Ox.nozzleLength, Ox.chamLength, Ox.AEq, Ox.alphaEq*2);
    fprintf(reportWrite, '\nFUEL DIMENSIONS\nNozzle Radius: %f m\nChamber Radius: %f m\nInlet Radius: %f m\nInlet Tilt: %f deg\nNumber of Inlet Passes: %f\nInlet Position: %f m\nInlet Length: %f m\nNozzle Length: %f m\nChamber Length: %f m\nA Parameter: %f\nSpray Cone Angle: %f deg\n\nEND\n\n', ...
        Fuel.Rn, Fuel.Rs, Fuel.rIn, Fuel.tilt, Fuel.numInPass, Fuel.RIn, ...
        Fuel.inletLength, Fuel.nozzleLength, Fuel.chamLength, Fuel.A, Fuel.alphaCalc*2);
    disp(numValid);
end

function findValidGeometry(speciesOne, speciesTwo, Chart, report, Param)
    numValid = 0;
    reportWrite = fopen(report, 'w');
    %for (angle = 45:60)
        %speciesOne.sprayConeHalfAngle = angle;
        %speciesTwo.sprayConeHalfAngle = speciesOne.sprayConeHalfAngle - 7.5;
        for i = 0.7: 0.1: 3 % Sweep through stage one geometries
            for j = 2:1:6
                speciesOne.numInPass = j;
                speciesOne.RInFactor = i;
                speciesOne = stageOneParam(speciesOne, Chart);
                if speciesOne.alphaEq > (speciesOne.sprayConeHalfAngle - 0.5) && ...
                        speciesOne.alphaEq < (speciesOne.sprayConeHalfAngle + 0.5) &&...
                        speciesOne.isConverged && speciesOne.AEq >= 1 && ...
                        speciesOne.AEq <= 13

                    for i2 = 0.7: 0.1 :3 % Sweep through stage two geometries only on valid stage one geometries
                        for j2 = 2:1:6
                            speciesTwo.numInPass = j2;
                            speciesTwo.RInFactor = i2;
                            speciesTwo = stageTwoParam(speciesTwo, speciesOne, Chart, Param);

                            if speciesTwo.isConverged && speciesTwo.A >= 1 && speciesTwo.A <= 13
                                numValid = numValid + 1;
                                disp(numValid);
                                reportGeometry(speciesOne, speciesTwo, reportWrite, Param, numValid);
                            end
                        end
                    end
                end
            end
        end
    %end
end

% 2D Interpolations

function A = thirtySevenInterpA(Chart, RbarQuery, alpha)
    targetAngle = 2*alpha;

    AValues = [];
    anglesAtRbar = [];

    for i = 1:4 % Get each of the 4 angles corresponding with Rbar at each curve
        currentCurve = Chart.angleRbarA{i};
        Rbar = currentCurve(:, 1);
        if RbarQuery < min(Rbar) || RbarQuery > max(Rbar) % Check if RbarQuery is outside range
            continue
        end
        angle = currentCurve(:, 2);
        anglesAtRbar(end + 1) =  interp1(Rbar, angle, RbarQuery, 'linear');
        AValues(end + 1) = i;
    end

    if numel(anglesAtRbar) < 2 % Check if there are curves to interpolate
        A = NaN;
        return;
    end

    A = interp1(anglesAtRbar, AValues, targetAngle, 'linear');
end

function alpha = thirtySevenInterpAlpha(Chart, RbarQuery, AQuery)
    AValues = [];
    anglesAtRbar = [];

    for i = 1:4 % Get each of the 4 angles corresponding with Rbar at each curve
        currentCurve = Chart.angleRbarA{i};
        Rbar = currentCurve(:, 1);
        if RbarQuery < min(Rbar) || RbarQuery > max(Rbar) % Check if RbarQuery is outside range
            continue
        end
        angle = currentCurve(:, 2);
        anglesAtRbar(end + 1) =  interp1(Rbar, angle, RbarQuery, 'linear');
        AValues(end + 1) = i;
    end

    if numel(anglesAtRbar) < 2 % Check if there are curves to interpolate
        alpha = NaN;
        return;
    end

    alpha = (interp1(AValues, anglesAtRbar, AQuery, 'linear'))/2;
end

function mu = thirtyEightInterp(Chart, RbarQuery, AQuery)

    AValues = [];
    muAtRbar = [];

    for i = 1:5 
        currentCurve = Chart.muRbarA{i};
        Rbar = currentCurve(:, 1);
        if RbarQuery < min(Rbar) || RbarQuery > max(Rbar) 
            continue
        end
        mu = currentCurve(:, 2);
        muAtRbar(end + 1) =  interp1(Rbar, mu, RbarQuery, 'linear');

        if (i == 5)
            AValues(end + 1) = 6;
        else
            AValues(end + 1) = i;
        end
    end

    if numel(muAtRbar) < 2 
        mu = NaN;
        return;
    end

    mu = interp1(AValues, muAtRbar, AQuery, 'linear');
end