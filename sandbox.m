%**************************************************************************
%
% Angular Spectrum Field Calculator
%
%   This script is for testing effects of modifications of the angular
%   spectrum method.
%
%
%
%                 Scott Schoen Jr | Georgia Tech | 2016
%
%**************************************************************************

clear all;
close all;
clc;

%% This block defines sensor positions

zDistance = 20E-2; % Axial position of sensor array [m]
dx = 5E-3;         % Linear "size" of sensor (same for all dim.) [m]

% Set linear array dimensions
Dx = 50E-3; % [m]
Dy = 0; % [m]

% Spatial extent to 0-pad (on each side!)
xPaddingDistance = Dx./2;  % [m] 
yPaddingDistance = Dy./2;  % [m] 

% Calculate total number of array points (sensors), 
% ---including padding points
numSpacePointsX = round( ...
    (Dx + 2.*xPaddingDistance)./dx ...
    ); 
numSpacePointsY = round( ...
    (Dy + 2.*yPaddingDistance)./dx ...
    ); % Total number of array points (sensors), including padding points

% Specify linear array if needed
numSpacePointsX = max( numSpacePointsX, 1);
numSpacePointsY = max( numSpacePointsY, 1);

% Bubble position
x0 = 0;
y0 = 0;
z0 = 0;

%% This block computes the bubble wall velocity

% Start very basically with a bubble excited by a speficied pulse
addpath( '../bubble/' ); % Add path with bubble functions

% Set medium properties
medium.p0 = 1.0e5;
medium.c0 = 1485;
medium.rho = 998;
medium.k = 1.33;
medium.sigma = 0.0725;  
medium.mu = 0.001;

% Set bubble properties
bubble.R0 = 1.0E-6;   % Equilibrium radius [m]
bubble.dR0 = 0;       % Initial wall velocity [m/s]
bubble.Pvap = 2.33E3; % Vapor Pressure [Pa]
bubble.hasShell = 0;

% Set simulation properties
tMin = 0;
tMax = 5E-4;

% Natural frequency for Rayleigh-Plesset
f0_rp = 1./(2*pi*bubble.R0).*sqrt( ...
    (1./medium.rho).*(3.*medium.k*(medium.p0 ) ...
    + 2*medium.sigma./bubble.R0.*( 3*medium.k - 1) ) ...
    );

% Excitation function properties
pAmp = 0.1.*101E3;   % [Pa]
f0 = 0.5.*f0_rp;        % [Hz]
dt = 1./( 15.*f0 ); % Time step [s] (above Nyquist of f0 at least)
numTimePoints = round((tMax - tMin)./dt); % Number of time points
omega0 = 2.*pi.*f0;
tSim = linspace( tMin, tMax, numTimePoints);
t0 = 10E-5; % [s] Make sure not too large that ODE solver misses excitation
BW = 0.3;  % Fractional bandwidth

% Create signal
normSignal = excitationPulse( tSim, f0, BW, t0, 0 );
excitation.signal = pAmp.*normSignal;
excitation.tVector = tSim;

% Solve Rayleigh-Plesset Equation
tSpan = [ tMin, tMax ];
initialConditions = [1*bubble.R0; bubble.dR0]';

[timeVector_sol, solution] = ode15s(  ...
    @(t, y) RPEqn(t, y, medium, bubble, excitation), ...
    tSpan, initialConditions );

% Get solution vectors of interets
R_sol = solution(:, 1);    % [m]
Rdot_sol = solution(:, 2); % [m/s]

% Interpolate the solution vectors to have uniform spacing in time
timeVector = tSim; % The originally specified time vector
R = interp1( timeVector_sol, R_sol, timeVector );
Rdot = interp1( timeVector_sol, Rdot_sol, timeVector );

% Plot bubble pulsation for reference
figure();
subplot( 2, 1, 1)
plot( excitation.tVector.*1E6, excitation.signal, 'k' );
ylabel('$s(t)$ [Pa]' );
zoom xon;
subplot( 2, 1, 2)
plot( timeVector.*1E6, R./bubble.R0, 'k' );
xlabel('Time [$\mu$s]' );
ylabel('$R/R_{0}$ [m]' );
zoom xon;



%% This block computes the pressure at each sensor

% Create sensor array such that:
%   sensorData( 1, 1, : ) - Time series recorded by sensor (1, 1)
%   sensorData( :, :, 1 ) - Pressure field at t = t(1)
sensorData = zeros( numSpacePointsX, numSpacePointsY, numTimePoints );

% Assign x- and y- positions of the sensors
xStart = -Dx./2 - xPaddingDistance;
xEnd = Dx./2 + xPaddingDistance;
xPositions = linspace( xStart, xEnd, numSpacePointsX );

yStart = -Dy./2 - yPaddingDistance;
yEnd = Dy./2 + yPaddingDistance;
yPositions = linspace( yStart, yEnd, numSpacePointsY );

% Create frequency vector
dt = timeVector(2) - timeVector(1); % Time step [s]
Fs = 1./dt; % Sampling frequency [Hz]
numFreqPoints = numTimePoints; % For readability
fVector = linspace( 0, Fs, numFreqPoints );

% Now get FT of bubble wall velocity
RdotTilde = fft(Rdot);

% Now for each frequency, compute the contribution at the sensors
rho0 = medium.rho;
c0 = medium.c0;
zs = zDistance;
for fCount = 1:numFreqPoints
    
    % Get current wavenumber and frequency
    omega = 2.*pi.*fVector( fCount );
    k0 = omega./c0;
    
    % Get the velocity time series for this frequency
    u0 = RdotTilde(fCount).*exp( 1j.*omega.*timeVector ); 
    
    % Compute received pressure at each sensor from this frequency
    for sensorCountX = 1:numSpacePointsX
        
        % Get this sensor's position and skip if we're in the padding
        % region
        xs = xPositions( sensorCountX );
       
        for sensorCountY = 1:numSpacePointsY
            
            % Get this sensor's position and skip if we're in the padding
            % region
            ys = yPositions( sensorCountY );
            
            % Get the distance from the source to the sensor
            r = sqrt( (x0 - xs).^(2) + (x0 - xs).^(2) + (z0 - zs).^(2) );
            
            % Define the radiation impedance
            Z = rho0.*c0./( 1 + 1./(1j.*k0.*r) );
            
            % Add in the contribution from this freqyency to the received 
            % pressure at the sensor
            currentTimeSeries = squeeze( ...
                sensorData( sensorCountX, sensorCountY, : ) )';
            freqContribution = real( Z.*u0.*exp(-1j.*k0.*r) );
            sensorData( sensorCountX, sensorCountY, : ) = ...
                currentTimeSeries + freqContribution;
                
        end
    end
end

%% Plotting

% Set x- or y-slice positions to plot
xSlicePos = 0; % [m]
ySlicePos = 0; % [m]

% Get normalized data
sensorDataNorm = sensorData./(max(max(max(abs(sensorData)))) );

% FFT of (normalized) received data
sensorDataNormTilde = fft( sensorDataNorm, [], 3 );

% Get meshgrid for plotting
[tPlotVectorX, xPosPlotVectorTime] = meshgrid( timeVector, xPositions);
[fPlotVectorX, xPosPlotVectorFreq] = meshgrid( fVector, xPositions);
[tPlotVectorY, yPosPlotVectorTime] = meshgrid( timeVector, yPositions);
[fPlotVectorY, yPosPlotVectorFreq] = meshgrid( fVector, yPositions);

% 2D Colorplot

% Get matrix ( x-receiver by time ) to plot for that y value
yIndex = find( yPositions >= ySlicePos, 1 );
xPlotData = squeeze( sensorDataNorm( :, yIndex, : ) );
xPlotDataTilde = abs(squeeze( sensorDataNormTilde( :, yIndex, : ) ));

figure()
subplot( 2, 1, 1 );
hold all;
pcolor( tPlotVectorX.*1E6, xPosPlotVectorTime.*1E3, xPlotData );
shading interp
set( gca, 'clim', [-1, 1] );
ylabel('Distance [mm]');
xlabel('Time [$\mu$s]');

subplot( 2, 1, 2 );
hold all;
pcolor( fPlotVectorX.*1E-6, xPosPlotVectorFreq.*1E3, xPlotDataTilde );
shading interp
ylabel('Distance [mm]');
xlabel('Frequency [MHz]');
xlim( [0, (Fs.*1E-6)/2] );
set( gca, 'clim', [0, 8] );

zoom xon;

