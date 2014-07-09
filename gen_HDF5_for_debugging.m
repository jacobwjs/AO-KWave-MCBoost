
% Defining An Ultrasound Transducer Example
%
% In principle, simulations using ultrasound transducers can be run
% following the examples given under Time Varying Source Problems. However,
% assigning the grid points that belong to each transducer element, and
% then assigning the correctly delayed input signals to each point of each
% element can be a laborious task. For this purpose, a special input object
% created using makeTransducer can be substituted for either the source or
% sensor inputs (or both). This example illustrates how this object is
% created and can be used to simulate the field produced by an ultrasound
% transducer.
%
% Note, transducer inputs can only be used in 3D simulations and thus these
% examples are inherently memory and CPU intensive. Whilst the grid sizes
% and source frequencies used in the examples have been scaled down for the
% purpose of demonstrating the capabilities of the toolbox (the inputs do
% not necessarily represent realistic ultrasound settings), they still
% require a comparatively large amount of computational resources. To
% reduce this load, it is advised to run the simulations in single
% precision by setting the optional input 'DataCast' to 'single'.
% Similarly, if you have access to a recent model GPU and a GPU toolbox
% (for example, the MATLAB Parallel Computing Toolbox R2012a or later, or
% Accelereyes Jacket), the simulation times can be significantly reduced by
% setting 'DataCast' to 'gpuArray-single' or 'gsingle'. Alternatively, the
% simulations can be run using the optimised C++ code. See the k-Wave
% Manual for more information.   
%
% The creation of a kWaveTransducer object will only work in versions of
% MATLAB recent enough to support user defined classes. 
%
% author: Bradley Treeby
% date: 20th July 2011
% last update: 25th September 2012
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2012 Bradley Treeby and Ben Cox

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

clear all;


% Decide to run the simulation in matlab or save it to an h5 file for running 
% in the AO_sim for debugging.
SAVE_TO_DISK = false;
PLANAR_WAVE = true;

% simulation settings
DATA_CAST = 'single';

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 10;            % [grid points]
PML_Y_SIZE = 10;            % [grid points]
PML_Z_SIZE = 10;            % [grid points]

% set total number of grid points not including the PML
x_axis_num_points = 128;
y_axis_num_points = 128;
z_axis_num_points = 128;
Nx = x_axis_num_points; 
%Nx = (x_axis_num_points - 2*PML_X_SIZE);    % [grid points]
Ny = y_axis_num_points; 
%Ny = (y_axis_num_points - 2*PML_Y_SIZE);    % [grid points]
Nz = z_axis_num_points; 
%Nz = (z_axis_num_points - 2*PML_Z_SIZE);     % [grid points]

% set desired grid size in the x-direction not including the PML
x = 40e-3;                  % [m]

% calculate the spacing between the grid points
dx = x/Nx;                  % [m]
dx = dx/2;
dy = dx;                    % [m]
dz = dx;                    % [m]

% create the k-space grid
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================
c0 = 1500;
rho0 = 1000;
% define the properties of the propagation medium
%medium.sound_speed = 1500;      % [m/s]
%medium.density = 1000;          % [kg/m^3]
%medium.alpha_coeff = 0.0;      % [dB/(MHz^y cm)]
%medium.alpha_power = 0.0;
%medium.BonA = 0;



% Acoustically homogeneous medium
medium.sound_speed = c0;
medium.density     = rho0;

% Acoustically heterogeneous medium
% -------------------------------------------------------
% define a random distribution of ultrasound scatterers for the medium
% background_map_mean = 1;
% background_map_std = 0.008;
% background_map = background_map_mean + background_map_std*randn([Nx, Ny, Nz]);
% 
% % % define properties
% sound_speed_map = c0*ones(Nx, Ny, Nz).*background_map;
% density_map     = rho0*ones(Nx, Ny, Nz).*background_map;
% medium.sound_speed = sound_speed_map(:, :, :);  % [m/s]
% medium.density = density_map(:, :, :);          % [kg/m^3]

%Courant-Friedrichs-Lewy (CFL) stability level (k-Wave default is 0.3) 
cfl = 0.3; 

% Simulation runtime
% Only transmitting, so t_end is only gtbased on the time needed to reach the
% bottom of the medium (plus a little more, thus the 1.1 factor).
% Subtraction of 100 is just to reduce computation time. We don't need the
% full medium simulated when testing/debugging.
t_end = (Nx*dx)*1.1/c0;                % [s]
% Calculate time step.  Based on the CFL, max SOS and the minimum voxel
% size.
%dt = cfl*dx/medium.sound_speed;
US_freq = 1.0e6;
lambda = c0 / US_freq;
% Time it takes to propagate the US 180degrees. For the AO sim to be
% correct 'dt' must evenly divide this number.
pi_phase_shift = lambda/2 * (1/c0);
display('To meet criteria of the medium, max time step allowed is: ');
cfl*dx/c0
display('Setting time step to: ');
dt = (pi_phase_shift/16)
pause(2);
% Calculate the number of steps we must take to allow the ultrasound to
% reach the distance created by t_end/dt.
Nt = floor(t_end/dt);
% Form the time array from the above defined values.
kgrid.t_array = 0:dt:(Nt-1)*dt;


% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
%
% source_strength = 0.889825*1e6;   % [Pa] 11 cycles
% num_cycles = 11;    

% source_strength = 0.8612*1e6;  % [Pa] 10cycles
% num_cycles = 10;

% source_strength = 0.9038*1e6; % [Pa] 9cycles
% num_cycles = 9;

% source_strength = 0.8612*1e6; % [Pa] 8cycles
% num_cycles = 8;

% source_strength = 0.932*1e6;  % [Pa] 7cycles
% num_cycles = 7;

% source_strength = 0.8619*1e6; % [Pa] 6cycles
% num_cycles = 6;

source_strength = 1.0052*1e6;      % [Pa] 5cycles
num_cycles = 5;

tone_burst_freq = US_freq;     % [Hz]
tone_burst_cycles = num_cycles;

% create the input signal using toneBurst 
if (PLANAR_WAVE)
    tone_burst_cycles = 100;
    input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles, 'Envelope', 'Rectangular');
else
    input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
end

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength./(c0*rho0)).*input_signal;

% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% physical properties of the transducer
transducer.number_elements = 64;    % total number of transducer elements
transducer.element_width = 1;       % width of each element [grid points]
if (PLANAR_WAVE)
    transducer.element_length = 64;     % length of each element [grid points]
else
    transducer.element_length = 12;
end
transducer.element_spacing = 0;     % spacing (kerf width) between the elements [grid points]
transducer.radius = inf;            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements*transducer.element_width ...
    + (transducer.number_elements - 1)*transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
x_offset = PML_X_SIZE+5;
transducer.position = round([x_offset, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = 1500;              % sound speed [m/s]
if (PLANAR_WAVE)
    transducer.focus_distance = inf;
    transducer.elevation_focus_distance = inf;
else
    transducer.focus_distance = 20e-3;          % focus distance [m]
    transducer.elevation_focus_distance = 19e-3;% focus distance in the elevation plane [m]
end
transducer.steering_angle = 0;              % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Rectangular';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
transducer.active_elements = zeros(transducer.number_elements, 1);
transducer.active_elements(1:end) = 1;

% append input signal used to drive the transducer
transducer.input_signal = input_signal;

% create the transducer using the defined settings
transducer = makeTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties;

% =========================================================================
% DEFINE SENSOR MASK
% =========================================================================
% Define the region that data is recorded over the medium.
sensor_Nx = PML_X_SIZE*3:(Nx - PML_X_SIZE*3);
sensor_Ny = PML_Y_SIZE*3:(Ny - PML_Y_SIZE*3);
%sensor_Nz = PML_Z_SIZE*3:(Nz - PML_Z_SIZE*3);
sensor_Nz = Nz/2;
sensor_dims = [size(sensor_Nx,2), size(sensor_Ny,2), size(sensor_Nz,2)];
reshape_dims = sensor_dims(sensor_dims ~= 1);

% create a binary sensor mask with four detection positions
if (SAVE_TO_DISK)
    % Create sensor map of entire medium for use with AO simulation.
    sensor.mask = ones(Nx, Ny, Nz);
    
%     sensor.mask = zeros(Nx, Ny, Nz);
%     sensor.mask(sensor_Nx,...
%                 sensor_Ny,...
%                 sensor_Nz) = 1;
else
    % Otherwise only interested in a small portion.
    sensor.mask = zeros(Nx, Ny, Nz);
    sensor.mask(sensor_Nx,...
                sensor_Ny,...
                sensor_Nz) = 1;
            
    % Define what to save over the sensor mask.
    sensor.record = {'p_max', 'p', 'I'};
end
% =========================================================================
% RUN THE SIMULATION
% =========================================================================

% set the input settings
if (SAVE_TO_DISK)
    input_args = {...
    'PMLInside', true, 'PlotSim', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
    %'DisplayMask', transducer.mask,...
    %'PlotScale', [-(source_strength+source_strength/2), (source_strength+source_strength/2)],...
     };
else
    input_args = {'DisplayMask', transducer.all_elements_mask | sensor.mask, ...
        'PMLInside', true, 'PlotPML', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
        'DataCast', DATA_CAST, 'PlotScale', [-source_strength/2, source_strength/2],...
        'PerfectPlanar', 'x-axis'};
end

% run the simulation
if (SAVE_TO_DISK)
    filename = ['AO_sim_Debug_', num2str(tone_burst_cycles), 'cycles_INPUT.h5'];
    kspaceFirstOrder3D(kgrid, medium, transducer, sensor, 'SaveToDisk', filename, input_args{:});
else
    [sensor_data] = kspaceFirstOrder3D(kgrid, medium, transducer, sensor, input_args{:});
    % calculate the amplitude spectrum of the input signal and the signal
    % recorded each of the sensor positions
    [f_input, as_input] = spect([input_signal, zeros(1, 2*length(input_signal))], 1/kgrid.dt);
    [f, as_1] = spect(sensor_data(1, :), 1/kgrid.dt);
    [f, as_2] = spect(sensor_data(2, :), 1/kgrid.dt);
    [f, as_3] = spect(sensor_data(3, :), 1/kgrid.dt);
end


% =========================================================================
% VISUALISATION
% =========================================================================
if (~SAVE_TO_DISK)
    % plot the input signal and its frequency spectrum
    figure;
    subplot(2, 1, 1), plot((0:kgrid.dt:(length(input_signal)-1)*kgrid.dt)*1e6, input_signal, 'k-');
    xlabel('Time [\mus]');
    ylabel('Input Particle Velocity [m/s]');
    subplot(2, 1, 2), plot(f_input/1e6, as_input./max(as_input(:)), 'k-');
    hold on;
    line([tone_burst_freq, tone_burst_freq]/1e6, [0 1], 'Color', 'k', 'LineStyle', '--');
    xlabel('Frequency [MHz]');
    ylabel('Relative Amplitude Spectrum [au]');
    f_max = medium.sound_speed / (2*dx);
    set(gca, 'XLim', [0 f_max/1e6]);
    
    % plot the recorded time series
    figure;
    stackedPlot(kgrid.t_array*1e6, {'Sensor Position 1', 'Sensor Position 2', 'Sensor Position 3'}, sensor_data);
    xlabel('Time [\mus]');
    
    % plot the corresponding amplitude spectrums
    figure;
    plot(f/1e6, as_1./max(as_1(:)), 'k-', f/1e6, as_2./max(as_1(:)), 'b-', f/1e6, as_3./max(as_1(:)), 'r-');
    legend('Sensor Position 1', 'Sensor Position 2', 'Sensor Position 3');
    xlabel('Frequency [MHz]');
    ylabel('Normalised Amplitude Spectrum [au]');
    f_max = medium.sound_speed / (2*dx);
    set(gca, 'XLim', [0 f_max/1e6]);
end

