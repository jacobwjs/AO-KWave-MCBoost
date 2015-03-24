
clear all;

SAVE_TO_DISK  = true;
VISUALISATION = false;

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
x_axis_num_points = 1024;
y_axis_num_points = 432; 
z_axis_num_points = 432; 
% x_axis_num_points = 1024/2;
% y_axis_num_points = 648/2; 
% z_axis_num_points = 648/2; 

Nx = x_axis_num_points; 
Ny = y_axis_num_points; 
Nz = z_axis_num_points; 

% set desired grid size in the x-direction not including the PML
x = 60e-3;        % [m]

% calculate the spacing between the grid points
dx = (x/Nx);      % [m]
dy = dx;          % [m]
dz = dx;          % [m]

% create the k-space grid
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================
% Density changes as a function of temperature.
rho0_10 = 999.7;  % Density of water at 10C. [kg/m^3]
rho0_20 = 998.2;  % Density of water at 20C.
rho0_30 = 995.7;  % Density of water at 30C.
rho0_40 = 992.2;  % Density of water at 40C.
rho0_50 = 988.1;  % Density of water at 50C.
rho0_60 = 983.2;  % Density of water at 60C.
rho0_70 = 977.8;  % Density of water at 70C.

% SOS changes as a function of temperature. [m/s]
c0_heated_10 = speedSoundWater(10); % SOS of water at 10C. 
c0_heated_20 = speedSoundWater(20); % SOS of water at 20C.
c0_heated_30 = speedSoundWater(30); % SOS of water at 30C.
c0_heated_37 = speedSoundWater(37); % SOS of water at 37C.
c0_heated_40 = speedSoundWater(40); % SOS of water at 40C.
c0_heated_50 = speedSoundWater(50); % SOS of water at 50C.
c0_heated_60 = speedSoundWater(60); % SOS of water at 60C.
c0_heated_70 = speedSoundWater(70); % SOS of water at 70C.

% % Density, attenuation, SOS, etc. of breast tissue.
% % -------------------------------------------------
% % Ref: T. L. Szabo, Diagnostic Ultrasound Imaging (Elsevier, Burlington, 2004), pp. 4?6.
% rho0_breast             = 1020;                 % [kg/m^3]
% c0_breast               = 1510;                 % [m/s]
% alpha_atten_breast      = 0.75;                 % [dB/(MHz^y cm)]
% alpha_power_breast      = 1.5;
% BonA_breast             = 9.63;

% Density, attenuation, SOS, etc. of Agar
% ------------------------------------------------------
rho0_agar           = 1024;
c0_agar             = 1500;
alpha_atten_agar    = 0.7;
alpha_power_agar    = 1.5;
BonA_agar           = 6.0;

% Define the acoustic properties of the medium
% -------------------------------------------------------
c0 = c0_agar;
rho0 = rho0_agar;
% Non-linear effects
% -------------------------------------------------------
medium.alpha_coeff = alpha_atten_agar; 	
medium.alpha_power = alpha_power_agar;
medium.BonA = BonA_agar;
% Acoustically homogeneous medium
% -------------------------------------------------------
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


% =========================================================================
% DEFINE THE HIFU PROBE
% =========================================================================
% The probe ships with the diameter and radius of curvature in the
% documentation. We need to find the height (h) of the curvature from that, and
% translate it to the discrete grid in kWave.
US_freq = 3.5e6;
source_strength = 1.8e6;    	% [Pa] 
num_cycles = 5;

radius_of_curvature = 0.035; % [m]
diameter = 0.023; % [m]
h = radius_of_curvature - sqrt(radius_of_curvature^2 - (diameter/2)^2); % [m]

radius_of_curv_voxel_cnt = round(radius_of_curvature / dy); % [voxels]
height_voxel_cnt         = round(h / dx);                   % [voxels]

[ss,~] = makeSphericalSection(radius_of_curv_voxel_cnt,...  % radius of curvature [voxels]
                              height_voxel_cnt,...               % height [voxels]
                              [],...                             % width [voxels]. Empty so no scaling.
                              false);                            % boolean to plot.

% Assign the HIFU probe to the source mask.
% -------------------------------------------------------
Probe_centered_in_medium_YZaxis = zeros(Nx, Ny, Nz);
Probe_centered_in_medium_YZaxis((PML_X_SIZE + 5): (PML_X_SIZE + 5) + (size(ss,1)-1),...
             Ny/2+1 - round(size(ss,2)/2) : Ny/2-1 + round(size(ss,2)/2),...
             Ny/2+1 - round(size(ss,2)/2) : Ny/2-1 + round(size(ss,2)/2)) = ...
             ss;
         
source.p_mask = Probe_centered_in_medium_YZaxis;                         
      


% =========================================================================
% DEFINE THE SIMULATION RUNTIME PARAMETERS
% =========================================================================
% Simulation runtime
% Only transmitting, so t_end is only based on the time needed to reach the
% bottom of the medium (plus a little more, thus the 1.1 factor).
t_end = (Nx*dx)*1.1/c0;                % [s]
lambda = c0 / US_freq;
% Time it takes to propagate the US 180degrees. For the AO sim to be
% correct 'dt' must evenly divide this number.
pi_phase_shift = lambda/2 * (1/c0);
fprintf('To meet criteria of the medium, max time step allowed is: %.12f [secs]\n', cfl*dx/c0);
% Calculate time step.  Based on the CFL, max SOS and the minimum voxel
% size.
dt = cfl*dx/medium.sound_speed;
%dt = (pi_phase_shift/32)
%dt = 10e-9;
fprintf('Setting time step to: %.12f [secs]\n', dt);
pause(1);
% Calculate the number of steps we must take to allow the ultrasound to
% reach the distance created by t_end/dt.
Nt = floor(t_end/dt);
% Form the time array from the above defined values.
kgrid.t_array = 0:dt:(Nt-1)*dt;


% =========================================================================
% DEFINE THE SOURCE SIGNAL APPLIED TO THE HIFU TRANSDUCER
% =========================================================================
% Create the input signal using a tone burst
tone_burst_freq = US_freq;     % [Hz]
tone_burst_cycles = num_cycles;
input_signal = source_strength*toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
%input_signal = repmat(input_signal, length(find(source.p_mask == 1)), 1);
% Filter the source to remove any high frequencies not supported by the
% grid.
input_signal = filterTimeSeries(kgrid, medium, input_signal);
% Scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
%input_signal = (source_strength./(c0*rho0)).*input_signal;
source.p = input_signal;




% =========================================================================
% DEFINE THE LOCATION OF THE SENSOR MASK
% =========================================================================
%sensor.mask = ones(Nx, Ny, Nz); % The full medium
sensor.mask = zeros(Nx, Ny, Nz);

% Single plane through the medium along the ultrasound axis.
sensor.mask(:, :,Nz/2) = 1; 
% sensor.mask(Nx/4:Nx-Nx/4,...
%             Ny/4:Ny-Ny/4,...
%             Nz/4:Nz-Nz/4) = 1;

%sensor.record = {'p_final', 'p_max'};


% % =========================================================================
% % RUN THE SIMULATION
% % ========================================================================= 
% set the input settings
input_args = {...
    'PMLInside', true, 'PlotSim', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
    'DisplayMask', source.p_mask,...
    %'PlotScale', [-(source_strength+source_strength/2), (source_strength+source_strength/2)],...
     };
 
if (SAVE_TO_DISK)
    filename = ['AO_sim_HIFU_3p50MHz_', num2str(tone_burst_cycles),...
                'cycles_', num2str(source_strength/1e6), 'MPa_INPUT.h5'];
    kspaceFirstOrder3D(kgrid, medium, source, sensor, 'SaveToDisk', filename, input_args{:});
else
    [sensor_data] = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
end


% =========================================================================
% VISUALISATION
% =========================================================================
if (VISUALISATION)
    
    % sensor_data = h5read('filename.h5', '/p');
    time_steps_recorded = size(sensor_data, 2);
    sensor_data = reshape(sensor_data, [Nx, Ny, time_steps_recorded]);
    
    % Compute the amplitude spectrum
    [freq, amp_spect] = spect(sensor_data, 1/kgrid.dt, 'Dim', 3);
    
    % Compute the index at which the source frequency and its harmonics
    % occur
    [f1_value, f1_index] = findClosest(freq, tone_burst_freq);
    [f2_value, f2_index] = findClosest(freq, tone_burst_freq*2);
    [f3_value, f3_index] = findClosest(freq, tone_burst_freq*3);
    
    % Extract the amplitude at the source frequency and store
    beam_pattern_f1 = amp_spect(:,:,f1_index);
    
    % Extract the amplitude at the 2nd harmonic and store
    beam_pattern_f2 = amp_spect(:,:,f2_index);
    
    % Extract the amplitude at the 3rd harmonic and store
    beam_pattern_f3 = amp_spect(:,:,f3_index);
    
    % Extract the integral of the total amplitude spectrum
    beam_pattern_total = sum(amp_spect, 3);
    
    % Plot the fundemental, 2nd harmonic, and 3rd harmonic beam patterns
    figure; imagesc(beam_pattern_f1);
    figure; imagesc(beam_pattern_f2);
    figure; imagesc(beam_pattern_f3);
    
    
    % calculate the amplitude spectrum of the input signal and the signal
    % recorded each of the sensor positions
    [f_input, as_input] = spect([input_signal, zeros(1, 2*length(input_signal))], 1/kgrid.dt);
    
    [f, as_1] = spect(sensor_data(587, 216, :), 1/kgrid.dt);
    [f, as_2] = spect(sensor_data(554, 216, :), 1/kgrid.dt);
    [f, as_3] = spect(sensor_data(531, 216, :), 1/kgrid.dt);
    f = squeeze(f);
    as_1 = squeeze(as_1);
    as_2 = squeeze(as_2);
    as_3 = squeeze(as_3);
    
    
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
    
    % plot the corresponding amplitude spectrum
    figure;
    plot(f/1e6, as_1./max(as_1(:)), 'k-', ...
         f/1e6, as_2./max(as_1(:)), 'b-', ...
         f/1e6, as_3./max(as_1(:)), 'r-');
    legend('Sensor Position 1', 'Sensor Position 2', 'Sensor Position 3');
    xlabel('Frequency [MHz]');
    ylabel('Normalised Amplitude Spectrum [au]');
    f_max = medium.sound_speed / (2*dx);
    set(gca, 'XLim', [0 f_max/1e6]);
end

