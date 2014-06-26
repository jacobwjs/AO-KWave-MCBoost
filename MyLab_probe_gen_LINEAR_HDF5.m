clear all;


PA_GUIDED_FOCUS = false;
MY_FIXED_FOCUS     = true;

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 20;            % [grid points]
PML_Y_SIZE = 20;            % [grid points]
PML_Z_SIZE = 20;            % [grid points]

% set total number of grid points not including the PML

Nx = 1296; % dx*Nx = 39.7 [mm]
Ny = 768;  % dy*Ny = 23.5 [mm]
Nz = 512;  % dz*Nz = 15.7 [mm]

%Nx = 1296; % dx*Nx = 39.7 [mm]
%Ny = 864;  % dy*Ny = 26.5 [mm]
%Nz = 512;  % dz*Nz = 15.7 [mm]



% calculate the spacing between the grid points.
% To better approximate the US transducer simulated, voxel size is based
% upon pitch size of the transducer.
% Definitions to match the SL3323 MyLAB probe
elevation_height = 5e-3;
pitch            = 0.245e-3;
SL3323_active_elements  = 64;
kerf = 1;                           % Assume kerf of 1.

% dx = x/Nx                  % [m]
dx = pitch/8;                % [m]
dy = dx;                     % [m]
dz = dx;                     % [m]

% create the k-space grid
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium

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


% Density, attenuation, SOS, etc. of breast tissue.
% -------------------------------------------------
% Ref: T. L. Szabo, Diagnostic Ultrasound Imaging (Elsevier, Burlington, 2004), pp. 4?6.
rho0_breast             = 1020;                 % [kg/m^3]
c0_breast               = 1510;                 % [m/s]
alpha_atten_breast      = 0.75;                 % [dB/(MHz^y cm)]
alpha_power_breast      = 1.5;
BonA_breast             = 9.63;

% Density, attenuation, SOS, etc. of Agar
% ---------------------------------------
rho0_agar           = 1024;
c0_agar             = 1500;
alpha_atten_agar    = 0.7;
alpha_power_agar    = 1.5;
BonA_agar           = 6.0;


rho0 = 1000;
c0 = 1500;

Nx_tot = Nx;
Ny_tot = Ny;
Nz_tot = Nz;

% define a random distribution of scatterers for the medium
%background_map_mean = 1;
%background_map_std = 0.008;
%background_map = background_map_mean + background_map_std*randn([Nx, Ny, Nz]);
%medium.sound_speed = sound_speed_map(:, :, :);  % [m/s]
%medium.density = density_map(:, :, :);          % [kg/m^3]


medium.sound_speed = c0;  % [m/s]
medium.density     = rho0;  % [kg/m^3]
medium.alpha_coeff = 0.0;   % [dB/(MHz^y cm)] 	
medium.alpha_power = 0.0;
medium.BonA = 0.0;



% =========================================================================
% DEFINE THE SENSOR PARAMETERS
% =========================================================================
% Place sensors over a region of interest in the medium.
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(PML_X_SIZE*2:(Nx-2*PML_X_SIZE),...
            PML_Y_SIZE*2:(Ny-2*PML_Y_SIZE),...
            PML_Z_SIZE*2:(Nz-2*PML_Z_SIZE)) = 1;


% create the time array
% -------------------------------------------------------------------------
%Courant-Friedrichs-Lewy (CFL) stability level (k-Wave default is 0.3) 
cfl = 0.3; 

% Simulation runtime
% Only transmitting, so t_end is only based on the time needed to reach the
% bottom of the medium (plus a little more, thus the 1.3 factor).
t_end = (Nx*dx)*1.1/c0;  

% Calculate time step.  Based on the CFL, max SOS and the minimum voxel
% size.
dt = cfl*dx/c0;
% Setting dt to 10 ns to make things easy in AO_sim, but make sure that it
% can be done.  That is if the CFL safe dt is larger than 10 ns, we can
% reduce it without worrying about accuracy (in fact it becomes more
% accurate).
if (dt > 5e-9)
    display('Setting dt = ');
    dt = 5e-9
else
    display('Cannot set dt=5 ns');
    display('Using default: ')
    dt
end


% Calculate the number of steps we must take to allow the ultrasound to
% reach the distance created by t_end/dt.
Nt = t_end/dt;

% !!!!!!!!!!!!!!!!!!!!!!!!!! FOR DEBUGGING AND TESTING !!!!!!!!!!!!!!!!!!!!
% For testing purposes (i.e. looking at particle displacement) we don't
% need to simulate the full medium, only to the point at which the wave has
% passed the sensor location.  Here we've set the focus to be less than
% half of the available medium, so cut simulation time in half to speed
% this process up.
%Nt = 40;


% Form the time array from the above defined values.
kgrid.t_array = 0:dt:(Nt-1)*dt;

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
% Measured peak-to-peak
source_strength = 1.0e6;    	% [Pa] 
tone_burst_freq = 5.0e6;        % [Hz]
source_freq = tone_burst_freq;
tone_burst_cycles = 5;


% create the input signal using toneBurst
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
%input_signal = sin(2*pi*source_freq*kgrid.t_array);

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity).  This is due
% to the fact that the ultrasound array is directional, so the input
% pressure is converted to a directional velocity.
input_signal = (source_strength./(c0*rho0)).*input_signal;



% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% physical properties of the transducer
transducer.number_elements = SL3323_active_elements;        % total number of transducer elements
transducer.element_width   = round(pitch/dy);               % width of each element [grid points]
transducer.element_length  = round(elevation_height/dz);  	% length of each element [grid points]
transducer.element_spacing = kerf;  	% spacing (kerf  width) between the elements [grid points]
transducer.radius = inf;                % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements*transducer.element_width ...
    + (transducer.number_elements - 1)*transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
% Note the placement below the PML with some cushion.
MIDDLE_of_medium = round([PML_X_SIZE+5,...
                          Ny/2 - transducer_width/2,...
                          Nz/2 - transducer.element_length/2]);
                      
FRONT_EDGE_of_medium_5mmOpticalDepth_0mmShift = round([PML_X_SIZE+5,...
                                                    Ny/2 - transducer_width/2,...
                                                    Nz/3 - transducer.element_length/2]);
% FIXME:
% - Why does the extra 64 pixels need to be included.
azimuthal_4mmShift = round(0.004/dy); 
FRONT_EDGE_of_medium_5mmOpticalDepth_4mmShift = round([PML_X_SIZE+5,...
                                                       Ny/2 - transducer_width/2 - azimuthal_4mmShift,...
                                                       Nz/3 - transducer.element_length/2]);
                                                
% Bead is ~6 mm deep (optical axis), so center the US probe that depth along
% z-axis (taking into account PML)
BEAD_A_centered = round([(PML_X_SIZE+5), ...                                        % x-axis
                         (Ny/2 - transducer_width/2), ...                           % y-axis
                         (PML_Z_SIZE + 6e-3/dz  + transducer.element_length/2)]);   % z-axis

% Assign the transducer position.
transducer.position = MIDDLE_of_medium;


% Define the focal depth
focal_depth = 21.5e-3;


% Create the steering angle for moving 4mm along the azimuthal plane.
horizontal_steering_angle_4mm = 180/pi*atan(4e-3/focal_depth);

% The steering angle along the azimuthal plane for focusing.
steering_angle = 0;

% properties used to derive the beamforming delays
transducer.sound_speed = c0;                % sound speed [m/s]

% Provide delays from an experiment.  Convert back from MyLAB scaled ticks into
% usecs.
if (PA_GUIDED_FOCUS)
    
    info.filename = '/Volumes/TJS CRUZER/Tagging Volume Experiments/20-5-2013/PA_signals_20-5-2013/PA_bead-A_6mmshift.rfe';
    info.middle_channel = 31;
    info.dt = 1/50e6;       % Set MyLAB acquisition time-step (based on 50 MHz sampling).
    delays = crossCorrelateTimeOfArrivals(info, []);
    smooth_delays = smooth(delays, 'rlowess');  % Remove any outliers by smoothing the delay curve.
    
    figure; plot(smooth_delays);
    pause(0.5);
    
    % k-Wave does a lot of manipulation with the beamforming delays it
    % calculates.  In order to minimize modifications to k-Wave, we put the
    % beamforming delays calculated above into a form k-Wave expects.
    smooth_delays = smooth_delays - max(smooth_delays);
    
    % Assign the delays to the transducer object.
    transducer.supplied_transmit_delays = smooth_delays;
elseif (MY_FIXED_FOCUS)
    F = focal_depth;  % Focal length [m]
	c = c0;   % speed-of-sound [m]
	N = SL3323_active_elements;     % Number of elements
	N_avg = (N-1)/2;
	d = pitch;    % The center-to-center spacing of elements on the probe (i.e. the pitch).
	
	%theta = horizontal_steering_angle_4mm*(pi/180);    % Steering angle from center of array.
    theta = steering_angle;
    
	n = 0:1:N-1;

	% Take into account the negative steering angle.  Just want the regular
	% case reflected across the y-axis, so use matlab 'fliplr' with positive
	% angle calculation.
	if (theta >= 0)
    		t_n = F/c * (sqrt(1 + (N_avg*d/F).^2 + ((2*N_avg*d)/F)*sin(theta)) - ...
        	sqrt(1+((n-N_avg)*d/F).^2 - (2*(n-N_avg)*d/F)*sin(theta)));
	else
    		theta = abs(theta);
    		t_n = F/c * (sqrt(1 + (N_avg*d/F).^2 + ((2*N_avg*d)/F)*sin(theta)) - ...
        	sqrt(1+((n-N_avg)*d/F).^2 - (2*(n-N_avg)*d/F)*sin(theta)));
    		% Reflect delays across the y-axis.
    		t_n = fliplr(t_n);
	end
	
	% k-Wave does a lot of manipulation with the beamforming delays it
    	% calculates.  In order to minimize modifications to k-Wave, we put the
    	% beamforming delays calculated above into a form k-Wave expects.
   	FF_delays = t_n - max(t_n);
   	
   	% Assign the delays to the transducer object.
   	transducer.supplied_transmit_delays = FF_delays;
else
    transducer.focus_distance = focal_depth;          % focus distance [m]
    %transducer.steering_angle = horizontal_steering_angle_4mm; % steering angle [degrees]
    transducer.steering_angle = steering_angle; % [degrees]
end
transducer.elevation_focus_distance = 21.5e-3;


% apodization
transducer.transmit_apodization = 'Rectangular';
transducer.receive_apodization = 'hanning';

% define the transducer elements that are currently active
number_active_elements = transducer.number_elements;
transducer.active_elements = ones(transducer.number_elements, 1);

% append input signal used to drive the transducer
transducer.input_signal = input_signal;

% create the transducer using the defined settings
transducer = makeTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties
% transducer.plot




 
% % =========================================================================
% % RUN THE SIMULATION
% % =========================================================================
% 
% 
% % Set whether or not to return velocity fields.
info.return_velocity = false;
% 
% set the input settings
input_args = {...
    'PMLInside', true, 'PlotSim', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
    %'DisplayMask', transducer.mask,...
    %'PlotScale', [-(source_strength+source_strength/2), (source_strength+source_strength/2)],...
     };

 
% save input files to disk
if (PA_GUIDED_FOCUS)
    PA_file = strsplit(info.filename, '/');             % Split up string
    PA_file = char(PA_file(1,end));                     % Convert last cell to char array
    filename = ['MyLAB_', PA_file(1:end-4),
                num2str(source_strength), 'Pa', '_',...
                num2str(tone_burst_freq), 'Hz','LINEAR_INPUT.h5'];     % Form new file name
else

    filename = ['MyLAB_FF', num2str(transducer.focus_distance), '_',...
                num2str(source_strength), 'Pa', '_',...
                num2str(tone_burst_freq), 'Hz', 'LINEAR_INPUT.h5'];
end
kspaceFirstOrder3D(kgrid, medium, transducer, sensor, 'SaveToDisk', filename, input_args{:});
