% Creates phasor diagram
% ----------------------------------

% Example of taking AO_sim exit data and forming electric field. Need to
% get here before phasor plot can take place.
% --------------------------------------------------------------
% tX = dlmread('trans_tX_exit_data.dat');
% phasor_tX = sqrt(tX(:,1).*exp(1i*k*tX(:,4));  % combined OPL is position 4.


DISPLACEMENT = true;
REFRACTIVE   = true;
COMBINED     = true;

k = 2*pi/532e-9;

extension = '.dat';

% we want all files in the current directory.
extension = ['*', extension];

%list all the files in current folder
filelist = dir(extension);

% Perform a natural sort on the filelist names since the threaded nature of
% the AO_simulation might write data at different times, but the timestamp
% creation of the file is in the name, which is what should be sorted by.
[files, ~] = sort_nat({filelist.name});

% How many files are there that matched the extension in this directory.
num_files = length(filelist(:));

for i=1:num_files
    % Pre-allocate the array that holds the image data.
    temp = dlmread(char(files(i)));
    
    % Displacement
    if (DISPLACEMENT)
        exit_data(i).phasor_D = sqrt(temp(:,1)).*exp(1i*k*temp(:,2));
    end
    
    % Refractive index
    if (REFRACTIVE)
        exit_data(i).phasor_N = sqrt(temp(:,1)).*exp(1i*k*temp(:,3));
    end
    
    % Combination (displacement + refractive)
    if (COMBINED)
        exit_data(i).phasor_C = sqrt(temp(:,1)).*exp(1i*k*temp(:,4));   
    end
end

% Find which paths in the exit data file that were tagged by looking at
% differences in OPLs between the modulated and unmodulated case. Note that
% it's assumed the first frame is the static (i.e. unperturbed medium, no US)
% case.
unmodulated = dlmread(char(files(1)));
modulated   = dlmread(char(files(2)));
% Find the indices where the OPLs have changed between the two files.
modulation_indices_Dmap = find(unmodulated(:,2) ~= modulated(:,2));
modulation_indices_Nmap = find(unmodulated(:,3) ~= modulated(:,3));
modulation_indices_Cmap = find(unmodulated(:,4) ~= modulated(:,4));
% Assume the combination is turned on, which means the indices of tagging
% are the same for all mechanisms.
tagged_index = modulation_indices_Cmap;

% All phasors begin at the origin.
%origin = zeros(num_files, 2);

% Which photon path are we interested in.
path = 1944;

% The first file in the exit-data for the 'AO_sim_loadData' case is the
% unmodulated paths. We grab those for reference.
if (DISPLACEMENT)
    Pd = [real(exit_data(1).phasor_D(tagged_index(path))) ...
          imag(exit_data(1).phasor_D(tagged_index(path)))];
end
if (REFRACTIVE)
    Pn = [real(exit_data(1).phasor_N(tagged_index(path))) ...
          imag(exit_data(1).phasor_N(tagged_index(path)))];
end
if (COMBINED)
    Pc = [real(exit_data(1).phasor_C(tagged_index(path))) ...
          imag(exit_data(1).phasor_C(tagged_index(path)))];
end
 


% Build up the coordinates for the arrows in the phasor diagram.
for i=2:num_files
    if (DISPLACEMENT)
        Pd = [Pd; real(exit_data(i).phasor_D(tagged_index(path))) ...
                  imag(exit_data(i).phasor_D(tagged_index(path)))];
    end
    if (REFRACTIVE)
        Pn = [Pn; real(exit_data(i).phasor_N(tagged_index(path))) ...
                  imag(exit_data(i).phasor_N(tagged_index(path)))];
    end
    if (COMBINED)
        Pc = [Pc; real(exit_data(i).phasor_C(tagged_index(path))) ...
                  imag(exit_data(i).phasor_C(tagged_index(path)))];
    end
end



if (DISPLACEMENT)
    figure;
    
    hs_polar = subplot(1,2,2);
    polar(0,1),
    hold on;
    
    hs_complex = subplot(1,2,1);
    hold on;
    line([0 0], [-1 1], 'LineStyle', '--', 'Color', 'k');
    line([-1 1], [0 0], 'LineStyle', '--', 'Color', 'k');
    axis([-1, 1, -1, 1])
    
    subplot(hs_polar);
    arrow3([0 0], Pd(1,:), 'r');
    
    subplot(hs_complex);
    arrow3([0 0], Pd(1,:), 'r');
    for i = 2:num_files
        pause(0.25);
        
        subplot(hs_polar);
        arrow3([0 0], Pd(i,:), 'colors');
        
        subplot(hs_complex);
        arrow3([0 0], Pd(i,:), 'colors');
    end
    hold off;
end

if (REFRACTIVE)
    figure;
    
    hs_polar = subplot(1,2,2);
    polar(0,1),
    hold on;
    
    hs_complex = subplot(1,2,1);
    hold on;
    line([0 0], [-1 1], 'LineStyle', '--', 'Color', 'k');
    line([-1 1], [0 0], 'LineStyle', '--', 'Color', 'k');
    axis([-1, 1, -1, 1])
    
    subplot(hs_polar);
    arrow3([0 0], Pn(1,:), 'r');
    
    subplot(hs_complex);
    arrow3([0 0], Pn(1,:), 'r');
    for i = 2:num_files
        pause(0.25);
        
        subplot(hs_polar);
        arrow3([0 0], Pn(i,:), 'colors');
        
        subplot(hs_complex);
        arrow3([0 0], Pn(i,:), 'colors');
    end
    hold off;
end

if (COMBINED)
    figure;
    
    hs_polar = subplot(1,2,2);
    polar(0,1),
    hold on;
    
    hs_complex = subplot(1,2,1);
    hold on;
    line([0 0], [-1 1], 'LineStyle', '--', 'Color', 'k');
    line([-1 1], [0 0], 'LineStyle', '--', 'Color', 'k');
    axis([-1, 1, -1, 1])
    
    subplot(hs_polar);
    arrow3([0 0], Pc(1,:), 'r');
    
    subplot(hs_complex);
    arrow3([0 0], Pc(1,:), 'r');
    for i = 2:num_files
        pause(0.25);
        
        subplot(hs_polar);
        arrow3([0 0], Pc(i,:), 'colors');
        
        subplot(hs_complex);
        arrow3([0 0], Pc(i,:), 'colors');
    end
    hold off;
end




p0 = mod(angle(exit_data(1).phasor_C(tagged_index(path)))*180/pi+360, 360);
p1 = mod(angle(exit_data(8).phasor_C(tagged_index(path)))*180/pi+360, 360);

min(abs(p1-p0))
max(abs(p1-p0))


% origin = zeros(3,2);
% P = [real(phasor_t1(tagged_index(1))) imag(phasor_t1(tagged_index(1)))];
% P = [P; real(phasor_t2(tagged_index(1))) imag(phasor_t2(tagged_index(1)))];
% P = [P; real(phasor_t13(tagged_index(1))) imag(phasor_t13(tagged_index(1)))]
% 
% figure;
% arrow(origin, P);

