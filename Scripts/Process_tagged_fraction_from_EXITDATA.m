% This file loads and processes exit data from an AO simulation to form the
% tagged fraction of light.

function [tagged_frac] = Process_tagged_fraction_from_EXITDATA(extension, scaling)


% Wavenumber (assuming 532nm wavelength)
k = 2*pi/532e-9;

% we want all files in the current directory that match the extension (e.g. '.dat').
extension = ['*', extension];

%list all the files in current folder
filelist = dir(extension);

% Perform a natural sort on the filelist names since the threaded nature of
% the AO_simulation might write data at different times, but the timestamp
% creation of the file is in the name, which is what should be sorted by.
[files, ~] = sort_nat({filelist.name});

% How many files are there that matched the extension in this directory.
num_files = length(filelist(:));

% First file is the unmodulated (i.e. no ultrasound), while the second to
% the last all undergo modulation.
Ea = dlmread(char(files(2)));
tagged_frac = [];
for i=2:num_files
    
    Eb = dlmread(char(files(i)));
    
    % Get the tagged fraction of displacement and refractive index
    % contributions.
    alpha_displaced  = mean(exp(1i*k*scaling*(Eb(:,2)-Ea(:,2))));
    alpha_refractive = mean(exp(1i*k*scaling*(Eb(:,3)-Ea(:,3)))); 
    alpha_combined   = mean(exp(1i*k*scaling*(Eb(:,4)-Ea(:,4))));
    
    tagged_frac.displaced(i-1)  = (1-(abs(alpha_displaced)^2))/2;
    tagged_frac.refractive(i-1) = (1-(abs(alpha_refractive)^2))/2;
    tagged_frac.combined(i-1)   = (1-(abs(alpha_combined)^2))/2;
end

figure; hold on;
plot(tagged_frac.displaced, '-b');
plot(tagged_frac.refractive, '-g');
plot(tagged_frac.combined, '-k');
hold off;

end % end file