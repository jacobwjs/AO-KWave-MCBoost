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
MOVIE        = false;

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
        %exit_data(i).phasor_D = 1.*exp(1i*k*temp(:,2));
    end
    
    % Refractive index
    if (REFRACTIVE)
        exit_data(i).phasor_N = sqrt(temp(:,1)).*exp(1i*k*temp(:,3));
        %exit_data(i).phasor_N = 1.*exp(1i*k*temp(:,3));
    end
    
    % Combination (displacement + refractive)
    if (COMBINED)
        exit_data(i).phasor_C = sqrt(temp(:,1)).*exp(1i*k*temp(:,4));  
        %exit_data(i).phasor_C = 1.*exp(1i*k*temp(:,4));
    end
end

% % Find which paths in the exit data file that were tagged by looking at
% % differences in OPLs between the modulated and unmodulated case. Note that
% % it's assumed the first frame is the static (i.e. unperturbed medium, no US)
% % case.
% unmodulated = dlmread(char(files(1)));
% modulated   = dlmread(char(files(2)));
% % Find the indices where the OPLs have changed between the two files.
% modulation_indices_Dmap = find(unmodulated(:,2) ~= modulated(:,2));
% modulation_indices_Nmap = find(unmodulated(:,3) ~= modulated(:,3));
% modulation_indices_Cmap = find(unmodulated(:,4) ~= modulated(:,4));
% % Assume the combination is turned on, which means the indices of tagging
% % are the same for all mechanisms.
% tagged_index = modulation_indices_Cmap;

% All phasors begin at the origin.
%origin = zeros(num_files, 2);

% % Which photon path are we interested in.
% a = dlmread('trans-251_11_38_50_exit_data.dat');
% b = dlmread('trans-251_11_39_6_exit_data.dat');
% find(abs(b(:,4)-a(:,4)) == max(abs(b(:,4)-a(:,4)));
% i = find((b(:,4)-a(:,4)).^2 == min((b(:,4)-a(:,4)).^2)); %find zeros
% t = 1:1:150e3;
% t(i) = []; % remove zeros
% find(abs(b(t,4)-a(t,4)) == min(abs(b(t,4)-a(t,4))))
max_path_1p25mm = 85976;
min_path_1p25mm = 607;

max_path_5p00mm = 25789;
min_path_5p00mm = 53394;

% This is a common modulated path (i.e. it is present in small spheres and
% larger).
path = 1292;
%path = max_path_1p25mm;

% The first file in the exit-data for the 'AO_sim_loadData' case is the
% unmodulated paths. We grab those for reference.
if (DISPLACEMENT)
    Pd = [real(exit_data(1).phasor_D(path))...
          imag(exit_data(1).phasor_D(path))];
end
if (REFRACTIVE)
    Pn = [real(exit_data(1).phasor_N(path)) ...
          imag(exit_data(1).phasor_N(path))];
end
if (COMBINED)
    Pc = [real(exit_data(1).phasor_C(path)) ...
          imag(exit_data(1).phasor_C(path))];
end
 


% Build up the coordinates for the arrows in the phasor diagram for a single path.
for i=2:num_files
    if (DISPLACEMENT)
        Pd = [Pd; real(exit_data(i).phasor_D(path)) ...
                  imag(exit_data(i).phasor_D(path))];
    end
    if (REFRACTIVE)
        Pn = [Pn; real(exit_data(i).phasor_N(path)) ...
                  imag(exit_data(i).phasor_N(path))];
    end
    if (COMBINED)
        Pc = [Pc; real(exit_data(i).phasor_C(path)) ...
                  imag(exit_data(i).phasor_C(path))];
    end
end

% Build up the phasor for the sum of all modulated and unmodulated light
% paths.
sum_Pd = [];
sum_Pn = [];
sum_Pc = [];
% Add in the unmodulated case (i.e. no ultrasound in the medium).
sum_Pd = sum(exit_data(1).phasor_D);
sum_Pd = sum_Pd / abs(sum_Pd); % normalize the phasor to have unity magnitude.
sum_Pd = [real(sum_Pd) imag(sum_Pd)];
    
sum_Pn = sum(exit_data(1).phasor_N);
sum_Pn = sum_Pn / abs(sum_Pn); 
sum_Pn = [real(sum_Pn) imag(sum_Pn)];

sum_Pc = sum(exit_data(1).phasor_C);
sum_Pc = sum_Pc / abs(sum_Pc);
sum_Pc = [real(sum_Pc) imag(sum_Pc)];
for i=2:num_files
    if (DISPLACEMENT)
        temp = sum(exit_data(i).phasor_D);
        temp = temp / abs(temp); % normalize the phasor to have unity magnitude.
        sum_Pd = [sum_Pd; real(temp) imag(temp)];
    end
    if (REFRACTIVE)
        temp = sum(exit_data(i).phasor_N);
        temp = temp / abs(temp); % normalize the phasor to have unity magnitude.
        sum_Pn = [sum_Pn; real(temp) imag(temp)];
    end
    if (COMBINED)
        temp = sum(exit_data(i).phasor_C);
        temp = temp / abs(temp); % normalize the phasor to have unity magnitude.
        sum_Pc = [sum_Pc; real(temp) imag(temp)];
    end
end

display('Assuming phasor plot to be sum of all components');
pause(2);
Pd = sum_Pd;
Pn = sum_Pn;
Pc = sum_Pc;

if (DISPLACEMENT)
    figure;
    
    hs_polar = subplot(1,2,2);
    polar(0,1),
    hold on;
    
    hs_complex = subplot(1,2,1);
    title('Displacement');
    hold on;
    line([0 0], [-1 1], 'LineStyle', '--', 'Color', 'k');
    line([-1 1], [0 0], 'LineStyle', '--', 'Color', 'k');
    axis([-1, 1, -1, 1])
    
    subplot(hs_polar);
    arrow3([0 0], Pd(1,:), 'r');
    
    subplot(hs_complex);
    arrow3([0 0], Pd(1,:), 'r');
    if (MOVIE)
        writerObj_disp = VideoWriter('disp_phasor.avi');
        writerObj_disp.FrameRate = 2;
        writerObj_disp.Quality   = 95;
        open(writerObj_disp);
        %set(gca, 'nextplot', 'replacechildren');
        set(gcf, 'Renderer', 'zbuffer');
        
        frame = getframe(gcf);
        writeVideo(writerObj_disp, frame);
    end
    for i = 2:num_files
        pause(0.25);
        
        subplot(hs_polar);
        arrow3([0 0], Pd(i,:), 'colors');
        
        subplot(hs_complex);
        arrow3([0 0], Pd(i,:), 'colors');
        axis([-1, 1, -1, 1]);
        
        if (MOVIE)
            frame = getframe(gcf);
            writeVideo(writerObj_disp, frame);
        end
    end
    hold off;
    
    if (MOVIE)
        close(writerObj_disp);
    end
end

if (REFRACTIVE)
    figure;
    
    hs_polar = subplot(1,2,2);
    polar(0,1),
    hold on;
    
    hs_complex = subplot(1,2,1);
    title('Refractive index');
    hold on;
    line([0 0], [-1 1], 'LineStyle', '--', 'Color', 'k');
    line([-1 1], [0 0], 'LineStyle', '--', 'Color', 'k');
    axis([-1, 1, -1, 1])
    
    subplot(hs_polar);
    arrow3([0 0], Pn(1,:), 'r');
    
    subplot(hs_complex);
    arrow3([0 0], Pn(1,:), 'r');
    
    if (MOVIE)
        writerObj_refract = VideoWriter('refractive_phasor.avi');
        writerObj_refract.FrameRate = 2;
        writerObj_refract.Quality   = 95;
        open(writerObj_refract);
        %set(gca, 'nextplot', 'replacechildren');
        set(gcf, 'Renderer', 'zbuffer');
        
        frame = getframe(gcf);
        writeVideo(writerObj_refract, frame);
    end
    for i = 2:num_files
        pause(0.25);
        
        subplot(hs_polar);
        arrow3([0 0], Pn(i,:), 'colors');
        
        subplot(hs_complex);
        arrow3([0 0], Pn(i,:), 'colors');
        
        if (MOVIE)
            frame = getframe(gcf);
            writeVideo(writerObj_refract, frame);
        end
    end
    hold off;
    
    if (MOVIE)
        close(writerObj_refract);
    end
end

if (COMBINED)
    figure;
    
    hs_polar = subplot(1,2,2);
    polar(0,1),
    hold on;
    
    hs_complex = subplot(1,2,1);
    title('Combined (displacement + refractive index)');
    hold on;
    line([0 0], [-1 1], 'LineStyle', '--', 'Color', 'k');
    line([-1 1], [0 0], 'LineStyle', '--', 'Color', 'k');
    axis([-1, 1, -1, 1])
    
    subplot(hs_polar);
    arrow3([0 0], Pc(1,:), 'r');
    
    subplot(hs_complex);
    arrow3([0 0], Pc(1,:), 'r');
    
    if (MOVIE)
        writerObj_combined = VideoWriter('combined_phasor.avi');
        writerObj_combined.FrameRate = 2;
        writerObj_combined.Quality   = 95;
        open(writerObj_combined);
        %set(gca, 'nextplot', 'replacechildren');
        set(gcf, 'Renderer', 'zbuffer');
        
        frame = getframe(gcf);
        writeVideo(writerObj_combined, frame);
    end
    for i = 2:num_files
        pause(0.25);
        
        subplot(hs_polar);
        arrow3([0 0], Pc(i,:), 'colors');
        
        subplot(hs_complex);
        arrow3([0 0], Pc(i,:), 'colors');
        
        if (MOVIE)
            frame = getframe(gcf);
            writeVideo(writerObj_combined, frame);
        end
    end
    hold off;
    
    if (MOVIE)
        close(writerObj_combined);
    end
end


for i=1:length(exit_data)
    if (DISPLACEMENT)
        modulation_index.D(i) = abs(abs(angle(exit_data(1).phasor_D(path))) - ...
                                    abs(angle(exit_data(i).phasor_D(path))));
    end
    if (REFRACTIVE)
        modulation_index.N(i) = abs(abs(angle(exit_data(1).phasor_N(path))) - ...
                                    abs(angle(exit_data(i).phasor_N(path))));
    end
    if (COMBINED)
        modulation_index.C(i) = abs(abs(angle(exit_data(1).phasor_C(path))) - ...
                                    abs(angle(exit_data(i).phasor_C(path))));
    end
end

p0 = mod(angle(exit_data(1).phasor_C(path)*180/pi+360), 360);
p1 = mod(angle(exit_data(8).phasor_C(path)*180/pi+360), 360);

min(abs(p1-p0))
max(abs(p1-p0))


% origin = zeros(3,2);
% P = [real(phasor_t1(tagged_index(1))) imag(phasor_t1(tagged_index(1)))];
% P = [P; real(phasor_t2(tagged_index(1))) imag(phasor_t2(tagged_index(1)))];
% P = [P; real(phasor_t13(tagged_index(1))) imag(phasor_t13(tagged_index(1)))]
% 
% figure;
% arrow(origin, P);




