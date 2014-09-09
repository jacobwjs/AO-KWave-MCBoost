
% 'data'   => Processed data containing OPLs and phase shifts
% 'info'   => Information about the settings used and others pertaining to
%             the calculations.
function [data, info] = Process_AO_sim_OPL_files(filename, MEAS_TYPE)


% Save which data type we are processing.
info.meas_type = MEAS_TYPE;

% Calculate the wavenumber
lambda = 532e-9;
k = 2*pi/lambda;

% Store the wavenumber so it's known at a later time how these phase shifts
% were calculated.
info.wavenumber = k;

% Array that contains the seeds that generated the paths of the photon(s)
% that were modulated and unmodulated.
info.modulated_seeds   = [];
info.unmodulated_seeds = [];

% Array that contains the indices of the modulated photons (i.e. OPL's that
% have changed between runs).
info.modulation_indices_Nmap = [];
info.modulation_indices_Dmap = [];
info.modulation_indices_Cmap = [];

% Holds the OPLs after being read in from disk.
data.displacement_OPLs = [];
data.refractive_OPLs   = [];
data.combined_OPLs     = [];

% Holds the tagged and untagged electric field for each mechanism
data.E_tagged_Dmap = [];
data.E_tagged_Nmap = [];
data.E_tagged_Cmap = [];

data.E_total_Dmap = [];
data.E_total_Nmap = [];
data.E_total_Cmap = [];




temp = [];


if (strcmp(MEAS_TYPE, 'OPL_data'))
    if filename == 0
        [filename,path] = uigetfile('*.dat','Open optical path length *.dat file');
        if filename==0 ,return ,end
        filename = strcat(path,filename);
    end
    
    % Open the file pertaining to the optical path lengths for the simulation.
    fid = fopen(filename);
    
    % Prime the loop.
    photon = fgetl(fid);
    photon = textscan(photon,'%f %f','delimiter',',');
    temp = photon{1};
    data.displacement_OPLs = [data.displacement_OPLs, temp];
    
    cnt = 0;
    while (~feof(fid))
        cnt = cnt+1;
        
        photon = fgetl(fid);
        photon = textscan(photon,'%f %f','delimiter',',');
        
        temp = photon{1};
        if (size(temp,1) ~= size(data.displacement_OPLs,1))
            display('Skipping data: ')
            cnt
            continue;
        end
        data.displacement_OPLs = [data.displacement_OPLs, temp];
        
        temp = photon{2};
        data.refractive_OPLs = [data.refractive_OPLs, temp];
        
        temp = photon{3};
        data.combined_OPLs = [data.combined_OPLs, temp];
    end
    
    fclose(fid);
    
    data.displacement_OPLs = data.displacement_OPLs';
    data.refractive_OPLs   = data.refractive_OPLs';
    data.combined_OPLs     = data.combined_OPLs';
    
    
    info.modulation_indices_Nmap = find(data.refractive_OPLs(:,1) ~=...
                                    data.refractive_OPLs(:,2));
    info.modulation_indices_Dmap = find(data.displacement_OPLs(:,1) ~=...
                                    data.displacement_OPLs(:,2));
    info.modulation_indices_Cmap = find(data.combined_OPLs(:,1) ~=...
                                    data.combined_OPLs(:,2)); 
                                
    
  
    % Calculate the phase shift (in degrees) for all photons that were modulated.
    % The OPLs, after multiplication with the wavenumber, are in radians. Here
    % we just finish the computation of phase shifts by putting them in
    % degrees.
    time_steps = size(data.refractive_OPLs,2);
    data.displacement_phase_shifts = k*180/pi*(data.displacement_OPLs(:,2:time_steps) - ...
                                 data.displacement_OPLs(:,ones(1,time_steps-1))); 
    data.refractive_phase_shifts = k*180/pi*(data.refractive_OPLs(:,2:time_steps) - ...
                                 data.refractive_OPLs(:,ones(1,time_steps-1)));                           
    data.combined_phase_shifts = k*180/pi*(data.combined_OPLs(:,2:time_steps) - ...
                                 data.combined_OPLs(:,ones(1,time_steps-1)));                         
                             
elseif (strcmp(MEAS_TYPE, 'EXIT_data'))
    % We need to open two files for comparing the exit data.
    if filename == 0
        [filename,path] = uigetfile('*.dat','Open optical path length *.dat file');
        if filename==0 ,return ,end
        filename = strcat(path,filename);
    end
    t0 = dlmread(filename);
    
    filename = 0;
    if filename == 0
        [filename,path] = uigetfile('*.dat','Open optical path length *.dat file');
        if filename==0 ,return ,end
        filename = strcat(path,filename);
    end
    t1 = dlmread(filename);
    
    % The exit data file is written out to disk as follows.
    % column(1)  => weight
    % column(2)  => displaced_OPL
    % column(3)  => refractive_OPL
    % column(4)  => combined_OPL
    % column(5)  => x-axis exit location
    % column(6)  => y-axis exit location
    % column(7)  => z-axis exit location
    % column(8)  => seeds.s1a
    % column(9)  => seeds.s2
    % column(10) => seeds.s3
    % column(11) => seeds.s4
    
    
    % Store the OPLs incase later processing is needed.
    data.displacement_OPLs(:,1) = t0(:,2);
    data.displacement_OPLs(:,2) = t1(:,2);
    data.refractive_OPLs(:,1) = t0(:,3);
    data.refractive_OPLs(:,2) = t1(:,3);
    data.combined_OPLs(:,1) = t0(:,4);
    data.combined_OPLs(:,2) = t1(:,4);
    
    % Find the indices where the OPLs have changed between the two files.
    info.modulation_indices_Dmap = find(t0(:,2) ~= t1(:,2));
    info.modulation_indices_Nmap = find(t0(:,3) ~= t1(:,3));
    info.modulation_indices_Cmap = find(t0(:,4) ~= t1(:,4));
    
    
    
    
    
    % Store the seeds that created photon paths that underwent modulation
    % and those that didn't (i.e. tagged and untagged photons).
    % NOTE:
    % - Not all mechanisms are turned on, so we check which does not have
    %   empty 'modulation_indices' and store those seeds.
    % - Seeds in 't1' and 't0' are the same, so just save one.
    %
    modulation_indices = [];
    if (~isempty(info.modulation_indices_Dmap))
        modulation_indices = info.modulation_indices_Dmap;
        data.modulated_exit_photons_Dmap_t0 = t0(modulation_indices,:);
        data.modulated_exit_photons_Dmap_t1 = t1(modulation_indices,:);
    end
    if (~isempty(info.modulation_indices_Nmap))
        modulation_indices = info.modulation_indices_Nmap;
        data.modulated_exit_photons_Nmap_t0 = t0(modulation_indices,:);
        data.modulated_exit_photons_Nmap_t1 = t1(modulation_indices,:);
    end
    if (~isempty(info.modulation_indices_Cmap))
        modulation_indices = info.modulation_indices_Cmap;
        data.modulated_exit_photons_Cmap_t0 = t0(modulation_indices,:);
        data.modulated_exit_photons_Cmap_t1 = t1(modulation_indices,:);
    end
    % Update the arrays that contain seeds based on the indices found
    % above.
    info.modulated_seeds = t1(modulation_indices, 8:11);
    % Grab all the seeds and remove the modulated path producing seeds.
    info.unmodulated_seeds = t1(:,8:11);
    info.unmodulated_seeds(modulation_indices,:) = [];
        
                                
    % Calculate the phase shift due to change in OPL for each mechanism.
    data.displacement_phase_shifts = k*180/pi*(t1(:,2) - t0(:,2));
    data.refractive_phase_shifts = k*180/pi*(t1(:,3) - t0(:,3));
    data.combined_phase_shifts = k*180/pi*(t1(:,4) - t0(:,4));
else
    display('Error: Measurement type not specified');
end



                             
               
% Specify the bin 'width' in degrees for the absolute value case.
magnitude_bin_width = 2:2:450;
% Specify the bin 'width' in degrees for the total value case.
total_bin_width = -450:2:450;

% Save these values for later.
info.magnitude_bin_width = magnitude_bin_width;
info.total_bin_width = total_bin_width;




% % Plot the histogram of the absolute phase shift in degrees.
if (~isempty(info.modulation_indices_Cmap))
    
    data.E_tagged_Cmap = sqrt(t1(modulation_indices,1)).* ...
                         exp(1i*k*(t1(modulation_indices,4) - t0(modulation_indices,4)));
    data.E_total_Cmap  = sqrt(t0(:,1)).* ...
                         exp(1i*k*t0(:,4));
    data.tagged_fraction_Cmap = mean(abs(data.E_tagged_Cmap).^2)/...
                                mean(abs(data.E_total_Cmap).^2);
    
    max_phase_shift = 0;
    if (max(data.combined_phase_shifts) >...
        abs(min(data.combined_phase_shifts)))
        max_phase_shift = max(data.combined_phase_shifts);
    else
        max_phase_shift = abs(min(data.combined_phase_shifts));
    end
    
    
%     figure;
%     %subplot(2,2,1);
%     hist(abs(data.combined_phase_shifts),...
%          2:2:max_phase_shift);
%     title('Magnitude of phase shifts (\Delta{n}+\Delta{d})');
%     xlabel('\Phi_j (degrees)');
%     ylabel('Number of photons');
%     set(gca,'FontSize',14,'fontWeight','normal')
%     set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
    
    figure;
    %subplot(2,2,2);
    hist(abs(data.combined_phase_shifts(info.modulation_indices_Cmap)),...
         2:2:max_phase_shift);
    title('Magnitude of phase shifts (\Delta{n}+\Delta{d})');
    xlabel('\Phi_j (degrees)');
    ylabel('Number of photons');
    set(gca,'FontSize',14,'fontWeight','normal')
    set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
    
    
%     subplot(2,2,3)
%     hist(data.combined_phase_shifts,...
%         -max_phase_shift:2:max_phase_shift);
%     title('Absolute phase shifts (\Delta{n}+\Delta{d})');
%     xlabel('\Phi_j (degrees)');
%     ylabel('Number of photons');
%     set(gca,'FontSize',14,'fontWeight','normal')
%     set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
    
    figure;
    %subplot(2,2,4)
    hist(data.combined_phase_shifts(info.modulation_indices_Cmap),...
         -max_phase_shift:2:max_phase_shift);
    title('Absolute phase shifts (\Delta{n}+\Delta{d})');
    xlabel('\Phi_j (degrees)');
    ylabel('Number of photons');
    set(gca,'FontSize',14,'fontWeight','normal')
    set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
end
% else

if (~isempty(info.modulation_indices_Nmap))
    
    
    data.E_tagged_Nmap = sqrt(t1(modulation_indices,1)).* ...
                         exp(1i*k*(t1(modulation_indices,3) - t0(modulation_indices,3)));
    data.E_total_Nmap  = sqrt(t0(:,1)).* ...
                         exp(1i*k*t0(:,3));
    data.tagged_fraction_Nmap = mean(abs(data.E_tagged_Cmap).^2)/...
                                mean(abs(data.E_total_Cmap).^2);
    
    
    max_phase_shift = 0;
    if (max(data.refractive_phase_shifts) >...
            abs(min(data.refractive_phase_shifts)))
        max_phase_shift = max(data.refractive_phase_shifts);
    else
        max_phase_shift = abs(min(data.refractive_phase_shifts));
    end
    
%     figure;
%     subplot(2,2,1);
%     hist(abs(data.refractive_phase_shifts),...
%         2:2:max_phase_shift);
%     title('Magnitude of phase shifts (\Delta{n})');
%     xlabel('\Phi_j (degrees)');
%     ylabel('Number of photons');
%     set(gca,'FontSize',14,'fontWeight','normal')
%     set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
    
    figure;
    %subplot(2,2,2);
    hist(abs(data.refractive_phase_shifts(info.modulation_indices_Nmap)),...
        2:2:max_phase_shift);
    title('Magnitude of phase shifts (\Delta{n})');
    xlabel('\Phi_j (degrees)');
    ylabel('Number of photons');
    set(gca,'FontSize',14,'fontWeight','normal')
    set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
    
%     subplot(2,2,3)
%     hist(data.refractive_phase_shifts,...
%         -max_phase_shift:2:max_phase_shift);
%     title('Absolute phase shifts (\Delta{n})');
%     xlabel('\Phi_j (degrees)');
%     ylabel('Number of photons');
%     set(gca,'FontSize',14,'fontWeight','normal')
%     set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
    
    figure;
    %subplot(2,2,4)
    hist(data.refractive_phase_shifts(info.modulation_indices_Nmap),...
        -max_phase_shift:2:max_phase_shift);
    title('Absolute phase shifts (\Delta{n})');
    xlabel('\Phi_j (degrees)');
    ylabel('Number of photons');
    set(gca,'FontSize',14,'fontWeight','normal')
    set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
end


if (~isempty(info.modulation_indices_Dmap))
    
    
    data.E_tagged_Dmap = sqrt(t1(modulation_indices,1)).* ...
                         exp(1i*k*(t1(modulation_indices,2) - t0(modulation_indices,2)));
    data.E_total_Dmap  = sqrt(t0(:,1)).* ...
                         exp(1i*k*t0(:,2));
    data.tagged_fraction_Dmap = mean(abs(data.E_tagged_Cmap).^2)/...
                                mean(abs(data.E_total_Cmap).^2);
    
    
    max_phase_shift = 0;
    if (max(data.displacement_phase_shifts) >...
            abs(min(data.displacement_phase_shifts)))
        max_phase_shift = round(max(data.displacement_phase_shifts));
    else
        max_phase_shift = round(abs(min(data.displacement_phase_shifts)));
    end
    
%     figure;
%     subplot(2,2,1);
%     hist(abs(data.displacement_phase_shifts),...
%         2:2:max_phase_shift);
%     title('Magnitude of phase shifts (\Delta{d})');
%     xlabel('\Phi_j (degrees)');
%     ylabel('Number of photons');
%     set(gca,'FontSize',14,'fontWeight','normal')
%     set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
    
    figure;
    %subplot(2,2,2);
    hist(abs(data.displacement_phase_shifts(info.modulation_indices_Dmap)),...
        2:2:max_phase_shift);
    title('Magnitude phase shifts (\Delta{d})');
    xlabel('\Phi_j (degrees)');
    ylabel('Number of photons');
    set(gca,'FontSize',14,'fontWeight','normal')
    set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
    
%     subplot(2,2,3);
%     hist(data.displacement_phase_shifts,...
%         -max_phase_shift:2:max_phase_shift);
%     title('Magnitude of phase shifts (\Delta{d})');
%     xlabel('\Phi_j (degrees)');
%     ylabel('Number of photons');
%     set(gca,'FontSize',14,'fontWeight','normal')
%     set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
    
    figure;
    %subplot(2,2,4);
    hist(data.displacement_phase_shifts(info.modulation_indices_Dmap),...
        -max_phase_shift:2:max_phase_shift);
    title('Absolute phase shifts (\Delta{d})');
    xlabel('\Phi_j (degrees)');
    ylabel('Number of photons');
    set(gca,'FontSize',14,'fontWeight','normal')
    set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','normal')
end


% figure; 
% plot(n_data.refractive_OPLs(:,1)./3e8);
% title('Time-of-flight');
% xlabel('Photon_j');
% ylabel('Time [sec]');
% set(gca,'FontSize',16,'fontWeight','normal')
% set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')


% % Example of calculating tagged fraction without creating speckle
% % patterns using 'Speckle-boost'
% E_t0 = (sqrt(t0(:,1)).*exp(1i*t0(:,4)*2*pi/532e-9));
% E_t1 = (sqrt(t1(:,1)).*exp(1i*t1(:,4)*2*pi/532e-9));
% E_tagged = E_t1 - E_t0;
% mean(abs(E_tagged).^2)/mean(abs(E_t1).^2)


% % Example of writing seeds out to disk after processing.
% %-------------------------------------------------------------------
% dlmwrite('trans_modulated_seeds_10659photons.dat', info.modulated_seeds,...
%          'delimiter', ' ',...
%          'precision', 10,...
%          'newline', 'pc');

% % Example of transparent, colored histograms
% %-------------------------------------------------------------------
% figure;
% hold on;
% hist(n_3p75mm.phase_shifts(info_3p75mm.modulation_indices), 10:10:4000);
% hist(n_2p5mm.phase_shifts(info_2p5mm.modulation_indices), 10:10:4000);
% hist(n_1p25mm.phase_shifts(info_1p25mm.modulation_indices), 10:10:4000);
% h = findobj(gca,'Type','patch');
% set(h(1), 'FaceColor', 'r', 'EdgeColor', 'k', 'FaceAlpha', 0.95)
% set(h(2), 'FaceColor', 'y', 'EdgeColor', 'k', 'FaceAlpha', 0.50)
% set(h(3), 'FaceColor', 'g', 'EdgeColor', 'k', 'FaceAlpha', 0.70)

end % end function 
