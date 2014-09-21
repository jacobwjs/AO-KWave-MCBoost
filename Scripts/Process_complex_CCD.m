% Processes complex to form speckle pattern from a virtual CCD detection.
%
% 'x_pixels' => Number of pixels of the virtual detector in the simulation
%               for the x-dimension. Typically 512.
% 'y_pixels' => Number of pixels of the virtual detector in the simulation
%               for the y-dimension. Typically 512.
% 'extension' => If we have multiple files to load (such as data collected
%                over an US period) we simply pass in the extension of the
%                file (e.g. '.dat') and load all files in the current
%                directory based on the extension.
function [CCD_complex_image, tagged_fraction, SPD_tagged_fraction] = Process_complex_CCD(x_pixels, y_pixels, extension)

% Speckle Pattern Difference method of tagged fraction.
SPD_tagged_fraction = [];

% Our theory.
tagged_fraction = [];

if (isempty(extension))
    [filename,path] = uigetfile('*.dat','Open optical path length *.dat file');
    if filename==0 ,return ,end
    filename = strcat(path,filename);
    
    fid = fopen(filename);
    
    CCD_complex_image = zeros(x_pixels, y_pixels);
    i = 1;
    tline = fgetl(fid);
    while ischar(tline)
        CCD_complex_image(i,:) = str2num(tline);
        i = i+1;
        tline = fgetl(fid);
    end
    
    fclose(fid);
    figure; imagesc(abs(CCD_complex_image).^2);
else
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
    
    % Allocate the object that will hold all complex images.
    CCD_complex_image = zeros(num_files, x_pixels, y_pixels);
    
    figure;
    for i=1:num_files
        fid = fopen(char(files(i)));
        
        j = 1;
        tline = fgetl(fid);
        while ischar(tline)
            CCD_complex_image(i, j,:) = str2num(tline);
            j = j+1;
            tline = fgetl(fid);
        end
        fclose(fid);
        
        imagesc(abs(squeeze(CCD_complex_image(i,:,:))).^2);
        drawnow
    end
    
%     
%     % Computation of the tagged fraction. This assumes data was collected
%     % over an US period, and the first image is the unmodulated (i.e. no
%     % pertubation to the medium) case.
    frame1 = squeeze(CCD_complex_image(2,:,:));
    I1 = abs(frame1).^2;
    for i=2:num_files
        frameN = squeeze(CCD_complex_image(i,:,:));
        I_N = abs(frame1 - frameN).^2;
        tagged_fraction(i) = 1/4*mean2(I_N)/mean2(I1);
        
%         tagged_fraction(i) = 2*mean2(abs(frame1 - frameN).^2) / ...
%                              (8*mean2(abs(frame1).^2));
    end
    figure; plot(tagged_fraction);
    

    % Computation of the tagged fraction using Steffen's speckle pattern
    % difference method.
    frame1 = squeeze(CCD_complex_image(2,:,:));
    I1 = abs(frame1).^2;
    for i=2:num_files
        frameN = squeeze(CCD_complex_image(i,:,:));
        In = abs(frameN).^2;
        SPD_tagged_fraction(i) = (mean2((In/mean2(In) - I1/mean2(I1)).^2)) /...
                                 (8);
    end
    figure; plot(SPD_tagged_fraction); 
        
end




