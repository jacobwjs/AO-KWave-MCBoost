% Processes complex to form speckle pattern from a virtual CCD detection.
function [CCD] = Process_complex_CCD(x_pixels, y_pixels)

[filename,path] = uigetfile('*.dat','Open optical path length *.dat file');
if filename==0 ,return ,end
filename = strcat(path,filename);

fid = fopen(filename);

CCD = zeros(x_pixels, y_pixels);
i = 1;
tline = fgetl(fid);
while ischar(tline)
    CCD(i,:) = str2num(tline);
    i = i+1;
    tline = fgetl(fid);
end

fclose(fid);

figure; imagesc(abs(CCD).^2);