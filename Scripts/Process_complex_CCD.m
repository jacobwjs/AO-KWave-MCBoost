% Processes complex to form speckle pattern from a virtual CCD detection.

fid = fopen('SPECKLE_trans-100photons 3.dat');

CCD = zeros(512, 512);
i = 1;
tline = fgetl(fid);
while ischar(tline)
    CCD(i,:) = str2num(tline);
    i = i+1;
    tline = fgetl(fid);
end

fclose(fid);

figure; imagesc(abs(CCD).^2);