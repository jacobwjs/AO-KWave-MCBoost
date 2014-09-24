% This file plots the magnitude of the sidebands for various modulation
% indices using Bessel functions of the 1st kind for integer orders.
function [stem_data] = sideband_stem_plot(mod_index, order)

stem_data = [];

figure; 

if (isempty(mod_index))
    modulation_index = 0.1;
    
    for i=1:20
        
        % Bessel function orders 0-to-'order'
        pos_temp = besselj(0:order, modulation_index);
        neg_temp = fliplr(pos_temp(2:end)); % Don't include the carrier.
        stem_data(i,:) = [neg_temp, pos_temp];  % Positive sidebands, negative sidebands, and the carrier.
        
        stem(stem_data(i,:), 'LineWidth', 2);
        axis([0 (order*2 + 2) 0 1]);
        
        pause(0.25);
        modulation_index = modulation_index + 0.1;
    end
   
else
    modulation_index = mod_index;
    % Bessel function orders 0-to-'order'
    pos_temp = besselj(0:order, modulation_index);
    neg_temp = fliplr(pos_temp(2:end)); % Don't include the carrier.
    stem_data = [neg_temp, pos_temp];   % Positive sidebands, negative sidebands, and the carrier.
    
    stem(stem_data, 'LineWidth', 2);
    axis([0 (order*2 + 2) 0 1]);
end

