% This file plots the magnitude of the sidebands for various modulation
% indices using Bessel functions of the 1st kind for integer orders.
function [stem_data] = sideband_stem_plot(mod_index, order)


MOVIE = true;


stem_data = [];


% Create the labels for the plot
bands = [-order:order];
clear xlabels;
xlabels = [];
xlabels = {' '};
for i=1:length(bands)
    if (bands(i) == 0)
        xlabels = [xlabels, {'Wc'}];
    else
        xlabels = [xlabels, {strcat(num2str(bands(i)), 'Wm')}];
    end
end
xlabels = [xlabels, {' '}];


figure('position', [100, 50, 1050, 600])  % create new figure with specified size
%figure; 

% Coord=get(gca,'Position');
% set(gca,'Position', [0,0,700,700]);

if (isempty(mod_index))
    modulation_index = 0.1;
    
     if (MOVIE)
        writerObj_disp = VideoWriter('sidebands.avi');
        writerObj_disp.FrameRate = 2;
        writerObj_disp.Quality   = 95;
        open(writerObj_disp);
        %set(gca, 'nextplot', 'replacechildren');
        set(gcf, 'Renderer', 'zbuffer');
        
        frame = getframe(gcf);
        writeVideo(writerObj_disp, frame);
    end
    
    for i=1:60
        
        % Bessel function orders 0-to-'order'
        %         pos_temp = besselj(0:order, modulation_index);
        %         neg_temp = fliplr(pos_temp(2:end)); % Don't include the carrier.
        %         stem_data(i,:) = [neg_temp, pos_temp];  % Positive sidebands, negative sidebands, and the carrier.
        %stem(stem_data(i,:), 'LineWidth', 2);
        
        stem_data = abs(besselj(-order:order, modulation_index));
        
        clf;
        hold on;
        stem(stem_data(1:order), 'or', 'LineWidth', 2);
        stem(order+1, stem_data(order+1), 'ob', 'LineWidth', 2);
        stem([order+2:order*2+1], stem_data(order+2:end), 'or', 'LineWidth', 2);
        hold off;
        axis([0 (order*2 + 2) 0 1]);
        set(gca, 'Xtick', 0:length(xlabels)); % Change x-axis ticks
        set(gca, 'XtickLabel', xlabels);
        legend(strcat('m=',num2str(modulation_index)));
        
        if (MOVIE)
            frame = getframe(gcf);
            writeVideo(writerObj_disp, frame);
        end
        
        pause(0.25);
        modulation_index = modulation_index + 0.1;
    end
    
    if (MOVIE)
        close(writerObj_disp);
    end
   
else
    modulation_index = mod_index;
    % Bessel function orders 0-to-'order'
    %pos_temp = besselj(0:order, modulation_index);
    %neg_temp = fliplr(pos_temp(2:end)); % Don't include the carrier.
    %stem_data = [neg_temp, pos_temp];   % Positive sidebands, negative sidebands, and the carrier.
    
    stem_data = abs(besselj(-order:order, modulation_index));
    
    hold on;
    stem(stem_data(1:order), 'or', 'LineWidth', 2);
    stem(order+1, stem_data(order+1), 'ob', 'LineWidth', 2);
    stem([order+2:order*2+1], stem_data(order+2:end), 'or', 'LineWidth', 2);
    hold off;
    
    axis([0 (order*2 + 2) 0 1]);
    set(gca, 'Xtick', 0:length(xlabels)); % Change x-axis ticks
    set(gca, 'XtickLabel', xlabels);
    
    legend(strcat('m=',num2str(modulation_index)));
end

