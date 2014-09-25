% This file plots the magnitude of the sidebands for various modulation
% indices using Bessel functions of the 1st kind for integer orders.
function [stem_data] = sideband_stem_plot(mod_index, order, fig_title)


MOVIE = true;


stem_data = [];

figure_title = fig_title;

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

% create new figure with specified size
figure('position', [100, 50, 1050, 600])  


if (MOVIE)
        writerObj = VideoWriter(strcat('sidebands_', figure_title, '.avi'));
        writerObj.FrameRate = 2;
        writerObj.Quality   = 100;
        open(writerObj);
        %set(gca, 'nextplot', 'replacechildren');
        set(gcf, 'Renderer', 'zbuffer');
        
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
    end


if (isempty(mod_index))
    modulation_index = 0.1;
       
    for i=1:60
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
        legend_h = legend(strcat('m=',num2str(modulation_index(i))));
        set(legend_h, 'FontSize', 14);
        
        
        if (MOVIE)
            frame = getframe(gcf);
            writeVideo(writerObj, frame);
        end
        
        pause(0.25);
        modulation_index = modulation_index + 0.1;
    end
else
    iterations = length(mod_index);
    modulation_index = mod_index;
    
    for i=1:iterations
        stem_data = abs(besselj(-order:order, modulation_index(i)));
        
        clf;
        hold on;
        stem(stem_data(1:order), 'or', 'LineWidth', 2);
        stem(order+1, stem_data(order+1), 'ob', 'LineWidth', 2);
        stem([order+2:order*2+1], stem_data(order+2:end), 'or', 'LineWidth', 2);
        hold off;
        
        axis([0 (order*2 + 2) 0 1]);
        set(gca, 'Xtick', 0:length(xlabels)); % Change x-axis ticks
        set(gca, 'XtickLabel', xlabels);
    
        legend_h = legend(strcat('m=',num2str(modulation_index(i))));
        set(legend_h, 'FontSize', 14);
        title(figure_title);
        
        if (MOVIE)
            frame = getframe(gcf);
            writeVideo(writerObj, frame);
        end
        
        pause(0.25);
    end
end

if (MOVIE)
    close(writerObj);
end



