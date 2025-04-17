close all;
clc;
clear;
load genus.mat;

% Example usage:
data = microbiome_data; % Example data (ensure no zeros)
microbiome_data = clr_transform(data);


save('C:\Users\AL-GHADIR\Desktop\table6\MEDA-Toolbox-1.6\MEDA-Toolbox-1.6\GenusData\data\CLR_Genus.mat', 'microbiome_data');


function clr_data = clr_transform(data)
    % Check if the data contains zeros (replace or handle appropriately)
    if any(data(:) == 0)
        warning('Data contains zeros. Consider using a pseudocount.');
        data(data == 0) =  1; % Small pseudocount to avoid log(0)
    end
    
    % Compute the geometric mean for each row
    geo_mean = exp(mean(log(data), 2));
    
    % Apply the CLR transformation
    clr_data = log(data ./ geo_mean);
end


