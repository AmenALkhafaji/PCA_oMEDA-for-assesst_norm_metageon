%% Please, read Guide line in the descrption below 

%  Code:  Valdtion of normalztion methods trhough PCA and Omaeda packges                                                           
%  Author: Amen Adnan Khabeer; Jose Comacho; Crolina Comez                                    
%  Date: 19/02/2024                                         
%  Email: amen.a.khabeer@uotechnology.edu.iq ;josecamacho@ugr.es;gomezll@ugr.es                

%% Data Preprocessing :   data Raw
% # Load Raw Data

close all;
clc;
%clear all;
load phylum.mat;
% ## Sorting data according to their depths
[~, sorting_index] = sort(sum(microbiome_data, 2));
sorted_transformed_data = microbiome_data(sorting_index, :);
sorted_infection_status= infection_status(sorting_index, :);


% ## Ploting Omeda for low  depths with oMEDA and ground truth
dummy = sorted_infection_status(1:100, :);
dummy(dummy == 2) = -1;
rj = omedaPca(sorted_transformed_data(1:100, :),[], sorted_transformed_data(1:100, :), dummy, 'VarsLabel', lab);
tot1 = rj / sqrt(sum((rj).^2));
gt1 = ground_truth / sqrt(sum((ground_truth).^2));
e1 = sum((tot1 - gt1).^2);
fprintf('Total Error at low dataset:%.4f\n', e1);
%fprintf(sprintf('%.4f\n', e1));
% ## Ploting Omeda for high depth and comparing to ground trutho
dummy = sorted_infection_status(200:300, :);
dummy(dummy == 2) = -1;
rj = omedaPca(sorted_transformed_data(200:300, :), [], sorted_transformed_data(200:300, :), dummy, 'VarsLabel', lab);

tot2 = rj / sqrt(sum((rj).^2));
gt2 = ground_truth / sqrt(sum((ground_truth).^2));
e2 = sum((tot2 - gt2).^2);
fprintf('Total Error at high dataset:%.4f\n', e2);
% fprintf('%.4f\n', e2);
% ## Ploting Omeda for entire dataset depth and comparing to ground truth
dummy = sorted_infection_status;
dummy(dummy == 2) = -1;

rj = omedaPca(sorted_transformed_data, [], sorted_transformed_data, dummy, 'VarsLabel', lab);
tot3 = rj / sqrt(sum((rj).^2));
gt3 = ground_truth / sqrt(sum((ground_truth).^2));
e3 = sum((tot3 - gt3).^2);
fprintf('Total Error at entire dataset:%.4f\n', e3);
%fprintf('%.4f\n', e3);
% ## Plotting all depths 
s = cellstr(categorical(infection_status, [1, 2], ["Healthy ecosystems", "Dybiosis ecosystems"]));
myemptylabel(1:300,1)="";
scoresPca(microbiome_data,'PCs',1:2,"ObsTest",microbiome_data,'ObsClass',s,'ObsLabel',myemptylabel,'PlotCal', false);
scoresPca(microbiome_data,"PCs",1:2,"ObsTest",microbiome_data,"ObsLabel",sum(microbiome_data,2),'ObsLabel',myemptylabel,'ObsClass',sum(microbiome_data,2),'PlotCal', false);
 loadingsPca(microbiome_data,"PCs",1:2,'VarsLabel',lab','VarsClass',lab');


