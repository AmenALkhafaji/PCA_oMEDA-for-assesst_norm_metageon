% Read the CSV file
data = readtable('CSS_Rarefaction_Genus.csv');
microbiome_data = table2array(data(1:end, :));
save('C:\Users\AL-GHADIR\Desktop\table6\MEDA-Toolbox-1.6\MEDA-Toolbox-1.6\GenusData\data\CSS_Rarefaction_Genus.mat', 'microbiome_data');
