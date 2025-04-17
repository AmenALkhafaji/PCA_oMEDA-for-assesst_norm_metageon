% Read the CSV file
data = readtable('Genus.csv');
microbiome_data = table2array(data(1:end, 1:36));
infection_status = table2array(data(1:end, 37));

G= readtable('GGT.xlsx');
ground_truth = table2array(G(1:end, 4));
lab = table2array(G(1:end, 1));

save('C:\Users\AL-GHADIR\Desktop\table6\MEDA-Toolbox-1.6\MEDA-Toolbox-1.6\GenusData\data\genus.mat', 'microbiome_data', 'infection_status', 'lab','ground_truth');
