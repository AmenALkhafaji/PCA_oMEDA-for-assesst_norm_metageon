% Read the CSV file
data = readtable('Phylum.csv');
microbiome_data = table2array(data(1:end, 1:5));
infection_status = table2array(data(1:end, 6));

G= readtable('GTP.xlsx');
ground_truth = table2array(G(1:end, 4));
lab = table2array(G(1:end, 1));

save('C:\Users\AL-GHADIR\Desktop\table6\MEDA-Toolbox-1.6\MEDA-Toolbox-1.6\PhylumData\data\phylum.mat', 'microbiome_data', 'infection_status', 'lab','ground_truth');
