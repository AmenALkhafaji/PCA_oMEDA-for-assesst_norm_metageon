% Read the CSV file
data = readtable('Deseq_rarefaction.csv');
microbiome_data = table2array(data(1:end,: ));


save('C:\Users\AL-GHADIR\Desktop\table6\MEDA-Toolbox-1.6\MEDA-Toolbox-1.6\PhylumData\data\Deseq_rarefaction.mat', 'microbiome_data');
