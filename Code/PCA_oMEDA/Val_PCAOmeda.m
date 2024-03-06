function   Val_PCAOmeda(data_frame,ground_truth)


%% Please, read Guide line in the descrption below 
% 
%
%
%  Code:  Valdtion of normalztion methods trhough PCA and Omaeda packges                                                           
%  Author: Amen Adnan Khabeer; Jose Comacho; Crolina Comez                                    
%  Date: 19/02/2024                                         
%  Email: amen.a.khabeer@uotechnology.edu.iq ;josecamacho@ugr.es;gomezll@ugr.es                
%                                      
%    Val_PCAOmeda( data_frame as (2D string) , ground_truth as ( vector))                                                                  #
%  
%    1 : Both arguements are reqiured, no optional arguements in this function.                                                                    #
%    2 : Please load your Simulation or data   
%   
%        Data_frame: is your normalized data ( 2D matrix with string values) ,where, row=no of
%                samples, clo represents no taxan, last col represent controls( healthy/unhealthy),  the header is the name of bactria
%    
%        Ground_truth: a vector with double values, holding ground truth values  
%
% 
%% A: Data_frame                                                       
%     1 :Load your data (as CSV)  in your workspace
%     2 :From the new window that open, choose output type : array string
%     3 :Before press import, please select all data include header
% 
%% B: ground_truth
%    1 : Load your ground truth (as CSV) in your workspce
%    2 : From the open window, output type : Numric matrix
%    3 : Press import 
%
%% LaTeX Inline Expression Example
% This is an expression: $e = (sum_{i=1}^{n} (ngt_i - {tot}_i))^2$ 
% 
% \[\text{tot1} = \frac{\sqrt{\lvert r_j \rvert} \cdot \text{sign}(r_j)}{\sum \sqrt{\lvert r_j \rvert}}\]
%% 
%
%

%% Assguining the arguments and split the depth for two area

 


% Exracuted the name of Classes ( Taxa) from simalution which is thie
% frist raw with out last col (Tags)
mystr=data_frame;
class_name=mystr(1,1:width(mystr)-1); 
ind1 = str2double(mystr);

% extacted the data only from our dataset without header and last col(Tags)
xx=ind1(2:height(ind1),1:width(ind1)-1);

% extracted the ecyosystems form our data (healthy-unhealthy ecosystems)
l1=xx(:,width(xx));

% split points according to depth to high and low

low=100;
high=200;


%% Ploting omeda for high depth and comparing to groundtruth

lh=(l1(high:size(xx,1)));
xh=xx(high:size(xx,1),:);
rj=omeda_pca(xh,1:2,[],lh);
tot= sqrt( abs(rj)).* sign(rj)/sum(sqrt( abs(rj)));
 
gt=(normalize(transpose(ground_truth),"norm",1));
ngt= ( abs(gt)).* sign(gt)/sum(( abs(gt)));
e=(sum(ngt'-tot)).^2;
fprintf('Total Error in high depth : %s\n', e);



%% Ploting omeda for low depth and comparing to groundtruth

ll=(l1(1:low));
xl=xx(1:low,:);
rj=omeda_pca(xl,1:2,[],ll);
tot= sqrt( abs(rj)).* sign(rj)/sum(sqrt( abs(rj)));
gt=(normalize(transpose(ground_truth),"norm",1));
ngt= ( abs(gt)).* sign(gt)/sum(( abs(gt)));
e1=((ngt'-tot).^2);
e=(sum(ngt'-tot)).^2;
fprintf('Total Error low depth : %s\n', e);






%% Ploting omeda for whole depth and comparing to groundtruth
% Find the spscies that contrubate this diffrence
% MATLAB code

rj=omeda_pca(xx,1:2,[],l1);
tot1= sqrt( abs(rj)).* sign(rj)/sum(sqrt( abs(rj)));
gt1=(normalize(transpose(ground_truth),"norm",1));
ngt1= ( abs(gt1)).* sign(gt1)/sum(( abs(gt1)));
e12=(sum(ngt1'-tot1)).^2;
disp(sprintf('Total Error whole depth: %s', e12));
end