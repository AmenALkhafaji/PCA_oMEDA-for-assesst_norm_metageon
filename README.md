Our proposed method to assess the normalization methods for metagenomics data.

Contact person: Amen A. Khabeer (amen.a.khabeer@uotechnology.edu.iq)
Last modification of this document: 06/03/24


Instercations: 

	- Download MEDA packge from github repostery through the following: https://github.com/josecamachop/MEDA-Toolbox. Please
          read the readme file to extracte it in your matalbworkspace
	- Download the Simalution & Normalzition folder and place it in R 
        - run simulacion_Filos.R and save the data in your path up to your drive
	- run the other codes for normalztion and save the dataset in the same path the you seclected before
        - pleace your dataset without processing firstly in the worksapce of matlab, please selcet all columan as arrary string
        - pleace your ground truth "Code\Data\Gund Truth\phylum" in your workspace, select all columan as double
        - extracut MEDA in your current folder
        - load Val_PCAOmeda function in your edtior. you can find it in the current folder: \Code\PCA_oMEDA 
        - set both phylum dataset and phylum ground truth to your function example:Val_PCAOmeda( phylums,phylumg) : where phylums denoted as simaltuion and phylumg denoted as ground truth
        - finally run the code. 
        - repeat same process for other dataset that resulting in R only replace the dataset  phylums with your current preprocessing dataset.
        - for Genus_Simaltion reapte same process, with conserdtion of replace ground truth and new dataset.
        - Here we go ,You did it! well done!!
        
       
Please, acknowledge the use of this software by refercing it: "Camacho, J., Pérez, A., Rodríguez, R., Jiménez-Mañas, E. Multivariate 
Exploratory Data Analysis (MEDA) Toolbox. Chemometrics and Intelligent Laboratory Systems, 2015, 143: 49-57, available at 
https://github.com/josecamachop/MEDA-Toolbox" Also, please check the documentation of the routines used for more related references. 

Please, do not forget that you are welcome to contact us through the email for more support. As well as, if there are any iusses related with MEDA please contact Dr Jose camacho through email (josecamacho@ugr.es)
We would like to thanks the direct or indirect contribution of colleague:

- Alejandro Garc´ ıa V´azquez for his contriubtion in the created simalution. 

Guied detials about the folder conatants: 
- Data : involves simaltuion datasets folder that you can run without R file and ground truth of each simaltuion
- PCA_oMEDA: involves the functions that recall pca_omeda to assesst your normalatzion method
-Result : involves three folders: high, low and whole depth with results plots (loading,distance omeda,score) in each depth
-Simalution & Normalzition: R based code for creat the simalution and normalzition.



Copyright (C) 2024  Universidad de Granada
Copyright (C) 2024  Amen Adnan,José Camacho Páez
 
This program is free scrpits: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
