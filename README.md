# Digital-twin-predicting-diet-response-before-and-after-long-term-fasting
Oscar Silfvergren, Christian Simonsson, Mattias Ekstedt, Peter Lundberg, Peter Gennemark, Gunnar Cedersund

Mathematical model that predicts diets consiting of protein, carbohydrates and fasting interventions through 
integrating partial insight from several clinical studies (estimation data) into an interconnective and big 
picture (predictions).

All scripts, functions, models, and data used are attached in this folder. Please run scripts initially through 
the script called "main" to ensure that all data and variables are loaded properly. 
Scritps have been created to make it easy to plot all figures used in article, see "Figure2.m", "Figure3.m" etc.

To be able to create a mex file you need a comparable complier (for example MinGW64 Compiler) and IQM toolbox installed. The script "LoadEverything.m" will help you install IQMtools if you unzip folder and place the single 
"IQMtools" folder, that contains all IQM files, in your directory. Otherwise, please see attached map "IQMtools" and read 
the instructions listed to manually install.

If needed, please contact gunnar.cedersund@liu.se for any questions.

Best regards,
Oscar Silfvergren



%--------- Update 2022-12-05 ---------%

* Fixed inconsistencies and problems regarding the simulation of figure 6. Decleration of body compostions in both optimisation and plotfigure are now correct and differences in simulation outside and inide costfunction are now removed. 

* A new map is added called "New improved simulation method" is added where several implementation and parameter estimation is improved. Please see map "new improved simulation method" if you wish to implement model
