Overall:
The main model is run with 
SWEDEN_FIRST_IO_Model.m
but requires and input data to be updated first.
Documentation of procedure at https://zenodo.org/record/1489943


Step 1:
INPUT DATA:
SWEDEN_FIRST_IO_Model is based on the data prepared in 
1. InputData, 
	- if you update exiobase data, then go to Step 5.
	- if you update swedish data, then go to Step 4.	
2. and the preparation of meta data/linking indicators to stressors in \Sweden Model\MetaData\
	



***********

Any change in stressors/indicators!
You need to make a link between the chosen indicators, the Swedish stressors and the exiobase stressors (you can just extend what is already there):
e.g. in:
MetaData\MatchStressors_SWE.xlsx

And then you need to re-run the update scripts for exiobase:
InputData\prep_exio_data.m
Note - parameters for dataset and year at the top
Note - make sure you have a copy of the exio data from: https://ntnu.box.com/s/75v6quxvuxlpw7fu4roa7lnrfc09jotr



********
Step 4:
If the Sweden data is updated,
you would need to change the data in the Swedish input data:
InputData\IO_SE_08_14_170607.xlsx
And adjust the reading script for the env extensions or other in:
InputData\prep_swe_data.m


***********
Step 5:
If EXIOBASE data is updated
run
InputData\prep_exio_data.m
Note - parameters for dataset and year at the top
Note - make sure you have a copy of the exio data from: https://ntnu.box.com/s/75v6quxvuxlpw7fu4roa7lnrfc09jotr


**********
Re-running results:
And then the model:
SWEDEN_FIRST_IO_Model.m

