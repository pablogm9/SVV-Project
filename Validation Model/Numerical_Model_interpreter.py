from SVV_input_interpreter import aileron_coords,deformation_bending,Jam_Bent_deformation, Jam_Straight_deformation,\
Bending_VMstresses_riba,Bending_VMstresses_ribb,Bending_VMstresses_ribc,Bending_VMstresses_ribd,\
JamBent_VMstresses_riba,JamBent_VMstresses_ribb,JamBent_VMstresses_ribc,JamBent_VMstresses_ribd,\
JamStraight_VMstresses_riba,JamStraight_VMstresses_ribb,JamStraight_VMstresses_ribc,JamStraight_VMstresses_ribd
import pandas as pd
'''IMPORTS:
    The pandas dataframes imported at the top of this file are all imported from the input interpreter file. The imports are
    made up of one set of deformed displacement points for each load case provided. Also imported are dataframes that include
    the stresses at each rib of each load case provided.
    '''
#Find the error between two sets of dispplacement data, the aileron coord displacements
#for ribs for each case do that for LE AND TE sets
#compare stresses at the ribs between the numerical values
#This should be the node positions AFTER the bending force is applied
Numer_model_def_bending = {'x':[0,0,0,0,0], 'y':[0,0,0,0,0], 'z': [0,0,0,0]}
Numer_model_def_JamBent = {'x':[0,0,0,0,0], 'y':[0,0,0,0,0], 'z': [0,0,0,0]}
Numer_model_def_JamStraight = {'x':[0,0,0,0,0], 'y':[0,0,0,0,0], 'z': [0,0,0,0]}




