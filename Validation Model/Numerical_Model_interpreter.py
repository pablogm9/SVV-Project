from SVV_input_interpreter import aileron_coords,deformation_bending,Jam_Bent_deformation, Jam_Straight_deformation,\
Bending_VMstresses_riba,Bending_VMstresses_ribb,Bending_VMstresses_ribc,Bending_VMstresses_ribd,\
JamBent_VMstresses_riba,JamBent_VMstresses_ribb,JamBent_VMstresses_ribc,JamBent_VMstresses_ribd,\
JamStraight_VMstresses_riba,JamStraight_VMstresses_ribb,JamStraight_VMstresses_ribc,JamStraight_VMstresses_ribd,\
coordinates_riba, coordinates_ribb, coordinates_ribc, coordinates_ribd,LE_coords, TE_coords
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

#For the Bending case ribs
Numer_model_bending_riba_stress = {'Von_Mises_riba' :[0,0,0,0]}
Numer_model_bending_riba_stress_df = pd.DataFrame(Numer_model_bending_riba_stress, columns = ['Von_Mises_riba'])
Numer_model_bending_ribb_stress = {'Von_Mises_ribb' :[0,0,0,0]}
Numer_model_bending_ribb_stress_df = pd.DataFrame(Numer_model_bending_ribb_stress, columns = ['Von_Mises_ribb'])
Numer_model_bending_ribc_stress = {'Von_Mises_ribc' :[0,0,0,0]}
Numer_model_bending_ribc_stress_df = pd.DataFrame(Numer_model_bending_ribc_stress, columns = ['Von_Mises_ribc'])
Numer_model_bending_ribd_stress = {'Von_Mises_ribd' :[0,0,0,0]}
Numer_model_bending_ribd_stress_df = pd.DataFrame(Numer_model_bending_ribd_stress, columns = ['Von_Mises_ribd'])
#For the JamBent case ribs
Numer_model_JamBent_riba_stress = {'Von_Mises_riba' :[0,0,0,0]}
Numer_model_JamBent_riba_stress_df = pd.DataFrame(Numer_model_JamBent_riba_stress, columns = ['Von_Mises_riba'])
Numer_model_JamBent_ribb_stress = {'Von_Mises_ribb' :[0,0,0,0]}
Numer_model_JamBent_ribb_stress_df = pd.DataFrame(Numer_model_JamBent_ribb_stress, columns = ['Von_Mises_ribb'])
Numer_model_JamBent_ribc_stress = {'Von_Mises_ribc' :[0,0,0,0]}
Numer_model_JamBent_ribc_stress_df = pd.DataFrame(Numer_model_JamBent_ribc_stress, columns = ['Von_Mises_ribc'])
Numer_model_JamBent_ribd_stress = {'Von_Mises_ribd' :[0,0,0,0]}
Numer_model_JamBent_ribd_stress_df = pd.DataFrame(Numer_model_JamBent_ribd_stress, columns = ['Von_Mises_ribd'])

#For the JamStraight case ribs
Numer_model_JamStraight_riba_stress = {'Von_Mises_riba' :[0,0,0,0]}
Numer_model_JamStraight_riba_stress_df = pd.DataFrame(Numer_model_JamStraight_riba_stress, columns = ['Von_Mises_riba'])
Numer_model_JamStraight_ribb_stress = {'Von_Mises_ribb' :[0,0,0,0]}
Numer_model_JamStraight_ribb_stress_df = pd.DataFrame(Numer_model_JamStraight_ribb_stress, columns = ['Von_Mises_ribb'])
Numer_model_JamStraight_ribc_stress = {'Von_Mises_ribc' :[0,0,0,0]}
Numer_model_JamStraight_ribc_stress_df = pd.DataFrame(Numer_model_JamStraight_ribc_stress, columns = ['Von_Mises_ribc'])
Numer_model_JamStraight_ribd_stress = {'Von_Mises_ribd' :[0,0,0,0]}
Numer_model_JamStraight_ribd_stress_df = pd.DataFrame(Numer_model_JamStraight_ribd_stress, columns = ['Von_Mises_ribd'])

#For error calculations
#Bending
deformation_bending['Deformation Error x'] = (deformation_bending['x']-Numer_model_def_bending['x'])/(deformation_bending['x'])
deformation_bending['Deformation Error y'] = (deformation_bending['y']-Numer_model_def_bending['y'])/(deformation_bending['y'])
deformation_bending['Deformation Error z'] = (deformation_bending['z']-Numer_model_def_bending['z'])/(deformation_bending['z'])




#JamBent
Jam_Bent_deformation['Deformation Error x'] = (Jam_Bent_deformation['x']-Numer_model_def_JamBent['x'])/(Jam_Bent_deformation['x'])
Jam_Bent_deformation['Deformation Error y'] = (Jam_Bent_deformation['y']-Numer_model_def_JamBent['y'])/(Jam_Bent_deformation['y'])
Jam_Bent_deformation['Deformation Error z'] = (Jam_Bent_deformation['z']-Numer_model_def_JamBent['z'])/(Jam_Bent_deformation['z'])



#JamStraight
Jam_Straight_deformation['Deformation Error x'] = (Jam_Straight_deformation['x']-Numer_model_def_JamStraight['x'])/(Jam_Straight_deformation['x'])
Jam_Straight_deformation['Deformation Error y'] = (Jam_Straight_deformation['y']-Numer_model_def_JamStraight['y'])/(Jam_Straight_deformation['y'])
Jam_Straight_deformation['Deformation Error z'] = (Jam_Straight_deformation['z']-Numer_model_def_JamStraight['z'])/(Jam_Straight_deformation['z'])

aileron_coords = aileron_coords.set_index('label')
riba_coords = aileron_coords.iloc(coordinates_riba)
ribb_coords = aileron_coords.iloc(coordinates_ribb)
ribc_coords = aileron_coords.iloc(coordinates_ribc)
ribd_coords = aileron_coords.iloc(coordinates_ribd)

print(riba_coords)

