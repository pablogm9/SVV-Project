import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
'''
CONTENTS OF THIS FILE: This file is used to interpret the B737.rpt and B737.inp 
'''
#TODO: Combine the stresses from the two regions into one pandas df
#      Average the inside and outside stresses (Mises1 and Mises2)
#      Plot the average stresses according to the x and z values
df = pd.read_table('B737.inp',skiprows=8)
aileron_coords = df['*Node'][:6588].str.split(',',expand = True)
aileron_coords=aileron_coords.rename(columns={0:'label',1:'x', 2:'y', 3:'z'})
aileron_coords['x'] = aileron_coords['x'].astype('float')
aileron_coords['y'] = aileron_coords['y'].astype('float')
aileron_coords['z'] = aileron_coords['z'].astype('float')
coordinates_riba=[ 12,  13,  15,  36,  37,  38, 196, 197, 198, 199, 200, 246, 247, 248, 249, 250,\
823, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887, 888, 889, 890,\
891, 892, 893, 894, 895, 896, 897, 898, 899, 900, 901, 902, 903, 904, 905, 906,\
907, 908, 909, 910, 911, 948, 956, 957, 958, 959, 960, 961, 962]

coordinates_ribb=[9,  10,  32,  33,  34,  35, 179, 180, 181, 182, 183, 697, 698, 699, 700, 701,\
702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 722, 730, 731,\
732, 733, 734, 735, 736, 749, 750, 751, 752, 753, 788, 824, 825, 826, 827, 828,\
829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841]

coordinates_ribc=[5,   6,  21,  22,  25,  28, 104, 105, 106, 107, 108, 390, 391, 392, 393, 394,\
395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 460, 478, 479,\
480, 481, 482, 483, 484, 542, 543, 544, 545, 546, 558, 584, 585, 586, 587, 588,\
589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601]

coordinates_ribd=[2,   3,  16,  17,  20,  23,  45,  46,  47,  48,  49, 285, 286, 287, 288, 289,\
290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 325, 326, 327,\
328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 455,\
459, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477]
LE_coords = [1, 2, 5, 7, 9, 11, 12, 43, 44, 57, 58, 59, 60, 61, 62, 63,\
64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,\
80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,\
96, 97, 98, 99, 100, 101, 102, 103, 156, 157, 158, 159, 160, 161, 173, 174,\
175, 176, 177, 178, 190, 191, 192, 193, 194, 195, 212, 213, 214, 215, 216, 217,\
218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233,\
234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245]
TE_coords = [
16,  19,  22,  30,  33,  37,  40, 323, 324, 408, 409, 410, 411, 412, 413, 414,\
415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430,\
431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446,\
447, 448, 449, 450, 451, 452, 453, 454, 602, 603, 604, 605, 606, 607, 715, 716,\
717, 718, 719, 720, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853,\
854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869,\
870, 871, 872, 873, 874, 875, 936, 937, 938, 939, 940, 941]

def deformed_coords(displacement):
    displacement = displacement.rename(columns = {0:'label',1:'magnitude',2:'x',3:'y',4:'z'})
    displacement['x'] = displacement['x'].astype('float')
    displacement['y'] = displacement['y'].astype('float')
    displacement['z'] = displacement['z'].astype('float')
    deformation = displacement[['x', 'y', 'z']] + aileron_coords[['x', 'y', 'z']]
    return deformation
#Bending displacemnt
displacement_bending = pd.read_csv('B737.rpt',skiprows=20073,nrows=6588).iloc[:,0][:6588].str.split(expand=True)[[0,1,2,3,4]]
deformation_bending = deformed_coords(displacement_bending)
#Jam Bent Displacement
Jam_Bent_displ = pd.read_csv('B737.rpt',skiprows=26723,nrows=6588).iloc[:,0][:6588].str.split(expand=True)[[0,1,2,3,4]]
Jam_Bent_deformation = deformed_coords(Jam_Bent_displ)
#Jam Straight Displacment
Jam_Straight_displ = pd.read_csv('B737.rpt',skiprows=33373,nrows=6588 ).iloc[:,0][:6588].str.split(expand=True)[[0,1,2,3,4]]
Jam_Straight_deformation = deformed_coords(Jam_Straight_displ)

#--------------------------------------STRESSES------------------------------------------
def average_stress_per_rib(stress):
    #For the average stresses
    stress = stress.set_index('label')
    #Average for points at rib a
    VMstresses_riba = stress.iloc[coordinates_riba]
    VMstresses_riba['VM_AV'] = (VMstresses_riba['Mises1']+VMstresses_riba['Mises2'])/2
    VMstresses_riba['S12_AV'] = (VMstresses_riba['S121']+VMstresses_riba['S122'])/2
    #Average for points at rib b
    VMstresses_ribb = stress.iloc[coordinates_ribb]
    VMstresses_ribb['VM_AV'] = (VMstresses_ribb['Mises1']+VMstresses_ribb['Mises2'])/2
    VMstresses_ribb['S12_AV'] = (VMstresses_ribb['S121']+VMstresses_ribb['S122'])/2
    #average for points at rib c
    VMstresses_ribc = stress.iloc[coordinates_ribc]
    VMstresses_ribc['VM_AV'] = (VMstresses_ribc['Mises1']+VMstresses_ribc['Mises2'])/2
    VMstresses_ribc['S12_AV'] = (VMstresses_ribc['S121']+VMstresses_ribc['S122'])/2
    #Average for points at rib d
    VMstresses_ribd = stress.iloc[coordinates_ribc]
    VMstresses_ribd['VM_AV'] = (VMstresses_ribd['Mises1']+VMstresses_ribd['Mises2'])/2
    VMstresses_ribd['S12_AV'] = (VMstresses_ribd['S121']+VMstresses_ribd['S122'])/2
    return VMstresses_riba,VMstresses_ribb,VMstresses_ribc,VMstresses_ribd
def VM_stress_whole(R1stresses,R2stresses):
    stress = pd.concat([R1stresses, R2stresses], ignore_index=True)
    stress = stress.rename(columns={0: 'label', 1: 'Int', 2: 'Mises1', 3: 'Mises2', 4: 'S121', 5: 'S122'})
    stress['label'] = stress['label'].astype('int')
    stress['Mises1'] = stress['Mises1'].astype('float')
    stress['Mises2'] = stress['Mises2'].astype('float')
    stress['S121'] = stress['S121'].astype('float')
    stress['S122'] = stress['S122'].astype('float')
    return stress
#--------------Bending-----------------

#Bending Von Mises Stresses
Bending_R1stresses = pd.read_csv('B737.rpt', skiprows=19,nrows=5778).iloc[:,0][:5578].str.split(expand=True)[[0,1,2,3,4,5]]
Bending_R2stresses = pd.read_csv('B737.rpt', skiprows=5815,nrows=856).iloc[:,0][:856].str.split(expand=True)[[0,1,2,3,4,5]]
Bending_stress = VM_stress_whole(Bending_R1stresses, Bending_R2stresses)
Bending_VMstresses_riba,Bending_VMstresses_ribb,Bending_VMstresses_ribc,Bending_VMstresses_ribd = average_stress_per_rib(Bending_stress)

#--------------Jam-Bent----------------
#Jam-Bent Von Mises Stresses
JamBent_R1stresses = pd.read_csv('B737.rpt', skiprows=6704,nrows=5778).iloc[:,0][:5578].str.split(expand=True)[[0,1,2,3,4,5]]
JamBent_R2stresses = pd.read_csv('B737.rpt', skiprows=12500,nrows=856).iloc[:,0][:856].str.split(expand=True)[[0,1,2,3,4,5]]
JamBent_stress = VM_stress_whole(JamBent_R1stresses, JamBent_R2stresses)
JamBent_VMstresses_riba,JamBent_VMstresses_ribb,JamBent_VMstresses_ribc,JamBent_VMstresses_ribd = average_stress_per_rib(JamBent_stress)
#------------Jam-Straight-----------------------
#Jam-Straight Von Mises Stresses
JamStraight_R1stresses = pd.read_csv('B737.rpt', skiprows=13389,nrows=5778).iloc[:,0][:5578].str.split(expand=True)[[0,1,2,3,4,5]]
JamStraight_R2stresses = pd.read_csv('B737.rpt',skiprows=19187,nrows=856).iloc[:,0][:856].str.split(expand=True)[[0,1,2,3,4,5]]
JamStraight_stress =pd.concat([JamStraight_R1stresses,JamStraight_R2stresses], ignore_index= True)
JamStraight_stress = JamStraight_stress.rename(columns = {0:'label',1:'Int',2:'Mises1',3:'Mises2',4:'S121',5:'S122'})
#JamStraight_stress['label'] = JamStraight_stress['label'].astype('int')
JamStraight_stress['Mises1'] = JamStraight_stress['Mises1'].astype('float')
JamStraight_stress['Mises2'] = JamStraight_stress['Mises2'].astype('float')
JamStraight_stress['S121'] = JamStraight_stress['S121'].astype('float')
JamStraight_stress['S122'] = JamStraight_stress['S122'].astype('float')
JamStraight_VMstresses_riba,JamStraight_VMstresses_ribb,JamStraight_VMstresses_ribc,JamStraight_VMstresses_ribd = average_stress_per_rib(JamStraight_stress)

#fig = plt.figure()
#plt.plot(list(Jam_Straight_displ['label']),list(Jam_Straight_displ['Mises1']))

#plt.show()
#ax.scatter(Jam_Straight_displ['label'],Jam_Straight_displ['Mises1'], c = 'r', marker = 'o')

# fig = plt.figure(figsize=(16,8))
# ax = fig.add_subplot(111,projection='3d')
# ax.set_title('Aileron Deformation Jam Bent')
# ax.scatter3D(Jam_Bent_deformation['z'],Jam_Bent_deformation['x'],Jam_Bent_deformation['y'])
# ax.set_xlabel('Z axis')
# ax.set_ylabel('X axis')
# ax.set_zlabel('Y axis')
# plt.show()
# fig = plt.figure(figsize=(16,8))
# ax = fig.add_subplot(221)
# ax.set_title('Aileron LE deformation for Bending Case')
# ax.scatter(aileron_coords['x'].iloc[LE_coords],aileron_coords['y'].iloc[LE_coords], c='g')
# ax.scatter(deformation_bending['x'].iloc[LE_coords],deformation_bending['y'].iloc[LE_coords], c='r')
# ax.grid()
# ax = fig.add_subplot(222)
# ax.set_title('Aileron LE Deformation for Jam_Bent Case')
# ax.scatter(aileron_coords['x'].iloc[LE_coords],aileron_coords['y'].iloc[LE_coords], c='g')
# ax.scatter(Jam_Bent_deformation['x'].iloc[LE_coords],Jam_Bent_deformation['y'].iloc[LE_coords], c='r')
# ax.grid()
# ax = fig.add_subplot(223)
# ax.set_title('Aileron LE deformation Jam_Straight Case')
# ax.scatter(aileron_coords['x'].iloc[LE_coords],aileron_coords['y'].iloc[LE_coords], c='g')
# ax.scatter(Jam_Straight_deformation['x'].iloc[LE_coords],Jam_Straight_deformation['y'].iloc[LE_coords], c='r')
# ax.grid()
# ax = fig.add_subplot(224)
# ax.set_title('Aileron deformation at rib D for Jam_Straight Case')
# ax.scatter(aileron_coords['x'].iloc[LE_coords],aileron_coords['y'].iloc[LE_coords], c='g')
# ax.scatter(Jam_Straight_deformation['x'].iloc[LE_coords],Jam_Straight_deformation['y'].iloc[LE_coords], c='r')
# ax.grid()
#plt.show()




'''IMPORTS:
    The pandas dataframes imported at the top of this file are all imported from the input interpreter file. The imports are
    made up of one set of deformed displacement points for each load case provided. Also imported are dataframes that include
    the stresses at each rib of each load case provided.
    '''
#Find the error between two sets of dispplacement data, the aileron coord displacements
#for ribs for each case do that for LE AND TE sets
#compare stresses at the ribs between the numerical values
#This should be the node positions AFTER the bending force is applied

aileron_coords = aileron_coords.set_index('label')
#Coordinates used by the numerical that are also used by the validation model
riba_coords = aileron_coords.iloc[coordinates_riba]
ribb_coords = aileron_coords.iloc[coordinates_ribb]
ribc_coords = aileron_coords.iloc[coordinates_ribc]
ribd_coords = aileron_coords.iloc[coordinates_ribd]


#Here is where you input the data the numerical model creates
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
#------------------------------------DEFORMATION ERROR--------------------------------

def error_calculator(Numer_model, FEM_model):
    return (FEM_model-Numer_model)/FEM_model

#Bending Deformation Error
deformation_bending['Deformation Error x'] = error_calculator(Numer_model_def_bending['x'], deformation_bending['x'])
deformation_bending['Deformation Error y'] = error_calculator(Numer_model_def_bending['y'], deformation_bending['y'])
deformation_bending['Deformation Error z'] = error_calculator(Numer_model_def_bending['z'], deformation_bending['z'])

#JamBent Deformation Error
Jam_Bent_deformation['Deformation Error x'] = error_calculator(Numer_model_def_JamBent['x'], Jam_Bent_deformation['x'])
Jam_Bent_deformation['Deformation Error y'] = error_calculator(Numer_model_def_JamBent['y'], Jam_Bent_deformation['y'])
Jam_Bent_deformation['Deformation Error z'] = error_calculator(Numer_model_def_JamBent['z'], Jam_Bent_deformation['z'])

#JamStraight Deformation Error
Jam_Straight_deformation['Deformation Error x'] = error_calculator(Numer_model_def_JamStraight['x'], Jam_Straight_deformation['x'])
Jam_Straight_deformation['Deformation Error y'] = error_calculator(Numer_model_def_JamStraight['y'], Jam_Straight_deformation['y'])
Jam_Straight_deformation['Deformation Error z'] = error_calculator(Numer_model_def_JamStraight['z'], Jam_Straight_deformation['z'])


#--------------------------------------STRESS ERROR---------------------------------------
#Bending Stress Error
Bending_VMstresses_riba['VM Stress Error A'] = error_calculator(Numer_model_bending_riba_stress_df['Von_Mises_riba'], Bending_VMstresses_riba['VM_AV'])
Bending_VMstresses_ribb['VM Stress Error B'] = error_calculator(Numer_model_bending_ribb_stress_df['Von_Mises_ribb'], Bending_VMstresses_ribb['VM_AV'])
Bending_VMstresses_ribc['VM Stress Error C'] = error_calculator(Numer_model_bending_ribc_stress_df['Von_Mises_ribc'], Bending_VMstresses_ribc['VM_AV'])
Bending_VMstresses_ribd['VM Stress Error D'] = error_calculator(Numer_model_bending_ribd_stress_df['Von_Mises_ribd'], Bending_VMstresses_ribd['VM_AV'])

#JamBent Error
JamBent_VMstresses_riba['VM Stress Error A'] = error_calculator(Numer_model_JamBent_riba_stress_df['Von_Mises_riba'], JamBent_VMstresses_riba['VM_AV'])
JamBent_VMstresses_ribb['VM Stress Error B'] = error_calculator(Numer_model_JamBent_ribb_stress_df['Von_Mises_ribb'], JamBent_VMstresses_ribb['VM_AV'])
JamBent_VMstresses_ribc['VM Stress Error C'] = error_calculator(Numer_model_JamBent_ribc_stress_df['Von_Mises_ribc'], JamBent_VMstresses_ribc['VM_AV'])
JamBent_VMstresses_ribd['VM Stress Error D'] = error_calculator(Numer_model_JamBent_ribd_stress_df['Von_Mises_ribd'], JamBent_VMstresses_ribd['VM_AV'])

#JamStraight Error
JamStraight_VMstresses_riba['VM Stress Error A'] = error_calculator(Numer_model_JamStraight_riba_stress_df['Von_Mises_riba'], JamStraight_VMstresses_riba['VM_AV'])
JamStraight_VMstresses_ribb['VM Stress Error B'] = error_calculator(Numer_model_JamStraight_ribb_stress_df['Von_Mises_ribb'], JamStraight_VMstresses_ribb['VM_AV'])
JamStraight_VMstresses_ribc['VM Stress Error C'] = error_calculator(Numer_model_JamStraight_ribc_stress_df['Von_Mises_ribc'], JamStraight_VMstresses_ribc['VM_AV'])
JamStraight_VMstresses_ribd['VM Stress Error D'] = error_calculator(Numer_model_JamStraight_ribd_stress_df['Von_Mises_ribd'], JamStraight_VMstresses_ribd['VM_AV'])



fig=plt.figure(figsize=(16,8))
ax = fig.add_subplot(221)
colorset=JamStraight_VMstresses_riba['VM_AV']
ax.set_title('Stresses at rib A for the JamStraight Deformation')

rib_a=ax.scatter(aileron_coords['z'].iloc[coordinates_riba], aileron_coords['y'].iloc[coordinates_riba], c=colorset)
colorplot = plt.colorbar(rib_a)
colorplot.ax.set_ylabel('avg Von Mises Stress (kN/$mm^2$)', rotation=90)
ax.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

ax = fig.add_subplot(222)
ax.set_title('Aileron stresses at rib B for the JamStraight Deformation')
colorset=JamStraight_VMstresses_ribb['VM_AV']
rib_b=ax.scatter(aileron_coords['z'].iloc[coordinates_ribb], aileron_coords['y'].iloc[coordinates_ribb], c=colorset)
colorplot = plt.colorbar(rib_b)
colorplot.ax.set_ylabel('avg Von Mises Stress (kN/$mm^2$)', rotation=90)
ax.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

ax = fig.add_subplot(223)
ax.set_title('Aileron stresses at rib C for the JamStraight Deformation')
colorset=JamStraight_VMstresses_ribc['VM_AV']
rib_c=ax.scatter(aileron_coords['z'].iloc[coordinates_ribc], aileron_coords['y'].iloc[coordinates_ribc], c=colorset)
colorplot = plt.colorbar(rib_c)
colorplot.ax.set_ylabel('avg Von Mises Stress (kN/$mm^2$)', rotation=90)
ax.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

ax = fig.add_subplot(224)
ax.set_title('Aileron stresses at rib D for the JamStraight Deformation')
colorset=JamStraight_VMstresses_ribd['VM_AV']
rib_d=ax.scatter(aileron_coords['z'].iloc[coordinates_ribd], aileron_coords['y'].iloc[coordinates_ribd], c=colorset)
colorplot = plt.colorbar(rib_d)
colorplot.ax.set_ylabel('avg Von Mises Stress (kN/$mm^2$)', rotation=90)
ax.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
plt.plot()
plt.show()

