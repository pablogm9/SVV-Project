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
coordinates_a=[12,  13,  15,  36,  37,  38, 196, 197, 198, 199, 200, 246, 247, 248, 249, 250,\
 823, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887, 888, 889, 890,\
 891, 892, 893, 894, 895, 896, 897, 898, 899, 900, 901, 902, 903, 904, 905, 906,\
 907, 908, 909, 910, 911, 948, 956, 957, 958, 959, 960, 961, 962]

coordinates_b=[9,  10,  32,  33,  34,  35, 179, 180, 181, 182, 183, 697, 698, 699, 700, 701,\
 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 722, 730, 731,\
 732, 733, 734, 735, 736, 749, 750, 751, 752, 753, 788, 824, 825, 826, 827, 828,\
 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841]

coordinates_c=[5,   6,  21,  22,  25,  28, 104, 105, 106, 107, 108, 390, 391, 392, 393, 394,\
 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 460, 478, 479,\
 480, 481, 482, 483, 484, 542, 543, 544, 545, 546, 558, 584, 585, 586, 587, 588,\
 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601]

coordinates_d=[2,   3,  16,  17,  20,  23,  45,  46,  47,  48,  49, 285, 286, 287, 288, 289,\
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


#Bending displacemnt
displacement_2 = pd.read_csv('B737.rpt',skiprows=20073,nrows=6588).iloc[:,0][:6588].str.split(expand=True)[[0,1,2,3,4]]
displacement_2 = displacement_2.rename(columns = {0:'label',1:'magnitude',2:'x',3:'y',4:'z'})
displacement_2['x'] = displacement_2['x'].astype('float')
displacement_2['y'] = displacement_2['y'].astype('float')
displacement_2['z'] = displacement_2['z'].astype('float')
deformation_bending = displacement_2[['x','y','z']] + aileron_coords[['x','y','z']]

#Jam Bent Displacement
Jam_Bent_displ = pd.read_csv('B737.rpt',skiprows=26723,nrows=6588).iloc[:,0][:6588].str.split(expand=True)[[0,1,2,3,4]]

Jam_Bent_displ = Jam_Bent_displ.rename(columns = {0:'label',1:'magnitude',2:'x',3:'y',4:'z'})
Jam_Bent_displ['x'] = Jam_Bent_displ['x'].astype('float')
Jam_Bent_displ['y'] = Jam_Bent_displ['y'].astype('float')
Jam_Bent_displ['z'] = Jam_Bent_displ['z'].astype('float')
Jam_Bent_deformation = Jam_Bent_displ[['x','y','z']] + aileron_coords[['x','y','z']]

#Jam Straight Displacment
Jam_Straight_displ = pd.read_csv('B737.rpt',skiprows=33373,nrows=6588 ).iloc[:,0][:6588].str.split(expand=True)[[0,1,2,3,4]]

Jam_Straight_displ = Jam_Straight_displ.rename(columns = {0:'label',1:'magnitude',2:'x',3:'y',4:'z'})
Jam_Straight_displ['x'] = Jam_Straight_displ['x'].astype('float')
Jam_Straight_displ['y'] = Jam_Straight_displ['y'].astype('float')
Jam_Straight_displ['z'] = Jam_Straight_displ['z'].astype('float')
Jam_Straight_deformation = Jam_Straight_displ[['x','y','z']] + aileron_coords[['x','y','z']]
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
#
# ax.scatter(deformations['z'],deformations['x'],deformations['y'], c = 'r', marker = 'o')





#Bending Von Mises Stresses
Bending_R1stresses = pd.read_csv('B737.rpt', skiprows=19,nrows=5778).iloc[:,0][:5578].str.split(expand=True)[[0,1,2,3,4,5]]
Bending_R2stresses = pd.read_csv('B737.rpt', skiprows=5815,nrows=856).iloc[:,0][:856].str.split(expand=True)[[0,1,2,3,4,5]]
Jam_Straight_stress =Bending_R1stresses.append(Bending_R2stresses)
Jam_Straight_stress = Jam_Straight_stress.rename(columns = {0:'label',1:'Int',2:'Mises1',3:'Mises2',4:'S121',5:'S122'})
Jam_Straight_stress['label'] = Jam_Straight_stress['label'].astype('int')
Jam_Straight_stress['Mises1'] = Jam_Straight_stress['Mises1'].astype('float')
Jam_Straight_stress['Mises2'] = Jam_Straight_stress['Mises2'].astype('float')
Jam_Straight_stress['S121'] = Jam_Straight_stress['S121'].astype('float')
Jam_Straight_stress['S122'] = Jam_Straight_stress['S122'].astype('float')



#Jam-Bent Von Mises Stresses
JamBent_R1stresses = pd.read_csv('B737.rpt', skiprows=6704,nrows=5778).iloc[:,0][:5578].str.split(expand=True)[[0,1,2,3,4,5]]
JamBent_R2stresses = pd.read_csv('B737.rpt', skiprows=12500,nrows=856).iloc[:,0][:856].str.split(expand=True)[[0,1,2,3,4,5]]
#Jam-Straight Von Mises Stresses
JamStraight_R1stresses = pd.read_csv('B737.rpt', skiprows=13389,nrows=5778).iloc[:,0][:5578].str.split(expand=True)[[0,1,2,3,4,5]]
JamStraight_R2stresses = pd.read_csv('B737.rpt',skiprows=19187,nrows=856).iloc[:,0][:856].str.split(expand=True)[[0,1,2,3,4,5]]


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
fig = plt.figure(figsize=(16,8))
ax = fig.add_subplot(221)
ax.set_title('Aileron LE deformation for Bending Case')
ax.scatter(aileron_coords['x'].iloc[LE_coords],aileron_coords['y'].iloc[LE_coords], c='g')
ax.scatter(deformation_bending['x'].iloc[LE_coords],deformation_bending['y'].iloc[LE_coords], c='r')
ax.grid()
ax = fig.add_subplot(222)
ax.set_title('Aileron LE Deformation for Jam_Bent Case')
ax.scatter(aileron_coords['x'].iloc[LE_coords],aileron_coords['y'].iloc[LE_coords], c='g')
ax.scatter(Jam_Bent_deformation['x'].iloc[LE_coords],Jam_Bent_deformation['y'].iloc[LE_coords], c='r')
ax.grid()
ax = fig.add_subplot(223)
ax.set_title('Aileron LE deformation Jam_Straight Case')
ax.scatter(aileron_coords['x'].iloc[LE_coords],aileron_coords['y'].iloc[LE_coords], c='g')
ax.scatter(Jam_Straight_deformation['x'].iloc[LE_coords],Jam_Straight_deformation['y'].iloc[LE_coords], c='r')
ax.grid()
# ax = fig.add_subplot(224)
# ax.set_title('Aileron deformation at rib D for Jam_Straight Case')
# ax.scatter(aileron_coords['x'].iloc[LE_coords],aileron_coords['y'].iloc[LE_coords], c='g')
# ax.scatter(Jam_Straight_deformation['x'].iloc[LE_coords],Jam_Straight_deformation['y'].iloc[LE_coords], c='r')
# ax.grid()
plt.show()
#fig.savefig('Rib_Deformation_Bending.eps')

print(Jam_Straight_stress.tail())
Jam_Straight_stress.sort_values(by= ['label'], inplace=True)

print(Jam_Straight_stress.tail())
