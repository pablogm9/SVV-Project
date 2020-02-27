import matplotlib.pyplot as plt

def summation_error(aero_loads, resultant_forces):
    to_return = (sum(resultant_forces)-np.sum(aero_loads))/np.sum(aero_loads)
    return to_return

def magnitude_check(new_aerodata,aero_loads):
    return (1.611*0.505*np.mean(new_aerodata)/np.sum(aero_loads)) < 1


def summed_charts(aero_loads):

    AAsummed = []
    for i in aero_loads:
        AAsummed.append(sum(i))
    aero_loads1 = zip(*reversed(aero_loads))
    AAsummedspan = []
    for k in aero_loads1:
        AAsummedspan.append(sum(k))

    print('this is the verification for the sum of the span wise points')
    plt.figure()
    plt.title('Chord wise summation of forces')
    plt.ylabel('Newtons')
    plt.xlabel('Meters')
    plt.plot(new_nodes_z,AAsummed)
    plt.show()
    print('same as above but for the chord')
    plt.figure()
    plt.title('Span wise summation of forces')
    plt.ylabel('Newtons')
    plt.xlabel('Meters')
    plt.plot(new_nodes_x,AAsummedspan)
    plt.show()