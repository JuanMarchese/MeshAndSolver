
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import math 

def calculate_modulus(element):
    return math.sqrt((float(element[0])**2.0) + (float(element[1])**2.0))

def plot_force_defor(simulado_defor,simulado_force,z_axis = True):
    NCURVES = len(simulado_defor.keys()) + 1
    values = range(NCURVES)
    
    plt.figure(figsize=(20, 10), dpi=50)

    jet = cm = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    i = 0
    colorVal = scalarMap.to_rgba(values[i])
    
    

    for simulado_key in sorted(simulado_defor.keys()):
        i += 1
        deform = simulado_defor[simulado_key][1]
        force = simulado_force[simulado_key][1]
        v = simulado_defor[simulado_key][0]

        x_axis = [abs(element[1]) for element in force]
        
        if z_axis:
            y_axis = [element[0] for element in deform]
        else:
            y_axis = [element[1] for element in deform]
        
        colorVal = scalarMap.to_rgba(values[i])
        plt.plot(x_axis,y_axis, marker="o", color=colorVal, label=simulado_key)

        for k in range(0,len(v)):

            plt.annotate(
                str(v[k]) + " kV/m", 
                xy = (x_axis[k], y_axis[k]), xytext = (-10, 10),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5))
        

    plt.ylabel('Radio nuevo / Radio original')
    plt.xlabel('[dyn]')

    plt.show()
    
import csv
import math

def get_paper_force_deformation():
    z = [1, 1.1559594097, 1.2175455646, 1.2738143261, 1.4195829382]
    r = [1, 0.9732103237, 0.8212460524, 0.7866042061, 0.7625862052]
    force = [0, 7.00E-07, 1.40E-06, 2.50E-06, 5.60E-06]

    return [force, z, r]

def calculate_alpha(x_med, y_med):
    if y_med >= 0:
        if x_med <= 0:
            alpha = math.pi * 0.5
        else:
            alpha = math.atan(y_med / x_med)
    else:
        if x_med <= 0:
            alpha = -1.0 * math.pi * 0.5
        else:
            alpha = math.atan(y_med / x_med)

    return alpha

def force_equation(row):
    result = [0, 0]

    z_med = float(row[1])
    r_med = float(row[2])

    radio = math.sqrt((z_med - 50.0)**2 + r_med**2)
    alpha = calculate_alpha(r_med, z_med - 50.0)
    t_zz = float(row[3])
    t_rr = float(row[4])

    lado = float(row[7])

    # force r
    result[0] = ((t_zz * math.cos(alpha)) - (t_rr * math.sin(alpha))) * 2 * math.pi * r_med * lado
    
    # force z
    result[1] = ((t_zz * math.sin(alpha)) + (t_rr * math.sin(alpha))) * 2 * math.pi * r_med * lado
    
    
    
    return result

def read_force(file_path):
    reader_force = open(file_path, 'r')

    position_force = 1
    fuerza = 0
    i = 0
    
    for line in reader_force:
        if i == position_force:
            row = line.rstrip('\n').split("\t")
            row = [float(element) for element in row]
            fuerza = force_equation(row)
            break
        i += 1

        
    reader_force.close()
    
    return fuerza

def file_force(file_path,out_file_path):
    reader_force = open(file_path, 'r')
    out_file = open(out_file_path, 'w')
    
    fuerza = 0
   
    for line in reader_force:
        row = line.rstrip('\n').split("\t")
        row = [float(element) for element in row]
        fuerza = force_equation(row)
        out_file.write(str(row[0]) + "\t" + str(fuerza[1]) + "\t" + str(fuerza[0]) + "\n")
        


    out_file.close()
    reader_force.close()
    
    return fuerza



def parse_file_name(file_name):
    parameters = file_name.split("--")
    dict_key = parameters[3] + "__" + parameters[4] + "__" + parameters[5]
    dict_key = dict_key
    density = float(parameters[0].replace("V", ""))
    
    return [dict_key,density]

def sort_dictionaries(dictionary, start_x, start_y):
    dictionary_result = {}
    
    for key in dictionary.keys():
        result = sorted(dictionary[key], key=lambda point: point[0])
        x = []
        y = []

        x.append(start_x)
        y.append(start_y)
    
        for point in result:
            x.append(point[0])
            y.append(point[1])
            
        dictionary_result[key] = [x, y]
        
    return dictionary_result
    

def get_simulation_force(path, values_file, max_density = -1.0):
    dictionary = {}
   

    f = open(path + values_file, 'r')

    for value in f:
        
        force_file_path = path + "Demo-" + value.rstrip() + "-001.tcsv"

        dict_key,density = parse_file_name(value)
        
        if not (dict_key in dictionary):
            dictionary[dict_key] = []

        if density < max_density or max_density < 0:
            dictionary[dict_key].append([density, read_force(force_file_path)])
                        
   
    f.close()

    dictionary_result = sort_dictionaries(dictionary, 0, [0,0])
    f.close()

    return dictionary_result

def get_simulation_defor(path, values_file, max_density = -1.0):
    dictionary = {}
    dictionary_result = {}
    

    f = open(path + values_file, 'r')

    for value in f:
        
        defor_file_path = path + "Demo-" + value.rstrip() + "-001.dcsv"
        
        dict_key,density = parse_file_name(value)
        
        defor_file = open(defor_file_path, 'r')

        line = defor_file.readline()
        line = defor_file.readline()
        line = defor_file.readline()
        line = defor_file.readline()
        
        columns = line.split('\t')
        defor_z = (float(columns[3])-50.0) / (float(columns[1])-50.0)
        
        for i in range(0,7):
            line = defor_file.readline()

        columns = line.split('\t')
        defor_r = float(columns[4]) / float(columns[2])
        
        if not (dict_key in dictionary):
            dictionary[dict_key] = []

        if density < max_density or max_density < 0:
            dictionary[dict_key].append([density, [defor_z,defor_r]])

        defor_file.close()
    
    
    dictionary_result = sort_dictionaries(dictionary, 0, [1,1])
    f.close()

    return dictionary_result
    
