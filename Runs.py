#!/usr/bin/python
import shlex
import sys
from subprocess import call
from joblib import Parallel, delayed
import multiprocessing


class CellParameters:
	def __init__(self):
		self.radius_pos = 0
		self.potential_pos = 1
		self.poisson_pos = 2
		self.young_pos = 3
		self.output_path_pos = 4
		self.radius_difference_pos = 5
		self.cell_iterations_pos = 6
		self.sigma_pos = 7
		self.membrane_layers_pos = 8
		self.sufix_func_pos = 9
		self.param_count = 10
		self.params = []

		for i in range(0,self.param_count):
			self.params.append([])

class CellConfiguration:

	def __init__(self):
		self.radius = 0
		self.potential = 0
		self.poisson_in = 0
		self.poisson_mem = 0
		self.poisson_out = 0
		self.young_in = 0
		self.young_mem = 0
		self.young_out = 0
		self.output_path = ''

		self.name = ''
		
		self.radius_difference = 0
		self.cell_iterations = 0
		self.sigma_inner_cell = 0
		self.sigma_outer_cell = 0
		self.sigma_membrane = 0

		self.tolerance = 0

		self.membrane_layers = 0

		self.sufix_func = 0

def config_cell(parameters, iteration_vector):

	cell = CellConfiguration()

	cell.radius = parameters.params[parameters.radius_pos][iteration_vector[parameters.radius_pos]]
	cell.potential = parameters.params[parameters.potential_pos][iteration_vector[parameters.potential_pos]]

	cell.output_path = parameters.params[parameters.output_path_pos][iteration_vector[parameters.output_path_pos]]

	cell.radius_difference = parameters.params[parameters.radius_difference_pos][iteration_vector[parameters.radius_difference_pos]]
	cell.cell_iterations = parameters.params[parameters.cell_iterations_pos][iteration_vector[parameters.cell_iterations_pos]]

	cell.membrane_layers = parameters.params[parameters.membrane_layers_pos][iteration_vector[parameters.membrane_layers_pos]]
	cell.sufix_func = parameters.params[parameters.sufix_func_pos][iteration_vector[parameters.sufix_func_pos]]

	cell.poisson_in = parameters.params[parameters.poisson_pos][iteration_vector[parameters.poisson_pos]][0]
	cell.poisson_mem = parameters.params[parameters.poisson_pos][iteration_vector[parameters.poisson_pos]][1]
	cell.poisson_out = parameters.params[parameters.poisson_pos][iteration_vector[parameters.poisson_pos]][2]

	cell.young_in = parameters.params[parameters.young_pos][iteration_vector[parameters.young_pos]][0]
	cell.young_mem = parameters.params[parameters.young_pos][iteration_vector[parameters.young_pos]][1]
	cell.young_out = parameters.params[parameters.young_pos][iteration_vector[parameters.young_pos]][2]

	cell.sigma_inner_cell = parameters.params[parameters.sigma_pos][iteration_vector[parameters.sigma_pos]][0]
	cell.sigma_membrane = parameters.params[parameters.sigma_pos][iteration_vector[parameters.sigma_pos]][1]
	cell.sigma_outer_cell = parameters.params[parameters.sigma_pos][iteration_vector[parameters.sigma_pos]][2]
	

	return cell

def pinzas_file_sigma(cell):
	sufix = "Pinzas-%sV-%s-%s-%s" % (
		"{:.1f}".format(cell.potential), "{:.15f}".format(cell.sigma_inner_cell), "{:.15f}".format(cell.sigma_outer_cell),
		"{:.15f}".format(cell.sigma_membrane))
		
	cell.name = "%s-%s-%s" % (
		"{:.15f}".format(cell.sigma_inner_cell), "{:.15f}".format(cell.sigma_outer_cell),
		"{:.15f}".format(cell.sigma_membrane))

	return sufix

def pinzas_file_poisson(cell):
	sufix = "Pinzas-%sV-%s-%s-%s" % (
		"{:.1f}".format(cell.potential), "{:.25f}".format(cell.poisson_in), "{:.25f}".format(cell.poisson_mem),
		"{:.25f}".format(cell.poisson_out))

	cell.name = "%s-%s-%s" % (
		"{:.25f}".format(cell.poisson_in), "{:.25f}".format(cell.poisson_mem),
		"{:.25f}".format(cell.poisson_out))

	return sufix

def pinzas_file_young(cell):
	sufix = "Pinzas-%sV-%s" % (
		"{:.1f}".format(cell.potential),
		"{:.25f}".format(cell.young_mem))
		
	cell.name = "%sV-%s" % (
		"{:.1f}".format(cell.potential),
                "{:.25f}".format(cell.young_mem))

	return sufix

def itv_ancho_potential(cell):
	sufix = "ITV-Ancho-%sum" % ("{:.5f}".format(cell.radius_difference))
	cell.name = "%s" % ("{:.5f}".format(cell.radius_difference))

	return sufix

def itv_radio_potential(cell):
	sufix = "ITV-Radio-%sum" % ("{:.1f}".format(cell.radius))
	cell.name = "%s" % ("{:.1f}".format(cell.radius))

	return sufix

def itv_v_potential(cell):
	sufix = "ITV-V-%s" % ("{:.1f}".format(cell.potential))
	cell.name = "%s" % ("{:.1f}".format(cell.potential))

	return sufix



def pinzas_file(cell):
	sufix = "Pinzas-%sV-%s-%s-%s-%s-%s-%s" % (
	"{:.1f}".format(cell.potential), "{:.10f}".format(cell.poisson_in), "{:.10f}".format(cell.poisson_mem),
	"{:.10f}".format(cell.poisson_out), "{:.10f}".format(cell.young_in), "{:.10f}".format(cell.young_mem),
	"{:.10f}".format(cell.young_out))

	return sufix

def old_file_sufix(cell):
	sufix = "P%s-R%s-PO%s-PM%s-PI%s-YO%s-YM%s-YI%s-1000" % (
	"{:.1f}".format(cell.potential), "{:.1f}".format(cell.radius), "{:.6f}".format(cell.poisson_in), "{:.6f}".format(cell.poisson_mem),
	"{:.6f}".format(cell.poisson_out), "{:.6f}".format(cell.young_in), "{:.6f}".format(cell.young_mem),
	"{:.6f}".format(cell.young_out))

	return sufix


def recursive_initializing(iteration, cells, parameters, iteration_vector):

	if iteration < parameters.param_count:
		for i in range(0, len(parameters.params[iteration])):
			iteration_vector[iteration] = i
			recursive_initializing(iteration + 1, cells, parameters, iteration_vector)

	else:
		cells.append(config_cell(parameters, iteration_vector))


def call_program(cell_config, deformation_rep):
	radius_difference = cell_config.radius_difference
	cell_iterations = cell_config.cell_iterations

	sigma_inner_cell =	cell_config.sigma_inner_cell
	sigma_membrane = cell_config.sigma_membrane
	sigma_outer_cell = cell_config.sigma_outer_cell

	ground = 0
	max_iterations = 1000000
	tolerance = 0.000000000001

	output_path = cell_config.output_path

	membrane_layers = cell_config.membrane_layers

	funyoung_in = cell_config.young_in
	funyoung_membrane = cell_config.young_mem
	funyoung_out = cell_config.young_out

	funposs_in = cell_config.poisson_in
	funposs_membrane = cell_config.poisson_mem
	funposs_out = cell_config.poisson_out

	radius = cell_config.radius
	potential = cell_config.potential
	# variacion de potencial

	# print "-----------------------Iterando-----------------------"
	inner_radius = radius - radius_difference;
	outer_radius = radius;
	fileSufix = cell_config.sufix_func(cell_config)
	command = "/home/juan/workspaceMars/MeshAndSolver/MeshAndSolver -p -ir %f -or %f -ci %d -ml %d -is %.25f -os %.25f -ms %.25f -v %f -g %f -i %d -t %.25f -dr %d -path %s -sufix %s -iy %.25f -my %.25f -oy %.25f -ip %.25f -mp %.25f -op %.25f " %	(inner_radius ,outer_radius ,cell_iterations, membrane_layers ,sigma_inner_cell ,sigma_outer_cell ,sigma_membrane ,potential ,ground ,max_iterations , tolerance, deformation_rep, output_path ,fileSufix,funyoung_in,funyoung_membrane,funyoung_out,funposs_in,funposs_membrane,funposs_out)
	
	call(shlex.split(command))
	#print(command)
	# print "---------------------------------------------------------------\n"

	
# 0 in
# 1 out
# 2 mem

def file_fuerza(cell):
        sufix = "F-V%s-kVm" % ("{:.1f}".format(cell.potential))
        cell.name = "V%s-kVm" % ("{:.1f}".format(cell.potential))

        return sufix

def file_alta_frecuencia(cell):
	sufix = "Pins-V%s--R%s--M%s--S%s-%s-%s--P%s-%s-%s--Y%s-%s-%s" % (
		"{:.1f}".format(cell.potential),"{:.1f}".format(cell.radius),"{:.5f}".format(cell.radius_difference),
		"{:.9f}".format(cell.sigma_inner_cell), "{:.9f}".format(cell.sigma_outer_cell),"{:.15f}".format(cell.sigma_membrane),
		"{:.7f}".format(cell.poisson_in), "{:.7f}".format(cell.poisson_out),"{:.7f}".format(cell.poisson_mem),
		"{:.9f}".format(cell.young_in), "{:.9f}".format(cell.young_out),"{:.9f}".format(cell.young_mem))
		
	cell.name = "V%s--R%s--M%s--S%s-%s-%s--P%s-%s-%s--Y%s-%s-%s" % (
		"{:.1f}".format(cell.potential),"{:.1f}".format(cell.radius),"{:.5f}".format(cell.radius_difference),
		"{:.9f}".format(cell.sigma_inner_cell), "{:.9f}".format(cell.sigma_outer_cell),"{:.15f}".format(cell.sigma_membrane),
		"{:.7f}".format(cell.poisson_in), "{:.7f}".format(cell.poisson_out),"{:.7f}".format(cell.poisson_mem),
		"{:.9f}".format(cell.young_in), "{:.9f}".format(cell.young_out),"{:.9f}".format(cell.young_mem))

	return sufix

def file_pinzas_opticas(cell):
        sufix = "Pins-V%s--R%s--M%s--Y%s-%s-%s" % (
                "{:.1f}".format(cell.potential),"{:.1f}".format(cell.radius),"{:.5f}".format(cell.radius_difference),
                "{:.9f}".format(cell.young_in), "{:.9f}".format(cell.young_out),"{:.9f}".format(cell.young_mem))

        cell.name = "V%s--R%s--M%s--Y%s-%s-%s" % (
                "{:.1f}".format(cell.potential),"{:.1f}".format(cell.radius),"{:.5f}".format(cell.radius_difference),
                "{:.9f}".format(cell.young_in), "{:.9f}".format(cell.young_out),"{:.9f}".format(cell.young_mem))

        return sufix

def file_pinzas_opticasPosion(cell):
        sufix = "Pins-V%s--R%s--M%s--P%s-%s-%s" % (
                "{:.1f}".format(cell.potential),"{:.1f}".format(cell.radius),"{:.5f}".format(cell.radius_difference),
                "{:.7f}".format(cell.poisson_in), "{:.7f}".format(cell.poisson_out),"{:.7f}".format(cell.poisson_mem))

        cell.name = "V%s--R%s--M%s--P%s-%s-%s" % (
                "{:.1f}".format(cell.potential),"{:.1f}".format(cell.radius),"{:.5f}".format(cell.radius_difference),
                "{:.7f}".format(cell.poisson_in), "{:.7f}".format(cell.poisson_out),"{:.7f}".format(cell.poisson_mem))

        return sufix

def file_pinzas_opticasSigma(cell):
        sufix = "Pins-V%s--R%s--M%s--S%s-%s-%s" % (
                "{:.1f}".format(cell.potential),"{:.1f}".format(cell.radius),"{:.5f}".format(cell.radius_difference),
                "{:.9f}".format(cell.sigma_inner_cell), "{:.9f}".format(cell.sigma_outer_cell),"{:.15f}".format(cell.sigma_membrane))

        cell.name = "V%s--R%s--M%s--S%s-%s-%s" % (
                "{:.1f}".format(cell.potential),"{:.1f}".format(cell.radius),"{:.5f}".format(cell.radius_difference),
                "{:.9f}".format(cell.sigma_inner_cell), "{:.9f}".format(cell.sigma_outer_cell),"{:.15f}".format(cell.sigma_membrane))

        return sufix


def read_parameters(cell_conf):
	
	for parameter in sys.argv:
		str_values = str(parameter).split(';')
		
		if str_values[0] == "v":		#potencial
			cell_conf.params[cell_conf.potential_pos].append(float(str_values[1]))
		
		elif str_values[0] == "r":		#radio celular
			cell_conf.params[cell_conf.radius_pos].append(float(str_values[1]))
			
		elif str_values[0] == "m":	#acho membrana
			cell_conf.params[cell_conf.radius_difference_pos].append(float(str_values[1]))
			
		elif str_values[0] == "p":		#poisson
			cell_conf.params[cell_conf.poisson_pos].append([float(str_values[1]), float(str_values[2]), float(str_values[3])])
		
		elif str_values[0] == "y":		#young
			cell_conf.params[cell_conf.young_pos].append([float(str_values[1]), float(str_values[2]), float(str_values[3])])
		
		elif str_values[0] == "s":		#sigma
			cell_conf.params[cell_conf.sigma_pos].append([float(str_values[1]), float(str_values[2]), float(str_values[3])])
	
		elif str_values[0] == "path":		#ruta
			cell_conf.params[cell_conf.output_path_pos].append(str_values[1])
	
if __name__ == '__main__':
	print("Comenzando procesamiento batch paralelo\n")

	cell_conf = CellParameters()

	cell_conf.params[cell_conf.cell_iterations_pos].append(400)
	cell_conf.params[cell_conf.membrane_layers_pos].append(1)

	cell_conf.params[cell_conf.sufix_func_pos].append(pinzas_file_young)
	read_parameters(cell_conf);

	iteration_vector = []
	cells = []

	for i in range(0,cell_conf.param_count):
		iteration_vector.append(0)

	recursive_initializing(0,cells,cell_conf,iteration_vector)

	for cell in cells:
		cell.sufix_func(cell)
	
	with open(cell_conf.params[cell_conf.output_path_pos][0] + "Values.txt", 'a') as outfile:
		for cell in cells:
			outfile.write(cell.name + '\n')
		
	
	num_cores = 8
	num_defor_iterations = 2
	results = Parallel(n_jobs=num_cores)(delayed(call_program)(current_cell,num_defor_iterations) for current_cell in cells)

	print("Fianlizado procesamiento")
