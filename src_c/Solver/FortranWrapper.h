/*
 * FortranWrapper.h
 *
 *  Created on: 16/11/2014
 *      Author: juan
 */

#ifndef FORTRANWRAPPER_H_
#define FORTRANWRAPPER_H_

#include "../MeshLogic/Data.h"

extern"C" {
	/*Fortran functions*/
	void allocate_memory_(int* number_nodes,int* nodes_per_elem,int* number_elements);
	void de_allocate_memory_();
	void de_allocate_external_memory_();
	void poisson_();
	void print_memory_();
	void deforma_();


	/*Fortran modules variables*/
	extern double* __def_solver_MOD_solucion;

	extern int __def_variables_MOD_nnodes;
	extern int __def_variables_MOD_nelements;
	extern int __def_variables_MOD_nodpel;
	extern int __def_variables_MOD_nodpel2;

	extern int* __def_variables_MOD_conect;
	extern int* __def_variables_MOD_vec_poten;
	extern int* __def_variables_MOD_vec_tierra;
	extern int* __def_variables_MOD_vec_quieto;
	extern int* __def_variables_MOD_material;
	extern double* __def_variables_MOD_coor_x;
	extern double* __def_variables_MOD_coor_y;
	extern double* __def_variables_MOD_tension;

	extern double __def_variables_MOD_potencial;
	extern double __def_variables_MOD_sigmaext;
	extern double __def_variables_MOD_sigmamem;
	extern double __def_variables_MOD_sigmaint;
	extern double __def_variables_MOD_tierra;

	extern double __def_variables_MOD_funyoung_in;
	extern double __def_variables_MOD_funyoung_membrane;
	extern double __def_variables_MOD_funyoung_out;

	extern double __def_variables_MOD_funposs_in;
	extern double __def_variables_MOD_funposs_membrane;
	extern double __def_variables_MOD_funposs_out;


	extern double __def_constantes_MOD_toler;
	extern int __def_constantes_MOD_itermax;
}


class FortranWrapper {
private:
	Data* data;
	ids_t maxIterations;
	coordinates_t tolerance;
	ids_t cellLayers;


public:
	FortranWrapper();
	virtual ~FortranWrapper();

	void setData(Data* data);
	void setParameters(ids_t maxIterations,coordinates_t tolerance,ids_t cellLayers);
	void initSolver();
	void solveElectromagnecticField();
	void solveDeformation();
	void destroySolver();

	void solveProblemsIteratively();



private:
	void copyCoordinatesData();
	void copyQuadsData();
	void copyBoundaryData();
	void copyMaterialData();
};

#endif /* FORTRANWRAPPER_H_ */
