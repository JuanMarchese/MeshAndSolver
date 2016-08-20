#include "MeshLogic/Data.h"
#include "MeshLogic/CellGenerator.h"
#include "MeshLogic/MeshGenerator.h"
#include "Solver/FortranWrapper.h"
#include "Problem/ConfigurationReader.h"
#include "Problem/ParameterReader.h"
#include <iostream>


const index_t membraneLayers = 2;

const std::string inner_radius = "inner_radius";
const std::string outer_radius = "outer_radius";
const std::string cell_iterations = "cell_iterations";
const std::string sigma_inner_cell = "sigma_inner_cell";
const std::string sigma_outer_cell = "sigma_outer_cell";
const std::string sigma_membrane = "sigma_membrane";

const std::string membrane_layers = "membrane_layers";
const std::string potential_density = "density";
const std::string max_iterations = "max_iterations";
const std::string tolerance = "tolerance";
const std::string output_path = "output_path";

const std::string funyoung_in = "funyoung_in";
const std::string funyoung_membrane = "funyoung_membrane";
const std::string funyoung_out = "funyoung_out";

const std::string funposs_in = "funposs_in";
const std::string funposs_membrane = "funposs_membrane";
const std::string funposs_out = "funposs_out";

const std::string deformation_rep = "deformation_rep";
const std::string file_sufix = "file_sufix";

/* Parameters position:
 *
 * 0 - executable
 * 1 - command
 *
 * 2 - inner radius (float)
 * 3 - outer radius (float)
 * 4 - cell iterations (int) must be mult of 4
 * 5 - sigma inner (float)
 * 6 - sigma outer (float)
 * 7 - sigma membrane (float)
 * 8 - potential density (float)
 * 9 - iterations (int)
 * 10- tolerance (float)
 * 11- output path (string)
 * 12- file prefix (string - optional)
 * */


/* Parameters position:
 *
 * 0 - executable
 * 1 - command
 *
 * 2 - path
 */



void readFromFile(char *argv[],
				  Data& data,
				  CellGenerator& generator,
				  MeshGenerator& generatorSquare,
				  FortranWrapper& solver){

	ConfigurationReader config;
	config.init("config.txt",'#','=');
	config.process();

    data.setPotentialDensity(std::stod(config.getParameter(potential_density,"1.0")));

    data.setSigmas(std::stod(config.getParameter(sigma_inner_cell,"1.0")),
    		 	   std::stod(config.getParameter(sigma_outer_cell,"1.0")),
    		 	   std::stod(config.getParameter(sigma_membrane,"0.0001")));

    data.setPath(config.getParameter(output_path,"./output/"));

    generator.initCellGeometrics(std::stod(config.getParameter(inner_radius,"10.0")),
    							 std::stod(config.getParameter(outer_radius,"10.1")));


    generator.initCellParameters(std::stoi(config.getParameter(membrane_layers,"2")),
    							 std::stoi(config.getParameter(cell_iterations,"40")));

    solver.setParameters(std::stoi(config.getParameter(max_iterations,"10000")),
    					 std::stod(config.getParameter(tolerance,"0.000001")),
    					 std::stoi(config.getParameter(membrane_layers,"2")));


    data.initDeformationIterations(std::stoi(config.getParameter(deformation_rep,"10")));


    data.initDeformationPoss(std::stod(config.getParameter(funyoung_in,"0.0001")),
    						 std::stod(config.getParameter(funyoung_membrane,"0.1")),
    						 std::stod(config.getParameter(funyoung_out,"0.0001")));

    data.initDeformationYoung(std::stod(config.getParameter(funposs_in,"0.49")),
    						  std::stod(config.getParameter(funposs_membrane,"0.49")),
    						  std::stod(config.getParameter(funposs_out,"0.49")));


    data.setFileSufix(config.getParameter(file_sufix,""));

    /*------------------------*/

}


void howToUse(char *argv[],ParameterReader& paramReader){

	std::cout<< "Para utilizar" << std::endl<< std::endl;
	std::cout<< "Parametros" << std::endl;
	std::cout<< argv[0] << " -p [-XX valor1 -YY valor2 -ZZ valor3 ...] "<<std::endl;

	paramReader.printUsage();

	std::cout << std::endl;

	std::cout<< "Archivo configuracion" << std::endl;
	std::cout<< argv[0] << " -f [ruta al archivo] "<<std::endl;

	std::cout << std::endl;

	std::cout<< "Automatico" << std::endl;
	std::cout<< argv[0] << " -a"<<std::endl;
}

int main(int argc, char *argv[])
{
    Data data;
    CellGenerator generator;
    MeshGenerator generatorSquare;
    FortranWrapper solver;
    ParameterReader paramReader;

    int dataMode;

    dataMode = paramReader.process(argc,argv);


	if((dataMode == ParameterReader::error) ||(dataMode == ParameterReader::helpMode)){
		howToUse(argv,paramReader);
		return 1;
	}

	if(dataMode == ParameterReader::parameterMode){
		paramReader.initData(data,generator,generatorSquare,solver);
		std::cout << "La configuracion a usar es: " << std::endl;
		paramReader.printParametersInUse();
	}else{
		if(dataMode == ParameterReader::parameterMode){
			readFromFile(argv,data,generator,generatorSquare,solver);
		}else{
			howToUse(argv,paramReader);
			return 1;
		}
	}


    coordinates_t maxX = 100.0;
    coordinates_t minX = 0.0;
    coordinates_t maxY = 50.0;
    coordinates_t minY = 0.0;

    index_t outerLayers = 10;
    coordinates_t outerLayersSeparation = 1.0;

    index_t innerLayers = 10;
    coordinates_t innerLayersSeparation = (generator.getInnerRadius()/2.25)/innerLayers;

    /*-----------------------------------------------*/
    coordinates_t min_length_x = (generator.getOuterRadius() + (((coordinates_t)outerLayers + 15) * outerLayersSeparation)) * 2.0;
    coordinates_t min_length_y = generator.getOuterRadius() + (((coordinates_t)outerLayers + 15) * outerLayersSeparation);

    if(maxX < min_length_x){
    	maxX = min_length_x;
    }

    if(maxY < min_length_y){
    	maxY = min_length_y;
    }
    /*-----------------------------------------------*/

    data.setBoundaries(minX,maxX,minY,maxY);

    index_t numberOfLayers = 5;
    index_t numberOfFixingLayers = 5;


    ids_t innerCellType = 1;
    ids_t inCellType = 2;
    ids_t outerCellType = 3;


    coordinates_t blockSize = 2.50;

    index_t blockXIterations = generator.getCellIterations() / 2 + 1;
    index_t blockYIterations = generator.getCellIterations() / 4 + 1;

    data.initTypes(innerCellType,outerCellType,inCellType);

    generator.initBoundaryParameters(outerLayers,outerLayersSeparation);
    generator.initInnerParameters(innerLayers,innerLayersSeparation);
    generator.initTypes(inCellType,outerCellType,innerCellType);


    generator.initData(&data);

    std::cout << "Creando la maya de la celula" << std::endl;
    generator.process();

    generatorSquare.initData(&data);
    generatorSquare.initTransitionData(numberOfFixingLayers);
    generatorSquare.setInnerVector(generator.getOuterPoints());
    generatorSquare.setInCellVector(generator.getInnerPoints());
    generatorSquare.initCenterPosition(data.getCellXCenter(),data.getCellYCenter());
    generatorSquare.initBlocksSizes(blockSize,blockXIterations,blockYIterations,numberOfLayers);
    generatorSquare.initTypes(outerCellType,innerCellType);

    std::cout << "Creando el resto de la maya" << std::endl;
    generatorSquare.process();

    std::cout << "Inicializando los modulos en Fortran" << std::endl;
    solver.setData(&data);


    solver.solveProblemsIteratively();

    return 0;
}
