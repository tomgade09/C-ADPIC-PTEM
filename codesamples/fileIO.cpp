#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

#define DEBUGSTDOUTPUT
#define DLLFILE

#ifdef DLLFILE
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT
#endif

//Eigen::VectorXd readbin(const int& dimensions, const std::vector<double>& dimensionArray);
DLLEXPORT double* readDblBin(const std::string& filename, unsigned int numOfDblsToRead)
{
	std::ifstream binFile;
	binFile.open(filename, std::ios::in | std::ios::binary);

	if (!binFile.is_open())
	{
		std::cout << "Warning: Could not open file " << filename << "\n";
		return nullptr;
	}

	double* dataArray = new double[numOfDblsToRead];
	binFile.read(reinterpret_cast<char*>(dataArray), std::streamsize(numOfDblsToRead * sizeof(double)));

	binFile.close();

	return dataArray;
}

DLLEXPORT void clrDataMemory(double* dataArray)
{
	delete[] dataArray;
}

DLLEXPORT unsigned long indexLambdaFcn(int xsteps, int ysteps, int zsteps, int xsize=(1+10/0.04), int ysize=(1+10/0.04), int zsize=(1+10/0.04)){
	return 3 * (zsteps * xsize * ysize + ysteps * xsize + xsteps); };

//Eigen::Vector3d interpolateB3D( const Eigen::Vector3d& pos, const Eigen::Vector2d& xrng, const Eigen::Vector2d& yrng,
//								const Eigen::Vector2d& zrng, const double& binSize, const Eigen::VectorXd& dataArray );
Eigen::Vector3d interpolateB3D( const Eigen::Vector3d& pos, const Eigen::Vector2d& xrng, const Eigen::Vector2d& yrng,
								const Eigen::Vector2d& zrng, const double& binSize, const double* dataArray )
{
	int xsize, ysize, zsize;
	xsize = static_cast<int>(1 + (xrng[1] - xrng[0]) / binSize);
	ysize = static_cast<int>(1 + (yrng[1] - yrng[0]) / binSize);
	zsize = static_cast<int>(1 + (zrng[1] - zrng[0]) / binSize);

	unsigned long nmDbls{ (static_cast<unsigned long>(xsize) * ysize * zsize * 3) };

	Eigen::Vector3i stepsFromMin{ static_cast<int>(floor((pos[0] - xrng[0]) / binSize)),
								  static_cast<int>(floor((pos[1] - yrng[0]) / binSize)), 
								  static_cast<int>(floor((pos[2] - zrng[0]) / binSize)) };
	Eigen::Vector3d binLessThanAbsPos{ xrng[0] + stepsFromMin[0] * binSize,
									   yrng[0] + stepsFromMin[1] * binSize,
									   zrng[0] + stepsFromMin[2] * binSize };
#ifdef DEBUGSTDOUTPUT
	std::cout << "In inner interpolate.  pos, steps from Min, binLTAbsPos : \n" << pos[0] << ", " << pos[1] << ", " << pos[2] << "\n";
	std::cout << stepsFromMin[0] << ", " << stepsFromMin[1] << ", " << stepsFromMin[2] << "\n";
	std::cout << binLessThanAbsPos[0] << ", " << binLessThanAbsPos[1] << ", " << binLessThanAbsPos[2] << "\n";
#endif
	
	if (stepsFromMin[0] > xsize || stepsFromMin[0] < 0)
	{
		std::cout << "X position is outside of range of calculated values.  Returning B = 0.\n";
		return{ 0.0, 0.0, 0.0 };
	}
	if (stepsFromMin[1] > ysize || stepsFromMin[1] < 0)
	{
		std::cout << "Y position is outside of range of calculated values.  Returning B = 0.\n";
		return{ 0.0, 0.0, 0.0 };
	}
	if (stepsFromMin[2] > zsize || stepsFromMin[2] < 0)
	{
		std::cout << "Z position is outside of range of calculated values.  Returning B = 0.\n";
		return{ 0.0, 0.0, 0.0 };
	}

	/*auto indexLambdaFcn = [=](int xsteps, int ysteps, int zsteps) -> unsigned int {
		return 3 * (zsteps * xsize * ysize + ysteps * xsize + xsteps); };*/

	std::vector<unsigned int> indicies;
/*	indicies = { //calculate indicies of B points located all around test charge (8 total)
		indexLambdaFcn(stepsFromMin[0], stepsFromMin[1], stepsFromMin[2]), //x, y, z
		indexLambdaFcn(stepsFromMin[0] + 1, stepsFromMin[1], stepsFromMin[2]), //x+1, y, z
		indexLambdaFcn(stepsFromMin[0], stepsFromMin[1] + 1, stepsFromMin[2]), //x, y+1, z
		indexLambdaFcn(stepsFromMin[0], stepsFromMin[1], stepsFromMin[2] + 1), //x, y, z+1
		indexLambdaFcn(stepsFromMin[0] + 1, stepsFromMin[1] + 1, stepsFromMin[2]), //x+1, y+1, z
		indexLambdaFcn(stepsFromMin[0] + 1, stepsFromMin[1], stepsFromMin[2] + 1), //x+1, y, z+1
		indexLambdaFcn(stepsFromMin[0], stepsFromMin[1] + 1, stepsFromMin[2] + 1), //x, y+1, z+1			//don't know if I like the below change yet
		indexLambdaFcn(stepsFromMin[0] + 1, stepsFromMin[1] + 1, stepsFromMin[2] + 1) }; //x+1, y+1, z+1*/  //but should fix issues with around zero being way off

	indicies = { //calculate indicies of B points located all around test charge (8 total)
		indexLambdaFcn(stepsFromMin[0] - 1, stepsFromMin[1] - 1, stepsFromMin[2] - 1), //x-1, y-1, z-1
		indexLambdaFcn(stepsFromMin[0] + 1, stepsFromMin[1] - 1, stepsFromMin[2] - 1), //x+1, y-1, z-1
		indexLambdaFcn(stepsFromMin[0] - 1, stepsFromMin[1] + 1, stepsFromMin[2] - 1), //x-1, y+1, z-1
		indexLambdaFcn(stepsFromMin[0] - 1, stepsFromMin[1] - 1, stepsFromMin[2] + 1), //x-1, y-1, z+1
		indexLambdaFcn(stepsFromMin[0] + 1, stepsFromMin[1] + 1, stepsFromMin[2] - 1), //x+1, y+1, z-1
		indexLambdaFcn(stepsFromMin[0] + 1, stepsFromMin[1] - 1, stepsFromMin[2] + 1), //x+1, y-1, z+1
		indexLambdaFcn(stepsFromMin[0] - 1, stepsFromMin[1] + 1, stepsFromMin[2] + 1), //x-1, y+1, z+1
		indexLambdaFcn(stepsFromMin[0] + 1, stepsFromMin[1] + 1, stepsFromMin[2] + 1) }; //x+1, y+1, z+1
	
#ifdef DEBUGSTDOUTPUT
	std::cout << "indicies calculated: ";
	for (int iii = 0; iii < 8; ++iii)
		std::cout << indicies[iii] << ", ";
	std::cout << "\n";
#endif

	Eigen::VectorXd distFromPos(8);
/*	distFromPos << //calculate scalar distance of each B point from position of test charge
		std::sqrt(pow((binLessThanAbsPos[0] - pos[0]), 2) + pow((binLessThanAbsPos[1] - pos[1]), 2) + pow((binLessThanAbsPos[2] - pos[2]), 2)), //x, y, z
		std::sqrt(pow((binSize + binLessThanAbsPos[0] - pos[0]), 2) + pow((binLessThanAbsPos[1] - pos[1]), 2) + pow((binLessThanAbsPos[2] - pos[2]), 2)), //x+1, y, z
		std::sqrt(pow((binLessThanAbsPos[0] - pos[0]), 2) + pow((binSize + binLessThanAbsPos[1] - pos[1]), 2) + pow((binLessThanAbsPos[2] - pos[2]), 2)), //x, y+1, z
		std::sqrt(pow((binLessThanAbsPos[0] - pos[0]), 2) + pow((binLessThanAbsPos[1] - pos[1]), 2) + pow((binSize + binLessThanAbsPos[2] - pos[2]), 2)), //x, y, z+1
		std::sqrt(pow((binSize + binLessThanAbsPos[0] - pos[0]), 2) + pow((binSize + binLessThanAbsPos[1] - pos[1]), 2) + pow((binLessThanAbsPos[2] - pos[2]), 2)), //x+1, y+1, z
		std::sqrt(pow((binSize + binLessThanAbsPos[0] - pos[0]), 2) + pow((binLessThanAbsPos[1] - pos[1]), 2) + pow((binSize + binLessThanAbsPos[2] - pos[2]), 2)), //x+1, y, z+1
		std::sqrt(pow((binLessThanAbsPos[0] - pos[0]), 2) + pow((binSize + binLessThanAbsPos[1] - pos[1]), 2) + pow((binSize + binLessThanAbsPos[2] - pos[2]), 2)), //x, y+1, z+1
		std::sqrt(pow((binSize + binLessThanAbsPos[0] - pos[0]), 2) + pow((binSize + binLessThanAbsPos[1] - pos[1]), 2) + pow((binSize + binLessThanAbsPos[2] - pos[2]), 2)); */ //x+1, y+1, z+1

	distFromPos << //calculate scalar distance of each B point from position of test charge
		std::sqrt(pow((binLessThanAbsPos[0] - binSize - pos[0]), 2) + pow((binLessThanAbsPos[1] - binSize - pos[1]), 2) + pow((binLessThanAbsPos[2] - binSize - pos[2]), 2)), //x-1, y-1, z-1
		std::sqrt(pow((binSize + binLessThanAbsPos[0] - pos[0]), 2) + pow((binLessThanAbsPos[1] - binSize - pos[1]), 2) + pow((binLessThanAbsPos[2] - binSize - pos[2]), 2)), //x+1, y-1, z-1
		std::sqrt(pow((binLessThanAbsPos[0] - binSize - pos[0]), 2) + pow((binSize + binLessThanAbsPos[1] - pos[1]), 2) + pow((binLessThanAbsPos[2] - binSize - pos[2]), 2)), //x-1, y+1, z-1
		std::sqrt(pow((binLessThanAbsPos[0] - binSize - pos[0]), 2) + pow((binLessThanAbsPos[1] - binSize - pos[1]), 2) + pow((binSize + binLessThanAbsPos[2] - pos[2]), 2)), //x-1, y-1, z+1
		std::sqrt(pow((binSize + binLessThanAbsPos[0] - pos[0]), 2) + pow((binSize + binLessThanAbsPos[1] - pos[1]), 2) + pow((binLessThanAbsPos[2] - binSize - pos[2]), 2)), //x+1, y+1, z-1
		std::sqrt(pow((binSize + binLessThanAbsPos[0] - pos[0]), 2) + pow((binLessThanAbsPos[1] - binSize - pos[1]), 2) + pow((binSize + binLessThanAbsPos[2] - pos[2]), 2)), //x+1, y-1, z+1
		std::sqrt(pow((binLessThanAbsPos[0] - binSize - pos[0]), 2) + pow((binSize + binLessThanAbsPos[1] - pos[1]), 2) + pow((binSize + binLessThanAbsPos[2] - pos[2]), 2)), //x-1, y+1, z+1
		std::sqrt(pow((binSize + binLessThanAbsPos[0] - pos[0]), 2) + pow((binSize + binLessThanAbsPos[1] - pos[1]), 2) + pow((binSize + binLessThanAbsPos[2] - pos[2]), 2)); //x+1, y+1, z+1

#ifdef DEBUGSTDOUTPUT
	std::cout << "dist (r vectors) calculated: ";
	for (int iii = 0; iii < 8; ++iii)
		std::cout << distFromPos[iii] << ", ";
	std::cout << "\n";
#endif
	
	double distSum = distFromPos.sum();
	Eigen::VectorXd bxVec(8);
	Eigen::VectorXd byVec(8);
	Eigen::VectorXd bzVec(8);

	for (int iii = 0; iii < 8; iii++)
	{
		if (indicies[iii] > nmDbls)
		{
			std::cout << "Error: Index out of range!!  Writing zeros for B corresponding to indicies that are out of range: " << indicies[iii] << "\n";
			indicies[iii] = 0;
			distFromPos[iii] = 0.0;
			bxVec[iii] = 0.0;
			byVec[iii] = 0.0;
			bzVec[iii] = 0.0;
		}
		else if (indicies[iii] < 0)
		{
			std::cout << "Error: Index less than 0!!  Writing zeros for B corresponding to indicies that are out of range: " << indicies[iii] << "\n";
			indicies[iii] = 0;
			distFromPos[iii] = 0.0;
			bxVec[iii] = 0.0;
			byVec[iii] = 0.0;
			bzVec[iii] = 0.0;
		}
		else
		{
			bxVec[iii] = dataArray[indicies[iii]];
			byVec[iii] = dataArray[indicies[iii] + 1];
			bzVec[iii] = dataArray[indicies[iii] + 2];
		}
	}

	return{ bxVec.sum() - bxVec.dot(distFromPos) / distSum,
			byVec.sum() - byVec.dot(distFromPos) / distSum,
			bzVec.sum() - bzVec.dot(distFromPos) / distSum };
}

/*Eigen::Vector3d interpolateB3D(const Eigen::Vector3d& pos, const Eigen::Vector2d& xrng, const Eigen::Vector2d& yrng,
	const Eigen::Vector2d& zrng, const double& binSize, const Eigen::VectorXd& dataArray)*/
DLLEXPORT double* interpolateB3DPyWrapper(double* pos, double* xrng, double* yrng, double* zrng, double binSize, double* dataArray)
{
	Eigen::Vector3d p;
	p = { pos[0], pos[1], pos[2] };

#ifdef DEBUGSTDOUTPUT
	std::cout << "On way in\n x,x  y,y  z,z: ";
	for (int iii = 0; iii < 3; iii++)
		std::cout << p[iii] << "  " << pos[iii] << "\n";
#endif
	
	Eigen::Vector2d x, y, z;
	x = { xrng[0], xrng[1] };
	y = { yrng[0], yrng[1] };
	z = { zrng[0], zrng[1] };

	Eigen::Vector3d result;
	result = interpolateB3D(p, x, y, z, binSize, dataArray);
	
#ifdef DEBUGSTDOUTPUT
	std::cout << "interpolate called (on way out).  B value returned: " << "\n" << result << "\n\n\n";
#endif

	double* res = new double[3];

	for (int iii = 0; iii < 3; ++iii)
		res[iii] = result[iii];

	return res;
}