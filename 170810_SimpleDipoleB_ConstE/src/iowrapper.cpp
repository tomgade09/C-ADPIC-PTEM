#include "include\iowrapper.h"

void writeParticlesToBin(double*** particles, std::string directory)
{
	std::string file1, file2, file3, file4, file5, file6;
	file1 = directory + "/e_vpara.bin";
	file2 = directory + "/e_vperp.bin";
	file3 = directory + "/e_z.bin";
	file4 = directory + "/i_vpara.bin";
	file5 = directory + "/i_vperp.bin";
	file6 = directory + "/i_z.bin";
	
	fileIO::writeDblBin(file1.c_str(), particles[0][0], NUMPARTICLES);
	fileIO::writeDblBin(file2.c_str(), particles[0][1], NUMPARTICLES);
	fileIO::writeDblBin(file3.c_str(), particles[0][2], NUMPARTICLES);
	fileIO::writeDblBin(file4.c_str(), particles[1][0], NUMPARTICLES);
	fileIO::writeDblBin(file5.c_str(), particles[1][1], NUMPARTICLES);
	fileIO::writeDblBin(file6.c_str(), particles[1][2], NUMPARTICLES);
}

double*** readParticlesFromBin(std::string directory)
{
	std::string file1, file2, file3, file4, file5, file6;
	file1 = directory + "/e_vpara.bin";
	file2 = directory + "/e_vperp.bin";
	file3 = directory + "/e_z.bin";
	file4 = directory + "/i_vpara.bin";
	file5 = directory + "/i_vperp.bin";
	file6 = directory + "/i_z.bin";

	double* e_vpara = fileIO::readDblBin(file1.c_str(), NUMPARTICLES);
	double* e_vperp = fileIO::readDblBin(file2.c_str(), NUMPARTICLES);
	double* e_z = fileIO::readDblBin(file3.c_str(), NUMPARTICLES);
	double* i_vpara = fileIO::readDblBin(file4.c_str(), NUMPARTICLES);
	double* i_vperp = fileIO::readDblBin(file5.c_str(), NUMPARTICLES);
	double* i_z = fileIO::readDblBin(file6.c_str(), NUMPARTICLES);

	double** electrons = new double*[4];
	double** ions = new double*[4];

	double*** particles = new double**[2];

	electrons[0] = e_vpara;
	electrons[1] = e_vperp;
	electrons[2] = e_z;
	electrons[3] = nullptr;

	ions[0] = i_vpara;
	ions[1] = i_vperp;
	ions[2] = i_z;
	ions[3] = nullptr;

	particles[0] = electrons;
	particles[1] = ions;

	return particles;
}