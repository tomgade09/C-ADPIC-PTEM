#include "include\Simulation170925.h"

void Simulation170925::loadDataFilesIntoParticleArray()
{//not super generic, but ok - the function is virtual
	std::string importdir{ rootdir_m };
	importdir = importdir + "\\in\\data\\";
	std::vector<std::vector<std::string>> files;
	files.resize(2);
	files[0].resize(3);
	files[1].resize(3);
	files[0][0] = "e_vpara.bin";
	files[0][1] = "e_vperp.bin";
	files[0][2] = "e_z.bin";
	files[1][0] = "i_vpara.bin";
	files[1][1] = "i_vperp.bin";
	files[1][2] = "i_z.bin";

	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		for (int jjj = 0; jjj < numberOfAttributesTracked_m; jjj++)
			loadFileIntoParticleAttribute(iii, jjj, importdir.c_str(), files[iii][jjj].c_str());
	}
}

void Simulation170925::convertVPerpToMu()
{
	if (mu_m)
	{
		std::cout << "v_perp has already been converted to mu.  Returning with no changes.\n";
		return;
	}

	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		for (int jjj = 0; jjj < numberOfParticlesPerType_m; jjj++)
			particles_m[iii][1][jjj] = 0.5 * mass_m[iii] * particles_m[iii][1][jjj] * particles_m[iii][1][jjj] / BFieldatZ(particles_m[iii][2][jjj], simTime_m);
	}
	mu_m = true;
}

void Simulation170925::convertMuToVPerp()
{
	if (!mu_m)
	{
		std::cout << "Quantity is v_perp (not mu).  Returning with no changes.\n";
		return;
	}

	for (int iii = 0; iii < numberOfParticleTypes_m; iii++)
	{
		for (int jjj = 0; jjj < numberOfParticlesPerType_m; jjj++)
			particles_m[iii][1][jjj] = sqrt(2 * particles_m[iii][1][jjj] * BFieldatZ(particles_m[iii][2][jjj], simTime_m) / mass_m[iii]);
	}
	
	mu_m = false;
}

void Simulation170925::prepareResults()
{
	convertMuToVPerp();
	serializeParticleArray();
	resultsPrepared_m = true;
}