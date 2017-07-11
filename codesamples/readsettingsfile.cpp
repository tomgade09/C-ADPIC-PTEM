#include <fstream>
#include <iostream>
#include <vector>
#include <string>

int main(int argc, char *argv[])
{
	//ensure program is being used right
	if (argc != 2)
	{
		std::cout << "Error: Wrong number of arguments!\n";
		std::cout << "Usage: " << argv[0] << " <csv config file>";
		return 0;
	}
	
	//open settings csv file
	std::ifstream settings(argv[1]);
	if (!settings.is_open())
	{
		std::cout << "Could not open file: " << argv[1];
		return 0;
	}
	settings.ignore(200, '\n'); //ignore characters until first line break - leaves room for title row

	//setup variables here
	double simTime{ 0.0 };
	double electronsInSim{ 0.0 };
	double a{ 0.0 };
	double timeLimit{ 0.0 };
	double dt{ 0.0 };

	int index{ 0 };
	std::vector<double*> parameters{&electronsInSim, &a, &timeLimit, &dt};//add pointers to variables

	while (true)
	{//reads settings from a csv settings file, seems to work well
		errno = 0;
		char in[50];
		settings.getline(in, 50, ',');
		if (settings.gcount() == 0) //if no more settings, break loop
			break;
		double d{ strtod(in, nullptr) };
		if (errno != 0)
		{
			std::cout << "Error reading settings file.  Setting is not a double, contains invalid characters, or is out of range for a double: " << in;
			return 0;
		}
		*parameters[index] = d;
		++index;
		//check for too many settings/not enough?? or just trust the user/config file?
	}

	std::cout << electronsInSim << "\n" << a << "\n" << dt << "\n" << timeLimit;

	/*
	std::vector<double> particles{ generateParticleDistribution(argv[1]) }; //need to write this function
	
	while (simtime < timelimit)
	{
		for (unsigned long part = 0, part < (sizeof(particles), ++part) //up to 4.29 billion particles, should be good
			update_r_v(particle[part]); //need to write this function
		simtime += dt;
	}
	*/
	return 0;
}