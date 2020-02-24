#include "utils/validation/PTEM.h"
#include "ionosphere/ionosphere.h"
#include "Simulation/Simulation.h"
#include "utils/numerical.h"

using ionosphere::dEFlux::newPA;

namespace validation
{
	bool PTEM_dipB_noE(double maxError, bool printError)
	{
		//create sim with small particle set from 0-180 deg, 10^(0.5-4.5)eV, iterate particles for a bit
		
		//compare distribution with "ideal" using newPA and generate error

		//print error if desired

		//return whether or not maximum error of particle set is greater than maxError

		//potentially do this for both DipoleB and DipoleBLUT

		return false;
	}

	bool PTEM_dipB_QSPS(double maxError, bool printError)
	{
		//create sim with QSPS and small particle set from 0-180 deg, 10^(0.5-4.5)eV, iterate particles for a bit
		
		//compare distribution with "ideal" using newPA to QSPS top, adding E over QSPS, then newPA from bottom and generate error
		
		//print error if desired
		
		//return whether or not maximum error of particle set is greater than maxError

		//potentially do this for both DipoleB and DipoleBLUT

		return false;
	}

	bool PTEM_dist_noE(double maxError, string simDataDir, bool printError)
	{
		Simulation sim(simDataDir);
		//if (sim.Emodel()->size() != 0) throw logic_error("validation::PTEM_dist_noE: sim loaded has a non-zero E Field.  Results will not be accurate.");

		Particle*  particle{ sim.particle(0) };
		BModel*    bmodel{ sim.Bmodel() };

		ratio maxPAerr{ 0.0 };
		ratio maxEerr{ 0.0 };

		const vector<vector<double>>& orig{ particle->data(true) };
		const vector<vector<double>>& curr{ particle->data(false) };
		vector<vector<double>> origEPA(2);
		vector<vector<double>> currEPA(2);
		utils::numerical::v2DtoEPitch(orig.at(0), orig.at(1), MASS_ELECTRON, origEPA.at(0), origEPA.at(1));
		utils::numerical::v2DtoEPitch(curr.at(0), curr.at(1), MASS_ELECTRON, currEPA.at(0), currEPA.at(1));

		auto err = [](double base, double compare)
		{
			return abs((base - compare) / base);
		};
		
		for (int part = 0; part < orig.at(0).size(); part++)
		{
			degrees newPAideal{ ionosphere::dEFlux::newPA(
						origEPA.at(1).at(part),
						bmodel->getBFieldAtS(orig.at(2).at(part), 0.0),
						bmodel->getBFieldAtS(curr.at(2).at(part), 0.0)) };

			ratio Eerr{ err(origEPA.at(0).at(part), currEPA.at(0).at(part)) };
			ratio PAerr{ err(newPAideal, currEPA.at(1).at(part)) };

			if (newPAideal < 0.0)
			{//if sim particle lower than physically possible, find s_reflect, use mirror force to calc PA of sim particle at that point
				meters sideal{ particleIdealMirrorAltitude(bmodel, orig.at(2).at(part), origEPA.at(1).at(part)) };
				tesla  Bideal{ bmodel->getBFieldAtS(sideal, 0.0) };
				degrees PAideal{ ionosphere::dEFlux::newPA(origEPA.at(1).at(part), bmodel->getBFieldAtS(orig.at(2).at(part), 0.0), Bideal) };
				degrees PAsim{ ionosphere::dEFlux::newPA(currEPA.at(1).at(part), bmodel->getBFieldAtS(curr.at(2).at(part), 0.0), Bideal) };

				PAerr = err(PAideal, PAsim);
			}
			if (PAerr > 1.0)
			{
				if (err(orig.at(2).at(part), curr.at(2).at(part)) < 0.01)
				{
					PAerr = err(180.0 - origEPA.at(1).at(part), currEPA.at(1).at(part));
				}
			}

			if (Eerr > maxEerr) maxEerr = Eerr;
			if (PAerr > maxPAerr) maxPAerr = PAerr;
		}

		if (printError) cout << "maxEerr:  " << maxEerr << "\n";
		if (printError) cout << "maxPAerr: " << maxPAerr << "\n";

		return (maxEerr < maxError && maxPAerr < maxError);
	}

	meters particleIdealMirrorAltitude(BModel* bmodel, meters sinit, degrees PAinit)
	{
		ratio   errorTolerance{ FLT_EPSILON };
		tesla   Binit{ bmodel->getBFieldAtS(sinit, 0.0) };
		meters  ds{ sinit / 50 }; //start iterating by 4%
		meters  scurr{ sinit };
		ratio   error{ 1.0 };
		degrees PAcurr{ PAinit };
		bool    up{ false }; //which direction are we moving s? up or down
		
		do //below code brackets s until we are very close to reflecting altitude (PA == 90.0)
		{			
			up ? scurr += ds : scurr -= ds; //if moving upward, add ds, else sub ds
			
			if (scurr < 0.0)
			{
				cout << "Particle s is < 0.0.  Bracketing failed.  Returning 0.\n";
				return 0.0;
			}

			PAcurr = ionosphere::dEFlux::newPA(PAinit, Binit, bmodel->getBFieldAtS(scurr, 0.0));
			
			if (PAcurr > 0.0 && up)
			{//upward moving bracketing, PA is valid - so cut down bracket step and move downward
				ds /= 10;
				up = false;
			}
			
			if (PAcurr < 0.0)
			{
				if (!up) ds /= 10; //cut step down in size if moving downward and...
				up = true;         //move upward
				continue;
			}

			error = abs((90.0 - PAcurr) / 90.0);
		} while (error > errorTolerance);

		return scurr;
	}
}