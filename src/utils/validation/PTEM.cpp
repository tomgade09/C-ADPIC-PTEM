#include "utils/validation/PTEM.h"
#include "ionosphere/ionosphere.h"
#include "Simulation/Simulation.h"
#include "utils/numerical.h"

using std::runtime_error;
using ionosphere::ParticleList;
using ionosphere::dEFlux::newPA;
using utils::numerical::EPitchTov2D;
using utils::numerical::v2DtoEPitch;

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

		Particles* particles{ sim.particles(0) };
		BModel* bmodel{ sim.Bmodel() };

		ratio maxPAerr{ 0.0 };
		ratio maxEerr{ 0.0 };

		const vector<vector<double>>& orig{ particles->data(true) };
		const vector<vector<double>>& curr{ particles->data(false) };
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

		return (maxEerr < maxError&& maxPAerr < maxError);
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

	void generateIdealDistributionAtSat(Simulation* simulation, string upwardSatName, string dnwardSatName, string bottomSatName)
	{
		cout << "Generating Ideal Distribution\n";

		Satellite* upwardsat{ nullptr };
		Satellite* dnwardsat{ nullptr };
		Satellite* bottomsat{ nullptr };
		Particles* particles{ nullptr };
		ParticleList orig;

		try
		{
			upwardsat = simulation->satellite(upwardSatName);
			dnwardsat = simulation->satellite(dnwardSatName);
			bottomsat = simulation->satellite(bottomSatName);
			particles = simulation->particles(upwardSatName);

			if (particles != simulation->particles(dnwardSatName) || particles != simulation->particles(bottomSatName))
				throw logic_error("validation::generateIdealDistributionAtSat: upward, downward, and bottom satellites"
					+ string(" do not refer to the same Particles class instance."));

			if (upwardsat->altitude() != dnwardsat->altitude())
				throw logic_error("validation::generateIdealDistributionAtSat: upward and downward satellites"
					+ string(" are not at the same altitude."));

			upwardsat->__data() = vector<vector<double>>(upwardsat->data().size(), vector<double>(upwardsat->data().at(0).size(), 0.0));
			dnwardsat->__data() = vector<vector<double>>(upwardsat->data().size(), vector<double>(upwardsat->data().at(0).size(), 0.0));
			bottomsat->__data() = vector<vector<double>>(upwardsat->data().size(), vector<double>(upwardsat->data().at(0).size(), 0.0));

			for (size_t part = 0; part < upwardsat->data().at(0).size(); part++)
			{
				upwardsat->__data().at(3).at(part) = -1.0;
				upwardsat->__data().at(4).at(part) = -1.0;
				dnwardsat->__data().at(3).at(part) = -1.0;
				dnwardsat->__data().at(4).at(part) = -1.0;
				bottomsat->__data().at(3).at(part) = -1.0;
				bottomsat->__data().at(4).at(part) = -1.0;
			}

			vector<vector<double>> data{ particles->data(true) };
			orig = ParticleList(data.at(0), data.at(1), particles->mass());
		}
		catch (const std::exception& e)
		{
			cout << e.what() << "\nReturning without modifying satellite distributions.\n";
			return;
		}

		bool qsps{ false };

		//guards that are necessary for now, perhaps in the future, a more robust generator will be written
		{
			for (int num = 0; num < simulation->Efield()->size(); num++)
			{
				if (simulation->Efield()->emodel(num)->type_m == EModel::Type::QSPS)
				{
					cout << "EField contains QSPS.\n";
					qsps = true;
				}
				else
				{
					cout << "Contains Alfven waves or type 'Other'.  No practical way to do this without time stepping."
						<< "  Returning without modifying satellite distributions.\n";
					return;
				}
			}

			if (simulation->Efield()->size() > 1)
			{
				cout << "EField contains more than 1 QSPS.  Right now this isn't supported."
					<< "  Returning without modifying satellite distributions.\n";
				return;
			}

			if (qsps)
			{ //guards that are needed for now to enforce simple QSPS's
				for (int num = 0; num < simulation->Efield()->size(); num++)
				{
					QSPS* q{ dynamic_cast<QSPS*>(simulation->Efield()->emodel(num)) };

					if (q->magnitude().size() > 1)
					{
						cout << "QSPS contains more than 1 region.  At this point, this isn't supported by this function."
							<< "  Returning without modifying satellite distributions.\n";
						return;
					}

					for (int region = 0; region < q->altMin().size(); region++)
					{
						if (q->altMin().at(region) < upwardsat->altitude())
						{
							cout << "Satellite is within or above QSPS which is not supported for ideal dist at this time."
								<< "Returning without modifying satellite distribution.\n";
							return;
						}
					}
				}
			}
		}

		//Now, finally, generate the distribution
		auto sumDeltaKMirror = [&](meters qspstop, meters qspsbtm, double mu, meters delta_s)
		{
			double sum{ 0.0 };

			for (double s = qspstop; s > qspsbtm; s -= delta_s)
				sum += -mu * simulation->Bmodel()->getGradBAtS(s, 0.0) * delta_s;

			return sum;
		};

		auto EPATov_Mirror = [&](eV E, degrees PA, meters s_init, Satellite* sat, size_t part, bool invertPA = false)
		{
			vector<vector<double>>& write{ sat->__data() };

			degrees PA_new{ newPA(PA, simulation->getBFieldAtS(s_init, 0.0), simulation->getBFieldAtS(sat->altitude(), 0.0)) };
			if (PA_new < 0.0)
			{
				write.at(0).at(part) = 0.0;
				write.at(1).at(part) = 0.0;
				write.at(2).at(part) = 0.0;
				write.at(3).at(part) = -1.0;
				write.at(4).at(part) = -1.0;

				return;
			}
			vector<double> vpara(1, 0.0);
			vector<double> vperp(1, 0.0);
			EPitchTov2D({ E }, { (invertPA ? 180.0 - PA_new : PA_new) }, MASS_ELECTRON, vpara, vperp);

			write.at(0).at(part) = vpara.at(0);
			write.at(1).at(part) = vperp.at(0);
			write.at(2).at(part) = sat->altitude();
			write.at(3).at(part) = 0.0;
			write.at(4).at(part) = (double)part;
		};

		auto modifyEPAwithAddedEpara = [&](eV* E, degrees* PA, eV Epara_add, Particles* part)
		{
			vector<double> vpara(1, 0.0);
			vector<double> vperp(1, 0.0);
			vector<double> vpara_add(1, 0.0);
			vector<double> Etot(1, 0.0);
			vector<double> PAtot(1, 0.0);
			vector<double> discard(1, 0.0);

			EPitchTov2D({ *E }, { *PA }, part->mass(), vpara, vperp);
			EPitchTov2D({ Epara_add }, { 0.0 }, part->mass(), vpara_add, discard);

			vpara.at(0) += vpara_add.at(0);

			v2DtoEPitch(vpara, vperp, part->mass(), Etot, PAtot);
			*E = Etot.at(0);
			*PA = PAtot.at(0);
		};

		//upgoing
		for (size_t part = 0; part < particles->data(true).size(); part++)
		{
			if (particles->data(true).at(2).at(part) < simulation->simMin() * 1.001)
			{ //ionospheric source
				EPATov_Mirror(orig.energy.at(part), orig.pitch.at(part), orig.s_pos.at(part), dnwardsat, part); //from ionosphere to sat downward detector
				if (qsps)
				{
					QSPS* q{ dynamic_cast<QSPS*>(simulation->Efield()->emodel(0)) };
					degrees PA_QSPS_bottom{ newPA(orig.pitch.at(part), simulation->getBFieldAtS(orig.s_pos.at(part), 0.0), q->altMin().at(0)) };

					double K_QSPS_Total{ particles->charge() * q->magnitude().at(0) * (q->altMax().at(0) - q->altMin().at(0)) }; //q * E * delta s
					double vpara_QSPS_bottom{ -sqrt(2 * orig.energy.at(part) * JOULE_PER_EV / particles->mass()) * cos(PA_QSPS_bottom * RADS_PER_DEG) };

					if (abs(K_QSPS_Total) > abs(0.5 * particles->mass() * vpara_QSPS_bottom * vpara_QSPS_bottom)) //if reflect
					{
						EPATov_Mirror(orig.energy.at(part), orig.pitch.at(part), orig.s_pos.at(part), upwardsat, true); //from ionosphere to QSPS, reflect, back down to satellite
						EPATov_Mirror(orig.energy.at(part), orig.pitch.at(part), orig.s_pos.at(part), bottomsat, part, true); //from iono to QSPS, reflect, back to bottom of sim
					}
				}
			}
			else if (particles->data(true).at(2).at(part) > simulation->simMax() * 0.999)
			{ //magnetospheric source
				if (qsps)
				{
					QSPS* q{ dynamic_cast<QSPS*>(simulation->Efield()->emodel(0)) };
					degrees PA_QSPS_bottom{ newPA(orig.pitch.at(part), simulation->getBFieldAtS(orig.s_pos.at(part), 0.0), q->altMin().at(0)) };
					degrees PA_QSPS_top{ newPA(orig.pitch.at(part), simulation->getBFieldAtS(orig.s_pos.at(part), 0.0), q->altMax().at(0)) };
					if (PA_QSPS_top < 0.0)
					{
						orig.pitch.at(part) = 180.0 - orig.pitch.at(part);
						continue;
					}

					double K_QSPS_Total{ particles->charge() * q->magnitude().at(0) * (q->altMax().at(0) - q->altMin().at(0)) }; //q * E * delta s
					double vpara_QSPS_top{ -sqrt(2 * orig.energy.at(part) * JOULE_PER_EV / particles->mass()) * cos(PA_QSPS_top * RADS_PER_DEG) };
					double Kpara_QSPS_top{ 0.5 * particles->mass() * vpara_QSPS_top * vpara_QSPS_top };
					if (PA_QSPS_bottom < 0.0)
					{
						double mu{ 0.5 * particles->mass() * orig.vperp.at(part) * orig.vperp.at(part) / simulation->Bmodel()->getBFieldAtS(orig.s_pos.at(part), 0.0) };
						double K_mirror_over_QSPS_region{ abs(sumDeltaKMirror(q->altMax().at(0), q->altMin().at(0), mu, 0.1)) };
						if (abs(Kpara_QSPS_top) + abs(K_QSPS_Total) < K_mirror_over_QSPS_region) //reflects
						{
							orig.pitch.at(part) = 180.0 - orig.pitch.at(part);
							continue;
						}
					}

					//if the code makes it here, the particle didn't reflect by the bottom of the QSPS
					orig.pitch.at(part) = newPA(orig.pitch.at(part), simulation->getBFieldAtS(orig.s_pos.at(part), 0.0), simulation->getBFieldAtS(q->altMax().at(0), 0.0));
					modifyEPAwithAddedEpara(&(orig.energy.at(part)), &(orig.pitch.at(part)), K_QSPS_Total, particles);
					orig.pitch.at(part) = newPA(orig.pitch.at(part), simulation->getBFieldAtS(q->altMax().at(0), 0.0), simulation->getBFieldAtS(q->altMin().at(0), 0.0));
					orig.s_pos.at(part) = q->altMin().at(0);
				}

				EPATov_Mirror(orig.energy.at(part), orig.pitch.at(part), orig.s_pos.at(part), upwardsat, part); //upward facing detector
				EPATov_Mirror(orig.energy.at(part), orig.pitch.at(part), orig.s_pos.at(part), bottomsat, part);
				if (newPA(orig.pitch.at(part), simulation->getBFieldAtS(orig.s_pos.at(part), 0.0), simulation->getBFieldAtS(bottomsat->altitude(), 0.0)) < 0.0)
					EPATov_Mirror(orig.energy.at(part), orig.pitch.at(part), orig.s_pos.at(part), dnwardsat, part, true);//need to remove those that make it to bottom
			}
			else
			{
				throw runtime_error("validation::generateIdealDistributionAtSat: particle is not able to be assigned to either ionospheric or magnetospheric source.");
			}
		}
	}
}
