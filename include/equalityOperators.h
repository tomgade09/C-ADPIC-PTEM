#ifndef EQUALITYOPERATORS_H
#define EQUALITYOPERATORS_H

#include "Simulation/Simulation.h"

inline bool operator==(const Simulation& x, const Simulation& y)
{
	try
	{
		return (x.dt() == y.dt()) &&
			(x.simMin() == y.simMin()) &&
			(x.simMax() == y.simMax()) &&
			(x.getNumberOfParticleTypes() == y.getNumberOfParticleTypes() &&
			(x.getNumberOfSatellites() == y.getNumberOfSatellites()));
	}
	catch (...)
	{
		return false;
	}
}

inline bool operator==(const Particle& x, const Particle& y)
{
	try
	{
		bool same{ (x.name() == y.name()) &&
			(x.mass() == y.mass()) &&
			(x.charge() == y.charge()) &&
			(x.getNumberOfParticles() == y.getNumberOfParticles()) };

		for (int iii = 0; iii < x.getNumberOfAttributes(); iii++)
			same &= (x.getAttrNameByInd(iii) == y.getAttrNameByInd(iii));

		same &= (x.data(true ) == y.data(true ));
		same &= (x.data(false) == x.data(false));

		return same;
	}
	catch (...)
	{
		return false;
	}
}

inline bool operator==(const Satellite& x, const Satellite& y)
{
	try
	{
		return (x.name() == y.name()) &&
			   (x.altitude() == y.altitude()) &&
		       (x.upward() == y.upward()) &&
		       (x.data() == y.data());
	}
	catch (...)
	{
		return false;
	}
}

inline bool operator==(const BField& x, const BField& y)
{
	try
	{
		if (!(
			(x.name() == y.name()) &&
			(x.getBFieldAtS(4.0e6, 0.0) == y.getBFieldAtS(4.0e6, 0.0)) &&
			(x.getGradBAtS(4.0e6, 0.0) == y.getGradBAtS(4.0e6, 0.0))
			))
			throw std::invalid_argument("");

		if (x.name() == "DipoleB")
		{
			if (((DipoleB*)&x)->getErrTol() != ((DipoleB*)&y)->getErrTol())
				throw std::invalid_argument("");
			if (((DipoleB*)&x)->getds() != ((DipoleB*)&y)->getds())
				throw std::invalid_argument("");
		}
		else if (x.name() == "DipoleBLUT")
		{
			if (((DipoleBLUT*)&x)->getErrTol() != ((DipoleBLUT*)&y)->getErrTol())
				throw std::invalid_argument("");
			if (((DipoleBLUT*)&x)->getds() != ((DipoleBLUT*)&y)->getds())
				throw std::invalid_argument("");
		}
		else
			std::invalid_argument("");

		return true;
	}
	catch (...)
	{
		std::cout << "caught something: BField\n";
		return false;
	}
}

inline bool operator==(const EElem& x, const EElem& y)
{
	try
	{
		if (x.name() != y.name())
			throw std::invalid_argument(""); //just throwing to return false

		if (x.name() == "QSPS")
		{
			if (((QSPS*)&x)->altMin() != ((QSPS*)&y)->altMin())
				throw std::invalid_argument("");

			if (((QSPS*)&x)->altMax() != ((QSPS*)&y)->altMax())
				throw std::invalid_argument("");

			if (((QSPS*)&x)->magnitude() != ((QSPS*)&y)->magnitude())
				throw std::invalid_argument("");
		}
		else
			throw std::invalid_argument(""); //later need to implement other classes (once added)

		return true;
	}
	catch (...)
	{
		return false;
	}
}

inline bool operator==(const EField& x, const EField& y)
{
	try
	{
		if (!(
			(x.capacity() == y.capacity()) &&
			(x.size() == y.size()) &&
			(x.getEFieldAtS(4.0e6, 0.0) == y.getEFieldAtS(4.0e6, 0.0))
			))
			throw std::invalid_argument("");

		for (int eelem = 0; eelem < x.size(); eelem++) //maybe shouldn't check elems here??
			if (!(*(x.element(eelem)) == *(y.element(eelem)))) //if the elements (references) are not equal - operator defined above
				throw std::invalid_argument("");

		return true;
	}
	catch (...)
	{
		return false;
	}
}

#endif