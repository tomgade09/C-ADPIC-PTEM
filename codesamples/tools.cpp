#include <Eigen/Dense>

inline Eigen::Vector3d lorentzForce(double q, const Eigen::Vector3d& E, const Eigen::Vector3d& v, const Eigen::Vector3d& B)
{ return q * (E + v.cross(B)); } //needed?

inline Eigen::Vector3d mirrorForce(double mu, Eigen::Vector3d dB_ds_var)
{ return -mu * dB_ds_var; } //needed?

inline Eigen::Vector3d dB_ds(const Eigen::Vector3d& p, double scaleLength)
{
	Eigen::Vector3d Bp = calcBatP(p); //what if Bp == 0.0 on some axis???
	Eigen::Vector3d ds = scaleLength * Bp / Bp.norm();
	
	return (calcBatP(p + ds) - calcBatP(p - ds)) / (2 * ds.norm());
}

inline double particleMu(double mass, double vPerpSq, double Bnorm)
{ return mass * vPerpSq / (2 * Bnorm); } //needed?

inline Eigen::Vector3d vParallel(Eigen::Vector3d& v, Eigen::Vector3d& B)
{ return v.dot(B) * B / B.norm(); } // needed?

inline Eigen::Vector3d vPerpendicular(Eigen::Vector3d& v, Eigen::Vector3d& vPara)
{ return (v - vPara); }

inline Eigen::Vector3d calcBatP(Eigen::Vector3d& p, double Bdata) //Probably need to change Bdata type
{
	//some method here - maybe read a file, or calculate in real time???  Either way, this will change based on the application
	return 0.0;
}

inline Eigen::Vector3d calcEatP(Eigen::Vector3d& p, double Edata) //Probably need to change Edata type
{
	//some method here - maybe read a file, or calculate in real time???  Either way, this will change based on the application
	return 0.0;
}

inline Eigen::Vector3d vPerpFromMu(Eigen::Vector3d& B, double mu)
{
	Eigen::Vector3d vPerpUV = ;
	
	return mu * 2 * B.norm() * vPerpUV;
} //not right - not a vector - how to find orthogonal vector???
												//F = v x B, v = B x F