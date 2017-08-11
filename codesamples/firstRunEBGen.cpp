#include <Eigen/Dense>
#include <cmath>
//B = (B0/r^3) * sqrt(1 + 3 * cos(theta)^2)
//E -> static

//				   x (p[0])
//				  __________
//				  |		   /
//				  |		  / <-- r
//		z (p[1])  |------/
//				  |theta/
//				_____  /		theta -> atan(x/z)
//			 /	  |	  / \
//		  /		  |	 /	   \
//	   /		  |	/	  	  \
//	 |	Earth 	  |/  		    |

double BFieldAtZ(const Eigen::Vector2d& p)
{//convert to polar, calc B, return B
	double r{ std::sqrt(p[0] * p[0] + p[1] * p[1]) };
	double theta{ std::atan2(p[0], p[1]) }; //atan(x / z)
	
	double B_r{ (B0 / pow(r, 3)) * std::sqrt(1 + 3 * pow(cos(theta), 2)) };
	//B is radially outward, so B_theta is theta
	double B_x{ B_r * sin(theta) };
	double B_z{ B_r * cos(theta) };

	return << B_x, B_z;
}

Eigen::Vector2d EFieldAtP(const Eigen::Vector2d& p)
{
	return << 2., 2.; //some realistic E value
}