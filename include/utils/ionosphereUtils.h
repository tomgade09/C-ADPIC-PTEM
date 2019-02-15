#ifndef IONOSPHERE_UTILS_H
#define IONOSPHERE_UTILS_H

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <functional>

using std::cout;
using std::string;
using std::vector;
using std::function;
using std::to_string;

//types to make function arguments more explicit
typedef double eV;
typedef double kg;
typedef double tesla;
typedef double meters;
typedef double dNflux;
typedef double dEflux;
typedef double degrees;
typedef double coulomb;
typedef vector<double> degrees_v1D;
typedef vector<double> dNflux_v1D;
typedef vector<double> dEflux_v1D;
typedef vector<vector<double>> dNflux_v2D;
typedef vector<vector<double>> dEflux_v2D;

typedef vector<double> double_v1D;
typedef vector<vector<double>> double_v2D;
//types to make function arguments more explicit

#define TESTVEC_ISZEROFRSTHALF(vec, name) vecTest(vec, [](double cnt) { return (cnt != 0.0); }, true, name, 0, (int)vec.size() / 2);
#define TESTVEC_ISZEROLASTHALF(vec, name) vecTest(vec, [](double cnt) { return (cnt != 0.0); }, true, name, (int)vec.size() / 2);
#define TESTVEC_ISZEROWHOLEVEC(vec, name) vecTest(vec, [](double cnt) { return (cnt != 0.0); }, true, name);
#define TESTVEC_NOTNEGWHOLEVEC(vec, name) vecTest(vec, [](double cnt) { return (cnt < 0.0); }, true, name);

inline bool vecTest(const double_v2D& vec, function<bool(double)> test, bool throwOnTrue = false, string label = "",
	int outStart = 0, int outStop = 0, int inStart = 0, int inStop = 0)
{
	if (outStop == 0) outStop = (int)vec.size();

	for (int out = outStart; out < outStop; out++)
	{
		if (inStop == 0) inStop = (int)vec.at(out).size();
		for (int in = 0; in < inStop; in++)
			if (test(vec.at(out).at(in)))
			{
				if (throwOnTrue)
					throw std::logic_error("vecTest: " + label + " condition met.  Throwing - out, in, data: " +
						to_string(out) + ", " + to_string(in) + ", " + to_string(vec.at(out).at(in)));

				return true;
			}
	}

	return false;
}

inline bool vecTest(const double_v1D& vec, function<bool(double)> test, bool throwOnTrue = false, string label = "",
	int start = 0, int stop = 0)
{
	double_v2D tmp;
	tmp.push_back(vec);

	return vecTest(tmp, test, throwOnTrue, label, 0, 0, start, stop);
}

inline double_v1D serialize2DVec(const double_v2D& in)
{
	double_v1D out;

	for (auto& v : in)
		out.insert(std::end(out), std::begin(v), std::end(v));

	return out;
}

inline void printVec2D(const double_v2D& prt, string name)
{
	cout << name << ":\n";
	cout << std::setprecision(10);
	for (auto& dblVec : prt)
	{
		for (auto& elem : dblVec)
			cout << elem << ",";
		cout << "\n";
	}
	cout << "\n";
};

#endif /* !IONOSPHERE_UTILS_H */