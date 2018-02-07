#ifndef SIMEXCEPHANDLER_H
#define SIMEXCEPHANDLER_H

#include <stdexcept> //necessary?
#include "ErrorHandling\SimFatalException.h"

//I think I want to throw a CUDA specific exception - if anything CUDA fails, might as well just exit(1/some other int)
#define SIM_API_EXCEP_CHECK(x) \
	try{ x; } \
	catch(const SimFatalException& fatal) { std::cerr << fatal.where() << " : " << "Simulation exception: " << fatal.what() << std::endl \
										<< "     >> " << fatal.args() << std::endl; exit(EXIT_FAILURE); } \
	catch(const SimException& simExp)     { std::cerr << simExp.where() << " : " << "Simulation exception: " << simExp.what() << std::endl \
										<< "     >> " << simExp.args() << std::endl; } \
	catch(const std::invalid_argument& a){std::cerr << __FILE__ << ":" << __LINE__ << " : " << "Invalid argument error: " << a.what() << std::endl; } \
	catch(const std::out_of_range& oor) { std::cerr << __FILE__ << ":" << __LINE__ << " : " << "Vector out of range exception: " << oor.what() << std::endl; } \
	catch(const std::logic_error& log)  { std::cerr << __FILE__ << ":" << __LINE__ << " : " << "Logic error exception: " << log.what() << std::endl; } \
	catch(const std::runtime_error& rte){ std::cerr << __FILE__ << ":" << __LINE__ << " : " << "Generic runtime error: " << rte.what() << std::endl; } \
	catch(const std::exception& exp)    { std::cerr << __FILE__ << ":" << __LINE__ << " : " << "Unhandled specific std::exception: " << exp.what() << std::endl; } \
	catch(...)                          { std::cerr << __FILE__ << ":" << __LINE__ << " : " << "Other unhandled exception - exiting out of precaution"  << std::endl; exit(EXIT_FAILURE);}

#endif /* SIMEXCEPHANDLER_H */