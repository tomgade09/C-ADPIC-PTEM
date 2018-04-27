#ifndef TEST_HELPERUTILS_H
#define TEST_HELPERUTILS_H

#ifdef _WIN32
#include "Windows.h"
#define GREEN_INTENSE_WITH_BLACK_TEXT (BACKGROUND_GREEN | BACKGROUND_INTENSITY)
#define RED_INTENSE_WITH_WHITE_TEXT 198
#else
//#ifdef LINUX something
//#include some linux header
//#define RED_INTENSE_WITH_WHITE_TEXT whatever this is in linux
//#define GREEN_INTENSE_WITH_BLACK_TEXT whatever in linux
//#endif /* LINUX something */
#endif

#include <string>
#include <iostream>

#include "test.h"

//Test Macros
#define TEST_RESULTS(name, pass) { std::string test{ name }; ((pass) ? coutColor("PASS", GREEN_INTENSE_WITH_BLACK_TEXT) : coutColor("FAIL", RED_INTENSE_WITH_WHITE_TEXT)); std::cout << ": " << test << "\n"; }
#define TEST_EXCEP_CHECK(x) try{x;}catch(std::exception& e){std::cout << __FILE__ << ":" << __LINE__ << " : " << e.what() << std::endl;}

namespace test
{
	//CPP functions (in .cpp file)
	void coutColor(std::string str, int color);
	void coutTextOptions();
	
	//CUDA functions (in .cu file)
	void checkGPUMemory(size_t& free, size_t& total);
}

#endif /* !TEST_HELPERUTILS_H */