#include "testHelperUtils.h"

#include <iostream>

namespace test
{
	//Pretty Console Colors
	#ifdef _WIN32
	void coutColor(std::string str, int color) //Windows way of changing console color
	{
		HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get handle to console
		
		CONSOLE_SCREEN_BUFFER_INFO origInfo;
		GetConsoleScreenBufferInfo(hConsole, &origInfo); //save original attributes

		SetConsoleTextAttribute(hConsole, color); //set color
		std::cout << str; //print string

		SetConsoleTextAttribute(hConsole, origInfo.wAttributes); //set back to normal
	}

	void coutTextOptions()
	{
		for (int iii = 0; iii < 512; iii++)
			coutColor(std::to_string(iii) + " : Test text\n", iii);
	}
	#else
	void coutColor(std::string str, int color)
	{
		std::cout << str; //for right now, color not defined on other OSes...build *nix version eventually
	}
	void printTextOptions()
	{
		coutColor("All you got for now\n", 0);
	}
	#endif
}