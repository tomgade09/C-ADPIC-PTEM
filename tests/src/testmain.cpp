#include "test.h"
#include "testHelperUtils.h"

int main()
{
	bool pass{ true };
	std::cout << "=================== Memory Management Tests ===================" << std::endl;
	TEST_EXCEP_CHECK(pass &= test::memLeakGPU());
	TEST_EXCEP_CHECK(pass &= test::simAttributeSaving());
	std::cout << std::endl;

	std::cout << "======================= Accuracy  Tests =======================" << std::endl;
	TEST_EXCEP_CHECK(pass &= test::dipoleAccuracy());
	TEST_EXCEP_CHECK(pass &= test::QSPSAccuracy());

	std::cout << std::endl;
	std::cout << "Summary of Tests:" << std::endl;
	if (pass)
	{
		test::coutColor("                \n", GREEN_INTENSE_WITH_BLACK_TEXT);
		test::coutColor("  >>> PASS <<<  \n", GREEN_INTENSE_WITH_BLACK_TEXT);
		test::coutColor("                \n", GREEN_INTENSE_WITH_BLACK_TEXT);
	}
	else
	{
		test::coutColor("                \n", RED_INTENSE_WITH_WHITE_TEXT);
		test::coutColor("  >>> FAIL <<<  \n", RED_INTENSE_WITH_WHITE_TEXT);
		test::coutColor("                \n", RED_INTENSE_WITH_WHITE_TEXT);
	}

	return 0;
}