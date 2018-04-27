#include "test.h"
#include "testHelperUtils.h"

int main()
{
	std::cout << "=================== Memory Management Tests ===================" << std::endl;
	TEST_EXCEP_CHECK(test::memLeakGPU());
	TEST_EXCEP_CHECK(test::simAttributeSaving());
	std::cout << std::endl;

	std::cout << "======================= Accuracy  Tests =======================" << std::endl;
	TEST_EXCEP_CHECK(test::dipoleAccuracy());
	TEST_EXCEP_CHECK(test::QSPSAccuracy());

	return 0;
}