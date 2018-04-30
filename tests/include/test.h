#ifndef TEST_ALLGPSTESTS_H
#define TEST_ALLGPSTESTS_H

//#define TESTS_VERBOSE //prints other useful info, often on failure, but sometimes regardless

namespace test
{
	bool memLeakGPU();
	bool simAttributeSaving(int runs = 500);
	bool dipoleAccuracy(bool save = false);
	bool QSPSAccuracy();
}

#endif /* !TEST_ALLGPSTESTS_H */