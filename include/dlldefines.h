#ifndef DLLIMPEXP_DEFINES_H
#define DLLIMPEXP_DEFINES_H

//Windows defines
#ifdef DLLFILE //DLLFILE is defined for GPS, not TESTS
#define DLLEXP_EXTC extern "C" __declspec(dllexport) //for python use of the library
#define DLLEXP __declspec(dllexport)
#define DLLCLEXP __declspec(dllexport)
#else
#define DLLEXP_EXTC
#define DLLEXP
#define DLLCLEXP __declspec(dllimport)
#endif /* DLLFILE */

#endif /* !DLLIMPEXP_DEFINES_H */