#ifndef DLLIMPEXP_DEFINES_H
#define DLLIMPEXP_DEFINES_H

//Windows defines
#ifdef _WIN32
#ifdef DLLFILE //DLLFILE is defined for GPS, not TESTS
#define DLLEXP_EXTC extern "C" __declspec(dllexport) //for python use of the library
#define DLLEXP __declspec(dllexport)
#define DLLCLEXP __declspec(dllexport)
#else
#define DLLEXP_EXTC
#define DLLEXP
#define DLLCLEXP __declspec(dllimport)
#endif /* DLLFILE */
#else /* !_WIN32 */
#define DLLEXP_EXTC
#define DLLEXP
#define DLLCLEXP
#define FLT_EPSILON 1.192092896e-7F
#endif /* _WIN32 */

using std::log10;
using std::pow;
using std::exp;
using std::sqrt;
using std::sin;
using std::cos;

#endif /* !DLLIMPEXP_DEFINES_H */
