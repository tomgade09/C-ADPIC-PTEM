#ifndef DLLIMPEXP_DEFINES_H
#define DLLIMPEXP_DEFINES_H

#ifdef DLLFILE
#define DLLEXP_EXTC extern "C" __declspec(dllexport)
#define DLLEXP __declspec(dllexport)
#define DLLIMP __declspec(dllimport)
#else
#define DLLEXP_EXTC
#define DLLEXP
#define DLLIMP
#endif /* DLLFILE */

#endif /* !DLLIMPEXP_DEFINES_H */