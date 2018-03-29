#ifndef DLLEXPORT_H
#define DLLEXPORT_H

#define DLLFILE
#ifdef DLLFILE
#define DLLEXPORT extern "C" __declspec(dllexport)
#define DLLEXP_NOEXTC __declspec(dllexport)
#else
#define DLLEXPORT
#define DLLEXP_NOEXTC
#endif /* DLLFILE */

#endif /* !DLLEXPORT_H */