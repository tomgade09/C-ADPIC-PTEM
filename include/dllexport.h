#ifndef DLLEXPORT_H
#define DLLEXPORT_H

#define DLLFILE
#ifdef DLLFILE
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT
#endif /* DLLFILE */

#endif /* !DLLEXPORT_H */