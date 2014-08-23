#if !defined(__EXPORT_TYPE)
 #if defined (__BUILDING_THE_DLL)
 #define __EXPORT_TYPE __export
 #elif defined (__BUILDING_GNU) && !defined (__BUILDING_LINUX)
 #define __EXPORT_TYPE __attribute__((dllexport))
 #elif defined (__BUILDING_GNU) && defined (__BUILDING_LINUX)
 #define __EXPORT_TYPE
 #else
 #define __EXPORT_TYPE __import
 #endif
#endif

#if defined(__BUILDING_GNU)
#define VEC_ITER_TYPE typename
#else
#define VEC_ITER_TYPE 
#endif
