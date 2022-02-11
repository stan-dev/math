
#ifndef SUNDIALS_EXPORT_H
#define SUNDIALS_EXPORT_H

#ifdef SUNDIALS_STATIC_DEFINE
#  define SUNDIALS_EXPORT
#  define SUNDIALS_NO_EXPORT
#else
#  ifndef SUNDIALS_EXPORT
#    ifdef sundials_generic_EXPORTS
        /* We are building this library */
#      define SUNDIALS_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define SUNDIALS_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef SUNDIALS_NO_EXPORT
#    define SUNDIALS_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef SUNDIALS_DEPRECATED
#  define SUNDIALS_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef SUNDIALS_DEPRECATED_EXPORT
#  define SUNDIALS_DEPRECATED_EXPORT SUNDIALS_EXPORT SUNDIALS_DEPRECATED
#endif

#ifndef SUNDIALS_DEPRECATED_NO_EXPORT
#  define SUNDIALS_DEPRECATED_NO_EXPORT SUNDIALS_NO_EXPORT SUNDIALS_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef SUNDIALS_NO_DEPRECATED
#    define SUNDIALS_NO_DEPRECATED
#  endif
#endif

#endif /* SUNDIALS_EXPORT_H */
