#ifndef CLASSMACROS_H
#define CLASSMACROS_H

/**
 * Convenience macro to declare that a class is not copyable.
 */

// (jenya) 2011/12/11: Undefing another definition of UNCOPYABLE_CLASS from toolbox.
// Proper fix is to remove Uncopable.h from toolbox.
#undef UNCOPYABLE_CLASS

#if __cplusplus >= 201103L
    #define UNCOPYABLE_CLASS(CLASS_NAME) \
        private:\
            CLASS_NAME(CLASS_NAME const&) = delete;\
            CLASS_NAME& operator=(CLASS_NAME const&) = delete;
#else
    #define UNCOPYABLE_CLASS(CLASS_NAME) \
        private:\
            CLASS_NAME(CLASS_NAME const&);\
            CLASS_NAME& operator=(CLASS_NAME const&)
#endif

#endif
