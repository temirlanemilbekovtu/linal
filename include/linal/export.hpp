#pragma once

#ifndef LINAL_EXPORT_HPP
#define LINAL_EXPORT_HPP

#if defined(_WIN32) || defined(__CYGWIN__)
    #ifdef LINAL_EXPORTS
        #define LINAL_API __declspec(dllexport)
    #else
        #define LINAL_API __declspec(dllimport)
    #endif
#elif defined(__GNUC__) && __GNUC__ >= 4
    #ifdef LINAL_EXPORTS
        #define LINAL_API __attribute__((visibility("default")))
    #else
        #define LINAL_API
    #endif
#else
    #define LINAL_API
#endif

#endif // LINAL_EXPORT_HPP