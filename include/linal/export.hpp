#pragma once

#ifndef LINAL_EXPORT_HPP
#define LINAL_EXPORT_HPP

#if defined(_WIN32) || defined(_WIN64)
    #if defined(LINAL_EXPORTS)
        #define LINAL_API __declspec(dllexport)
    #else
        #define LINAL_API __declspec(dllimport)
    #endif
#else
    #define LINAL_API
#endif

#endif // LINAL_EXPORT_HPP
