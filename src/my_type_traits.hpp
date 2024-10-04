#pragma once

#ifndef LINAL_TYPE_TRAITS_HPP
#define LINAL_TYPE_TRAITS_HPP

#include <type_traits>
#include <utility>

template <typename T, typename = void>
struct is_addable : std::false_type {};

template <typename T>
struct is_addable<T, std::void_t<decltype(std::declval<T>() + std::declval<T>())>> : std::true_type {};

template <typename T, typename = void>
struct is_multipliable : std::false_type {};

template <typename T>
struct is_multipliable<T, std::void_t<decltype(std::declval<T>() * std::declval<T>())>> : std::true_type {};

#endif //LINAL_TYPE_TRAITS_HPP
