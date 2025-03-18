
// #include <utility>
// #include <cstddef>

// template<std::size_t N>
// struct num { static constexpr std::size_t value = N; };

// template <std::size_t N, class F>
// void for_int(F func) {
//     // Lambda with templated parameter pack (C++20 feature)
//     ([]<std::size_t... Is>(F func, std::index_sequence<Is...>){
//         (func(num<Is>{}), ...);
//     })(func, std::make_index_sequence<N>{});
// }
#include <utility>
#include <cstddef>

template<std::size_t N>
struct num { static constexpr std::size_t value = N; };

namespace detail
{
    // Helper that unpacks the index_sequence into calls to func
    template <class F, std::size_t... Is>
    void for_int_impl(F func, std::index_sequence<Is...>)
    {
        (func(num<Is>{}), ...);
    }
} // namespace detail

template <std::size_t N, class F>
void for_int(F func)
{
    detail::for_int_impl(func, std::make_index_sequence<N>{});
}
