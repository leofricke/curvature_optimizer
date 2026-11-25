#pragma once

#include <memory>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>

namespace spmat
{

// Allocator that default-initializes instead of value-initializes, so resize() on trivial types
// leaves the new elements uninitialized instead of zero-filling them.
template <typename T, typename Alloc = std::allocator<T>>
class DefaultInitAllocator : public Alloc
{
    using Traits = std::allocator_traits<Alloc>;

public:
    template <typename U>
    struct rebind
    {
        using other = DefaultInitAllocator<U, typename Traits::template rebind_alloc<U>>;
    };

    using Alloc::Alloc;

    template <typename U>
    void construct(U* ptr) noexcept(std::is_nothrow_default_constructible_v<U>)
    {
        ::new (static_cast<void*>(ptr)) U;
    }

    template <typename U, typename... Args>
    void construct(U* ptr, Args&&... args)
    {
        Traits::construct(static_cast<Alloc&>(*this), ptr, std::forward<Args>(args)...);
    }
};

template <typename T>
using DefaultInitVector = std::vector<T, DefaultInitAllocator<T>>;

} // namespace spmat
