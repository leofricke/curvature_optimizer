#pragma once

#include <new>
#include <utility>
#include <cassert>
#include <algorithm>

template <typename T, std::size_t N>
class InlineVector
{
    static_assert(N >= 1, "InlineVector requires inline dimensions of at least 1.");

public:
    using iterator = T*;
    using const_iterator = const T*;

    InlineVector() noexcept
        : m_data(inline_data()), m_size(0), m_capacity(N)
    {
    }

    ~InlineVector() noexcept
    {
        clear();
        if (!using_inline())
        {
            operator delete[](m_data);
        }
    }

    InlineVector(const InlineVector& other) noexcept
        : m_data(inline_data()), m_size(0), m_capacity(N)
    {
        copy_from(other);
    }

    InlineVector& operator=(const InlineVector& other) noexcept
    {
        if (this != &other)
        {
            copy_from(other);
        }
        return *this;
    }

    InlineVector(InlineVector&& other) noexcept
        : m_data(inline_data()), m_size(0), m_capacity(N)
    {
        move_from(std::move(other));
    }

    InlineVector& operator=(InlineVector&& other) noexcept
    {
        if (this != &other)
        {
            clear();
            if (!using_inline())
            {
                operator delete[](m_data);
            }
            move_from(std::move(other));
        }
        return *this;
    }

    void push_back(const T& value)
    {
        emplace_back(value);
    }

    void push_back(T&& value)
    {
        emplace_back(std::move(value));
    }

    template <typename... Args>
    void emplace_back(Args&&... args)
    {
        if (m_size == m_capacity)
        {
            grow();
        }

        new (m_data + m_size) T(std::forward<Args>(args)...);
        ++m_size;
    }

    void pop_back()
    {
        assert(m_size > 0);
        --m_size;
        m_data[m_size].~T();
    }

    void clear()
    {
        for (std::size_t i = 0; i < m_size; ++i)
        {
            m_data[i].~T();
        }
        m_size = 0;
    }

    // clang-format off
    iterator begin() { return m_data; }
    iterator end() { return m_data + m_size; }

    const_iterator begin() const { return m_data; }
    const_iterator end() const { return m_data + m_size; }

    const_iterator cbegin() const { return m_data; }
    const_iterator cend() const { return m_data + m_size; }
    // clang-format on

    T& operator[](std::size_t i)
    {
        assert(i < m_size);
        return m_data[i];
    }

    const T& operator[](std::size_t i) const
    {
        assert(i < m_size);
        return m_data[i];
    }

    T* data()
    {
        return m_data;
    }

    const T* data() const
    {
        return m_data;
    }

    std::size_t size() const
    {
        return m_size;
    }

    std::size_t capacity() const
    {
        return m_capacity;
    }

    bool empty() const
    {
        return m_size == 0;
    }

private:
    alignas(T) unsigned char m_inline[sizeof(T) * N];
    T* m_data;
    std::size_t m_size;
    std::size_t m_capacity;

    T* inline_data()
    {
        return reinterpret_cast<T*>(m_inline);
    }

    const T* inline_data() const
    {
        return reinterpret_cast<const T*>(m_inline);
    }

    bool using_inline() const
    {
        return m_data == inline_data();
    }

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif
    void grow()
    {
        std::size_t new_capacity = m_capacity * 2;
        T* new_data = static_cast<T*>(operator new[](sizeof(T) * new_capacity));

        for (std::size_t i = 0; i < m_size; ++i)
        {
            new (new_data + i) T(std::move(m_data[i]));
        }

        for (std::size_t i = 0; i < m_size; ++i)
        {
            m_data[i].~T();
        }

        if (!using_inline())
        {
            operator delete[](m_data);
        }

        m_data = new_data;
        m_capacity = new_capacity;
    }

    void copy_from(const InlineVector& other)
    {
        m_size = other.m_size;
        m_capacity = other.m_capacity <= N ? N : other.m_capacity;

        if (!using_inline())
        {
            operator delete[](m_data);
        }
        if (other.using_inline())
        {
            m_data = inline_data();
        }
        else
        {
            m_data = static_cast<T*>(operator new[](sizeof(T) * m_capacity));
        }

        for (std::size_t i = 0; i < m_size; ++i)
        {
            new (m_data + i) T(other.m_data[i]);
        }
    }
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

    void move_from(InlineVector&& other)
    {
        m_size = other.m_size;
        m_capacity = other.m_capacity;

        if (other.using_inline())
        {
            m_data = inline_data();
            for (std::size_t i = 0; i < m_size; ++i)
            {
                new (m_data + i) T(std::move(other.m_data[i]));
            }
        }
        else
        {
            m_data = other.m_data;
            other.m_data = other.inline_data();
        }

        other.m_size = 0;
        other.m_capacity = N;
    }
};