#include <gtest/gtest.h>

#include "utils/InlineVector.h"
#include <utility>

template <typename T, std::size_t N>
bool usesInlineStorage(const InlineVector<T, N>& v)
{
    const auto* dataPtr = reinterpret_cast<const unsigned char*>(v.data());
    const auto* objStart = reinterpret_cast<const unsigned char*>(&v);
    const auto* objEnd = objStart + sizeof(InlineVector<T, N>);
    return dataPtr >= objStart && dataPtr < objEnd;
}

// Instrumentation type used for testing
// clang-format off
struct Tracker
{
    inline static int constructions = 0;
    inline static int destructions = 0;

    int value;

    Tracker(int v = 0) : value(v) { ++constructions; }
    Tracker(const Tracker& o) : value(o.value) { ++constructions; }
    Tracker(Tracker&& o) noexcept : value(o.value) { ++constructions; }
    ~Tracker() { ++destructions; }

    static void reset()
    {
        constructions = 0;
        destructions = 0;
    }
    static int alive() { return constructions - destructions; }
};
// clang-format on

TEST(InlineVectorTest, BasicTest)
{
    InlineVector<int, 4> v;
    EXPECT_TRUE(v.empty());
    EXPECT_EQ(v.size(), 0);
    EXPECT_EQ(v.capacity(), 4);

    v.push_back(10);
    v.push_back(20);
    v.emplace_back(30);

    ASSERT_EQ(v.size(), 3);
    EXPECT_FALSE(v.empty());
    EXPECT_EQ(v[0], 10);
    EXPECT_EQ(v[1], 20);
    EXPECT_EQ(v[2], 30);

    v[0] = 100;
    EXPECT_EQ(v[0], 100);

    int sum = 0;
    for (int x : v)
    {
        sum += x;
    }
    EXPECT_EQ(sum, 150);

    v.pop_back();
    EXPECT_EQ(v.size(), 2);

    v.clear();
    EXPECT_TRUE(v.empty());
    EXPECT_EQ(v.capacity(), 4);
}

TEST(InlineVectorTest, GrowsFromInlineToHeapPreservingValues)
{
    InlineVector<int, 4> v;
    for (int i = 0; i < 4; ++i)
    {
        v.push_back(i);
    }
    EXPECT_EQ(v.capacity(), 4);
    EXPECT_TRUE(usesInlineStorage(v));

    v.push_back(4); // triggers growth out of inline storage
    EXPECT_EQ(v.capacity(), 8);
    EXPECT_FALSE(usesInlineStorage(v));

    for (int i = 5; i < 20; ++i)
    {
        v.push_back(i);
    }

    ASSERT_EQ(v.size(), 20);
    for (int i = 0; i < 20; ++i)
    {
        EXPECT_EQ(v[i], i);
    }
}

TEST(InlineVectorTest, CopyIsIndependentDeepCopy)
{
    InlineVector<int, 4> a;
    for (int i = 0; i < 10; ++i)
    {
        a.push_back(i);
    }

    InlineVector<int, 4> b(a);
    InlineVector<int, 4> c;
    c = a;

    ASSERT_EQ(b.size(), 10);
    ASSERT_EQ(c.size(), 10);
    for (int i = 0; i < 10; ++i)
    {
        EXPECT_EQ(b[i], i);
        EXPECT_EQ(c[i], i);
    }

    // Mutating the copies must not affect the original.
    b[0] = 99;
    c[1] = 88;
    EXPECT_EQ(a[0], 0);
    EXPECT_EQ(a[1], 1);
}

TEST(InlineVectorTest, MoveTransfersContents)
{
    InlineVector<int, 4> a;
    for (int i = 0; i < 10; ++i)
    {
        a.push_back(i);
    }
    const int* originalData = a.data();

    InlineVector<int, 4> b(std::move(a));

    ASSERT_EQ(b.size(), 10);
    EXPECT_EQ(b.data(), originalData);
    for (int i = 0; i < 10; ++i)
    {
        EXPECT_EQ(b[i], i);
    }

    EXPECT_EQ(a.size(), 0);
    EXPECT_EQ(a.capacity(), 4);
}

TEST(InlineVectorTest, ManagesElementLifetimes)
{
    Tracker::reset();
    {
        InlineVector<Tracker, 4> v;
        for (int i = 0; i < 20; ++i)
        {
            v.emplace_back(i);
        }
        EXPECT_EQ(Tracker::alive(), 20);

        v.pop_back();
        EXPECT_EQ(Tracker::alive(), 19);

        v.clear();
        EXPECT_EQ(Tracker::alive(), 0);

        v.emplace_back(42);
    }
    EXPECT_EQ(Tracker::alive(), 0);
}