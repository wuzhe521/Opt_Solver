#include "../core/vector.h"
#include <gtest/gtest.h>

namespace gons {

/**
 * @brief Test case for the Vector constructor with initializer list.
 * This test verifies the correctness of the constructor when:
 * 1. The initializer list size matches the vector size.
 * 2. The initializer list size is smaller than the vector size.
 * 3. The initializer list size is larger than the vector size.
 */
TEST(VectorTest, InitializerListConstructor) {
    // Test case 1: Initializer list size matches vector size
    Vector<double, 3> vec1 = {1.0, 2.0, 3.0};
    EXPECT_DOUBLE_EQ(1.0, vec1(0, 0));
    EXPECT_DOUBLE_EQ(2.0, vec1(1, 0));
    EXPECT_DOUBLE_EQ(3.0, vec1(2, 0));

    // Test case 2: Initializer list size is smaller than vector size
    Vector<double, 4> vec2 = {1.0, 2.0};
    EXPECT_DOUBLE_EQ(1.0, vec2(0, 0));
    EXPECT_DOUBLE_EQ(2.0, vec2(1, 0));
    EXPECT_DOUBLE_EQ(0.0, vec2(2, 0)); // Default-initialized
    EXPECT_DOUBLE_EQ(0.0, vec2(3, 0)); // Default-initialized

    // Test case 3: Initializer list size is larger than vector size
    Vector<double, 2> vec3 = {1.0, 2.0, 3.0, 4.0};
    EXPECT_DOUBLE_EQ(1.0, vec3(0, 0));
    EXPECT_DOUBLE_EQ(2.0, vec3(1, 0));
    // Elements beyond R are ignored
}

/**
 * @brief Test case for the Vector constructor with default initialization.
 * This test verifies that the constructor initializes all elements to the default value of T.
 */
TEST(VectorTest, DefaultConstructor) {
    Vector<double, 3> vec;
    EXPECT_DOUBLE_EQ(0.0, vec(0, 0));
    EXPECT_DOUBLE_EQ(0.0, vec(1, 0));
    EXPECT_DOUBLE_EQ(0.0, vec(2, 0));
}

/**
 * @brief Test case for the Vector constructor with Matrix input.
 * This test verifies the correctness of the constructor when:
 * 1. The input Matrix is an lvalue.
 * 2. The input Matrix is an rvalue.
 */
TEST(VectorTest, MatrixConstructor) {
    // Test case 1: lvalue Matrix
    Matrix<double, 3, 1> mat1 = {{1.0}, {2.0}, {3.0}};
    Vector<double, 3> vec1(mat1);
    EXPECT_DOUBLE_EQ(1.0, vec1(0, 0));
    EXPECT_DOUBLE_EQ(2.0, vec1(1, 0));
    EXPECT_DOUBLE_EQ(3.0, vec1(2, 0));

    // Test case 2: rvalue Matrix
    Vector<double, 3> vec2(Matrix<double, 3, 1>{{4.0}, {5.0}, {6.0}});
    EXPECT_DOUBLE_EQ(4.0, vec2(0, 0));
    EXPECT_DOUBLE_EQ(5.0, vec2(1, 0));
    EXPECT_DOUBLE_EQ(6.0, vec2(2, 0));
}

} // namespace gons