#ifndef ARITHMETICS_H
#define ARITHMETICS_H

#include <concepts>

namespace util {

template <typename A>
concept Indexable = requires(const A& a, std::size_t idx) {
    a[idx];
};

template <typename A>
concept number = std::is_scalar_v<A> || std::same_as<A, std::complex<float>> ||
    std::same_as<A, std::complex<double>>;

// serve as a tag
struct ExpressionTemplate {};

template <typename T>
concept expression_template = std::is_base_of_v<ExpressionTemplate, T> &&
    requires(const T& t, std::size_t idx) {
    t[idx];
};

template <typename E1, typename E2, typename Op>
struct BinaryExpressionBase : ExpressionTemplate {
    const E1& e1;
    const E2& e2;

    BinaryExpressionBase(const E1& lhs, const E2& rhs) : e1(lhs), e2(rhs) {}
};

template <Indexable E1, Indexable E2, typename Op>
struct BinaryExpression : BinaryExpressionBase<E1, E2, Op> {
    using BinaryExpressionBase<E1, E2, Op>::BinaryExpressionBase;

    auto operator[](std::size_t idx) const {
        return Op{}(this->e1[idx], this->e2[idx]);
    }
};

template <Indexable E1, typename E2, typename Op>
struct BinaryExpressionA2 : BinaryExpressionBase<E1, E2, Op> {
    using BinaryExpressionBase<E1, E2, Op>::BinaryExpressionBase;
    auto operator[](std::size_t idx) const {
        return Op{}(this->e1[idx], this->e2);
    }
};
template <typename E1, Indexable E2, typename Op>
struct BinaryExpressionA1 : BinaryExpressionBase<E1, E2, Op> {
    using BinaryExpressionBase<E1, E2, Op>::BinaryExpressionBase;
    auto operator[](std::size_t idx) const {
        return Op{}(this->e1, this->e2[idx]);
    }
};

template <expression_template E1, expression_template E2>
auto operator+(const E1& e1, const E2& e2) {
    return BinaryExpression<E1, E2, std::plus<void>>{e1, e2};
}

template <expression_template E1, expression_template E2>
auto operator-(const E1& e1, const E2& e2) {
    return BinaryExpression<E1, E2, std::minus<void>>{e1, e2};
}

template <number E1, expression_template E2>
auto operator*(const E1& e1, const E2& e2) {
    return BinaryExpressionA1<E1, E2, std::multiplies<void>>{e1, e2};
}
template <expression_template E1, number E2>
auto operator*(const E1& e1, const E2& e2) {
    return BinaryExpressionA2<E1, E2, std::multiplies<void>>{e1, e2};
}

template <expression_template E1, number E2>
auto operator/(const E1& e1, const E2& e2) {
    return BinaryExpressionA2<E1, E2, std::divides<void>>{e1, e2};
}

}  // namespace util

#endif  // ARITHMETICS_H
