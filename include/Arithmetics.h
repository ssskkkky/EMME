#ifndef ARITHMETICS_H
#define ARITHMETICS_H

namespace util {

template <typename A>
concept Indexable = requires(const A& a, std::size_t idx) {
    a[idx];
};

template <std::size_t d>
struct Dimension {
    std::array<std::size_t, d> dim;
};

template <typename E1, typename E2, typename Op>
struct BinaryExpressionBase {
    const E1& e1;
    const E2& e2;
};

template <Indexable E1, Indexable E2, typename Op>
struct BinaryExpression : BinaryExpressionBase<E1, E2, Op> {
    auto operator[](std::size_t idx) const {
        return Op{}(this->e1[idx], this->e2[idx]);
    }
    auto get_dim() const { return this->e1.get_dim(); }
};

template <Indexable E1, typename E2, typename Op>
struct BinaryExpressionA2 : BinaryExpressionBase<E1, E2, Op> {
    auto operator[](std::size_t idx) const {
        return Op{}(this->e1[idx], this->e2);
    }
    auto get_dim() const { return this->e1.get_dim(); }
};
template <typename E1, Indexable E2, typename Op>
struct BinaryExpressionA1 : BinaryExpressionBase<E1, E2, Op> {
    auto operator[](std::size_t idx) const {
        return Op{}(this->e1, this->e2[idx]);
    }
    auto get_dim() const { return this->e2.get_dim(); }
};

}  // namespace util

#endif  // ARITHMETICS_H
