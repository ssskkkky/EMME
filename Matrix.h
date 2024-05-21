#ifndef MATRIX_H
#define MATRIX_H

#include <complex>
#include <stdexcept>
#include <tuple>
#include <vector>

#ifdef EMME_EXPRESSION_TEMPLATE
#include <functional>

#include "Arithmetics.h"

template <util::Indexable E1, util::Indexable E2>
auto operator-(const E1& e1, const E2& e2) {
    return util::BinaryExpression<E1, E2, std::minus<void>>{e1, e2};
}

template <util::Indexable E1, typename E2>
auto operator/(const E1& e1, const E2& e2) {
    return util::BinaryExpressionA2<E1, E2, std::divides<void>>{e1, e2};
}

#endif

template <typename T, typename A = std::allocator<T>>
class Matrix {
   public:
    using size_type = std::size_t;
    using value_type = T;
    using allocator_type = A;
    using matrix_type = Matrix<value_type, allocator_type>;

    // Constructor with dimensions
    Matrix(size_type rows, size_type cols)
        : rows_(rows), cols_(cols), data_(rows * cols) {}

    auto begin() const noexcept { return data_.begin(); }
    auto end() const noexcept { return data_.end(); }

    // Access element at (row, col)
    value_type& operator()(size_type row, size_type col) {
#ifdef EMME_DEBUG
        if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
            throw std::out_of_range("Matrix index out of bounds");
        }
#endif
        return data_[row * cols_ + col];
    }

    const value_type& operator()(size_type row, size_type col) const {
#ifdef EMME_DEBUG
        if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
            throw std::out_of_range("Matrix index out of bounds");
        }
#endif
        return data_[row * cols_ + col];
    }

    // Function to perform matrix subtraction
    matrix_type& operator-=(const matrix_type& other) {
#ifdef EMME_DEBUG
        // Check if dimensions are compatible
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument(
                "Matrices must have the same dimensions for subtraction.");
        }
#endif

        // Subtract corresponding elements from each matrix
        for (size_type i = 0; i < data_.size(); ++i) {
            data_[i] -= other.data_[i];
        }

        return *this;
    }

    // Function to perform division of the matrix by a scalar
    matrix_type& operator/=(const T& scalar) {
#ifdef EMME_DEBUG
        // Check for division by zero
        if (scalar.real() == 0 && scalar.imag() == 0) {
            throw std::invalid_argument("Division by zero is not allowed.");
        }
#endif
        // Divide each element by the scalar
        for (auto& v : data_) { v /= scalar; }
        return *this;
    }

#ifdef EMME_EXPRESSION_TEMPLATE
    explicit Matrix(util::Dimension<2> dim) : Matrix(dim.dim[0], dim.dim[1]) {}

    template <util::Indexable U>
    matrix_type& operator=(const U& other) {
        for (size_type i = 0; i < data_.size(); ++i) { data_[i] = other[i]; }
        return *this;
    }

    template <util::Indexable U>
    Matrix(const U& other) requires(
        !std::is_same_v<std::remove_cvref_t<U>, matrix_type>)
        : Matrix(other.get_dim()) {
        operator=(other);
    }

    auto get_dim() const {
        return util::Dimension<2>{{rows_, cols_}};
    }

    const value_type& operator[](size_type idx) const {
        return data_[idx];
    }

#else
    friend matrix_type operator-(matrix_type a, const matrix_type& b) {
        a -= b;
        return a;
    }
    friend matrix_type operator/(matrix_type m, const value_type& a) {
        m /= a;
        return m;
    }
#endif
    // Get number of rows
    size_type getRows() const {
        return rows_;
    }

    // Get number of columns
    size_type getCols() const {
        return cols_;
    }

    // Get a specific column
    const std::vector<T> getCol(size_type col) const {
#ifdef EMME_DEBUG
        if (col < 0 || col >= cols_) {
            throw std::out_of_range("Matrix column index out of bounds");
        }
#endif
        std::vector<T> col_view(rows_);
        for (size_type i = 0; i < rows_; ++i) { col_view[i] = (*this)(i, col); }
        return col_view;
    }

    // Get a specific row
    const std::vector<T> getRow(int row) const {
#ifdef EMME_DEBUG
        if (row < 0 || row >= rows_) {
            throw std::out_of_range("Matrix row index out of bounds");
        }
#endif
        return std::vector<T>(data_.begin() + row * cols_,
                              data_.begin() + (row + 1) * cols_);
    }

    // Assign values to a specific column
    void setCol(size_type col, const std::vector<value_type>& new_col) {
#ifdef EMME_DEBUG
        if (col < 0 || col >= cols_) {
            throw std::out_of_range("Matrix column index out of bounds");
        }
        if (new_col.size() != rows_) {
            throw std::invalid_argument("Column size mismatch");
        }
#endif
        for (size_type i = 0; i < rows_; ++i) {
            operator()(i, col) = new_col[i];
        }
    }

    // Assign values to a specific row
    void setRow(int row, const std::vector<value_type>& new_row) {
#ifdef EMME_DEBUG
        if (row < 0 || row >= rows_) {
            throw std::out_of_range("Matrix row index out of bounds");
        }
        if (new_row.size() != cols_) {
            throw std::invalid_argument("Row size mismatch");
        }
#endif
        std::copy(new_row.begin(), new_row.end(), data_.begin() + row * cols_);
    }

    // Calculate the trace of the matrix
    value_type trace() const {
#ifdef EMME_DEBUG
        if (rows_ != cols_) {
            throw std::invalid_argument(
                "Matrix must be square for trace calculation");
        }
#endif
        value_type sum{};
        for (size_type i = 0; i < rows_; ++i) {
            sum += operator()(i, i);  // Access diagonal element using (i, i)
        }
        return sum;
    }

    std::tuple<Matrix<T>, Matrix<T>, std::vector<int>> luDecomposition() const {
        int n = this->getRows();

        // Allocate memory for L, U, perm (assuming they have the same size as
        // J)
        Matrix<T> L(n, n);
        Matrix<T> U(n, n);
        std::vector<int> perm(n);
        auto Jv = (*this);

        for (int i = 0; i < n; ++i) { perm[i] = i; }

        // Perform LU decomposition (manual implementation)
        for (int i = 0; i < n; ++i) {
            // Partial pivoting (find the largest element in the current column
            // below the diagonal)
            int pivot = i;
            for (int j = i + 1; j < n; ++j) {
                if (std::abs(Jv(j, i)) > std::abs(Jv(pivot, i))) { pivot = j; }
            }

            // Swap rows if necessary (update permutation vector)
            if (pivot != i) {
                std::swap(perm[i], perm[pivot]);
                for (int j = 0; j < n; ++j) {
                    std::swap(Jv(i, j), Jv(pivot, j));
                }
            }

            // LU decomposition for the current column
            for (int j = i + 1; j < n; ++j) {
                std::complex<double> factor = Jv(j, i) / Jv(i, i);
                Jv(j, i) = factor;
                for (int k = i + 1; k < n; ++k) {
                    Jv(j, k) -= factor * Jv(i, k);
                }
            }
        }

        // Extract L and U from the modified Jv matrix
        for (int i = 0; i < n; ++i) {
            L(i, i) = 1;
            for (int j = 0; j < i; ++j) { L(i, j) = Jv(i, j); }
            for (int j = i; j < n; ++j) { U(i, j) = Jv(i, j); }
        }

        return std::make_tuple(L, U, perm);
    }

    auto data() {
        return data_.data();
    }

   private:
    size_type rows_;
    size_type cols_;
    std::vector<value_type, allocator_type> data_;
};

#endif
