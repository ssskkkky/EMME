#ifndef MATRIX_H  // Replace MATRIX_H with your unique guard macro name
#define MATRIX_H

#include <complex>
#include <stdexcept>
#include <tuple>
#include <vector>

template <typename T>
class Matrix {
   public:
    // Constructor with dimensions
    Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
        data_ = std::vector<T>(rows * cols);
    }

    const T* begin() const { return data_.data(); }
    const T* end() const { return data_.data() + data_.size(); }

    // Access element at (row, col)
    T& operator()(unsigned int row, unsigned int col) {
        if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
            throw std::out_of_range("Matrix index out of bounds");
        }
        return data_[row * cols_ + col];
    }

    const T& operator()(unsigned int row, unsigned int col) const {
        if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
            throw std::out_of_range("Matrix index out of bounds");
        }
        return data_[row * cols_ + col];
    }

    // Function to perform matrix subtraction
    Matrix<T> operator-(const Matrix<T>& other) const {
        // Check if dimensions are compatible
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument(
                "Matrices must have the same dimensions for subtraction.");
        }

        // Create a result matrix with the same dimensions
        Matrix<T> result(rows_, cols_);

        // Subtract corresponding elements from each matrix
        for (unsigned int i = 0; i < rows_; ++i) {
            for (unsigned int j = 0; j < cols_; ++j) {
                result(i, j) =
                    data_[i * cols_ + j] - other.data_[i * cols_ + j];
            }
        }

        return result;
    }

    // Function to perform division of the matrix by a scalar
    Matrix<T> operator/(const T& scalar) const {
        // Check for division by zero
        if (scalar.real() == 0 && scalar.imag() == 0) {
            throw std::invalid_argument("Division by zero is not allowed.");
        }

        // Create a result matrix with the same dimensions
        Matrix<T> result(rows_, cols_);

        // Divide each element by the scalar
        for (unsigned int i = 0; i < rows_; ++i) {
            for (unsigned int j = 0; j < cols_; ++j) {
                result(i, j) = data_[i * cols_ + j] / scalar;
            }
        }

        return result;
    }

    // Element-wise matrix multiplication (overloaded operator *)
    Matrix<T> operator*(const Matrix<T>& other) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument(
                "Matrices must have the same dimensions for element-wise "
                "multiplication");
        }
        Matrix<T> result(rows_, cols_);
        for (unsigned int i = 0; i < rows_; ++i) {
            for (unsigned int j = 0; j < cols_; ++j) {
                result(i, j) = (*this)(i, j) * other(i, j);
            }
        }
        return result;
    }

    // Get number of rows
    int getRows() const { return rows_; }

    // Get number of columns
    int getCols() const { return cols_; }

    // Get a reference to a specific column
    const std::vector<T> getCol(unsigned int col) const {
        if (col < 0 || col >= cols_) {
            throw std::out_of_range("Matrix column index out of bounds");
        }

        std::vector<T> col_view(rows_);
        for (unsigned int i = 0; i < rows_; ++i) {
            col_view[i] = (*this)(i, col);
        }
        return col_view;
    }

    // // Get a reference to a specific column
    // std::vector<T>& getCol(int col) {
    //     if (col < 0 || col >= cols_) {
    //         throw std::out_of_range("Matrix column index out of bounds");
    //     }

    //     std::vector<T> col_view(rows_);
    //     for (int i = 0; i < rows_; ++i) { col_view[i] = (*this)(i, col); }
    //     return col_view;
    // }

    // Get a reference to a specific row
    const std::vector<T>& getRow(int row) const {
        if (row < 0 || row >= rows_) {
            throw std::out_of_range("Matrix row index out of bounds");
        }

        // Return a reference to a sub-vector of the data_ vector starting at
        // row*cols_
        return std::vector<T>(data_.begin() + row * cols_,
                              data_.begin() + (row + 1) * cols_);
    }

    // Assign values to a specific column
    void setCol(unsigned int col, const std::vector<T>& new_col) {
        if (col < 0 || col >= cols_) {
            throw std::out_of_range("Matrix column index out of bounds");
        }
        if (new_col.size() != rows_) {
            throw std::invalid_argument("Column size mismatch");
        }

        for (unsigned int i = 0; i < rows_; ++i) {
            (*this)(i, col) = new_col[i];
        }
    }

    // Assign values to a specific row
    void setRow(int row, const std::vector<T>& new_row) {
        if (row < 0 || row >= rows_) {
            throw std::out_of_range("Matrix row index out of bounds");
        }
        if (new_row.size() != cols_) {
            throw std::invalid_argument("Row size mismatch");
        }

        std::copy(new_row.begin(), new_row.end(), data_.begin() + row * cols_);
    }

    // Calculate the trace of the matrix
    T trace() const {
        if (rows_ != cols_) {
            throw std::invalid_argument(
                "Matrix must be square for trace calculation");
        }

        T sum = 0.0;
        for (unsigned int i = 0; i < rows_; ++i) {
            sum += (*this)(i, i);  // Access diagonal element using (i, i)
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

   private:
    unsigned int rows_;
    unsigned int cols_;
    std::vector<T> data_;
};

#endif
