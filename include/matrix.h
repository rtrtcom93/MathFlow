#ifndef MATRIX_H
#define MATRIX_H

#include "pch.h"
#include "slice.h"
#include "vector.h"

template<typename T>
class Matrix {
protected:
    Vector<T> mat;
    size_t rows;
    size_t cols;

public :
    //Overloaded consturctors
    //Initializer Constructors
    Matrix();
    Matrix(size_t nrow);
    Matrix(size_t nrow, size_t ncol);
    Matrix(size_t nrow, size_t ncol, T init_val);
    Matrix(std::initializer_list<std::initializer_list<T>> init_val);
    Matrix(std::vector<std::vector<T>> init_vec);

    //Copy constructor
    Matrix(const Matrix<T> &source);

    //Move constructor
    Matrix(Matrix<T> &&source) noexcept;

    //Destructor
    ~Matrix();

    //Overloaded operators
    //Copy assignment
    Matrix<T> &operator=(const Matrix<T> &rhs);
    //Move assignment
    Matrix<T> &operator=(Matrix<T> &&rhs);
    Matrix<T> &operator=(T scalar);

    //In&Equality
    bool operator==(const Matrix<T> &rhs) const;
    bool operator!=(const Matrix<T> &rhs) const;

    //Accessor operators
    T& operator()(size_t row, size_t col);
    const T& operator()(size_t row, size_t col) const;

    class SliceProxy {
    private:
        Matrix<T>& mat_ref;
        Slice row_slc;
        Slice col_slc;
        bool is_row_fixed = false;
        bool is_col_fixed = false;
        size_t fixed_index = 0;  // 고정된 행 또는 열의 인덱스

    public:
        // 행과 열 슬라이스에 대한 일반 생성자
        SliceProxy(Matrix<T>& mat_ref, const Slice& row_slc, const Slice& col_slc)
            : mat_ref(mat_ref), row_slc(row_slc), col_slc(col_slc) {}

        // 행 고정에 대한 생성자
        SliceProxy(Matrix<T>& mat_ref, size_t row, const Slice& col_slc)
            : mat_ref(mat_ref), col_slc(col_slc), is_row_fixed(true), fixed_index(row) {}

        // 열 고정에 대한 생성자
        SliceProxy(Matrix<T>& mat_ref, const Slice& row_slc, size_t col)
            : mat_ref(mat_ref), row_slc(row_slc), is_col_fixed(true), fixed_index(col) {}

        // 단일 값을 할당
        SliceProxy& operator=(const T& value) {
            if (is_row_fixed) {
                for (size_t j = col_slc.start; j < col_slc.end; j += col_slc.step) {
                    mat_ref(fixed_index, j) = value;  // 고정된 행에서 열 슬라이스 할당
                }
            } else if (is_col_fixed) {
                for (size_t i = row_slc.start; i < row_slc.end; i += row_slc.step) {
                    mat_ref(i, fixed_index) = value;  // 고정된 열에서 행 슬라이스 할당
                }
            } else {
                for (size_t i = row_slc.start; i < row_slc.end; i += row_slc.step) {
                    for (size_t j = col_slc.start; j < col_slc.end; j += col_slc.step) {
                        mat_ref(i, j) = value;
                    }
                }
            }
            return *this;
        }

        // 고정된 행 또는 열에 대해 벡터 값을 할당
        SliceProxy& operator=(const Vector<T>& vec) {
            if (is_row_fixed) {
                // 크기 검증: 열 슬라이스의 크기와 벡터 크기가 일치해야 함
                size_t expected_size = (col_slc.end - col_slc.start + col_slc.step - 1) / col_slc.step;
                if (vec.size() != expected_size) {
                    throw std::invalid_argument("Size mismatch between the column slice and the assigned vector.");
                }

                // 값 할당
                size_t k = 0;
                for (size_t j = col_slc.start; j < col_slc.end; j += col_slc.step) {
                    mat_ref(fixed_index, j) = vec[k++];
                }
            } else if (is_col_fixed) {
                // 크기 검증: 행 슬라이스의 크기와 벡터 크기가 일치해야 함
                size_t expected_size = (row_slc.end - row_slc.start + row_slc.step - 1) / row_slc.step;
                if (vec.size() != expected_size) {
                    throw std::invalid_argument("Size mismatch between the row slice and the assigned vector.");
                }

                // 값 할당
                size_t k = 0;
                for (size_t i = row_slc.start; i < row_slc.end; i += row_slc.step) {
                    mat_ref(i, fixed_index) = vec[k++];
                }
            } else {
                throw std::logic_error("Vector assignment is only allowed when a row or column is fixed.");
            }

            return *this;
        }

        // 다른 행렬 값을 슬라이스에 할당
        SliceProxy& operator=(const Matrix<T>& rhs) {
            size_t row_count = (row_slc.end - row_slc.start + row_slc.step - 1) / row_slc.step;
            size_t col_count = (col_slc.end - col_slc.start + col_slc.step - 1) / col_slc.step;

            if (rhs.nrow() != row_count || rhs.ncol() != col_count) {
                throw std::invalid_argument("Size mismatch between the slice and the assigned matrix.");
            }

            if (is_row_fixed) {
                for (size_t j = 0, col = col_slc.start; col < col_slc.end; col += col_slc.step, ++j) {
                    mat_ref(fixed_index, col) = rhs(0, j);  // 고정된 행에 대해 복사
                }
            } else if (is_col_fixed) {
                for (size_t i = 0, row = row_slc.start; row < row_slc.end; row += row_slc.step, ++i) {
                    mat_ref(row, fixed_index) = rhs(i, 0);  // 고정된 열에 대해 복사
                }
            } else {
                for (size_t i = 0, row = row_slc.start; row < row_slc.end; row += row_slc.step, ++i) {
                    for (size_t j = 0, col = col_slc.start; col < col_slc.end; col += col_slc.step, ++j) {
                        mat_ref(row, col) = rhs(i, j);
                    }
                }
            }
            return *this;
        }
    };

    SliceProxy operator()(const Slice& row_slc, size_t col) {
    #ifdef DEBUG
        row_slc.check(this->rows);
        if (col >= this->cols)
            throw std::out_of_range("Column index out of range");
    #endif
        return SliceProxy(*this, row_slc, col);
    }
    
    SliceProxy operator()(size_t row, const Slice& col_slc) {
    #ifdef DEBUG
        col_slc.check(this->cols);
        if (row >= this->rows)
            throw std::out_of_range("Row index out of range");
    #endif
        return SliceProxy(*this, row, col_slc);
    }

    SliceProxy operator()(const Slice& row_slc, const Slice& col_slc) {
    #ifdef DEBUG
        row_slc.check(this->rows);
        col_slc.check(this->cols);
    #endif
        return SliceProxy(*this, row_slc, col_slc);
    }

    T* operator[](size_t row);
    const T* operator[](size_t row) const;
   
    //Overloaded insertion operator
    template<typename U>
    friend std::ostream &operator<<(std::ostream &os, const Matrix<U> &rhs);
    template<typename U>
    friend std::istream &operator>>(std::istream &is, Matrix<U> &rhs);

    //Arithmetic operators for a vector space
    //Matrix addition
    Matrix<T> operator-() const;
    Matrix<T> operator+(const Matrix<T> &rhs) const;
    Matrix<T> operator-(const Matrix<T> &rhs) const;

    //Scalar-Matrix addition (Special)
    Matrix<T> operator+(T Scalar) const;
    Matrix<T> operator-(T Scalar) const;
        //Commutativity 
    template<typename U>
    friend Matrix<U> operator+(U scalar, const Matrix<U> &rhs);
    template<typename U>
    friend Matrix<U> operator-(U scalar, const Matrix<U> &rhs);

    //Multiplication
    //Scalar muliplication
    Matrix<T> operator*(T scalar) const;
    Matrix<T> operator/(T scalar) const;
        //Commutativity
    template<typename U>
    friend Matrix<U> operator*(U scalar, const Matrix<U> &rhs);

    //Vector-Matrix multiplication
    Vector<T> operator*(const Vector<T> &vec) const;
        //Commutativity
    template<typename U>
    friend Vector<U> operator*(const Vector<U> &lhs, const Matrix<U> &rhs);
        //Method
    template<typename U>
    friend Vector<U> matmul(const Vector<U> &lhs, const Matrix<U> &rhs);
    
    //Matrix-Matrix multiplication
    Matrix<T> operator*(const Matrix<T> &rhs) const;
        //Method
    template<typename U>
    friend Matrix<U> matmul(const Matrix<U> &lhs, const Matrix<U> &rhs);

    //Display and Accessment options
    size_t size() const;
    size_t capacity() const;
    size_t nrow() const;
    size_t ncol() const;
    void shape() const;
    void print() const;

    //Accessment
    T& at(size_t row, size_t col);
    const T& at(size_t row, size_t col) const;
    const Matrix<T>& get_mat() const;

    typename std::vector<T>::iterator begin();
    typename std::vector<T>::iterator end();
    typename std::vector<T>::const_iterator begin() const;
    typename std::vector<T>::const_iterator end() const;

    typename std::vector<T>::iterator row_begin(size_t row);
    typename std::vector<T>::iterator row_end(size_t row);
    typename std::vector<T>::const_iterator row_begin(size_t row) const;
    typename std::vector<T>::const_iterator row_end(size_t row) const;

    class ColIterator {
    private:
        Matrix<T>& mat_ref;       // Referencing Matrix Object 
        size_t col;               // selected column
        size_t row_idx;           // row index

    public:
        // ColIterator constructor
        ColIterator(Matrix<T>& mat, size_t col, size_t row_idx)
            : mat_ref(mat), col(col), row_idx(row_idx) {}

        // Copy constructor
        ColIterator(const ColIterator& source)
            : mat_ref(source.mat_ref), col(source.col), row_idx(source.row_idx) {}
        
        // Move constructor
        ColIterator(const ColIterator&& source) noexcept
            : mat_ref(source.mat_ref), col(source.col), row_idx(source.row_idx) {}

        // Overloaded dereference operator(*)
        T& operator*() {
            return mat_ref.at(row_idx, col);
            // return mat_ref.at(row_idx, col);
        }

        // += operator : Move forward with n
        ColIterator& operator+=(size_t n) {
            row_idx += n;
            return *this;
        }

        // -= operator : Move backward with n
        ColIterator& operator-=(size_t n) {
            row_idx -= n;
            return *this;
        }

        // + 연산자
        ColIterator operator+(size_t n) const {
            ColIterator temp(*this);
            temp += n;
            return temp;
        }

        // - 연산자 (선택 사항)
        ColIterator operator-(size_t n) const {
            ColIterator temp(*this);
            temp -= n;
            return temp;
        }

        // Pre-increment 연산자 (++iterator)
        ColIterator& operator++() {
            ++row_idx;  // Move to next row
            return *this;
        }

        // Post-increment
        ColIterator operator++(int) {
            ColIterator temp(*this);
            ++row_idx;  // Move to next row
            return temp;
        }

        // Comparison operators (<, >, <=, >=)
        bool operator<(const ColIterator& other) const {
            return row_idx < other.row_idx;
        }

        bool operator>(const ColIterator& other) const {
            return row_idx > other.row_idx;
        }

        bool operator<=(const ColIterator& other) const {
            return row_idx <= other.row_idx;
        }

        bool operator>=(const ColIterator& other) const {
            return row_idx >= other.row_idx;
        }

        /*Overloaded equality operators*/
        bool operator!=(const ColIterator& rhs) const {
            return row_idx != rhs.row_idx || col != rhs.col;
        }

        bool operator==(const ColIterator& rhs) const {
            return !(*this != rhs);
        }
    };

    class ConstColIterator {
    private:
        const Matrix<T>& mat_ref;  // Matrix 객체에 대한 상수 참조
        size_t col;                // 현재 열 번호
        size_t row_idx;            // 현재 행 번호

    public:
        // Constructor
        ConstColIterator(const Matrix<T>& mat, size_t col, size_t row_idx)
            : mat_ref(mat), col(col), row_idx(row_idx) {}

        // Copy constructor
        ConstColIterator(const ConstColIterator& source)
            : mat_ref(source.mat_ref), col(source.col), row_idx(source.row_idx) {}

        // Move constructor
        ConstColIterator(ConstColIterator&& source) noexcept
            : mat_ref(source.mat_ref), col(source.col), row_idx(source.row_idx) {}

        // Dereference operator (*)
        const T operator*() const {
            return mat_ref.at(row_idx, col);
        }

        // += 연산자: 특정 거리만큼 이동
        ConstColIterator& operator+=(size_t n) {
            row_idx += n;
            return *this;
        }

        // -= 연산자: 특정 거리만큼 뒤로 이동
        ConstColIterator& operator-=(size_t n) {
            row_idx -= n;
            return *this;
        }

        // + 연산자
        ConstColIterator operator+(size_t n) const {
            ConstColIterator temp(*this);
            temp += n;
            return temp;
        }

        // - 연산자
        ConstColIterator operator-(size_t n) const {
            ConstColIterator temp(*this);
            temp -= n;
            return temp;
        }

        // Pre-increment 연산자 (++iterator)
        ConstColIterator& operator++() {
            ++row_idx;
            return *this;
        }

        // Post-increment 연산자 (iterator++)
        ConstColIterator operator++(int) {
            ConstColIterator temp(*this);
            ++row_idx;
            return temp;
        }

        // Comparison operators (<, >, <=, >=)
        bool operator<(const ConstColIterator& other) const {
            return row_idx < other.row_idx;
        }

        bool operator>(const ConstColIterator& other) const {
            return row_idx > other.row_idx;
        }

        bool operator<=(const ConstColIterator& other) const {
            return row_idx <= other.row_idx;
        }

        bool operator>=(const ConstColIterator& other) const {
            return row_idx >= other.row_idx;
        }

        // Equality operators
        bool operator!=(const ConstColIterator& rhs) const {
            return row_idx != rhs.row_idx || col != rhs.col;
        }

        bool operator==(const ConstColIterator& rhs) const {
            return !(*this != rhs);
        }
    };

    // Column-wise 순회를 위한 시작(begin) 및 끝(end) Iterator
    ColIterator col_begin(size_t col) {
        if (col >= this->ncol()) {
            throw std::out_of_range("Column index out of range");
        }
        return ColIterator(*this, col, 0);  // 첫 번째 행에서 시작
    }

    ColIterator col_end(size_t col) {
        if (col >= this->ncol()) {
            throw std::out_of_range("Column index out of range");
        }
        return ColIterator(*this, col, this->nrow());  // 마지막 행 다음
    }

    ConstColIterator col_begin(size_t col) const {
        if (col >= this->ncol()) {
            throw std::out_of_range("Column index out of range");
        }
        return ConstColIterator(*this, col, 0);  // 첫 번째 행에서 시작
    }

    ConstColIterator col_end(size_t col) const {
        if (col >= this->ncol()) {
            throw std::out_of_range("Column index out of range");
        }
        return ConstColIterator(*this, col, this->nrow());  // 마지막 행 다음
    }

    //Dynamical extension
    void push_back(const Vector<T>& vec, int axis);
    void push_back(const Matrix<T>& rhs, int axis);
    void pop_back(int axis);
    void clear();

    //Mathematical implementation
    Matrix<T> transpose() const;
    Matrix<T> diag() const;
    Matrix<T> upp_diag(size_t upper) const;
    Matrix<T> low_diag(size_t lower) const;
    T trace() const;

    T det() const;
};

#include "matrix.hpp"

#endif