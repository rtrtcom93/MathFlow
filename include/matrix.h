#ifndef MATRIX_H
#define MATRIX_H

#include "pch.h"
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
    const T& operator()(size_t r, size_t c) const;
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
    Matrix<T> operator*(const Vector<T> &vec) const;
        //Commutativity
    template<typename U>
    friend Matrix<U> operator*(const Vector<U> &lhs, const Matrix<U> &rhs);
        //Method
    template<typename U>
    friend Matrix<U> matmul(const Vector<U> &lhs, const Matrix<U> &rhs);
    
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
        Matrix<T>& mat_ref;  // Matrix 객체에 대한 참조
        size_t col;               // 현재 열 번호
        size_t row_idx;           // 현재 행 번호

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

        // Dereference 연산자 (*)
        T& operator*() {
            return mat_ref.at(row_idx, col);
            // return mat_ref.at(row_idx, col);
        }

        // += 연산자: 특정 거리만큼 이동
        ColIterator& operator+=(size_t n) {
            row_idx += n;
            return *this;
        }

        // -= 연산자: 특정 거리만큼 뒤로 이동 (선택 사항)
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
            ++row_idx;  // 다음 행으로 이동
            return *this;
        }

        // Post-increment
        ColIterator operator++(int) {
            ColIterator temp(*this);
            ++row_idx;  // 다음 행으로 이동
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

        // 동등성 비교 연산자
        bool operator!=(const ColIterator& rhs) const {
            return row_idx != rhs.row_idx || col != rhs.col;
        }

        // 동등성 확인 연산자
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
    void push_back(const Matrix<T>& vec, int axis);
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