#ifndef TENSOR_H
#define TENSOR_H

#include "pch.h"

template<typename T>
class Tensor {
protected:
    std::vector<T> vec;

public :
    //Overloaded constructors
    Tensor();
    Tensor(size_t size);
    Tensor(size_t size, T init_val);
    Tensor(std::initializer_list<T> init_val);
    Tensor(std::vector<T> init_vec);

    //Copy constructor
    Tensor(const Tensor<T> &v);

    //Move constructor
    Tensor(Tensor<T> &&source) noexcept;

    //Destructor
    ~Tensor();

    //Overloaded operators
    //Copy assignment
    Tensor<T> &operator=(const Tensor<T> &rhs);
    //Move assignment
    Tensor<T> &operator=(Tensor<T> &&rhs);
    Tensor<T> &operator=(T scalar);

    //Equality
    bool operator==(const Tensor<T> &rhs) const;

    //Accessor operators
    T& operator[](size_t index);
    const T& operator[](size_t index) const;

    //Extraction & Insertion operators 
    template<typename U>
    friend std::ostream &operator<<(std::ostream &os, const Tensor<U> &rhs);
    template<typename U>
    friend std::istream &operator>>(std::istream &is, Tensor<U> &rhs);

    //Arithmetic operators for a vector space
    //Vector addition
    Tensor<T> operator-() const;
    Tensor<T> operator+(const Tensor<T> &rhs) const;
    Tensor<T> operator-(const Tensor<T> &rhs) const;
    //Scalar-Vector addition (Special)
    Tensor<T> operator+(T Scalar) const;
    Tensor<T> operator-(T Scalar) const;
    //Commutativity for Scalar-Vector addition
    template<typename U>
    friend Tensor<U> operator+(U scalar, const Tensor<U> &rhs);
    template<typename U>
    friend Tensor<U> operator-(U scalar, const Tensor<U> &rhs);

    //Scalar muliplication
    Tensor<T> operator*(T scalar) const;
    Tensor<T> operator/(T scalar) const;
    //Commutativity for Scalar multiplication
    template<typename U>
    friend Tensor<U> operator*(U scalar, const Tensor<U> &rhs);

    //Product
    T dot(const Tensor<T> &rhs) const;
    template<typename U>
    friend double dot(const Tensor<U> &lhs, const Tensor<U> &rhs);

    //Display and Accessment options
    //Display
    size_t size();
    const size_t size() const ;
    void print() const;

    //Accessment
    T& at(size_t index);
    const T& at(size_t index) const;
    const std::vector<T>& get_vec() const;

    typename std::vector<T>::iterator begin();
    typename std::vector<T>::iterator end();
    typename std::vector<T>::const_iterator begin() const;
    typename std::vector<T>::const_iterator end() const;

    //Mathematical implementation
    T sum() const;
    T mag() const;
    Tensor<T> unit();
};

#include "tensor.hpp"

#endif