#ifndef VECTOR_SPACE_H
#define VECTOR_SPACE_H

#include <iostream>
#include <initializer_list>

template<typename T>
class VectorSpace 
{
private:
    size_t dim;

protected :
    T *vec {nullptr};

public :
    //Overloaded constructors
    VectorSpace();
    VectorSpace(size_t dim);
    VectorSpace(std::initializer_list<T> init_val);
    VectorSpace(size_t dim, std::initializer_list<T> init_val); 

    //Copy constructor
    VectorSpace(const VectorSpace &v );

    //Move constructor
    VectorSpace(VectorSpace<T> &&source) noexcept;

    //Destructor
    ~VectorSpace();

    //Display options
    size_t get_dim();
    T get_val() const;
    void print_dim();
    void print_val();
    void print_val(const size_t &d, const T &v);

    //Overloaded operators
    //Copy assignment
    VectorSpace<T> &operator=(const VectorSpace<T> &rhs);
    VectorSpace<T> &operator=(VectorSpace<T> &&rhs)

};

#include "vector_space.hpp"

#endif