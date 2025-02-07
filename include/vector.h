#ifndef VECTOR_H
#define VECTOR_H

#include "pch.h"
#include "slice.h"

template<typename T>
class Vector {
protected:
    std::vector<T> vec;
    size_t dim;
    
public :
    //Overloaded constructors
    Vector();
    Vector(size_t size);
    Vector(size_t size, T init_val);
    Vector(std::initializer_list<T> init_val);
    Vector(std::vector<T> init_vec);

    //Copy constructor
    Vector(const Vector<T> &v);

    //Move constructor
    Vector(Vector<T> &&source) noexcept;

    //Destructor
    ~Vector();

    //Overloaded operators
    //Copy assignment
    Vector<T> &operator=(const Vector<T> &rhs);
    //Move assignment
    Vector<T> &operator=(Vector<T> &&rhs);
    Vector<T> &operator=(T scalar);

    //In&Equality
    bool operator==(const Vector<T> &rhs) const;
    bool operator!=(const Vector<T> &rhs) const;

    //Accessor operators
    class SliceProxy {
    private:
        Vector<T>& vec_ref;
        Slice slc;

    public:
        SliceProxy(Vector<T>& vec_ref, const Slice& slc)
            : vec_ref(vec_ref), slc(slc) {}

        operator Vector<T>() const {
            size_t slice_size=0;
            Vector<T> result;
            size_t k = 0;
            for (size_t i = slc.start; i < slc.end; i += slc.step) {
                result[k++] =  vec_ref(i);
            }
            
            return result;
        }

        friend std::ostream& operator<<(std::ostream& os, const SliceProxy& proxy) {
            os << "[ ";
            for (size_t i = proxy.slc.start; i < proxy.slc.end; i+= proxy.slc.step){
                os << proxy.vec_ref(i) << " ";
            }
            os << "]";
            return os;
        }

        // 단일 값을 할당
        SliceProxy& operator=(const T& value) {
            for (size_t i = slc.start; i < slc.end; i += slc.step) {
                vec_ref[i] = value;
            }
            return *this;
        }

        // 다른 벡터 값을 슬라이스에 할당
        SliceProxy& operator=(const Vector<T>& rhs) {
            size_t k = 0;
            size_t expected_size = (slc.end - slc.start + slc.step - 1) / slc.step;
            if (rhs.size() != expected_size) {
                throw std::invalid_argument("Size mismatch between the slice and the assigned vector.");
            }

            for (size_t i = slc.start; i < slc.end; i += slc.step) {
                vec_ref[i] = rhs[k++];
            }
            return *this;
        }
    };

    T& operator()(size_t idx);
    const T& operator()(size_t idx) const;

    SliceProxy operator()(const Slice& slc) {
    #ifdef DEBUG
        slc.check(this->dim);
    #endif
        return SliceProxy(*this, slc);
    }
    
    T& operator[](size_t idx);
    const T& operator[](size_t idx) const;

    //Extraction & Insertion operators 
    template<typename U>
    friend std::ostream &operator<<(std::ostream &os, const Vector<U> &rhs);
    template<typename U>
    friend std::istream &operator>>(std::istream &is, Vector<U> &rhs);

    //Arithmetic operators for a vector space
    //Vector addition
    Vector<T> operator-() const;
    Vector<T> operator+(const Vector<T> &rhs) const;
    Vector<T> operator-(const Vector<T> &rhs) const;
    //Scalar-Vector addition (Special)
    Vector<T> operator+(T Scalar) const;
    Vector<T> operator-(T Scalar) const;
    //Commutativity for Scalar-Vector addition
    template<typename U>
    friend Vector<U> operator+(U scalar, const Vector<U> &rhs);
    template<typename U>
    friend Vector<U> operator-(U scalar, const Vector<U> &rhs);

    //Scalar muliplication
    Vector<T> operator*(T scalar) const;
    Vector<T> operator/(T scalar) const;
    //Commutativity for Scalar multiplication
    template<typename U>
    friend Vector<U> operator*(U scalar, const Vector<U> &rhs);

    //Product
    T dot(const Vector<T> &rhs) const;
    template<typename U>
    friend double dot(const Vector<U> &lhs, const Vector<U> &rhs);

    //Display and Accessment options
    //Display
    size_t size() const ;
    size_t capacity() const;
    void print() const;

    //Accessment
    T& at(size_t index);
    const T& at(size_t index) const;
    const std::vector<T>& get_vec() const;

    typename std::vector<T>::iterator begin();
    typename std::vector<T>::iterator end();
    typename std::vector<T>::const_iterator begin() const;
    typename std::vector<T>::const_iterator end() const;

    //Dynamic extension
    void push_back(const T& value);
    void pop_back();

    template<typename... Args>
    void emplace_back(Args&&... args) {
        vec.emplace_back(std::forward<Args>(args)...);
    }

    void resize(size_t new_size, T init_val = 0);
    void reserve(size_t new_size);
    void clear();
    typename std::vector<T>::iterator insert(size_t pos, const T& value);
    typename std::vector<T>::iterator erase(size_t pos);
    
    //Mathematical implementation
    T sum() const;
    T mag() const;
    Vector<T> unit();
};

#include "vector.hpp"

#endif