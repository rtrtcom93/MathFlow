#ifndef VECTOR_SPACE_HPP
#define VECTOR_SPACE_HPP

//Overloaded constructors
template<typename T>
VectorSpace<T>::VectorSpace() 
    : dim{1} {
        vec = new T[dim];
    }

template<typename T>
VectorSpace<T>::VectorSpace(size_t dim)
    : dim{dim} {
        vec = new T[dim];
        for (size_t i = 0; i < dim; ++i) {
            vec[i] = 0;
        }
    }

template<typename T>
VectorSpace<T>::VectorSpace(std::initializer_list<T> init_val) 
    : dim{init_val.size()} {
        vec = new T[dim];
        size_t i = 0;
        for (const T& val : init_val) {
            vec[i++] = val;
        }
    }

template<typename T>
VectorSpace<T>::VectorSpace(size_t dim, std::initializer_list<T> init_val) 
    : dim{dim} {
        vec = new T[dim];
        if (init_val.size() <= dim) {
            size_t i = 0;
            for (const T& val : init_val) {
                vec[i++] = val;
            }
        } else if (init_val.size() > dim) {
            throw std::invalid_argument("Initializer list size exceeds specified dimension.");
        }
    }

//Copy constructor
template<typename T>
VectorSpace<T>::VectorSpace(const VectorSpace<T> &source)
    : dim{source.dim} {
        vec = new T[dim];
        for (size_t i = 0; i < dim; ++i) {
            vec[i] = source.vec[i];
        }
        // std::cout << "Copy constructor - deep" << std::endl;
    }

//Move constructor
template<typename T>
VectorSpace<T>::VectorSpace(VectorSpace<T> &&source) noexcept
    : dim(source.dim), 
      vec(source.vec)  {
        source.dim=0;
        source.vec=nullptr;
        // std::cout << "Move constructor" << std::endl;
    }

//Destructor
template<typename T>
VectorSpace<T>::~VectorSpace() {
    // if (vec !=nullptr){
    //     std::cout << "Destructor freeing vector" << std::endl;
    // } else {
    //     std::cout << "Destructor freeing vector for nullptr" << std::endl;
    // }
    delete[] vec;
}

//Display options
template<typename T>
size_t VectorSpace<T>::get_dim() {
        return dim;
}

template<typename T>
T VectorSpace<T>::get_val() const
{
    return vec;
}

template<typename T>
void  VectorSpace<T>::print_dim()
{
    std::cout<< dim << std::endl;
}

template<typename T>
void VectorSpace<T>::print_val()
{
    for (size_t i=0; i < dim; ++i)
    {
        std::cout << vec[i] << " ";
    }
    std::cout << "\n";
}

template<typename T>
void VectorSpace<T>::print_val(const size_t &d, const T &v)
{
    for (size_t i=0; i < dim; ++i)
    {
        std::cout << v[i] << " ";
    }
    std::cout << "\n";
}

//Overloaded operators
//Copy assignment
template<typename T>
VectorSpace<T> &VectorSpace<T>::operator=(const VectorSpace<T> &rhs) {
    std::cout << "Copy assignment" << std::endl;
    if (this==&rhs) 
        return *this;
    delete[] this->vec;
    this->dim = rhs.dim;
    this->vec = new T[this->dim];
    for (size_t i = 0; i < dim; ++i) {
        vec[i] = rhs.vec[i];
    }
    return *this;
}

//Move assignement
template<typename T>
VectorSpace<T> &VectorSpace<T>::operator=(VectorSpace<T> &&rhs) {
    std::cout << "Move assignment" << std::endl;
    if(this==&rhs)
        return *this;
    delete[] this->vec;
    vec = rhs.vec;
    rhs.vec = nullptr;
    return *this;
}

#endif