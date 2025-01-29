#ifndef TENSOR_HPP
#define TENSOR_HPP

constexpr const char* DIMENSION_ERROR = "The dimensions of vector spaces are not same.";
constexpr const char* DIVISION_ERROR  = "Division by zero";

/******************Overloaded constructors*******************/
//Default constructors
template<typename T>
Tensor<T>::Tensor() 
    : vec(1, 0) {
        // std::cout << "Default initialization" << std::endl;
    }

template<typename T>
Tensor<T>::Tensor(size_t size)
    : vec(size, 0) {
        // std::cout << "Initialize with size" << std::endl;
    }

template<typename T>
Tensor<T>::Tensor(size_t size, T init_val)
    : vec(size, init_val) {
        // std::cout << "Initialize with size and list" << std::endl;
    }
    
template<typename T>
Tensor<T>::Tensor(std::initializer_list<T> init_val) 
    : vec(init_val) {
        // std::cout << "Initialize with list" << std::endl;
    }

template<typename T>
Tensor<T>::Tensor(std::vector<T> init_vec) 
    : vec(init_vec) {
        // std::cout << "Initialize with vector" << std::endl;
    }

//Copy constructor
template<typename T>
Tensor<T>::Tensor(const Tensor<T> &source)
    : vec(source.vec) {
        // std::cout << "Copy constructor - deep" << std::endl;
    }

//Move constructor
template<typename T>
Tensor<T>::Tensor(Tensor<T> &&source) noexcept
    : vec(std::move(source.vec))  {
        std::cout << "Move constructor" << std::endl;
    }

//Destructor
template<typename T>
Tensor<T>::~Tensor() {
    // std::cout << "Destroyed" << std::endl;
}

/******************Overloaded operators*******************/
//Copy assignment
template<typename T>
Tensor<T> &Tensor<T>::operator=(const Tensor<T> &rhs) {
    //std::cout << "Copy assignment" << std::endl;
    if (this==&rhs) 
        return *this;
    this->vec = rhs.vec;
    return *this;
}

//Move assignement
template<typename T>
Tensor<T> &Tensor<T>::operator=(Tensor<T> &&rhs) {
    // std::cout << "Move assignment" << std::endl;
    if(this==&rhs)
        return *this;
    vec = std::move(rhs.vec); //Copy reference value and Transfer of ownership
    return *this;
}

template<typename T>
Tensor<T> &Tensor<T>::operator=(T scalar) {
    // std::cout << "Move assignment" << std::endl;
    Tensor<T> temp(vec.size(), scalar);
    vec = std::move(temp.vec); //Copy reference value and Transfer of ownership
    return *this;
}

//Equality
template<typename T>
bool Tensor<T>::operator==(const Tensor<T> &rhs) const {
    return (vec == rhs.vec);
}

//Accessor operators
template<typename T>
T &Tensor<T>::operator[](size_t index) {
    return vec[index];
}

template<typename T>
const T &Tensor<T>::operator[](size_t index) const {
    return vec[index];
}

//Overloaded insertion operator
template<typename U>
std::ostream &operator<<(std::ostream &os, const Tensor<U> &rhs) {
    os << "[ ";
    for (const U& v : rhs.vec) {
        os << v << " ";
    }
    os << "]";
    return os;
}

//Overloaded extraction operator
template<typename U>
std::istream &operator>>(std::istream &is, Tensor<U> &rhs) {
    std::vector<U> buff;
    U input;
    std::cout << "Enter value(s) of a vector (press Enter to finish): ";
    std::string line;
    if (std::getline(is, line)) {
        std::istringstream iss(line);
        while (iss >> input) {
            buff.push_back(input);
        }
    }

    rhs = Tensor<U>{buff};
    return is;
}

//Arithmetic operators for a vector space
//Vector addition
template<typename T>
Tensor<T> Tensor<T>::operator-() const {
    Tensor<T> temp(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) 
        temp[i] = -vec[i];
    return temp;
}

template<typename T>
Tensor<T> Tensor<T>::operator+(const Tensor<T> &rhs) const {
#ifdef DEBUG
    if (this->size() != rhs.size()) {
        throw std::invalid_argument(DIMENSION_ERROR);
    } 
#endif
    Tensor<T> temp(vec.size());
    for (size_t i = 0; i < this->vec.size(); ++i)
        temp[i] = vec[i] + rhs[i];
    return temp;
}

template<typename T>
Tensor<T> Tensor<T>::operator-(const Tensor<T> &rhs) const {
#ifdef DEBUG
    if (this->size() != rhs.size()) {
        throw std::invalid_argument(DIMENSION_ERROR);
    }
#endif
    Tensor<T> temp(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        temp[i] = vec[i] - rhs[i];
    return temp;
}

//Scalar-Vector addition (Special)
template<typename T>
Tensor<T> Tensor<T>::operator+(T scalar) const {
    Tensor<T> temp(vec.size());
    for (size_t i = 0; i < this->vec.size(); ++i)
        temp[i] = vec[i] + scalar;
    return temp;
}

template<typename T>
Tensor<T> Tensor<T>::operator-(T scalar) const {
    Tensor<T> temp(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        temp[i] = vec[i] - scalar;
    return temp;
}

//Commutativity for Scalar-Vector addition
template<typename U>
Tensor<U> operator+(U scalar, const Tensor<U> &rhs) {
    Tensor<U> temp(rhs.size());
    for (size_t i = 0; i < rhs.size(); ++i)
        temp[i] = scalar + rhs[i];
    return temp;
}

template<typename U>
Tensor<U> operator-(U scalar, const Tensor<U> &rhs) {
    Tensor<U> temp(rhs.size());
    for (size_t i = 0; i < rhs.size(); ++i)
        temp[i] = scalar - rhs[i];
    return temp;
}

//Scalar muliplication
template<typename T>
Tensor<T> Tensor<T>::operator*(T scalar) const {
    Tensor<T> temp(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        temp[i] = vec[i]*scalar;
    return temp;
}

template<typename T>
Tensor<T> Tensor<T>::operator/(T scalar) const {
#ifdef DEBUG
    if (scalar == 0) {
        throw std::invalid_argument(DIVISION_ERROR);
    }
#endif    
    Tensor<T> temp(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        temp[i] = vec[i]/scalar;
    return temp;
}

//Commutativity for Scalar multiplication
template<typename U>
Tensor<U> operator*(U scalar, const Tensor<U> &rhs) {
    Tensor<U> temp(rhs.size());
    for (size_t i = 0; i < rhs.size(); ++i)
        temp[i] = scalar*rhs[i];
    return temp;
}

//Dot product of 1D
template<typename T>
T Tensor<T>::dot(const Tensor<T> &rhs) const {
#ifdef DEBUG
    if (this->size() != rhs.size()) {
        throw std::invalid_argument(DIMENSION_ERROR);
    }
#endif
    T temp{0};
    for (size_t i = 0; i < vec.size(); ++i)
        temp = temp + vec[i]*rhs[i];
    return temp; 
}

template<typename U>
double dot(const Tensor<U> &lhs, const Tensor<U> &rhs) {
#ifdef DEBUG
    if (lhs.vec.size() != rhs.vec.size()) {
        throw std::invalid_argument(DIMENSION_ERROR);
    }
#endif
    U temp{0};
    for (size_t i = 0; i < lhs.vec.size(); ++i)
        temp = temp + lhs[i]*rhs[i];
    return temp;   
}

/******************Display and Accessment options*******************/
//Display member methods
template<typename T>
size_t Tensor<T>::size() {
    return vec.size();
}

template<typename T>
const size_t Tensor<T>::size() const {
    return vec.size();
}

template<typename T>
void Tensor<T>::print() const {
    std::cout << "[ ";
    for (const T& v : this->vec) {
        std::cout << v << " ";
    }
    std::cout << "]";
    std::cout << "\n";
}

/********************Accesment on data***********************/
//Access member methods
template<typename T>
T& Tensor<T>::at(size_t index) {
    return vec.at(index);
}

template<typename T>
const T& Tensor<T>::at(size_t index) const{
    return vec.at(index);
}

template<typename T>
const std::vector<T>& Tensor<T>::get_vec() const {
    return vec;
}

template<typename T>
typename std::vector<T>::iterator Tensor<T>::begin() {
    return vec.begin();
}

template<typename T>
typename std::vector<T>::iterator Tensor<T>::end() {
    return vec.end();
}

template<typename T>
typename std::vector<T>::const_iterator Tensor<T>::begin() const {
    return vec.begin();
}

template<typename T>
typename std::vector<T>::const_iterator Tensor<T>::end() const {
    return vec.end();
}

/******************Mathematical implementation************************/
//Magnitude
template<typename T>
T Tensor<T>::sum() const {
    T temp{0};
    for (size_t i = 0; i < vec.size(); ++i)
        temp += vec[i];
    return std::sqrt(temp);
}

//Magnitude
template<typename T>
T Tensor<T>::mag() const {
    T temp{0};
    for (size_t i = 0; i < vec.size(); ++i)
        temp += std::pow(vec[i], 2);
    return std::sqrt(temp);
}

template<typename T>
Tensor<T> Tensor<T>::unit() {
#ifdef DEBUG
    if (this.mag() == 0) 
        throw std::invalid_argument(DIVISION_ERROR);
#endif
    Tensor<T> temp{*this};
    T length = this->mag();
    return temp/length;
}

#endif