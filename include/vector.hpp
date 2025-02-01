#ifndef VECTOR_HPP
#define VECTOR_HPP

constexpr const char* DIMENSION_ERROR = "Dimension mismatch in the shapes of vectors/matrices.";
constexpr const char* DIVISION_ERROR  = "Division by zero";

/******************Overloaded constructors*******************/
//Default constructors
template<typename T>
Vector<T>::Vector() 
    : vec(), dim(0) {
        // std::cout << "Default initialization" << std::endl;
    }

template<typename T>
Vector<T>::Vector(size_t size)
    : vec(size, 0), dim(size) {
        // std::cout << "Initialize with size" << std::endl;
    }

template<typename T>
Vector<T>::Vector(size_t size, T init_val)
    : vec(size, init_val), dim(size) {
        // std::cout << "Initialize with size and list" << std::endl;
    }
    
template<typename T>
Vector<T>::Vector(std::initializer_list<T> init_val) 
    : vec(init_val), dim(init_val.size()) {
        // std::cout << "Initialize with list" << std::endl;
    }

template<typename T>
Vector<T>::Vector(std::vector<T> init_vec) 
    : vec(init_vec), dim(init_vec.size()) {
        // std::cout << "Initialize with vector" << std::endl;
    }

//Copy constructor
template<typename T>
Vector<T>::Vector(const Vector<T> &source)
    : vec(source.vec), dim(source.dim) {
        // std::cout << "Copy constructor - deep" << std::endl;
    }

//Move constructor
template<typename T>
Vector<T>::Vector(Vector<T> &&source) noexcept
    : vec(std::move(source.vec)), dim(source.dim)  {
        // std::cout << "Move constructor" << std::endl;
    }

//Destructor
template<typename T>
Vector<T>::~Vector() {
    // std::cout << "Destroyed" << std::endl;
}

/******************Overloaded operators*******************/
//Copy assignment
template<typename T>
Vector<T> &Vector<T>::operator=(const Vector<T> &rhs) {
    //std::cout << "Copy assignment" << std::endl;
    if (this==&rhs) 
        return *this;
    this->vec = rhs.vec;
    return *this;
}

//Move assignement
template<typename T>
Vector<T> &Vector<T>::operator=(Vector<T> &&rhs) {
    // std::cout << "Move assignment" << std::endl;
    if(this==&rhs)
        return *this;
    vec = std::move(rhs.vec); //Copy reference value and Transfer of ownership
    dim = rhs.dim;
    return *this;
}

template<typename T>
Vector<T> &Vector<T>::operator=(T scalar) {
    // std::cout << "Move assignment" << std::endl;
    Vector<T> temp(vec.size(), scalar);
    vec = std::move(temp.vec); //Copy reference value and Transfer of ownership
    return *this;
}

//Equality
template<typename T>
bool Vector<T>::operator==(const Vector<T> &rhs) const {
    return (vec == rhs.vec);
}

template<typename T>
bool Vector<T>::operator!=(const Vector<T> &rhs) const {
    return (vec != rhs.vec);
}

//Accessor operators
template<typename T>
T &Vector<T>::operator[](size_t index) {
    return vec[index];
}

template<typename T>
const T &Vector<T>::operator[](size_t index) const {
    return vec[index];
}

//Overloaded insertion operator
template<typename U>
std::ostream &operator<<(std::ostream &os, const Vector<U> &rhs) {
    os << "[ ";
    for (const U& v : rhs.vec) {
        os << v << " ";
    }
    os << "]";
    return os;
}

//Overloaded extraction operator
template<typename U>
std::istream &operator>>(std::istream &is, Vector<U> &rhs) {
    std::vector<U> buff;
    U input;
    // std::cout << "Enter value(s) of a vector (press Enter to finish): "
    std::string line;
    if (std::getline(is, line)) {
        std::istringstream iss(line);
        while (iss >> input) {
            buff.push_back(input);
        }
    }
    
    rhs = Vector<U>{buff};
    return is;
}

//Arithmetic operators for a vector space
//Vector addition
template<typename T>
Vector<T> Vector<T>::operator-() const {
    Vector<T> temp(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) 
        temp[i] = -vec[i];
    return temp;
}

template<typename T>
Vector<T> Vector<T>::operator+(const Vector<T> &rhs) const {
#ifdef DEBUG
    if (this->size() != rhs.size()) {
        throw std::invalid_argument(DIMENSION_ERROR);
    } 
#endif
    Vector<T> temp(vec.size());
    for (size_t i = 0; i < this->vec.size(); ++i)
        temp[i] = vec[i] + rhs[i];
    return temp;
}

template<typename T>
Vector<T> Vector<T>::operator-(const Vector<T> &rhs) const {
#ifdef DEBUG
    if (this->size() != rhs.size()) {
        throw std::invalid_argument(DIMENSION_ERROR);
    }
#endif
    Vector<T> temp(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        temp[i] = vec[i] - rhs[i];
    return temp;
}

//Scalar-Vector addition (Special)
template<typename T>
Vector<T> Vector<T>::operator+(T scalar) const {
    Vector<T> temp(vec.size());
    for (size_t i = 0; i < this->vec.size(); ++i)
        temp[i] = vec[i] + scalar;
    return temp;
}

template<typename T>
Vector<T> Vector<T>::operator-(T scalar) const {
    Vector<T> temp(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        temp[i] = vec[i] - scalar;
    return temp;
}

//Commutativity for Scalar-Vector addition
template<typename U>
Vector<U> operator+(U scalar, const Vector<U> &rhs) {
    Vector<U> temp(rhs.size());
    for (size_t i = 0; i < rhs.size(); ++i)
        temp[i] = scalar + rhs[i];
    return temp;
}

template<typename U>
Vector<U> operator-(U scalar, const Vector<U> &rhs) {
    Vector<U> temp(rhs.size());
    for (size_t i = 0; i < rhs.size(); ++i)
        temp[i] = scalar - rhs[i];
    return temp;
}

//Scalar muliplication
template<typename T>
Vector<T> Vector<T>::operator*(T scalar) const {
    Vector<T> temp(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        temp[i] = vec[i]*scalar;
    return temp;
}

template<typename T>
Vector<T> Vector<T>::operator/(T scalar) const {
#ifdef DEBUG
    if (scalar == 0) {
        throw std::invalid_argument(DIVISION_ERROR);
    }
#endif    
    Vector<T> temp(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        temp[i] = vec[i]/scalar;
    return temp;
}

//Commutativity for Scalar multiplication
template<typename U>
Vector<U> operator*(U scalar, const Vector<U> &rhs) {
    Vector<U> temp(rhs.size());
    for (size_t i = 0; i < rhs.size(); ++i)
        temp[i] = scalar*rhs[i];
    return temp;
}

//Dot product of 1D
template<typename T>
T Vector<T>::dot(const Vector<T> &rhs) const {
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
double dot(const Vector<U> &lhs, const Vector<U> &rhs) {
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
size_t Vector<T>::size() {
    return dim;
}

template<typename T>
const size_t Vector<T>::size() const {
    return dim;
}

template<typename T>
void Vector<T>::print() const {
    std::cout << "[ ";
    for (const T& v : this->vec) {
        std::cout << v << " ";
    }
    std::cout << "]";
    std::cout << "\n";
}

/********************Accessment on data***********************/
//Access member methods
template<typename T>
T& Vector<T>::at(size_t index) {
    return vec.at(index);
}

template<typename T>
const T& Vector<T>::at(size_t index) const{
    return vec.at(index);
}

template<typename T>
const std::vector<T>& Vector<T>::get_vec() const {
    return vec;
}

template<typename T>
typename std::vector<T>::iterator Vector<T>::begin() {
    return vec.begin();
}

template<typename T>
typename std::vector<T>::iterator Vector<T>::end() {
    return vec.end();
}

template<typename T>
typename std::vector<T>::const_iterator Vector<T>::begin() const {
    return vec.begin();
}

template<typename T>
typename std::vector<T>::const_iterator Vector<T>::end() const {
    return vec.end();
}

//Dynamical extension
template<typename T>
void Vector<T>::push_back(const T& value) {
    vec.push_back(value);
    dim = vec.size();
}

template<typename T>
void Vector<T>::pop_back() {
    if (vec.empty()) {
        throw std::underflow_error("Vector is empty.");
    }
    vec.pop_back();
    dim = vec.size();
}

template<typename T>
void Vector<T>::resize(size_t new_size, T init_val) {
    vec.resize(new_size, init_val);
    dim = vec.size();
}

template<typename T>
void Vector<T>::reserve(size_t new_size) {
    vec.reserve(new_size);
    dim = vec.size();
}

template<typename T>
void Vector<T>::clear() {
    vec.clear();
    dim = vec.size();
}

template<typename T>
typename std::vector<T>::iterator Vector<T>::insert(size_t pos, const T& value) {
    if (pos > vec.size()) {
        throw std::out_of_range("Insert position out of bounds");
    }
    auto it = vec.insert(vec.begin() + pos, value);
    dim = vec.size(); 
    return it;
}

template<typename T>
typename std::vector<T>::iterator Vector<T>::erase(size_t pos) {
    if (pos >= vec.size()) {
        throw std::out_of_range("Erase position out of bounds");
    }
    auto it = vec.erase(vec.begin() + pos);
    dim = vec.size();
    return it;
}

/******************Mathematical implementation************************/
//Summation
template<typename T>
T Vector<T>::sum() const {
    T temp{0};
    for (size_t i = 0; i < vec.size(); ++i)
        temp += vec[i];
    return temp;
}

//Magnitude
template<typename T>
T Vector<T>::mag() const {
    T temp{0};
    for (size_t i = 0; i < vec.size(); ++i)
        temp += std::pow(vec[i], 2);
    return std::sqrt(temp);
}

template<typename T>
Vector<T> Vector<T>::unit() {
#ifdef DEBUG
    if (this.mag() == 0) 
        throw std::invalid_argument(DIVISION_ERROR);
#endif
    Vector<T> temp{*this};
    T length = this->mag();
    return temp/length;
}

#endif