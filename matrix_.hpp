#ifndef MATRIX_CA_HPP
#define MATRIX_CA_HPP

/******************Overloaded constructors*******************/
//Default constructors
template<typename T> 
Matrix<T>::Matrix() 
    : mat() {
        //std::cout << "Default constructors" << std::endl;
    }

template<typename T> 
Matrix<T>::Matrix(size_t nrow) 
    : mat(nrow, Vector<T>(1)) {
        //std::cout << "Initialize with size" << std::endl;
    } 

template<typename T> 
Matrix<T>::Matrix(size_t nrow, size_t ncol) 
    : mat(nrow, Vector<T>(ncol)) {
        //std::cout << "Initialize with size" << std::endl;
    } 

template<typename T> 
Matrix<T>::Matrix(size_t nrow, size_t ncol, T init_val) 
    : mat(nrow, Vector<T>(ncol, init_val)) {
        //std::cout << "Initialize with size and list" << std::endl;
    }

template<typename T> 
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> init_val) 
    : mat(init_val.size()) {
        size_t i{0};
        for (const auto& row : init_val) 
            mat[i++] = row;
        // std::cout << "Initialize with list" << std::endl;
    }
    
template<typename T> 
Matrix<T>::Matrix(std::vector<std::vector<T>> init_vec) 
    : mat(init_vec.size()) {
        size_t i{0};
        for (const auto& row : init_vec)
            mat[i++] = std::move(Vector<T>(row)) ;
        // std::cout << "Initialize with vector" << std::endl;
    }

//Copy constructor
template<typename T> 
Matrix<T>::Matrix(const Matrix<T> &source) 
    : mat(source.mat) {
        // std::cout << "Copy constructor - deep" << std::endl;
    }

//Move constructor
template<typename T> 
Matrix<T>::Matrix(Matrix<T> &&source) noexcept
    : mat(std::move(source.mat)) {}
    
//Destructor
template<typename T> 
Matrix<T>::~Matrix() {};
    
/******************Overloaded operators*******************/
//Copy assignment
template<typename T> 
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &rhs) {
    //std::cout << "Copy assignment" << std::endl;
    if (this==&rhs) 
        return *this;
    this->mat = rhs.mat;
    return *this;
}
    
//Move assignment
template<typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T> &&rhs) {
    // std::cout << "Move assignment" << std::endl;
    if(this==&rhs)
        return *this;
    mat = std::move(rhs.mat); //Copy reference value and Transfer of ownership
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(T scalar) {
    // std::cout << "Move assignment" << std::endl;
    Matrix<T> temp(this->nrow(), this->ncol(), scalar);
    mat = std::move(temp.mat); //Copy reference value and Transfer of ownership
    return *this;
}

template<typename T>
bool Matrix<T>::operator==(const Matrix<T> &rhs) const {
    return (mat == rhs.mat);
}

template<typename T>
bool Matrix<T>::operator!=(const Matrix<T> &rhs) const {
    return (mat != rhs.mat);
}

//Accessor operators
template<typename T>
Vector<T>& Matrix<T>::operator[](size_t index) {
    return mat[index];
}

template<typename T>
const Vector<T>& Matrix<T>::operator[](size_t index) const {
    return mat[index];
}

//Overloaded insertion operator
template<typename U>
std::ostream &operator<<(std::ostream &os, const Matrix<U> &rhs) {
    std::cout << "[";
    std::cout << std::endl;
    for (size_t i = 0; i < rhs.mat.size(); ++i)
            std::cout << " " << rhs.mat[i] << std::endl;
    std::cout << "]";
    return os;
}
    
template<typename U>
std::istream &operator>>(std::istream &is, Matrix<U> &rhs) {
    std::vector<std::vector<U>> buff;
    U input;
    // std::cout << "Enter value(s) of row(s) (press Enter to input new row): ";
    std::string line;
    while (std::getline(is, line)) {  // 여러 줄 입력 처리
        if (line.empty()) {
            break;  // 빈 줄이 입력되면 종료
        }
        std::istringstream iss(line);
        std::vector<U> row;  // 각 행에 해당하는 Vector<T>
        U value;
        while (iss >> value) {
            row.push_back(value);  // 각 값을 행에 추가
        }
        buff.push_back(row);  // 행을 mat에 추가
    }
    rhs = Matrix<U>{buff};
    return is;
}

//Arithmetic operators for a vector space
//Matrix addition
template<typename T>
Matrix<T> Matrix<T>::operator-() const {
    Matrix<T> temp(mat.size());
    for (size_t i = 0; i < mat.size(); ++i) 
        temp[i] = std::move(-mat[i]);
    return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &rhs) const {
// #ifdef DEBUG
    if (this->nrow() != rhs.nrow() || this->ncol() != rhs.ncol())
        throw std::invalid_argument(DIMENSION_ERROR);

// #endif
    Matrix<T> temp(this->nrow(), this->ncol());
    for (size_t i = 0; i < mat.size(); ++i) 
        temp[i] = mat[i] + rhs[i];
    return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &rhs) const {
// #ifdef DEBUG
    if (this->nrow() != rhs.nrow() || this->ncol() != rhs.ncol()) 
        throw std::invalid_argument(DIMENSION_ERROR);
// #endif
    Matrix<T> temp(this->nrow(), this->ncol());
    for (size_t i = 0; i < mat.size(); ++i)
        temp[i] = mat[i] - rhs[i];
    return temp;
}

//Scalar-Matrix addition (Special)
template<typename T>
Matrix<T> Matrix<T>::operator+(T scalar) const {
    Matrix<T> temp(this->nrow(), this->ncol());
    for (size_t i = 0; i < mat.size(); ++i)
        temp[i] = mat[i] + scalar;
    return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(T scalar) const {
    Matrix<T> temp(this->nrow(), this->ncol());
    for (size_t i = 0; i < mat.size(); ++i)
        temp[i] = mat[i] - scalar;
    return temp;
}

//Commutativity for Scalar-Matrix addition
template<typename U>
Matrix<U> operator+(U scalar, const Matrix<U> &rhs) {
    Matrix<U> temp(rhs.nrow(), rhs.ncol());
    for (size_t i = 0; i < rhs.nrow(); ++i)
        temp[i] = scalar + rhs.mat[i];
    return temp;
}

template<typename U>
Matrix<U> operator-(U scalar, const Matrix<U> &rhs) {
    Matrix<U> temp(rhs.nrow(), rhs.ncol());
    for (size_t i = 0; i < rhs.nrow(); ++i)
        temp[i] = scalar - rhs.mat[i];
    return temp;
}

//Scalar muliplication
template<typename T>
Matrix<T> Matrix<T>::operator*(T scalar) const {
    Matrix<T> temp(this->nrow());
    for (size_t i = 0; i < this->nrow(); ++i)
        temp[i] = mat[i]*scalar;
    return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator/(T scalar) const {
#ifdef DEBUG
    if (scalar == 0) {
        throw std::invalid_argument(DIVISION_ERROR);
    }
#endif    
    Matrix<T> temp(this->nrow());
    for (size_t i = 0; i < this->nrow(); ++i)
        temp[i] = mat[i]/scalar;
    return temp;
}

//Commutativity for Scalar multiplication
template<typename U>
Matrix<U> operator*(U scalar, const Matrix<U> &rhs) {
    Matrix<U> temp(rhs.nrow());
    for (size_t i = 0; i < rhs.nrow(); ++i)
        temp[i] = scalar*rhs[i];
    return temp;
}

//Product
template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &rhs) const {
// #ifdef DEBUG
    if (this->ncol() != rhs.nrow())
        throw std::invalid_argument("Matrix multiplication error: A.ncol() = " +
                                    std::to_string(this->ncol()) + ", B.nrow() = " +
                                    std::to_string(rhs.nrow()) +
                                    ". Ensure A.cols == B.rows.");
// #endif
    Matrix<T> temp(this->nrow(), rhs.ncol());
    for (size_t i = 0; i < this->nrow(); ++i) {
        for (size_t j = 0; j < rhs.ncol(); ++j) {
            for (size_t k = 0; k < this->ncol(); ++k) {
                temp[i][j] += mat[i][k]*rhs.mat[k][j];
            }
        }
    }   
    return temp; 
}

template<typename U>
Vector<U> operator*(const Vector<U> &lhs, const Matrix<U> &rhs) {
// #ifdef DEBUG 
    if (lhs.size() != rhs.nrow())
        throw std::invalid_argument("Matrix multiplication error: V.size() = " +
                                    std::to_string(lhs.size()) + ", M.nrow() = " +
                                    std::to_string(rhs.nrow()) +
                                    ". Ensure A.cols == B.rows.");
// #endif
    Matrix<U> temp(rhs.ncol());
    for (size_t i = 0; i < rhs.ncol(); ++i) {
        for (size_t j = 0; j < lhs.size(); ++j) {
            temp[i] += lhs[j]*rhs[j][i];
        }
    }   
    return temp;     
}

template<typename U>
Matrix<U> matmul(const Matrix<U> &lhs, const Matrix<U> &rhs) {
// #ifdef DEBUG 
    if (lhs.ncol() != rhs.nrow())
        throw std::invalid_argument("Matrix multiplication error: A.ncol() = " +
                                    std::to_string(lhs.ncol()) + ", B.nrow() = " +
                                    std::to_string(rhs.nrow()) +
                                    ". Ensure A.cols == B.rows.");
// #endif
    Matrix<U> temp(lhs.nrow(), rhs.ncol());
    for (size_t i = 0; i < lhs.nrow(); ++i) {
        for (size_t j = 0; j < rhs.ncol(); ++j) {
            for (size_t k = 0; k < lhs.ncol(); ++k) {
                temp[i][j] += lhs.mat[i][k] * rhs.mat[k][j];
            }
        }
    } 
    return temp; 
}

template<typename U>
Matrix<U> matmul(const Vector<U> &lhs, const Matrix<U> &rhs) {
// #ifdef DEBUG 
    if (lhs.size() != rhs.nrow()) 
        throw std::invalid_argument("Matrix multiplication error: V.size() = " +
                                    std::to_string(lhs.size()) + ", M.nrow() = " +
                                    std::to_string(rhs.nrow()) +
                                    ". Ensure A.cols == B.rows.");
// #endif
    Matrix<U> temp(rhs.ncol());
    for (size_t i = 0; i < rhs.ncol(); ++i) {
        for (size_t j = 0; j < lhs.size(); ++j) {
            temp[i] += lhs[j]*rhs[j][i];
        }
    }   
    return temp; 
}

//Display and Accessment options
//Display
template<typename T>
size_t Matrix<T>::nrow() {
    return mat.size();
}

template<typename T>
const size_t Matrix<T>::nrow() const {
    return mat.size();
}

template<typename T>
size_t Matrix<T>::ncol() {
    return mat[0].size();
}

template<typename T>
const size_t Matrix<T>::ncol() const {
    return mat[0].size();
}

template<typename T>
void Matrix<T>::print() const {
    std::cout << "[";
    std::cout << std::endl;
    for (size_t i = 0; i < mat.size(); ++i)
            std::cout << " " << mat[i] << std::endl;
    std::cout << "]";
    std::cout << "\n";
}

//Accessment
template<typename T>
T& Matrix<T>::at(size_t row, size_t col) {
    return mat.at(row).at(col);
}

template<typename T>
const T& Matrix<T>::at(size_t row, size_t col) const {
    return mat.at(row).at(col);
}

template<typename T>
const Matrix<T> &Matrix<T>::get_mat() const {
    return mat;
}

template<typename T>
typename std::vector<T>::iterator Matrix<T>::row_begin(size_t row) {
    if (row >= mat.size()) 
        throw std::out_of_range("Row index out of range");
    return mat[row].begin();
}

template<typename T>
typename std::vector<T>::iterator Matrix<T>::row_end(size_t row) {
    if (row >= mat.size()) 
        throw std::out_of_range("Row index out of range");
    return mat[row].end();
}

template<typename T>
typename std::vector<T>::const_iterator Matrix<T>::row_begin(size_t row) const {
    if (row >= mat.size()) 
        throw std::out_of_range("Row index out of range");
    return mat[row].begin();
}

template<typename T>
typename std::vector<T>::const_iterator Matrix<T>::row_end(size_t row) const {
    if (row >= mat.size()) 
        throw std::out_of_range("Row index out of range");
    return mat[row].end();
}

//Mathematical implementation
template<typename T>
Matrix<T> Matrix<T>::transpose() const{
    Matrix<T> temp(this->ncol(), this->nrow());
    for (size_t i = 0; i < this->nrow(); ++i) {
        for (size_t j = 0; j < this->ncol(); ++j) {
            temp[j][i] = mat[i][j];
        }
    }
    return temp;
}

template<typename T>
Vector<T> Matrix<T>::diag() const{
    size_t ndiag=
        (this->ncol() >= this->nrow()) ? 
        this->nrow() : this->ncol();
    Vector<T> temp(ndiag);
    for (size_t i = 0; i < ndiag; ++i)
        temp[i] = mat[i][i];
    return temp;
} 

template<typename T>
Vector<T> Matrix<T>::upp_diag(size_t upper) const{
    size_t ndiag=
        (this->ncol() >= this->nrow()) ? 
        this->nrow() : this->ncol();
    if (upper >= this->ncol())
        throw std::out_of_range("Upper index out of range");
    Vector<T> temp(ndiag);
    for (size_t i = 0; i < this->ncol()-upper; ++i) 
        temp[i] = mat[i][i+upper];
    return temp;
}

template<typename T>
Vector<T> Matrix<T>::low_diag(size_t lower) const{
    size_t ndiag=
        (this->ncol() >= this->nrow()) ? 
        this->nrow() : this->ncol();
    if (lower >= this->nrow())
        throw std::out_of_range("Lower index out of range");
    Vector<T> temp(ndiag);
    for (size_t i = 0; i < this->nrow()-lower; ++i) 
        temp[i] = mat[i+lower][i];
    return temp;
}

// template<typename T>
// Vector<T> Matrix<T>::low_diag(size_t lower) const{

// }

#endif