#ifndef MATRIX_HPP
#define MATRIX_HPP

/******************Overloaded constructors*******************/
//Default constructors
template<typename T> 
Matrix<T>::Matrix() 
    : mat(), rows(0), cols(0) {
        //std::cout << "Default constructors" << std::endl;
    }

template<typename T> 
Matrix<T>::Matrix(size_t nrow) 
    : mat(nrow*1), rows(nrow), cols(1) {
        //std::cout << "Initialize with size" << std::endl;
    } 

template<typename T> 
Matrix<T>::Matrix(size_t nrow, size_t ncol) 
    : mat(nrow*ncol), rows(nrow), cols(ncol) {
        //std::cout << "Initialize with size" << std::endl;
    } 

template<typename T> 
Matrix<T>::Matrix(size_t nrow, size_t ncol, T init_val) 
    : mat(nrow*ncol, init_val), rows(nrow), cols(ncol) {
        //std::cout << "Initialize with size and list" << std::endl;
    }

template<typename T> 
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> init_val) {
    rows = init_val.size();
    if (rows == 0) {
        cols = 0;
        return;
    }
    cols = init_val.begin()->size();
    for (auto &row : init_val) {
        if (row.size() != cols) {
            throw std::invalid_argument("Ragged initializer_list not allowed.");
        }
    }
    this->mat = Vector<T>(rows * cols);
    size_t idx = 0;
    for (auto &row_list : init_val) { 
        for (auto &val : row_list) {
            this->mat[idx++] = val;
        }
    }
    // std::cout << "Initialize with initializer_list" << std::endl;
}
    
template<typename T> 
Matrix<T>::Matrix(std::vector<std::vector<T>> init_vec) {
    rows == init_vec.size();
    if (rows == 0) {
        cols = 0;
        return;
    }
    cols = init_vec[0].size();
    for (auto &row : init_vec) {
        if (row.size() != cols) {
            throw std::invalid_argument("Ragged initializer_list not allowed.");
        }
    }
    this->mat = Vector<T>(rows * cols);
    size_t idx = 0;
    for (auto &row_list : init_vec) { 
        for (auto &val : row_list) {
            this->mat[idx++] = val;
        }
    }
    // std::cout << "Initialize with vector" << std::endl;
}

//Copy constructor
template<typename T> 
Matrix<T>::Matrix(const Matrix<T> &source) 
    : mat(source.mat), rows(source.rows), cols(source.cols) {
        // std::cout << "Copy constructor - deep" << std::endl;
    }

//Move constructor
template<typename T> 
Matrix<T>::Matrix(Matrix<T> &&source) noexcept
    : mat(std::move(source.mat)), rows(source.rows), cols(source.cols) {}
    
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
    this->mat  = rhs.mat;
    this->rows = rhs.rows;
    this->cols = rhs.cols;
    return *this;
}
    
//Move assignment
template<typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T> &&rhs) {
    // std::cout << "Move assignment" << std::endl;
    if(this==&rhs)
        return *this;
    this->mat  = std::move(rhs.mat); //Copy reference value and Transfer of ownership
    this->rows = rhs.rows;
    this->cols = rhs.cols;
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
T& Matrix<T>::operator()(size_t row, size_t col) {
#ifdef DEBUG
    if (r >= rows || c >= cols)
        throw std::out_of_range("Matrix index out of range");
#endif
    return mat[row * cols + col];
}

template<typename T>
const T& Matrix<T>::operator()(size_t row, size_t col) const {
#ifdef DEBUG
    if (r >= rows || c >= cols)
        throw std::out_of_range("Matrix index out of range");
#endif
    return mat[row * cols + col];
}

template<typename T>
T* Matrix<T>::operator[](size_t row) {
    return &mat[row * cols];  
}

template<typename T>
const T* Matrix<T>::operator[](size_t row) const {
    return &mat[row * cols];
}

//Overloaded insertion operator
template<typename U>
std::ostream &operator<<(std::ostream &os, const Matrix<U> &rhs) {
    os << "[\n";
    for (size_t i = 0; i < rhs.nrow(); ++i) {
        os << " [ ";  
        for (size_t j = 0; j < rhs.ncol(); ++j) {
            os << rhs(i, j) << " ";
        }
        os << "]\n";
    }
    os << "]";
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

//Arithmetic operators for matrices
//Matrix addition
template<typename T>
Matrix<T> Matrix<T>::operator-() const {
    Matrix<T> temp(this->rows, this->cols);
    temp.mat = -mat;
    return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &rhs) const {
// #ifdef DEBUG
    if (this->rows != rhs.rows || this->cols != rhs.cols)
        throw std::invalid_argument(DIMENSION_ERROR);

// #endif
    Matrix<T> temp(this->rows, this->cols);
    temp.mat = this->mat + rhs.mat;
    return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &rhs) const {
// #ifdef DEBUG
    if (this->rows != rhs.rows || this->cols != rhs.cols)
        throw std::invalid_argument(DIMENSION_ERROR);

// #endif
    Matrix<T> temp(this->rows, this->cols);
    temp.mat = this->mat - rhs.mat;
    return temp;
}

//Scalar-Matrix addition (Special)
template<typename T>
Matrix<T> Matrix<T>::operator+(T scalar) const {
    Matrix<T> temp(this->rows, this->cols);
    temp.mat = this->mat + scalar;
    return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(T scalar) const {
    Matrix<T> temp(this->rows, this->cols);
    temp.mat = this->mat - scalar;
    return temp;
}

//Commutativity for Scalar-Matrix addition
template<typename U>
Matrix<U> operator+(U scalar, const Matrix<U> &rhs) {
    Matrix<U> temp(rhs.rows, rhs.cols);
    temp.mat = scalar + rhs.mat;
    return temp;
}

template<typename U>
Matrix<U> operator-(U scalar, const Matrix<U> &rhs) {
    Matrix<U> temp(rhs.rows, rhs.cols);
    temp.mat = scalar - rhs.mat;
    return temp;
}

//Scalar muliplication
template<typename T>
Matrix<T> Matrix<T>::operator*(T scalar) const {
    Matrix<T> temp(this->rows, this->cols);
    temp.mat = this->mat*scalar;
    return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator/(T scalar) const {
#ifdef DEBUG
    if (scalar == 0) {
        throw std::invalid_argument(DIVISION_ERROR);
    }
#endif    
    Matrix<T> temp(this->rows, this->cols);
    temp.mat = mat/scalar;
    return temp;
}

//Commutativity for Scalar multiplication
template<typename U>
Matrix<U> operator*(U scalar, const Matrix<U> &rhs) {
    Matrix<U> temp(rhs.rows, rhs.cols);
    temp.mat = scalar*rhs.mat;
    return temp;
}

//Vector-Matrix multiplication
template<typename T>
Matrix<T> Matrix<T>::operator*(const Vector<T> &vector) const {
// #ifdef DEBUG 
    if (this->ncol() != vector.size())
        throw std::invalid_argument("Matrix multiplication error: M.ncol() = " +
                                    std::to_string(this->ncol()) + ", V.size() = " +
                                    std::to_string(vector.size()) +
                                    ". Ensure M.ncol() == V.size()");
// #endif
    Matrix<T> temp(this->nrow());
    for (size_t i = 0; i < this->nrow(); ++i) {
        for (size_t j = 0; j < vector.size(); ++j) {
            temp[i] += this[i][j]*vector[j];
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

//Matrix-Matrix multiplication
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
                temp[i][j] += mat[i*cols + k]*rhs[k][j];
            }
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

//Display and Accessment options
//Display
template<typename T>
size_t Matrix<T>::nrow() {
    return this->rows;
}

template<typename T>
const size_t Matrix<T>::nrow() const {
    return this->rows;
}

template<typename T>
size_t Matrix<T>::ncol() {
    return this->cols;
}

template<typename T>
const size_t Matrix<T>::ncol() const {
    return this->cols;
}

template<typename T>
void Matrix<T>::print() const {
    std::cout << "[\n";
    for (size_t i = 0; i < this->rows; ++i) {
        std::cout << " [ ";  
        for (size_t j = 0; j < this->cols; ++j) {
            std::cout << mat[i*cols + j] << " ";
        }
        std::cout << "]\n";
    }
    std::cout << "]\n";
}

//Accessment
template<typename T>
T& Matrix<T>::at(size_t row, size_t col) {
    return mat.at(row*this->cols + col);
}

template<typename T>
const T& Matrix<T>::at(size_t row, size_t col) const {
    return mat.at(row*this->cols + col);
}

template<typename T>
const Matrix<T> &Matrix<T>::get_mat() const {
    return mat;
}

template<typename T>
typename std::vector<T>::iterator Matrix<T>::begin() {
    return mat.begin();
}

template<typename T>
typename std::vector<T>::iterator Matrix<T>::end() {
    return mat.end();
}

template<typename T>
typename std::vector<T>::const_iterator Matrix<T>::begin() const {
    return mat.begin();
}

template<typename T>
typename std::vector<T>::const_iterator Matrix<T>::end() const {
    return mat.end();
}

template<typename T>
typename std::vector<T>::iterator Matrix<T>::row_begin(size_t row) {
    if (row >= rows)
        throw std::out_of_range("row_begin: row index out of range");
    return mat.begin() + row*cols;
}

template<typename T>
typename std::vector<T>::iterator Matrix<T>::row_end(size_t row) {
    if (row >= rows)
        throw std::out_of_range("row_end: row index out of range");
    return mat.begin() + row*cols + cols;
}

template<typename T>
typename std::vector<T>::const_iterator Matrix<T>::row_begin(size_t row) const {
    if (row >= rows)
        throw std::out_of_range("row_begin: row index out of range");
    return mat.begin() + row*cols;
}

template<typename T>
typename std::vector<T>::const_iterator Matrix<T>::row_end(size_t row) const {
    if (row >= rows)
        throw std::out_of_range("row_end: row index out of range");
    return mat.begin() + row*cols + rows;
}

//Mathematical implementation
template<typename T>
Matrix<T> Matrix<T>::transpose() const{
    Matrix<T> temp(this->ncol(), this->nrow());
    for (size_t i = 0; i < this->nrow(); ++i) {
        for (size_t j = 0; j < this->ncol(); ++j) {
            temp[j][i] = mat[i*cols + j];
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
        temp[i] = mat[i*cols + i];
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
        temp[i] = mat[i*cols + i + upper];
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
        temp[i] =mat[(i + lower)*cols + i];
    return temp;
}

// template<typename T>
// Vector<T> Matrix<T>::low_diag(size_t lower) const{

// }

#endif