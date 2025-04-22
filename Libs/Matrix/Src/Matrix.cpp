#include "Matrix.hpp"

Matrix::Matrix() noexcept
        : n_(0ull), m_(0ull), matrix_(nullptr)
{}

Matrix::Matrix(const std::size_t n, const std::size_t m)
        : n_(n), m_(m), matrix_(std::make_unique<double[]>(n * m))
{}

Matrix::Matrix(const std::size_t n, const std::size_t m, const double val)
        : Matrix(n, m)
{
    for (std::size_t e = 0ull; e < getRows() * getCols(); ++e) {
        matrix_[e] = val;
    }
}

Matrix::Matrix(const std::size_t n)
        : Matrix(n, n) {}

Matrix::Matrix(const std::size_t n, const double val)
        : Matrix(n, n, val)
{}

Matrix::Matrix(const double val)
        : Matrix(1ull, 1ull, double(val))
{}
template<std::size_t N, std::size_t M>
Matrix::Matrix(std::array<std::array<double, M>, N> arr)
:Matrix(N, M)
{
    for (std::size_t r = 0ull; r < getRows(); ++r)
    {
        for(std::size_t c = 0ull; c < getCols(); ++c)
        {
            matrix_[r * getRows() + c] = arr[r][c];
        }
    }

}
template<std::size_t M>
Matrix::Matrix(std::array<double, M> arr)
: Matrix(1ull, M)
{
    for(std::size_t c = 0ull; c < getCols(); ++c)
    {
        matrix_[c] = arr[c];
    }
}

double &Matrix::operator()(const std::size_t r, const std::size_t c) {
    return matrix_[r * getRows() + c];
}

const double &Matrix::operator()(const std::size_t r, const std::size_t c) const {
    return matrix_[r * getRows() + c];
}


std::size_t Matrix::getRows() const noexcept {return n_;}

std::size_t Matrix::getCols() const noexcept {return m_;}

Matrix Matrix::Identity() {
    return {1ull, 1ull, 1.0};
}
Matrix Matrix::Identity(const std::size_t n) {
    return Identity(n, n);
}

Matrix Matrix::Identity(const std::size_t n, const std::size_t m) {
    Matrix out(n, m, 0.0);
    for (std::size_t i = 0; i < std::min(out.getCols() , out.getCols()); ++i)
        out(i, i) = 1.0;
    return out;
}

Matrix Matrix::Zeros() {
    return {1ull, 1ull, 0.0};
}

Matrix Matrix::Zeros(const std::size_t n) {
    return {n, n, 0.0};
}

void Matrix::transpose() {
    std::unique_ptr<double[]> newMatrix(new double[n_ * m_]);
    for (std::size_t r = 0ull; r < getRows(); ++r)
    {
        for(std::size_t c = 0ull; c < getCols(); ++c)
        {
            newMatrix[c * getCols() + r] = matrix_[r * getRows() + c];
        }
    }
    matrix_ = std::move(newMatrix);
    const std::size_t m = m_;
    m_ = n_;
    n_ = m;
}

const double *Matrix::operator[](std::size_t i) const {
    return &matrix_[i * m_];
}

double *Matrix::operator[](std::size_t i) {
    return &matrix_[i * m_];
}

Matrix Matrix::operator+(const Matrix &other) const
{
    if (other.getCols() != this->getCols()) throw std::invalid_argument("Matrix::operator+ -> wrong number of rows");
    if (other.getRows() != this->getRows()) throw std::invalid_argument("Matrix::operator+ -> wrong number of cols");
    Matrix out(n_, m_);
    for (std::size_t i = 0ull; i < m_ * n_; ++i)
    {
        for (std::size_t j = 0ull; j<n_;++j)
        {
        out[i][j] = this->operator[](i)[j] + other[i][j];

        }
    }
    return out;
}

Matrix &Matrix::operator+=(const Matrix &other)
{
    if (other.getCols() != this->getCols()) throw std::invalid_argument("Matrix::operator+= -> wrong number of rows");
    if (other.getRows() != this->getRows()) throw std::invalid_argument("Matrix::operator+= -> wrong number of cols");
    for (std::size_t i = 0ull; i < m_; ++i)
    {
        for (std::size_t j = 0ull; j < n_; ++j)
        {
            this->operator[](i)[j] += other[i][j];
        }

    }
    return *this;
}

Matrix Matrix::operator-(const Matrix &other) const {
    if (other.getCols() != this->getCols()) throw std::invalid_argument("Matrix::operator- -> wrong number of rows");
    if (other.getRows() != this->getRows()) throw std::invalid_argument("Matrix::operator- -> wrong number of cols");
    Matrix out(n_, m_);
    for (std::size_t i = 0ull; i < m_; ++i)
    {
        for (std::size_t j = 0ull; j<n_;++j)
        {
            out[i][j] = this->operator[](i)[j] - other[i][j];

        }
    }
    return out;
}

Matrix &Matrix::operator-=(const Matrix &other)
{
    if (other.getCols() != this->getCols()) throw std::invalid_argument("Matrix::operator-= -> wrong number of rows");
    if (other.getRows() != this->getRows()) throw std::invalid_argument("Matrix::operator-= -> wrong number of cols");
    for (std::size_t i = 0ull; i < m_; ++i)
    {
        for (std::size_t j = 0ull; j < n_; ++j)
        {
            this->operator[](i)[j] -= other[i][j];
        }

    }
    return *this;
}

Matrix Matrix::operator*(double scalar) const
{
    Matrix out(n_, m_);
    for (std::size_t i = 0ull; i < m_; ++i)
    {
        for (std::size_t j = 0ull; j<n_;++j)
        {
            out[i][j] = this->operator[](i)[j] * scalar;

        }
    }
    return out;
}

Matrix &Matrix::operator*=(double scalar) {
    for (std::size_t i = 0ull; i < m_; ++i)
    {
        for (std::size_t j = 0ull; j < n_; ++j)
        {
            this->operator[](i)[j] *= scalar;
        }

    }
    return *this;
}

Matrix Matrix::operator/(double scalar) const {
    if (scalar == 0) throw std::invalid_argument("Matrix::operator/ -> division by 0");
    Matrix out(n_, m_);
    for (std::size_t i = 0ull; i < m_; ++i)
    {
        for (std::size_t j = 0ull; j<n_;++j)
        {
            out[i][j] = this->operator[](i)[j] / scalar;

        }
    }
    return out;
}

Matrix &Matrix::operator/=(double scalar) {
    if (scalar == 0) throw std::invalid_argument("Matrix::operator/= -> division by 0");
    for (std::size_t i = 0ull; i < m_; ++i)
    {
        for (std::size_t j = 0ull; j < n_; ++j)
        {
            this->operator[](i)[j] /= scalar;
        }
    }
    return *this;
}

Matrix Matrix::operator*(const Matrix &other) const {
    if (other.getRows() != getCols()) throw std::invalid_argument("Matrix::operator* -> rows of A != cols of B");
    Matrix out(getRows(), other.getCols());

    for (std::size_t r = 0ull; r < out.getRows(); ++r)
    {
        for (std::size_t c = 0ull; c < out.getCols(); ++c)
        {
            for (std::size_t e = 0ull; e < other.getRows(); ++e)
            {
                out[r][c] += this->operator[](r)[e] * other[e][c];
            }
        }
    }
    return out;
}

Matrix Matrix::operator*(const std::vector<double> &other) const {
    if (1ull != getCols()) throw std::invalid_argument("Matrix::operator* -> rows of A != cols of B");
    Matrix out(getRows(), other.size());

    for (std::size_t r = 0ull; r < out.getRows(); ++r)
    {
        for (std::size_t c = 0ull; c < out.getCols(); ++c)
        {
            for (std::size_t e = 0ull; e < other.size(); ++e)
            {
                out[r][c] += this->operator[](r)[e] * other[c];
            }
        }
    }
    return out;
}

Matrix Matrix::transposed() const
{
    Matrix out(getCols(), getRows());

    for (std::size_t r = 0ull; r < getRows(); ++r)
    {
        for(std::size_t c = 0ull; c < getCols(); ++c)
        {
            out[c][r] = operator[](r)[c];
        }
    }
    return out;
}

void Matrix::swapRows(std::size_t r1, std::size_t r2)
{
    if (r1 >= getRows()) throw std::invalid_argument("r1 is out of range");
    if (r2 >= getRows()) throw std::invalid_argument("r2 is out of range");
    if (r1 == r2) return;
    for(std::size_t c = 0ull; c < getCols(); ++c)
    {
        std::swap(matrix_[r1 * getRows() + c], matrix_[r2 * getRows() + c]);
    }
}

void Matrix::swapCols(std::size_t c1, std::size_t c2)
{
    if (c1 >= getCols()) throw std::invalid_argument("c1 is out of range");
    if (c2 >= getCols()) throw std::invalid_argument("r2 is out of range");
    if (c1 == c2) return;
    for(std::size_t r = 0; r < getRows(); ++r)
    {
        std::swap(matrix_[r * getRows() + c1], matrix_[r * getRows() + c2]);
    }
}

bool Matrix::isSquare() const
{
    return getCols() == getRows();
}





