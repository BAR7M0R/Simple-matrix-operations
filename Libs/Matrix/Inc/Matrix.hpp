#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cstddef>
#include <memory>
#include <vector>
class Matrix {
public:
    Matrix() noexcept;
    Matrix(std::size_t n,std::size_t m);
    Matrix(std::size_t n, std::size_t m, double val);
    explicit Matrix(std::size_t n);
    Matrix(std::size_t n, double val);
    explicit Matrix(double val);
    
    template<std::size_t N, std::size_t M>
    explicit Matrix(std::array<std::array<double, M>, N> arr);

    template<std::size_t N>
    explicit Matrix(std::array<double, N> arr);
    
    Matrix(const Matrix&) = delete;
    Matrix& operator=(const Matrix&) = delete;
    Matrix(Matrix&&) noexcept = default;
    Matrix& operator=(Matrix&&) noexcept = default;

    double& operator()(std::size_t r, std::size_t c);
    const double& operator()(std::size_t r, std::size_t c) const;
    double* operator[](std::size_t i);
    const double* operator[](std::size_t i) const;

    [[nodiscard]] std::size_t getRows() const noexcept;
    [[nodiscard]] std::size_t getCols() const noexcept;
    static Matrix Identity();
    static Matrix Identity(std::size_t n);
    static Matrix Identity(std::size_t n, std::size_t m);
    static Matrix Zeros();
    static Matrix Zeros(std::size_t n);

    Matrix operator+(const Matrix& other) const;
    Matrix& operator+=(const Matrix& other);
    Matrix operator-(const Matrix& other) const;
    Matrix& operator-=(const Matrix& other);
    Matrix operator*(double scalar) const;
    Matrix& operator*=(double scalar);
    Matrix operator/(double scalar) const;
    Matrix& operator/=(double scalar);
    Matrix operator*(const Matrix& other) const;
    Matrix operator*(const std::vector<double>& vec) const;

    void transpose();
    Matrix transposed() const;

    void swapRows(std::size_t r1, std::size_t r2);
    void swapCols(std::size_t c1, std::size_t c2);

    [[nodiscard]] bool isSquare() const;

private:
    std::unique_ptr<double[]> matrix_;
    std::size_t n_;
    std::size_t m_;
};

#endif /* MATRIX_HPP */

/*
class Matrix {
public:

    Matrix()
    :arr_()
    ,rows_(0)
    ,columns_(0)
    ,isSquare_(false)
    {}

    Matrix(std::vector<std::vector<T>>& Arr)
    :arr_(Arr)
    ,isSquare_(false)
    ,rows_(Arr.size())
    ,columns_(Arr[0].size())
    {
        if (rows_ == columns_)
        {
            
            isSquare_ = true;
        }
    }
    Matrix(std::vector<T>& Arr)
    :arr_()
    ,isSquare_(false)
    ,rows_(0)
    ,columns_(0)
    {
        arr_(Arr);
        rows_ = arr_.size();
        columns_ = arr_[0].size();
        if (rows_ == columns_)
        {
            isSquare_ = true;
        }
    }
    std::vector<std::vector<T>> getRawData()
    {
        return arr_;
    }
    
    void printRowData()
    {
        static size_t CTR= 0;
        for(size_t i = 0; i<rows_; i++)
        {
            for(size_t j=0; j<columns_; j++)
            {
                std::cout << arr_.at(i).at(j) << " ";
            }
            std::cout << "\n";
        }
        std::cout << "raw arr was printed " << "CTR:" << CTR++ <<"\n";
        return;
    }
    Matrix operator+(const Matrix& secondMatrix)
    {
        if ((rows_!= secondMatrix.getRows()) 
            ||
            (columns_ != secondMatrix.getColumns()))
        {
            throw std::invalid_argument("dim is not equal");;
        }
        auto second = secondMatrix.getArray();
        auto outArr = arr_;
        for (size_t i =0; i<rows_; ++i)
        {
            for(size_t j =0; j<columns_; ++j)
            {
                outArr[i][j] = arr_[i][j] + second[i][j];
            }
        }
        return outArr;
    }
    Matrix operator*(const T& scalar)
    {
        auto ret = this->getArray();
        for (size_t i =0; i<rows_; ++i)
        {
            for(size_t j =0; j<columns_; ++j)
            {
                ret[i][j] = ret[i][j]*scalar;
            }
        }
        return ret;

    }
    Matrix operator*(const Matrix& secondMatrix)
    {
        std::vector<std::vector<T>> firstArr = arr_;
        std::vector<std::vector<T>> secondArr = secondMatrix.getArray();
        std::vector<std::vector<T>> newMatrix;

        size_t rowsNew = firstArr.size();
        size_t colisNew = secondArr[0].size();

        if (firstArr.size() != secondArr.size())
        {
            std::cout << "size(A): " << firstArr.size() << " | size(B): " << secondArr[0].size();
            throw std::invalid_argument("matrix has wrong dimention");;
        }
        size_t commonDim = firstArr.size();

       
        std::vector<T> subVector;
        T subElement = 0;
        for (size_t i =0; i<rowsNew; ++i)
        {
            for(size_t j =0; j<colisNew; ++j)
            {
                subElement = 0;
                for(size_t k = 0; k < commonDim; ++k)
                {
                    subElement = subElement + firstArr[i][k]*secondArr[k][j];
                }
                subVector.insert(subVector.end(), subElement);
            }
            newMatrix.insert(newMatrix.end(), subVector);
            subVector.clear();
        }
        return newMatrix;
    }
    bool operator==(const Matrix& second)
    {
        for (size_t i =0; i<rows_; ++i)
        {
            for(size_t j = 0; j<columns_; ++j)
            {
                if( this->getElement(i+1,j+1) != second.getElement(i+1,j+1))
                {
                    return false;
                }
            }
        }
        return true;
    }
    bool getIsSquare()
    {
        return isSquare_;
    }
    std::vector<std::vector<T>> getArray() const
    {
        return arr_;
    }
    size_t getRows() const
    {
        return rows_;
    }
    size_t getColumns() const
    {
        return columns_;
    }
    Matrix& transpose()
    {
        std::vector<std::vector<T>> arr;
        std::vector<T> subVector;
        for(size_t i = 0; i<columns_; ++i)
        {
            for(size_t j = 0; j<rows_; ++j)
            {
                subVector.insert(subVector.end(), arr_[j][i]);
            }
            arr.insert(arr.end(),subVector);
            subVector.clear();
        }
        arr_ = arr;
        columns_ = arr[0].size();
        rows_ = arr.size();
        return *this;
    }
    Matrix& setElement(const size_t& row,const size_t& column, const T& value)
    {
        if((row < 1) || (row> rows_))
        {
            throw std::invalid_argument("wrong row index");;
        }
        if((column < 1) || (column> columns_))
        {
            throw std::invalid_argument("wrong column index");;
        }
        arr_[row-1][column-1] = value;
        return *this;
    }
    T getElement(const size_t& row,const size_t& column) const
    {
        if((row < 1) || (row> rows_))
        {
            throw std::invalid_argument("wrong row index");;
        }
        if((column < 1) || (column> columns_))
        {
            throw std::invalid_argument("wrong column index");;
        }

        return arr_[row-1][column-1];
    }
    Matrix& operationToElement(T (Matrix<T>::*op)(const T&, const T&) const, const T& value, const size_t& row,const size_t& column)
    {
        if((row < 1) || (row> rows_))
        {
            throw std::invalid_argument("wrong row index");;
        }
        if((column < 1) || (column> columns_))
        {
            throw std::invalid_argument("wrong column index");;
        }
        arr_[row-1][column-1] = (this->*op)(arr_[row-1][column-1], value);
        return *this;
    }
    Matrix& operationToRow(T (Matrix<T>::*op)(const T&, const T&) const, const T& value , const size_t& row)
    {
        if((row < 1) || (row > rows_))
        {
            throw std::invalid_argument("wrong row index");;
        }
        for(size_t i = 0; i<columns_; i++)
        {
            arr_[row-1][i] = (this->*op)(arr_[row-1][i],value);
        }
        return *this;
    }
    Matrix& operationToColumn(T (Matrix<T>::*op)(const T&, const T&) const, const T& value, const size_t& column)
    {
        if((column < 1) || (column> columns_))
        {
            throw std::invalid_argument("wrong column index");;
        }
        for(size_t i = 0; i<columns_; i++)
        {
            arr_[i][column-1] = (this->*op)(arr_[i][column-1],value);
        }
        return *this;
    }
    Matrix getRow(const size_t& row)
    {   
        if((row < 1) || (row> rows_))
        {
            throw std::invalid_argument("wrong row index");;
        }
        matrix vecMatr(arr_[row-1]);
        return vecMatr;
    }

    Matrix getColumn(const size_t column)
    {
        if((column < 1) || (column>= columns_))
        {
            throw std::invalid_argument("wrong column index");;
        }
        std::vector<T> temp;
        for(size_t i = 0; i<rows_; ++i)
        {
            temp.push_back(arr_[i][column-1]);
        }
        matrix vecMatr(temp);
        return vecMatr;
    }
    Matrix& zeros(const size_t& size)
    {
        arr_.clear();
        rows_ = size;
        columns_ = size;
        std::vector<std::vector<T>> subVector(size, std::vector<T>(size));
        arr_ =subVector;
        return *this;
    }
    Matrix& identity(const size_t& size)
    {
        this->zeros(size);
        for(size_t i = 0;i<size;++i)
        {
            arr_[i][i] = 1;
        }
        return *this;
    }

    Matrix triDagonalDecomp()
    {
        if (arr_[0][0] != 0 || arr_[rows_-1][columns_-1] != 0 || columns_ != 3)
        {
            throw std::invalid_argument("it's not triediagonal matrix");
        }
        size_t newRows = rows_;
        size_t newColumns = rows_;
        std::vector<std::vector<T>> temp(newRows, std::vector<T>(newColumns));
        size_t position = 0;
        size_t elementCounter = 0;
        for(size_t i=0; i<rows_;++i)
        { 
            for(size_t j=0; j<columns_;++j)
            {
                 
                if((elementCounter) > 0 && (elementCounter < rows_*columns_-1))
                {
                    
                    temp[i][position]=arr_[i][j];
                    ++position;
                }
                ++elementCounter;
                
            }
            position -= 2;
        }
        Matrix ret(temp);
        return ret;
    }
    Matrix& clear()
    {
        for(size_t i = 0; i < rows_; ++i)
        {
            for(size_t j=0; j< columns_; ++j)
            {
                arr_.at(i).at(j) = 0;
            }
        }
        return *this;
    }
    T minValue()
    {
        T minVal = arr_[0][0];
        for(size_t i = 0; i<rows_;++i)
        {
            for(size_t j = 0; j<columns_;++j)
            {
                if(minVal > arr_[i][j])
                {
                    minVal = arr_[i][j];
                }
            }
        }
    }
    T maxValue()
    {
        T maxVal = arr_[0][0];
        for(size_t i = 0; i<rows_;++i)
        {
            for(size_t j = 0; j<columns_;++j)
            {
                if(maxVal < arr_[i][j])
                {
                    maxVal = arr_[i][j];
                }
            }
        }
    }
    std::vector<T> vectorization() 
        {
        std::vector<T> retVec;
        retVec.reserve(rows_*columns_);
        for(size_t i = 0; i<columns_;++i)
        {
            for(size_t j = 0; j<rows_;++j)
            {
                retVec.push_back(this->arr_.at(i).at(j));
            }
        }
        return retVec;
    }
    T add(const T& addendA, const T& addendB) const
    {
        return(addendA + addendB);
    }
    T subtract(const T& minuend, const T& subtrahend) const
    {
        return(minuend - subtrahend);
    }
    T multiply(const T& multiplierA, const T& multiplierB) const
    {
        return(multiplierA * multiplierB);
    }
    T devide(const T& numerator,const T&  denominator) const
    {
        if(denominator == 0)
        {
            throw std::invalid_argument("denominator is 0");
        }
        return(numerator / denominator);
    }
private:
    bool isSquare_;
    size_t rows_;
    size_t columns_;
    std::vector<std::vector<T>> arr_; 
};*/
