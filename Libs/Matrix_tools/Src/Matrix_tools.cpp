/******************************************************************************
 * Copyright 2025, Bartłomiej Głodek
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/

#include "Matrix_tools.hpp"

std::tuple<Matrix, Matrix> Matrix_tools::LU(const Matrix &m) {
    if (!m.isSquare()) throw std::invalid_argument("matrix is not square");

    size_t n = m.getRows();
    Matrix L = Matrix::Identity(n);
    Matrix U = Matrix::Zeros(n);

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t k = i; k < n; ++k)
        {
            double sum = 0.0;
            for (size_t crossAxisIndex = 0ull; crossAxisIndex < i; ++crossAxisIndex)
                sum += L[i][crossAxisIndex] * U[crossAxisIndex][k];

            U[i][k] = m[i][k] - sum;
        }

        for (size_t k = i; k < n; ++k)
        {
            if (i != k)
            {
                double sum = 0.0;
                for (size_t crossAxisIndex = 0ull; crossAxisIndex < i; ++crossAxisIndex)
                    sum += L[k][crossAxisIndex] * U[crossAxisIndex][i];

                L[k][i] = (m[k][i] - sum) / U[i][i];
            }
        }
    }

    return std::make_tuple(std::move(L),std::move(U));
}

double Matrix_tools::trace(const Matrix &m) {
    if (!m.isSquare()) throw std::invalid_argument("matrix it is not square");

    double out = 0.0;

    for (std::size_t e = 0ull; e < m.getRows();++e)
    {
        out += m[e][e];
    }
    return out;
}

double Matrix_tools::frobeniusNorm(const Matrix &m) {
    if (!m.isSquare()) throw std::invalid_argument("matrix it is not square");
    double out = 0.0;
    for (std::size_t r = 0; r < m.getRows(); ++r)
    {
        for(std::size_t c = 0; c < m.getCols(); ++c)
        {
            out += pow(m[r][c],2);
        }
    }
    return sqrt(out);
}
