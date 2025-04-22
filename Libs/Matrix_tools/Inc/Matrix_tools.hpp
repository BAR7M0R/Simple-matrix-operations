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

#ifndef MATRIX_TOOLS_HPP
#define MATRIX_TOOLS_HPP

#include "Matrix.hpp"
#include <math.h>
#include <tuple>

namespace Matrix_tools
{
    double det(const Matrix& m);
    std::tuple<Matrix, Matrix> LU(const Matrix& m);
    double trace(const Matrix& m);
    double frobeniusNorm(const Matrix& m);
}

#endif //MATRIX_TOOLS_HPP
