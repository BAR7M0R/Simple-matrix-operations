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
 ******************************************************************************//*

#include <iostream>
#include <vector>
#include <include/matrix.hpp>
#include <include/myStatistic.hpp>
//cm#include <matplot/matplot.hpp>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <random>

template <typename T>
linaryAlgebra::matrix<T> matrixRead(const std::string& path);
template <typename T>
linaryAlgebra::matrix<T> matrixRead(const std::string& path, const std::string keyword);
template <typename T>
std::vector<T> generateRandomScope(T minValue, T maxValue, size_t elementsNumber);
int main(int, char**)
{
    std::string outputFileName = "console_output.txt";
    std::ofstream outputFile(outputFileName);
    std::streambuf *coutbuf = std::cout.rdbuf();
    std::cout.rdbuf(outputFile.rdbuf());
    */
/**
     * @brief task 1: Calculate L and U matrices for the A matrix without using Matlab functions.
     *  Write the L and U matrices to the report.
     *  Tips: The formulas for U and L matrices can be found on the first slide in the file
     *  matrix_lecture.pdf. Please note the correct use of indexes, e.g. ij or ji.
     *  Start from the following matrices: L=eye(N); U=zeros(N,N);
     *  Finally, check if L*U provides the original matrix.
     *//*

    std::cout << "begin task: 1" <<std::endl;
    std::cout << "entry data:" << std::endl;
    std::string pathToTask12("data\\MO-Task12-Data5.txt");
    linaryAlgebra::matrix dataTask12 = matrixRead<double>(pathToTask12);
    dataTask12.printRowData();
    //calculate LU
    std::vector<linaryAlgebra::matrix<double>> LUmatrix(dataTask12.LU());
    std::cout << "L matrix:" << std::endl;
    LUmatrix.at(0).printRowData();
    std::cout << "U matrix:" << std::endl;
    LUmatrix.at(1).printRowData();
    linaryAlgebra::matrix aftercalc = LUmatrix.at(0)*LUmatrix.at(1);
    if (aftercalc==dataTask12)
    {
        std::cout << "LU work function work" << std::endl;
    }
    for (size_t i = 0; i < 50; i++)
    {
        std::cout << "-";
    }
    std::cout << std::endl;
    */
/**
     * @brief task 3: Rewrite two arrays from sparse form into classical matrices. Multiply them and write
     * the resulting matrix to the report.
     * Tip: You should write a code with loops rewriting the values from a sparse form to classical
     * form of matrix.
     *//*


    std::cout << "begin task: 3" <<std::endl;
    std::string pathToTask3("data\\MO-Task3-Data5.txt");
    linaryAlgebra::matrix dataTask3A = matrixRead<double>(pathToTask3, std::string("sa"));
    linaryAlgebra::matrix dataTask3B = matrixRead<double>(pathToTask3, std::string("sb"));
    std::cout << "entry data set sa:" << std::endl;
    dataTask3A.printRowData();
    std::cout << "Unpack the matrix" << std::endl;
    dataTask3A.triDagonalDecomp().printRowData();
    std::cout << "entry data set sb:" << std::endl;
    dataTask3A.printRowData();
    std::cout << "Unpack the matrix" << std::endl;
    dataTask3A.triDagonalDecomp().printRowData();
    for (size_t i = 0; i < 50; i++)
    {
        std::cout << "-";
    }
    std::cout << std::endl;
    */
/**
     * @brief task 4: Randomly generate (use Matlab rand() function) a set of 100 data described by 5 features.
     * Ranges of values for individual features are provided in the text file. Normalize this data:
     * - with zero mean value and unit variance,
     * - to the range [-1, 1].
     * In case of normalization to the range [-1, 1], check its correctness by:
     * - finding minimum and maximum values,
     * - making charts (first 2 features).
     * Include unnormalized data and both types of normalized data in your report.
     * Data normalization (description based on [1])
     * Normalization is performed to "equalize" the impact of individual features in the feature
     * vector that differ significantly in size or range of values. In such a situation, features with
     * larger values have a greater impact in the various criteria used in classification algorithms
     * than features with relatively small values, although this does not necessarily reflect their
     * "importance". The data should then be normalized - linear or non-linear scaling of all data to
     * the appropriate range. Normalization is often used, which results in features having zero mean
     * value and unit variance, or normalization to a fixed range (both types of normalization will be
     * discussed in class).

     *//*

    std::cout << "begin task: 4" <<std::endl;
    std::cout << "generated random vector:" << std::endl;
    std::vector<double> randomVector = generateRandomScope<double>(-100., 100., 100);
    for (auto ele : randomVector)
    {
        std::cout << ele << " | ";
    }
    std::cout << std::endl;
    std::cout << "convert random vector with zero Mean Value and Unit Variance" << std::endl;
    std::vector<double>  zeroMeanValueUnitVariance = myStatistic::zeroMeanNormalization(randomVector);
    for (auto ele : zeroMeanValueUnitVariance)
    {
        std::cout << ele << " | ";
    }
    std::cout << std::endl;
    std::cout << "Normalize the previously generated vector to a specific range [-1;1]" << std::endl;
    std::vector<double>  normalizeVector = myStatistic::minMaxNormalization(zeroMeanValueUnitVariance);
    for (auto ele : normalizeVector)
    {
        std::cout << ele << " | ";
    }
    std::cout << std::endl;
    std::cout << "min value: "
              << *std::min_element(normalizeVector.begin(), normalizeVector.end())
              << " max value: "
              << *std::max_element(normalizeVector.begin(), normalizeVector.end())
              << std::endl;
    std::cout.rdbuf(coutbuf);
    return 0;
}

template <typename T>
linaryAlgebra::matrix<T> matrixRead(const std::string& path)
{


    std::fstream inputFile(path);
    std::vector<std::vector<T>> formatedData;
    std::string line;

    if (!inputFile.is_open()) {
        throw std::runtime_error("filed to read selected file");
    }
    while (std::getline(inputFile, line))
    {
        if (line.empty())
        {
            continue;
        }
        std::istringstream iss(line);
        std::vector<double> row;
        double value;
        while (iss >> value) {
            row.push_back(value);
        }
        formatedData.push_back(row);
    }
    return linaryAlgebra::matrix<T>(formatedData);
}
template <typename T>
linaryAlgebra::matrix<T> matrixRead(const std::string& path, const std::string keyword)
{
    std::fstream inputFile(path);
    std::vector<std::vector<T>> formatedData;
    std::string line;
    if (!inputFile.is_open())
    {
        throw std::runtime_error("filed to read selected file");
    }
    bool hasKeyword = false;
    while (std::getline(inputFile, line))
    {
        if (line.find(keyword) != std::string::npos)
        {
            hasKeyword = true;
            continue;
        }
        if (line.empty())
        {
            hasKeyword = false;
            continue;
        }
        if (hasKeyword)
        {
            std::istringstream iss(line);
            std::vector<double> row;
            double value;
            while (iss >> value)
            {
                row.push_back(value);
            }
            formatedData.push_back(row);
        }
    }
    return linaryAlgebra::matrix<T>(formatedData);
}
template <typename T>
std::vector<T> generateRandomScope(T minValue, T maxValue, size_t elementsNumber)
{
    std::random_device randomDevice;
    std::mt19937 randomGenerator(randomDevice());
    std::uniform_real_distribution<T> distribution(minValue, maxValue);
    std::vector<T> randomElements;
    for (int i = 0; i < elementsNumber; ++i)
    {
        randomElements.push_back(distribution(randomGenerator));
    }
    return randomElements;
}*/
