#ifndef HUNGARIAN_HPP
#define HUNGARIAN_HPP

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <list>
#include <string>
#include <type_traits>
#include <vector>

#include <iostream>

    

void step1(std::vector<std::vector<double>>& matrix, int& step);


void clear_covers(std::vector<int>& cover);


void step2(const std::vector<std::vector<double>>& matrix, 
           std::vector<std::vector<int>>& M, 
           std::vector<int>& RowCover,
           std::vector<int>& ColCover, 
           int& step);


void step3(const std::vector<std::vector<int>>& M, std::vector<int>& ColCover, int& step);


void find_a_zero(int& row, 
                 int& col,
                 const std::vector<std::vector<double>>& matrix,
                 const std::vector<int>& RowCover,
                 const std::vector<int>& ColCover);


bool star_in_row(int row, const std::vector<std::vector<int>>& M);


void find_star_in_row(int row, int& col, const std::vector<std::vector<int>>& M);


void step4(const std::vector<std::vector<double>>& matrix, 
           std::vector<std::vector<int>>& M, 
           std::vector<int>& RowCover,
           std::vector<int>& ColCover,
           int& path_row_0,
           int& path_col_0,
           int& step);


void find_star_in_col(int c, int& r, const std::vector<std::vector<int>>& M);


void find_prime_in_row(int r, int& c, const std::vector<std::vector<int>>& M);


void augment_path(std::vector<std::vector<int>>& path, int path_count, std::vector<std::vector<int>>& M);


void erase_primes(std::vector<std::vector<int>>& M);


void step5(std::vector<std::vector<int>>& path, 
           int path_row_0, 
           int path_col_0, 
           std::vector<std::vector<int>>& M, 
           std::vector<int>& RowCover,
           std::vector<int>& ColCover,
           int& step);


void find_smallest(double& minval, 
                   const std::vector<std::vector<double>>& matrix, 
                   const std::vector<int>& RowCover,
                   const std::vector<int>& ColCover);


void step6(std::vector<std::vector<double>>& matrix, 
           const std::vector<int>& RowCover,
           const std::vector<int>& ColCover,
           int& step);


std::vector<int> hungarian(std::vector<std::vector<double>>& original);



#endif
