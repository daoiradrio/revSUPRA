/* This is an implementation of the Hungarian algorithm in C++
 * The Hungarian algorithm, also know as Munkres or Kuhn-Munkres
 * algorithm is usefull for solving the assignment problem.
 *
 * Assignment problem: Let C be an n x n matrix 
 * representing the costs of each of n workers to perform any of n jobs.
 * The assignment problem is to assign jobs to workers so as to 
 * minimize the total cost. Since each worker can perform only one job and 
 * each job can be assigned to only one worker the assignments constitute 
 * an independent set of the matrix C.
 * 
 * It is a port heavily based on http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
 * 
 * This version is written by Fernando B. Giannasi */

#ifndef HUNGARIAN_HPP
#define HUNGARIAN_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <string>
#include <type_traits>
#include <vector>


    
/* Utility function to print Matrix */
template<template <typename, typename...> class Container,
                   typename T,
                   typename... Args>
//disable for string, which is std::basic_string<char>, a container itself
typename std::enable_if<!std::is_convertible<Container<T, Args...>, std::string>::value &&
                        !std::is_constructible<Container<T, Args...>, std::string>::value,
                            std::ostream&>::type
operator<<(std::ostream& os, const Container<T, Args...>& con)
{
    os << " ";
    for (auto& elem: con)
        os << elem << " ";

    os << "\n";
    return os;
}

/* Handle negative elements if present. If allowed = true, add abs(minval) to 
 * every element to create one zero. Else throw an exception */
template<typename T>
void handle_negatives(std::vector<std::vector<T>>& matrix, 
                      bool allowed = true){}

/* Ensure that the matrix is square by the addition of dummy rows/columns if necessary */
template<typename T>
void pad_matrix(std::vector<std::vector<T>>& matrix){}

/* For each row of the matrix, find the smallest element and subtract it from every 
 * element in its row.  
 * For each col of the matrix, find the smallest element and subtract it from every 
 * element in its col. Go to Step 2. */
template<typename T>
void step1(std::vector<std::vector<T>>& matrix, 
           int& step){}

/* helper to clear the temporary vectors */
inline void clear_covers(std::vector<int>& cover){}

/* Find a zero (Z) in the resulting matrix.  If there is no starred zero in its row or 
 * column, star Z. Repeat for each element in the matrix. Go to Step 3.  In this step, 
 * we introduce the mask matrix M, which in the same dimensions as the cost matrix and 
 * is used to star and prime zeros of the cost matrix.  If M(i,j)=1 then C(i,j) is a 
 * starred zero,  If M(i,j)=2 then C(i,j) is a primed zero.  We also define two vectors 
 * RowCover and ColCover that are used to "cover" the rows and columns of the cost matrix.
 * In the nested loop (over indices i and j) we check to see if C(i,j) is a zero value 
 * and if its column or row is not already covered.  If not then we star this zero 
 * (i.e. set M(i,j)=1) and cover its row and column (i.e. set R_cov(i)=1 and C_cov(j)=1).
 * Before we go on to Step 3, we uncover all rows and columns so that we can use the 
 * cover vectors to help us count the number of starred zeros. */
template<typename T>
void step2(const std::vector<std::vector<T>>& matrix, 
           std::vector<std::vector<int>>& M, 
           std::vector<int>& RowCover,
           std::vector<int>& ColCover, 
           int& step){}

/* Cover each column containing a starred zero.  If K columns are covered, the starred 
 * zeros describe a complete set of unique assignments.  In this case, Go to DONE, 
 * otherwise, Go to Step 4. Once we have searched the entire cost matrix, we count the 
 * number of independent zeros found.  If we have found (and starred) K independent zeros 
 * then we are done.  If not we procede to Step 4.*/
void step3(const std::vector<std::vector<int>>& M, 
           std::vector<int>& ColCover,
           int& step){}

// Following functions to support step 4
template<typename T>
void find_a_zero(int& row, 
                 int& col,
                 const std::vector<std::vector<T>>& matrix,
                 const std::vector<int>& RowCover,
                 const std::vector<int>& ColCover){}

bool star_in_row(int row, 
                 const std::vector<std::vector<int>>& M){}

void find_star_in_row(int row,
                      int& col, 
                      const std::vector<std::vector<int>>& M){}

/* Find a noncovered zero and prime it.  If there is no starred zero in the row containing
 * this primed zero, Go to Step 5.  Otherwise, cover this row and uncover the column 
 * containing the starred zero. Continue in this manner until there are no uncovered zeros
 * left. Save the smallest uncovered value and Go to Step 6. */
template<typename T>
void step4(const std::vector<std::vector<T>>& matrix, 
           std::vector<std::vector<int>>& M, 
           std::vector<int>& RowCover,
           std::vector<int>& ColCover,
           int& path_row_0,
           int& path_col_0,
           int& step){}

// Following functions to support step 5
void find_star_in_col(int c, 
                      int& r,
                      const std::vector<std::vector<int>>& M){}

void find_prime_in_row(int r, 
                       int& c, 
                       const std::vector<std::vector<int>>& M){}

void augment_path(std::vector<std::vector<int>>& path, 
                  int path_count, 
                  std::vector<std::vector<int>>& M){}

void erase_primes(std::vector<std::vector<int>>& M){}

/* Construct a series of alternating primed and starred zeros as follows.  
 * Let Z0 represent the uncovered primed zero found in Step 4.  Let Z1 denote the 
 * starred zero in the column of Z0 (if any). Let Z2 denote the primed zero in the 
 * row of Z1 (there will always be one).  Continue until the series terminates at a 
 * primed zero that has no starred zero in its column.  Unstar each starred zero of 
 * the series, star each primed zero of the series, erase all primes and uncover every 
 * line in the matrix.  Return to Step 3.  You may notice that Step 5 seems vaguely 
 * familiar.  It is a verbal description of the augmenting path algorithm (for solving
 * the maximal matching problem). */
void step5(std::vector<std::vector<int>>& path, 
           int path_row_0, 
           int path_col_0, 
           std::vector<std::vector<int>>& M, 
           std::vector<int>& RowCover,
           std::vector<int>& ColCover,
           int& step){}

// methods to support step 6
template<typename T>
void find_smallest(T& minval, 
                   const std::vector<std::vector<T>>& matrix, 
                   const std::vector<int>& RowCover,
                   const std::vector<int>& ColCover){}

/* Add the value found in Step 4 to every element of each covered row, and subtract it 
 * from every element of each uncovered column.  Return to Step 4 without altering any
 * stars, primes, or covered lines. Notice that this step uses the smallest uncovered 
 * value in the cost matrix to modify the matrix.  Even though this step refers to the
 * value being found in Step 4 it is more convenient to wait until you reach Step 6 
 * before searching for this value.  It may seem that since the values in the cost 
 * matrix are being altered, we would lose sight of the original problem.  
 * However, we are only changing certain values that have already been tested and 
 * found not to be elements of the minimal assignment.  Also we are only changing the 
 * values by an amount equal to the smallest value in the cost matrix, so we will not
 * jump over the optimal (i.e. minimal assignment) with this change. */
template<typename T>
void step6(std::vector<std::vector<T>>& matrix, 
           const std::vector<int>& RowCover,
           const std::vector<int>& ColCover,
           int& step){}

/* Calculates the optimal cost from mask matrix */
template<template <typename, typename...> class Container,
         typename T,
         typename... Args>
T output_solution(const Container<Container<T,Args...>>& original,
                  const std::vector<std::vector<int>>& M){}

/* Main function of the algorithm */
template<template <typename, typename...> class Container,
         typename T,
         typename... Args>
typename std::enable_if<std::is_integral<T>::value, T>::type // Work only on integral types
hungarian(const Container<Container<T,Args...>>& original,
          bool allow_negatives = true){}

#endif
