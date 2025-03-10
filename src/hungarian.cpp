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
 * This version is written by Fernando B. Giannasi
 *
 * Cloned repository from https://github.com/phoemur/hungarian_algorithm */

#include <hungarian.hpp>
    

/* For each row of the matrix, find the smallest element and subtract it from every 
 * element in its row.  
 * For each col of the matrix, find the smallest element and subtract it from every 
 * element in its col. Go to Step 2. */
void step1(std::vector<std::vector<double>>& matrix, int& step){
    // process rows
    for (auto& row: matrix) {
        auto smallest = *std::min_element(begin(row), end(row));
        if (smallest > 0)        
            for (auto& n: row)
                n -= smallest;
    }
    
    // process cols
    int sz = matrix.size(); // square matrix is granted
    for (int j=0; j<sz; ++j) {
        double minval = std::numeric_limits<double>::max();
        for (int i=0; i<sz; ++i) {
            minval = std::min(minval, matrix[i][j]);
        }
        
        if (minval > 0) {
            for (int i=0; i<sz; ++i) {
                matrix[i][j] -= minval;
            }
        }
    }
   
    step = 2;

    return;
}

/* helper to clear the temporary vectors */
void clear_covers(std::vector<int>& cover) 
{
    for (auto& n: cover) n = 0;
    return;
}

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
void step2(const std::vector<std::vector<double>>& matrix, 
           std::vector<std::vector<int>>& M, 
           std::vector<int>& RowCover,
           std::vector<int>& ColCover, 
           int& step)
{
    int sz = matrix.size();
    
    for (int r=0; r<sz; ++r) 
        for (int c=0; c<sz; ++c) 
            if (matrix[r][c] == 0)
                if (RowCover[r] == 0 && ColCover[c] == 0) {
                    M[r][c] = 1;
                    RowCover[r] = 1;
                    ColCover[c] = 1;
                }
            
    clear_covers(RowCover); // reset vectors for posterior using
    clear_covers(ColCover);
    
    step = 3;

    return;
}


/* Cover each column containing a starred zero.  If K columns are covered, the starred 
 * zeros describe a complete set of unique assignments.  In this case, Go to DONE, 
 * otherwise, Go to Step 4. Once we have searched the entire cost matrix, we count the 
 * number of independent zeros found.  If we have found (and starred) K independent zeros 
 * then we are done.  If not we procede to Step 4.*/
void step3(const std::vector<std::vector<int>>& M, std::vector<int>& ColCover, int& step){
    int sz = M.size();
    int colcount = 0;
    
    for (int r=0; r<sz; ++r)
        for (int c=0; c<sz; ++c)
            if (M[r][c] == 1)
                ColCover[c] = 1;
            
    for (auto& n: ColCover)
        if (n == 1)
            colcount++;
    
    if (colcount >= sz) {
        step = 7; // solution found
    }
    else {
        step = 4;
    }

    return;
}

// Following functions to support step 4
void find_a_zero(int& row, 
                 int& col,
                 const std::vector<std::vector<double>>& matrix,
                 const std::vector<int>& RowCover,
                 const std::vector<int>& ColCover)
{
    int r = 0;
    int c = 0;
    int sz = matrix.size();
    bool done = false;
    row = -1;
    col = -1;
    
    while (!done) {
        c = 0;
        while (true) {
            if (matrix[r][c] == 0 && RowCover[r] == 0 && ColCover[c] == 0) {
                row = r;
                col = c;
                done = true;
            }
            c += 1;
            if (c >= sz || done)
                break;
        }
        r += 1;
        if (r >= sz)
            done = true;
    }

    return;
}

bool star_in_row(int row, const std::vector<std::vector<int>>& M){
    bool tmp = false;
    for (unsigned c = 0; c < M.size(); c++)
        if (M[row][c] == 1)
            tmp = true;
    
    return tmp;
}


void find_star_in_row(int row, int& col, const std::vector<std::vector<int>>& M){
    col = -1;
    for (unsigned c = 0; c < M.size(); c++)
        if (M[row][c] == 1)
            col = c;
    return;
}


/* Find a noncovered zero and prime it.  If there is no starred zero in the row containing
 * this primed zero, Go to Step 5.  Otherwise, cover this row and uncover the column 
 * containing the starred zero. Continue in this manner until there are no uncovered zeros
 * left. Save the smallest uncovered value and Go to Step 6. */
void step4(const std::vector<std::vector<double>>& matrix, 
           std::vector<std::vector<int>>& M, 
           std::vector<int>& RowCover,
           std::vector<int>& ColCover,
           int& path_row_0,
           int& path_col_0,
           int& step)
{
    int row = -1;
    int col = -1;
    bool done = false;

    while (!done){
        find_a_zero(row, col, matrix, RowCover, ColCover);
        
        if (row == -1){
            done = true;
            step = 6;
        }
        else {
            M[row][col] = 2;
            if (star_in_row(row, M)) {
                find_star_in_row(row, col, M);
                RowCover[row] = 1;
                ColCover[col] = 0;
            }
            else {
                done = true;
                step = 5;
                path_row_0 = row;
                path_col_0 = col;
            }
        }
    }

    return;
}

// Following functions to support step 5
void find_star_in_col(int c, int& r, const std::vector<std::vector<int>>& M){
    r = -1;
    for (unsigned i = 0; i < M.size(); i++)
        if (M[i][c] == 1)
            r = i;
    return;
}

void find_prime_in_row(int r, int& c, const std::vector<std::vector<int>>& M){
    for (unsigned j = 0; j < M.size(); j++)
        if (M[r][j] == 2)
            c = j;
    return;
}

void augment_path(std::vector<std::vector<int>>& path, int path_count, std::vector<std::vector<int>>& M){
    for (int p = 0; p < path_count; p++)
        if (M[path[p][0]][path[p][1]] == 1)
            M[path[p][0]][path[p][1]] = 0;
        else
            M[path[p][0]][path[p][1]] = 1;
    return;
}

void erase_primes(std::vector<std::vector<int>>& M){
    for (auto& row: M)
        for (auto& val: row)
            if (val == 2)
                val = 0;
    return;
}


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
           int& step)
{
    int r = -1;
    int c = -1;
    int path_count = 1;
    
    path[path_count - 1][0] = path_row_0;
    path[path_count - 1][1] = path_col_0;
    
    bool done = false;
    while (!done) {
        find_star_in_col(path[path_count - 1][1], r, M);
        if (r > -1) {
            path_count += 1;
            path[path_count - 1][0] = r;
            path[path_count - 1][1] = path[path_count - 2][1];
        }
        else {done = true;}
        
        if (!done) {
            find_prime_in_row(path[path_count - 1][0], c, M);
            path_count += 1;
            path[path_count - 1][0] = path[path_count - 2][0];
            path[path_count - 1][1] = c;
        }
    }
    
    augment_path(path, path_count, M);
    clear_covers(RowCover);
    clear_covers(ColCover);
    erase_primes(M);
    
    step = 3;

    return;
}

// methods to support step 6
void find_smallest(double& minval, 
                   const std::vector<std::vector<double>>& matrix, 
                   const std::vector<int>& RowCover,
                   const std::vector<int>& ColCover)
{
    for (unsigned r = 0; r < matrix.size(); r++)
        for (unsigned c = 0; c < matrix.size(); c++)
            if (RowCover[r] == 0 && ColCover[c] == 0)
                if (minval > matrix[r][c])
                    minval = matrix[r][c];
    return;
}

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
void step6(std::vector<std::vector<double>>& matrix, 
           const std::vector<int>& RowCover,
           const std::vector<int>& ColCover,
           int& step)
{
    double minval = std::numeric_limits<double>::max();
    find_smallest(minval, matrix, RowCover, ColCover);
    
    int sz = matrix.size();
    for (int r = 0; r < sz; r++)
        for (int c = 0; c < sz; c++) {
            if (RowCover[r] == 1)
                matrix[r][c] += minval;
            if (ColCover[c] == 0)
                matrix[r][c] -= minval;
    }
    
    step = 4;

    return;
}


/* Main function of the algorithm */
std::vector<int> hungarian(std::vector<std::vector<double>>& original){  
    /* Initialize data structures */
    
    // Work on a vector copy to preserve original matrix
    // Didn't passed by value cause needed to access both
    std::vector<std::vector<double>> matrix(original.size(), std::vector<double>(original.begin()->size()));
    
    std::vector<std::vector<double>>::iterator it = original.begin();
    for (auto& vec: matrix){
        std::copy(it->begin(), it->end(), vec.begin());
        it = std::next(it);
    }
    
    std::size_t sz = matrix.size();
    
    /* The masked matrix M.  If M(i,j)=1 then C(i,j) is a starred zero,  
     * If M(i,j)=2 then C(i,j) is a primed zero. */
    std::vector<std::vector<int>> M (sz, std::vector<int>(sz, 0));
    
    /* We also define two vectors RowCover and ColCover that are used to "cover" 
     *the rows and columns of the cost matrix C*/
    std::vector<int> RowCover (sz, 0);
    std::vector<int> ColCover (sz, 0);
    
    int path_row_0, path_col_0; //temporary to hold the smallest uncovered value
    
    // Array for the augmenting path algorithm
    std::vector<std::vector<int>> path (sz+1, std::vector<int>(2, 0));
    
    /* Now Work The Steps */
    bool done = false;
    int step = 1;
    while (!done) {
        switch (step) {
            case 1:
                step1(matrix, step);
                break;
            case 2:
                step2(matrix, M, RowCover, ColCover, step);
                break;
            case 3:
                step3(M, ColCover, step);
                break;
            case 4:
                step4(matrix, M, RowCover, ColCover, path_row_0, path_col_0, step);
                break;
            case 5:
                step5(path, path_row_0, path_col_0, M, RowCover, ColCover, step);
                break;
            case 6:
                step6(matrix, RowCover, ColCover, step);
                break;
            case 7:
                for (auto& vec: M) {vec.resize(original.begin()->size());}
                M.resize(original.size());
                done = true;
                break;
            default:
                done = true;
                break;
        }
    }
    
    std::vector<int> assignment;
    for (int i = 0; i < M.size(); i++){
        for (int j = 0; j < M[i].size(); j++){
            if (M[i][j] == 1){
                assignment.push_back(j);
                break;
            }
        }
    }

    return assignment;
}
