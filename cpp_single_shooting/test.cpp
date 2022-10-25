#include "lseq_gaussian_solver.cpp"

int main()
{
    std::vector<std::vector<double>> mat
    {
        {1, 2},
        {4, 5}
    };
    std::vector<double> result = gaussianElimination(mat);
    int x = 5;
}