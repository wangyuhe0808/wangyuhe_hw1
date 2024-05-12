#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    if (a.rows != b.rows || a.cols != b.cols)
    {
        printf("Error: Matrix a and b must have the same rows and cols.");
        return create_matrix(0, 0);
    }
    else
    {
        Matrix m = create_matrix(a.rows, a.cols);
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < a.cols; j++)
                m.data[i][j] = a.data[i][j] + b.data[i][j];
        }
        return m;
    }
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    if (a.rows != b.rows || a.cols != b.cols)
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    else
    {
        Matrix m = create_matrix(a.rows, a.cols);
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < a.cols; j++)
                m.data[i][j] = a.data[i][j] - b.data[i][j];
        }
        return m;
    }
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    if (a.cols != b.rows)
    {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }
    Matrix m = create_matrix(a.rows, b.cols);
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < b.cols; j++)
        {
            double sum = 0.0;
            for (int k = 0; k < a.cols; k++)
            {
                sum += a.data[i][k] * b.data[k][j];
            }
            m.data[i][j] = sum;
        }
    }
    return m;
}

Matrix scale_matrix(Matrix a, double k)
{
    Matrix m = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
            m.data[i][j] = a.data[i][j] * k;
    }
    return m;
}

Matrix transpose_matrix(Matrix a)
{
    Matrix m = create_matrix(a.cols, a.rows);
    for (int i = 0; i < a.cols; i++)
    {
        for (int j = 0; j < a.cols; j++)
            m.data[i][j] = a.data[j][i];
    }
    return m;
}

double algebraic_cofactor(Matrix a, int c)
{
    int m = a.rows - 1;
    int q = 0;
    Matrix n = create_matrix(m, m);
    for (int i = 0; i < m; i++)
    {
        q = 0;
        for (int j = 0; j < a.rows; j++)
        {
            if (j == c)
            {
                continue;
            }
            n.data[i][q] = a.data[i + 1][j];
            q++;
        }
    }
    double result = det_matrix(n);
    return result;
}
int pow_(int a)
{
    int result = 1;
    for (int i = 1; i <= a; i++)
    {
        result *= -1;
    }
    return result;
}

double det_matrix(Matrix a)
{
    double sum = 0.0;
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    else
    {
        if (a.rows == 1)
        {
            return a.data[0][0];
        }
        else
        {
            for (int i = 0; i < a.rows; i++)
            {
                sum += pow_(i) * a.data[0][i] * algebraic_cofactor(a, i);
            }
        }
        return sum;
    }
}
Matrix adjoint_matrix(Matrix a)
{
    Matrix n = create_matrix(a.rows, a.rows);
    int b = 0, c = 0;
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.rows; j++)
        {
            Matrix m = create_matrix(a.rows - 1, a.rows - 1);
            b = 0;
            for (int p = 0; p < a.rows; p++)
            {
                if (p == j)
                {
                    continue;
                }
                else
                {
                    c = 0;
                    for (int q = 0; q < a.rows; q++)
                    {
                        if (q == i)
                        {
                            continue;
                        }
                        else
                        {
                            m.data[b][c] = a.data[p][q];
                        }
                        c++;
                    }
                    b++;
                }
            }
            n.data[i][j] = pow_(i + j) * det_matrix(m);
        }
    }
    return n;
}

Matrix inv_matrix(Matrix a)
{
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    else if (a.rows == 1 || det_matrix(a) == 0)
    {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }
    else
    {
        Matrix n = create_matrix(a.rows, a.rows);
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < a.rows; j++)
            {
                n.data[i][j] = adjoint_matrix(a).data[i][j] / det_matrix(a);
            }
        }
        return n;
    }
}
int rank_matrix(Matrix a)
{
    {
        int rank = 0;
        for (int col = 0, row = 0; col < a.cols && row < a.rows; ++col)
        {
            int max_row = row;
            for (int i = row + 1; i < a.rows; ++i)
            {
                if (fabs(a.data[i][col]) > fabs(a.data[max_row][col]))
                {
                    max_row = i;
                }
            }
            if (fabs(a.data[max_row][col]) < 1e-10)
                continue;
            if (max_row != row)
            {
                swap_rows(a, row, max_row);
            }
            for (int i = row + 1; i < a.rows; ++i)
            {
                double factor = a.data[i][col] / a.data[row][col];
                for (int j = col; j < a.cols; ++j)
                {
                    a.data[i][j] -= factor * a.data[row][j];
                }
            }
            ++rank;
            ++row;
        }
        return rank;
    }
}

Matrix swap_rows(Matrix a, int row1, int row2)
{
    for (int i = 0; i < a.cols; i++)
    {
        double temp = a.data[row1][i];
        a.data[row1][i] = a.data[row2][i];
        a.data[row2][i] = temp;
    }
    return a;
}

double trace_matrix(Matrix a)
{
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    else
    {
        int t = 0;
        for (int i = 0; i < a.rows; i++)
            t += a.data[i][i];
        return t;
    }
}
void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}