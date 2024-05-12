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

Matrix cofactor_matrix(Matrix a, int row, int col) // 计算代数余子式
{
    Matrix m = create_matrix(a.cols - 1, a.rows - 1);
    for (int i = 0, ti = 0; i < a.rows - 1; i++)
    {
        if (i == row)
            continue;
        for (int j = 0, tj = 0; j < a.cols - 1; j++)
        {
            if (j == col)
                continue;
            m.data[ti][tj++] = a.data[i][j];
        }
        ti++;
    }
    return m;
}

double det_matrix(Matrix a)
{
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    else if (a.rows == 1)
        return a.data[0][0];
    else if (a.rows == 2)
        return a.data[0][0] * a.data[1][1] - a.data[1][0] * a.data[0][1];
    else
    {
        int sum = 0, n = 1;
        for (int j = 0; j < a.cols; j++)
        {
            int n = (j % 2 == 0) ? 1 : -1;
            sum += n * a.data[0][j] * det_matrix(cofactor_matrix(a, 0, j));
        }
        return sum;
    }
}

Matrix inv_matrix(Matrix a)
{
    if (det_matrix(a) == 0)
    {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }
    else if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    else
    {
        Matrix m = create_matrix(a.rows, a.cols);
        double det = det_matrix(a);
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < a.cols; j++)
            {
                double n = ((i + j) % 2 == 0) ? 1.0 : -1.0;
                m = scale_matrix(cofactor_matrix(a, j, i), n);
            }
        }
        return scale_matrix(m, 1.0 / det);
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