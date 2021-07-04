//    MIT License
//
//    Copyright (c) 2021 Laxman Desai
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy
//            of this software and associated documentation files (the "Software"), to deal
//    in the Software without restriction, including without limitation the rights
//            to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//            copies of the Software, and to permit persons to whom the Software is
//    furnished to do so, subject to the following conditions:
//
//    The above copyright notice and this permission notice shall be included in all
//            copies or substantial portions of the Software.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//            AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//    SOFTWARE.

#ifndef _MATRIX_H_
#define _MATRIX_H_
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <iterator>
#include <cmath>
#include <climits>

#define _MN_ERR 1e-12

struct span {
    size_t l, r, len;
    span(size_t l, size_t r) : l(l), r(r), len(r-l+1) {
        assert(0<=l && l<=r);
    }
};

template <typename T>
class matrix {

public:
    // constructor with memory allocation
    matrix(size_t r, size_t c, T x=0) : R(r), C(c) {
        M.assign(R, std::vector<T>(C, x));
    }

    // default constructor
    matrix() : R(0), C(0) {}
    // constructor from initializer_list
    matrix(const std::initializer_list<std::vector<T>> m) : M(m) {
        R = M.size(), C = R > 0 ? M[0].size() : 0;
    }

    //  access size (rows, columns)
    bool empty() const {
        return R*C == 0;
    }
    //  access size (rows, columns)
    std::pair<size_t, size_t> size() const {
        return std::make_pair(R, C);
    }
    void resize(size_t r, size_t c, T x=0) {
        R = r, C = c;
        M.assign(R, std::vector<T>(C, x));
    }

    // Algebric Operations:
    // 1 matrix(A), 1 Scalar(b): A+I*b, A-I*b, A*b
    // 2 Matrices(A, B): A+B, A-B, A*B
    matrix operator+=(const matrix &rhs) {
        assert(R == rhs.R && C == rhs.C);
        for (size_t r = 0; r < R; r++) {
            for (size_t c = 0; c < C; c++) {
                M[r][c] += rhs(r, c);
            }
        }
        return *this;
    }
    matrix operator-=(const matrix &rhs) {
        assert(R == rhs.R && C == rhs.C);
        for (size_t r = 0; r < R; r++) {
            for (size_t c = 0; c < C; c++) {
                M[r][c] -= rhs(r, c);
            }
        }
        return *this;
    }
    matrix operator*=(const matrix &rhs) {
        assert(C == rhs.R);
        matrix res(R, rhs.C, 0);
        for (size_t r = 0; r < res.R; r++) {
            for (size_t c = 0; c<res.C; c++) {
                for (size_t i = 0; i < C; i++) {
                    res(r, c) += M[r][i] * rhs(i, c);
                }
            }
        }
        return res;
    }
    matrix operator+(const matrix &rhs) const {
        return matrix(*this) += rhs;
    }
    matrix operator-(const matrix &rhs) const {
        return matrix(*this) -= rhs;
    }
    matrix operator*(const matrix &rhs) const {
        return matrix(*this) *= rhs;
    }
    matrix operator*=(const T &rhs) {
        for (size_t r = 0; r < R; r++) {
            for (size_t c = 0; c < C; c++) {
                M[r][c] *= rhs;
            }
        }
        return *this;
    }
    matrix operator*(const T &rhs) const {
        return matrix(*this) *= rhs;
    }
    // To permit: 5*M, without this only M*5 is permitted
    friend matrix<T> operator*(const T &lhs, const matrix<T> &rhs) {
        return rhs * lhs;
    }
    matrix operator-() const {
        return matrix(R, C, 0) - *this;
    }

    // L2 norm of A-B (ie ||A-B||^2) must be lesser than _MN_ERR (=1e-12)
    bool operator==(const matrix &rhs) const {
        if (empty() || R != rhs.R || C != rhs.C) {
            return false;
        }
        double L2 = 0;
        for (size_t i = 0; i < R; i++) {
            for (size_t j = 0; j < C; j++) {
                L2 += (M[i][j] - rhs(i, j)) * (M[i][j] - rhs(i, j));
            }
        }
        L2 /= R*C;
        return L2 < _MN_ERR;
    }
    bool operator!=(const matrix &rhs) const {
        return !(*this == rhs);
    }

    // Element Numbering (0 based, Column  Major):
    //    C 0  1  2   3
    //  R
    //  0   0  3  6   9
    //  1   1  4  7  10
    //  2   2  5  8  11

    // access data operators
    T& operator[](size_t i) {
        return M[i % R][i / R];
    }
    T& operator()(size_t r, size_t c) {
        assert(0 <= r && r < R && 0 <= c && c < C);
        return M[r][c];
    }
    T operator()(size_t r, size_t c) const {
        assert(0 <= r && r < R && 0 <= c && c < C);
        return M[r][c];
    }
    matrix operator()(const span &row_range, const span &col_range) const {
        assert(row_range.r < R && col_range.r < C);
        matrix<T> res(row_range.len, col_range.len, 0);
        for (size_t r = 0; r < row_range.len; r++) {
            for (size_t c = 0; c < col_range.len; c++)
                res(r, c) = M[row_range.l + r][col_range.l + c];
        }
        return res;
    }

    // compute transpose
    matrix t() {
        matrix res(C, R);
        for (size_t r = 0; r < R; r++) {
            for (size_t c = 0; c < C; c++) {
                res(c, r) += M[r][c];
            }
        }
        return res;
    }

    // compute dth minor
    matrix compute_minor(int d) {
        matrix res(R, C);
        for (int i = 0; i < d; i++) {
            res(i, i) = 1;
        }
        for (int r = d; r < R; r++) {
            for (int c = d; c < C; c++) {
                res(r, c) = M[r][c];
            }
        }
        return res;
    }

    // return c-th column of m
    matrix column(int c) const {
        matrix<T> res(R, 1);
        for (size_t r = 0; r < R; r++) {
            res(r, 0) = M[r][c];
        }
        return res;
    }

    // return r-th row of m
    matrix row(int r) const {
        return {M[r]};
    }

    // display our matrix
    friend std::ostream& operator<<(std::ostream &cout, matrix<T> m) {
        size_t R = m.size().first, C = m.size().second;
        for (size_t r = 0; r < R; r++) {
            cout << "[";
            for (size_t c = 0; c < C; c++) {
                cout << m(r, c);
                if (c+1 != C) {
                    cout << " ";
                }
                else {
                    cout << "]" << std::endl;
                }
            }
        }
        return cout;
    }

private:
    size_t R, C;
    std::vector<std::vector<T>> M;

};

#endif