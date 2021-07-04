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

#ifndef _MATRIX_UTILITY_H_
#define _MATRIX_UTILITY_H_
#include "matrix.hpp"
#include <tuple>


// Zero/null matrix of size r*c
template <typename T>
matrix<T> zeros(size_t r, size_t c) {
    return matrix<T>(r, c, 0);
}

template <typename T>
matrix<T> zeros(size_t r) {
    return matrix<T>(r, r, 0);
}

// Identity matrix of size r*r
template<typename T>
matrix<T> eye(size_t r) {
    matrix<T> res(r, r, 0);
    for (size_t i = 0; i < r; i++) {
        res(i, i) = 1;
    }
    return res;
}


// matrix m to exponent n (Binary Exponentiation, O(logN))
template<typename T>
matrix<T> pow(matrix<T> m, int k) {
    size_t R = m.size().first, C = m.size().second;
    assert(R == C && k >= 0);
    if (k == 0) {
        return eye<T>(R);
    }
    matrix<T> res = m;
    for (int i = k; i > 0; i>>=1) {
        res = (i&1) ? res*m : res;
        m *= m;
    }
    return res;
}


// Norm of matrix m
// Options available:
// L1: Manhatten Norm
// L2: Eucledian Norm (under development)
// LF: Frobenius Norm
// Linf: Maximum Norm
template<typename T>
double norm(matrix<T> m, std::string NORM = "LF") {
    double res = 0;
    size_t R = m.size().first, C = m.size().second;
    if (NORM == "L1") {
        for (size_t c = 0; c < C; c++) {
            double cur = 0;
            for (size_t r = 0; r < R; r++) {
                cur += abs(m(r, c));
            }
            res = std::max(res, cur);
        }
    }
    else if (NORM == "LF"){
        for (size_t r = 0; r < R; r++) {
            for (size_t c = 0; c < C; c++) {
                res += m(r, c) * m(r, c);
            }
        }
        res = sqrt(res);
    }
    else if (NORM == "Linf") {
        for (size_t r = 0; r < R; r++) {
            double cur = 0;
            for (size_t c = 0; c < C; c++) {
                cur += abs(m(r, c));
            }
            res = std::max(res, cur);
        }
    }
    else {
        assert(false);
    }
    return res;
}


// MaxVal - MinVal over axis
// axis 1: range of each column
// axis 2: range of each row
// template<typename T>
// matrix<T> range(matrix<T> m, int axis = 1) {
//     assert(axis == 1 || axis == 2);
//     size_t R = m.size().first, C = m.size().second;
//     if (axis == 1) {
//         matrix<T> res(1, C);
//         for (size_t c = 0; c < C; c++) {
//             T mnVal = std::numeric_limits<T>::max();
//             T mxVal = std::numeric_limits<T>::lowest();
//             for (size_t r = 0; r < R; r++) {
//                 mnVal = std::min(mnVal, m(r, c));
//                 mxVal = std::max(mxVal, m(r, c));
//             }
//             res(0, c) = mxVal - mnVal;
//         }
//         return res;
//     }
//     else {
//         matrix<T> res(R, 1);
//         for (size_t r = 0; r < R; r++) {
//             T mnVal = std::numeric_limits<T>::max();
//             T mxVal = std::numeric_limits<T>::lowest();
//             for (size_t c = 0; c < C; c++) {
//                 mnVal = std::min(mnVal, m(r, c));
//                 mxVal = std::max(mxVal, m(r, c));
//             }
//             res(r, 0) = mxVal - mnVal;
//         }
//         return res;
//     }
// }


// uniform range (similar working to python)
template<typename T>
matrix<T> range(T start, T stop, T step = 1) {
    if (step == 0) {
    	return matrix<T>(1, 0);
    }
    int n = ceil(1.0 * (stop - start) / step);
    assert(n >= 0);
    matrix<T> res(1, n);
    for (size_t i = 0; i < n; i++) {
    	res(0, i) = start + i * step;
    }
    return res;
}


template<typename T>
matrix<T> range(T stop) {
	return range<T>(0, stop, 1);
}


// divide data by factor
template<typename T>
matrix<T> rescale(const matrix<T> &m, T factor=0) {
    if (factor == 0) {
        factor = norm(m);
        assert(factor != 0);
    }
    size_t R = m.size().first, C = m.size().second;
    matrix<T> res(R, C);
    for (size_t c = 0; c < C; c++) {
        for (size_t r = 0; r < R; r++) {
            res(r, c) /= factor;
        }
    }
    return res;
}


// Householder Method for QR decomposition
template<typename T>
std::tuple<matrix<T>, matrix<T>> qr_householder(matrix<T> m) {
    size_t R = m.size().first, C = m.size().second;
    // Final Q, R vertices
    matrix<T> mQ(R, C, 0), mR(C, C, 0);
    // Array of identity matrices
    std::vector<matrix<T>> qv(R, eye<T>(R));
    matrix<T> z(m), z1;

    for (size_t k = 0; k < C && k < R - 1; k++) {
        matrix<T> e(R, 1);
        z1 = z.compute_minor(k);
        matrix<T> x = z1.column(k);

        T a = norm(x);
        if (m(k, k) > 0) a *= -1;

        for (size_t r = 0; r < R; r++)
            e(r, 0) = (r == k);
        e = rescale(x + a*e);

        qv[k] -= 2 * e * e.t();
        z = qv[k] * z1;
    }

    mQ = qv[0];
    for (size_t i = 1; i < C && i < R - 1; i++) {
        z1 = qv[i] * mQ;
        mQ = z1;
    }

    mR = mQ * m;
    mQ = mQ.t();
    return std::make_tuple(mQ, mR);
}


// QR decomposition
template<typename T>
std::tuple<matrix<T>, matrix<T>> qr(matrix<T> A) {
    size_t m = A.size().first, n = A.size().second;
    matrix<T> Q(m, n, 0), R(n, n, 0);
    for (size_t k = 0; k < n; k++) {
        matrix<T> v = A({0, m-1}, {k, k});
        if (k >= 1) {
            matrix<T> temp = Q({0, m-1}, {0, k-1}).t() * A({0, m-1}, {k, k});
            for (size_t i = 0; i < k; i++) {
                R(i, k) = temp(i, 0);
            }
            v -= Q({0, m-1}, {0, k-1}) * R({0, k-1}, {k, k});
        }
        R(k, k) = norm(v);
        for (size_t i = 0; i < m; i++) {
            Q(i, k) = v(i, 0) / R(k, k);
        }
    }
    return std::make_tuple(Q, R);
}

#endif