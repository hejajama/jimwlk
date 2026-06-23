#include "Matrix.h"

#include <cmath>
#include <complex>
#include <sstream>
#include <string>
#include <iostream>

using namespace std;

#include <gtest/gtest.h>

namespace
{

bool approx(double a, double b, double eps = 1e-10)
{
    return std::abs(a - b) <= eps;
}

bool approxComplex(const std::complex<double>& a, const std::complex<double>& b, double eps = 1e-10)
{
    return std::abs(a - b) <= eps;
}

bool approxMatrix(const Matrix& a, const Matrix& b, double eps = 1e-10)
{
    if (a.getNDim() != b.getNDim())
    {
        return false;
    }
    for (int i = 0; i < a.getNN(); ++i)
    {
        if (!approxComplex(a(i), b(i), eps))
        {
            return false;
        }
    }
    return true;
}

Matrix make3x3(std::complex<double> a00, std::complex<double> a01, std::complex<double> a02,
               std::complex<double> a10, std::complex<double> a11, std::complex<double> a12,
               std::complex<double> a20, std::complex<double> a21, std::complex<double> a22)
{
    Matrix m(3);
    m.set(0, 0, a00);
    m.set(0, 1, a01);
    m.set(0, 2, a02);
    m.set(1, 0, a10);
    m.set(1, 1, a11);
    m.set(1, 2, a12);
    m.set(2, 0, a20);
    m.set(2, 1, a21);
    m.set(2, 2, a22);
    return m;
}

} // namespace

TEST(MatrixTest, ConstructorsAndBasicAccessors)
{
    Matrix z(3);
    EXPECT_EQ(z.getNDim(), 3);
    EXPECT_EQ(z.getNN(), 9);
    for (int i = 0; i < z.getNN(); ++i)
    {
        EXPECT_TRUE(approxComplex(z.get(i), std::complex<double>(0.0, 0.0)));
    }

    Matrix id(3, 1.0);
    EXPECT_TRUE(approxComplex(id.get(0, 0), {1.0, 0.0}));
    EXPECT_TRUE(approxComplex(id.get(0, 1), {0.0, 0.0}));

    Matrix m(3);
    m.set(0, std::complex<double>(1.0, 2.0));
    m.set(0, 1, std::complex<double>(3.0, 4.0));
    m.set(0, 2, std::complex<double>(-1.0, 0.0));
    m.set(1, 0, std::complex<double>(5.0, 6.0));
    m.set(1, 1, std::complex<double>(7.0, 8.0));
    m.set(1, 2, std::complex<double>(2.0, -1.0));
    m.set(2, 0, std::complex<double>(0.5, 1.5));
    m.set(2, 1, std::complex<double>(-2.0, 0.0));
    m.set(2, 2, std::complex<double>(1.0, 1.0));

    EXPECT_TRUE(approxComplex(m.get(0), {1.0, 2.0}));
    EXPECT_TRUE(approxComplex(m.get(0, 1), {3.0, 4.0}));
    EXPECT_TRUE(approxComplex(m.get(0, 2), {-1.0, 0.0}));

    m.setRe(0, 9.0);
    m.setIm(0, -2.0);
    EXPECT_TRUE(approx(m.getRe(0), 9.0));
    EXPECT_TRUE(approx(m.getIm(0), -2.0));

    m.setRe(2, 2, 10.0);
    m.setIm(2, 2, -3.0);
    EXPECT_TRUE(approxComplex(m.get(2, 2), {10.0, -3.0}));

    EXPECT_TRUE(approxComplex(m(1, 0), m.get(1, 0)));
    EXPECT_TRUE(approxComplex(m(4), m.get(4)));

    std::complex<double>* ptr = m.data();
    ptr[0] = std::complex<double>(11.0, 12.0);
    const Matrix& cm = m;
    const std::complex<double>* cptr = cm.data();
    EXPECT_TRUE(approxComplex(cptr[0], {11.0, 12.0}));

    m.setZero();
    for (int i = 0; i < m.getNN(); ++i)
    {
        EXPECT_TRUE(approxComplex(m(i), {0.0, 0.0}));
    }

    Matrix cpy(3);
    cpy = z;
    EXPECT_TRUE(cpy == z);

    Matrix mvSrc(3, 1.0);
    Matrix mvDst(3);
    mvDst = std::move(mvSrc);
    EXPECT_TRUE(approxComplex(mvDst.get(0, 0), {1.0, 0.0}));
}

TEST(MatrixTest, ArithmeticAndComparisonOperators)
{
    Matrix a = make3x3({1.0, 1.0}, {2.0, 0.0}, {0.0, 0.0}, 
                       {0.0, -1.0}, {3.0, 2.0}, {1.0, 0.0},
                       {0.5, 0.0}, {0.0, 0.5}, {2.0, 0.0});
    Matrix b = make3x3({-1.0, 0.0}, {0.5, 0.5}, {0.0, 1.0},
                       {4.0, -1.0}, {1.0, 0.0}, {0.0, 0.0},
                       {1.0, 1.0}, {0.0, -1.0}, {1.0, 0.0});

    Matrix c = a + b;
    EXPECT_TRUE(approxComplex(c.get(0, 0), {0.0, 1.0}));

    Matrix d = a - b;
    EXPECT_TRUE(approxComplex(d.get(0, 0), {2.0, 1.0}));

    Matrix e = a;
    e += b;
    EXPECT_TRUE(e == c);

    e -= b;
    EXPECT_TRUE(e == a);

    Matrix f = a * 2.0;
    EXPECT_TRUE(approxComplex(f.get(1, 1), {6.0, 4.0}));

    Matrix g = 2.0 * a;
    EXPECT_TRUE(g == f);

    
    std::complex<double> s(0.0, 1.0);
    Matrix h = s * a;
    EXPECT_TRUE(approxComplex(h.get(0, 0), {-1.0, 1.0}));

    Matrix m = a;
    m *= std::complex<double>(2.0, 0.0);
    EXPECT_TRUE(approxComplex(m.get(0, 0), {2.0, 2.0}));

    m /= std::complex<double>(2.0, 0.0);
    EXPECT_TRUE(m == a);

    m /= 2.0;
    EXPECT_TRUE(approxComplex(m.get(1, 1), {1.5, 1.0}));

    Matrix q = a / 2.0;
    EXPECT_TRUE(approxComplex(q.get(0, 0), {0.5, 0.5}));

    Matrix rr = (Matrix(3, 1.0) + a) - b;
    EXPECT_EQ(rr.getNDim(), 3);

    EXPECT_TRUE(a == a);
    EXPECT_TRUE(a != b);
}

TEST(MatrixTest, MultiplicationAndUtilityMethods)
{
    Matrix a = make3x3({1.0, 0.0}, {2.0, 0.0}, {0.0, 0.0},
                       {3.0, 0.0}, {4.0, 0.0}, {0.0, 0.0},
                       {0.0, 0.0}, {0.0, 0.0}, {5.0, 0.0});
    Matrix b = make3x3({5.0, 0.0}, {6.0, 0.0}, {0.0, 0.0},
                       {7.0, 0.0}, {8.0, 0.0}, {0.0, 0.0},
                       {0.0, 0.0}, {0.0, 0.0}, {9.0, 0.0});

    Matrix c = a * b;
    EXPECT_TRUE(approxComplex(c.get(0, 0), {19.0, 0.0}));
    EXPECT_TRUE(approxComplex(c.get(2, 2), {45.0, 0.0}));

    Matrix out(3);
    Matrix::mult(a, b, out);
    EXPECT_TRUE(approxMatrix(c, out));

    Matrix axpy(3);
    axpy.addMultiple(2.0, a);
    EXPECT_TRUE(approxComplex(axpy.get(0, 0), {2.0, 0.0}));
    axpy.addMultiple(std::complex<double>(0.0, 1.0), a);
    EXPECT_TRUE(approxComplex(axpy.get(0, 0), {2.0, 1.0}));

    EXPECT_TRUE(approxComplex(a.trace(), {10.0, 0.0}));
    EXPECT_TRUE(approx(a.FrobeniusNorm(), std::sqrt(55.0)));
    EXPECT_TRUE(approx(a.OneNorm(), 6.0));
    EXPECT_TRUE(approx(a.square(), 27.5));
}

TEST(MatrixTest, ConjugationImagInverseAndFunctions)
{
    Matrix a = make3x3({1.0, 2.0}, {3.0, 4.0}, {0.0, 0.0},
                       {5.0, -1.0}, {-2.0, 1.0}, {1.0, 0.0},
                       {0.0, 1.0}, {0.0, -1.0}, {2.0, 0.0});
    Matrix ah = a;
    ah.conjg();
    EXPECT_TRUE(approxComplex(ah.get(0, 1), std::conj(a.get(1, 0))));
    EXPECT_TRUE(approxComplex(ah.get(1, 0), std::conj(a.get(0, 1))));

    Matrix zero(3);
    zero.expm();
    EXPECT_TRUE(approxMatrix(zero, Matrix(3, 1.0)));

    Matrix d = make3x3({2.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
                       {0.0, 0.0}, {3.0, 0.0}, {0.0, 0.0},
                       {0.0, 0.0}, {0.0, 0.0}, {4.0, 0.0});
    Matrix dExp = d;
    dExp.expm();
    EXPECT_TRUE(approxComplex(dExp.get(0, 0), {std::exp(2.0), 0.0}, 1e-8));
    EXPECT_TRUE(approxComplex(dExp.get(1, 1), {std::exp(3.0), 0.0}, 1e-8));
    EXPECT_TRUE(approxComplex(dExp.get(2, 2), {std::exp(4.0), 0.0}, 1e-8));

    Matrix j = make3x3({1.0, 0.0}, {1.0, 0.0}, {0.0, 0.0},
                       {0.0, 0.0}, {1.0, 0.0}, {1.0, 0.0},
                       {0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0});
    Matrix jExp = j;
    jExp.expm();
    const double ee = std::exp(1.0);
    EXPECT_TRUE(approxComplex(jExp.get(0, 0), {ee, 0.0}, 1e-8));
    EXPECT_TRUE(approxComplex(jExp.get(0, 1), {ee, 0.0}, 1e-8));
    EXPECT_TRUE(approxComplex(jExp.get(0, 2), {0.5 * ee, 0.0}, 1e-8));
    EXPECT_TRUE(approxComplex(jExp.get(1, 1), {ee, 0.0}, 1e-8));
    EXPECT_TRUE(approxComplex(jExp.get(1, 2), {ee, 0.0}, 1e-8));
    EXPECT_TRUE(approxComplex(jExp.get(2, 2), {ee, 0.0}, 1e-8));

    Matrix iPlus = make3x3({1.1, 0.0}, {0.0, 0.0}, {0.0, 0.0},
                           {0.0, 0.0}, {1.2, 0.0}, {0.0, 0.0},
                           {0.0, 0.0}, {0.0, 0.0}, {1.3, 0.0});
    Matrix lp = iPlus;
    lp -= Matrix(3, 1.0);
    lp.logm_pade(6);
    EXPECT_TRUE(approxComplex(lp.get(0, 0), {std::log(1.1), 0.0}, 1e-8));
    EXPECT_TRUE(approxComplex(lp.get(1, 1), {std::log(1.2), 0.0}, 1e-8));
    EXPECT_TRUE(approxComplex(lp.get(2, 2), {std::log(1.3), 0.0}, 1e-8));

    Matrix l = iPlus;
    l.logm();
    EXPECT_TRUE(approxComplex(l.get(0, 0), {std::log(1.1), 0.0}, 1e-7));
    EXPECT_TRUE(approxComplex(l.get(1, 1), {std::log(1.2), 0.0}, 1e-7));
    EXPECT_TRUE(approxComplex(l.get(2, 2), {std::log(1.3), 0.0}, 1e-7));

    Matrix jLogInput = make3x3({ee, 0.0}, {ee, 0.0}, {0.5 * ee, 0.0},
                               {0.0, 0.0}, {ee, 0.0}, {ee, 0.0},
                               {0.0, 0.0}, {0.0, 0.0}, {ee, 0.0});
    Matrix jLog = jLogInput;
    jLog.logm();
    EXPECT_TRUE(approxComplex(jLog.get(0, 0), {1.0, 0.0}, 1e-7));
    EXPECT_TRUE(approxComplex(jLog.get(0, 1), {1.0, 0.0}, 1e-7));
    EXPECT_TRUE(approxComplex(jLog.get(0, 2), {0.0, 0.0}, 1e-7));
    EXPECT_TRUE(approxComplex(jLog.get(1, 1), {1.0, 0.0}, 1e-7));
    EXPECT_TRUE(approxComplex(jLog.get(1, 2), {1.0, 0.0}, 1e-7));
    EXPECT_TRUE(approxComplex(jLog.get(2, 2), {1.0, 0.0}, 1e-7));

    Matrix s = make3x3({4.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
                       {0.0, 0.0}, {9.0, 0.0}, {0.0, 0.0},
                       {0.0, 0.0}, {0.0, 0.0}, {16.0, 0.0});
    Matrix root = s;
    root.sqrtm();
    Matrix sq = root * root;
    EXPECT_TRUE(approxMatrix(sq, s, 1e-7));
}

TEST(MatrixTest, FormattingAndStreamOutput)
{
    Matrix a = make3x3({1.25, -2.5}, {3.75, 4.125}, {0.0, -1.0},
                       {2.0, 0.0}, {0.5, 0.5}, {1.0, -1.0},
                       {0.0, 0.0}, {1.5, 0.0}, {-1.0, 2.0});

    std::string text = a.getElementsText();
    EXPECT_FALSE(text.empty());
    EXPECT_NE(text.find("1.25"), std::string::npos);

    std::ostringstream os;
    os << a;
    const std::string printed = os.str();
    EXPECT_NE(printed.find("("), std::string::npos);
}

TEST(MatrixTest, MissingDeclaredButUndefinedMethods)
{
    SUCCEED() << "Matrix.h declares normAm(int), unary operator-, and Matrix/Matrix division, "
              << "but Matrix.cpp does not define them in this codebase.";
}
