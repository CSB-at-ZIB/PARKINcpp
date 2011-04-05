// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-11-15 td
// last changed:
//
#include <iostream>
#include <vector>

#include "linalg/Matrix.h"
#include "linalg/Vector.h"
// #include "linalg/QRDecomp.h"
#include "linalg/QRconDecomp.h"

using namespace PARKIN;

QRconDecomp doit(Matrix const& A)
{
    long rank = std::min( A.nr(), A.nc() );
    QRconDecomp qr;

    qr = A.factorQRcon(0,rank,1.0/sqrtEPMACH);

    return qr;
}

int testlinalg()
{
    Vector                  x0, x, b0, b;
    Matrix                  A;
    // QRDecomp                qrA;
    QRconDecomp             qrA;
    long                    m = 4, n = 7;
    std::vector<Vector>     sol;

    sol.clear();
    sol.resize(n);

    x0.ones(n);
    A.zeros(m,n);

    for (long j = 1; j <= A.nr(); ++j)
        for (long k = 1; k <= A.nc(); ++k)
            A(j,k) = /*exp*/( - 1.0 / (j+k-1) );

    // A.randm(m,n);
    b0 = A*x0;
    b = b0;

    std::cout << std::fixed;

    std::cout << std::endl;
    std::cout << " b0 = A*x0\n";
    std::cout << "    = " << b0.t(); // << std::endl;
    std::cout << " x0 = " << x0.t(); // << std::endl;

    std::cout << std::endl;
    std::cout << "pre : rk = " << qrA.getRank() <<
                    ", sc = " << qrA.getSubCond() <<
                    std::endl;

    // qrA = A.factorQR();
    // qrA = A.factorQRcon();
    qrA = doit(A);              // just to check cpy-ctor/assign

    // std::cout << std::endl;
    std::cout << "post: rk = " << qrA.getRank() <<
                    ", sc = " << qrA.getSubCond() <<
                    std::endl;
    std::cout << " A  = " << std::endl;
    std::cout << A;
    std::cout << std::endl;
    std::cout << " d  = " << qrA.getDiag().t();
    std::cout << std::endl;
    std::cout << "qrAH= " << std::endl;
    std::cout << qrA.getMatH();

    qrA.solve(b,x);

    std::cout << std::endl;
    std::cout << " x  = " << x.t();
    std::cout << " b  = " << b.t();
    std::cout << " b0 = " << b0.t();
    std::cout << std::endl;

    std::cout << "sol = " << sol.size() << "\n";

    return 0;
}
