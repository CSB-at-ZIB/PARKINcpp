// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// Original (Fortran/Matlab) code by
// P.Deuflhard, U.Nowak, L.Weimann, (c) 1993 - 2006
//
// first added : 16.03.2010 td
// last changed:
//

#include <cassert>
#include <cmath>
#include <algorithm>

#include "QRDecomp.h"

using namespace PARKIN;


//-----------------------------------------------------------------------------
void
QRDecomp::getFirstFactor(Matrix& Q) const
{
    Real   s;
    Vector b;
    long   m = _qrA.nr();
    // long   n = _qrA.nc();

    // if ( !isValid() ) return;

    Q.zeros(m,m);

    for (long k = 1; k <= m; ++k)
    {
        b.zeros(m);
        b(k) = 1.0;

        for (long j = 1; j <= _rank; ++j)
        {
            // s = _qrA(k,j);
            s = 0.0;
            for (long l = j; l <= m; ++l)
                s += _qrA(l,j) * b(l);

            s /= _diag(j)*_qrA(j,j);

            for (long l = j; l <= m; ++l)
                b(l) += _qrA(l,j) * s;
        }

        Q.set_rowm(k) = b.t();
    }

    return;
}
//-----------------------------------------------------------------------------
Matrix
QRDecomp::getFirstFactor() const
{
    Matrix Q;
    getFirstFactor(Q);
    return Q;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void
QRDecomp::getSecondFactor(Matrix& R) const
{
    long m = _qrA.nr();
    long n = _qrA.nc();
    long r = std::min(m,_rank);

    // if ( !isValid() ) return;

    R.zeros(m,n);

    for (long j = 1; j <= r; ++j)
    {
        R(j,j) = _diag(j);
        for (long k = j+1; k <= n; ++k) R(j,k) = _qrA(j,k);
    }

    return;
}
//-----------------------------------------------------------------------------
Matrix
QRDecomp::getSecondFactor() const
{
    Matrix R;
    getSecondFactor(R);
    return R;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void
QRDecomp::getThirdFactor(Matrix& P) const
{
    long n = _qrA.nc();

    // if ( !isValid() ) return;

    P.zeros(n,n);

    for (long k = 1; k <= n; ++k)
        P(k,_pivot(k)) = 1.0;

    return;
}
//-----------------------------------------------------------------------------
Matrix
QRDecomp::getThirdFactor() const
{
    Matrix P;
    getThirdFactor(P);
    return P;
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
ErrType const&
QRDecomp::decompose(Matrix const& A)
{
    long m = A.nr();
    long n = A.nc();

	_qrA = A;
	// _qrAH.zeros(n,n);
	_diag.zeros(n);
	_pivot.zeros(n);
    // for (long j = 1; j <= n; ++j) _pivot(j) = j;

	_rank = _rankMax;
	_cond = _condMax;

	_valid = computeHouseholderReflections(m,n);

	// _kerdim = (_rank < n) ? n - _rank : 0;

	return *_err;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
ErrType const&
QRDecomp::solve(Matrix const& A, Vector const& b, Vector& x)
{
    Vector c = b;
    decompose(A);

    if ( isValid() ) solve(c,x);

    return *_err;
}

//-----------------------------------------------------------------------------
ErrType const&
QRDecomp::solve(Vector& b, Vector& x) const
{
    Vector  v;
    Real    s;  // , sh;
    // long    rk1, j1;
    long    m = _qrA.nr();
    long    n = _qrA.nc();

    v.zeros(n);

	// --------------------------------------------------------------------
	// 1 Solution for pseudo-rank zero
	if ( _rank == 0 )
	{
	    x.zeros(n);
	    _err->setIerr(0);
	    return *_err;
	}

	// --------------------------------------------------------------------
	// 2 Householder transformation of right-hand side
    if ( (m != 1) || (n != 1) )
    {
        for (long j = 1; j <= _rank; ++j)
        {
            s = 0.0;
            for (long l = j; l <= m; ++l)
                s += _qrA(l,j) * b(l);
            s /= _diag(j)*_qrA(j,j);

            for (long l = j; l <= m; ++l)
                b(l) += _qrA(l,j)*s;
        }
    }

	// --------------------------------------------------------------------
	// 3 Solution for upper triangular system

	_pinv -> solveR(_qrA, _diag, _rank, b, v);

    /*

    rk1 = _rank+1;
    for (long jj = 1; jj <= _rank; ++jj)
    {
        long j;
        j = rk1 - jj;
        j1 = j+1;
        s = b(j);
        if (jj != 1)
        {
            sh = 0.0;
            for (long l = j1; l <= _rank; ++l)
                sh += _qrA(j,l) * v(l);
            s -= sh;
        }
        v(j) = s / _diag(j);
    }

    */

	// --------------------------------------------------------------------
	// 3.1 Computation of best least-squares solution
	if ( _rank != n )
	{
	    _pinv -> solvePInv(_rank, v);

        /*

	    for (long j = rk1; j <= n; ++j)
	    {
            s = 0.0;
            for (long l = 1; l < j; ++l)
                s += _qrAH(l,j) * v(l);
            v(j) = -s / _diag(j);
	    }
	    j1 = n+1;
	    for (long jj = 1; jj <= n; ++jj)
	    {
	        long j;
	        j = n - jj + 1;
	        s = 0.0;
	        if ( jj != 1 )
	        {
	            sh = 0.0;
                for (long l = j1; l <= n; ++l)
                    sh += _qrAH(j,l) * v(l);
                s = sh;
	        }
	        if ( jj != 1 && j <= _rank )
	        {
	            v(j) -= s;
	        }
	        else
	        {
	            j1 = j;
	            v(j) = -(s + v(j)) / _diag(j);
	        }
	    }

	    */
	}
	// --------------------------------------------------------------------
	// 4 Back-permutation of solution components
    x.zeros(n);
    for (long l = 1; l <= n; ++l)
        x(_pivot(l)) = v(l);

    _err->setIerr(0);
	return *_err;
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
ErrType const&
QRDecomp::solveR(Vector const& b, Vector& x) const
{
    Vector  v;
    // Real    s, sh;
    // long    rk1, j1;
    // long    m = _qrA.nr();
    long    n = _qrA.nc();

    v.zeros(n);

	// --------------------------------------------------------------------
	// 1 Solution for pseudo-rank zero
	if ( _rank == 0 )
	{
	    x.zeros(n);
	    _err->setIerr(0);
	    return *_err;
	}

	// --------------------------------------------------------------------
	// 2 Householder transformation of right-hand side
	/*
    if ( (kred >= 0) && (m != 1) || (n != 1) )
    {
        for (long j = 1; j <= _rank; ++j)
        {
            s = 0.0;
            for (long l = j; l <= m; ++l)
                s += _qrA(l,j) * b(l);
            s /= _diag(j)*_qrA(j,j);

            for (long l = j; j <= m; ++l)
                b(l) += _qrA(l,j)*s;
        }
    }
    */

	// --------------------------------------------------------------------
	// 3 Solution for upper triangular system

	_pinv -> solveR(_qrA, _diag, _rank, b, v);

    /*

    rk1 = _rank+1;
    for (long jj = 1; jj <= _rank; ++jj)
    {
        long j;
        j = rk1 - jj;
        j1 = j+1;
        s = b(j);
        if (jj != 1)
        {
            sh = 0.0;
            for (long l = j1; l <= _rank; ++l)
                sh += _qrA(j,l) * v(l);
            s -= sh;
        }
        v(j) = s / _diag(j);
    }

    */

	// --------------------------------------------------------------------
	// 3.1 Computation of best least-squares solution
	if ( _rank != n )
	{
	    _pinv -> solvePInv(_rank, v);

        /*

	    for (long j = rk1; j <= n; ++j)
	    {
            s = 0.0;
            for (long l = 1; l < j; ++l)
                s += _qrAH(l,j) * v(l);
            v(j) = -s / _diag(j);
	    }
	    j1 = n+1;
	    for (long jj = 1; jj <= n; ++jj)
	    {
	        long j;
	        j = n - jj + 1;
	        s = 0.0;
	        if ( jj != 1 )
	        {
	            sh = 0.0;
                for (long l = j1; l <= n; ++l)
                    sh += _qrAH(j,l) * v(l);
                s = sh;
	        }
	        if ( jj != 1 && j <= _rank )
	        {
	            v(j) -= s;
	        }
	        else
	        {
	            j1 = j;
	            v(j) = -(s + v(j)) / _diag(j);
	        }
	    }

	    */
	}
	// --------------------------------------------------------------------
	// 4 Back-permutation of solution components
    x.zeros(n);
    for (long l = 1; l <= n; ++l)
        x(_pivot(l)) = v(l);

    _err->setIerr(0);
	return *_err;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
ErrType const&
QRDecomp::setNewRank(long rank)
{
    long m = _qrA.nr();
    long n = _qrA.nc();

    if ( !isValid() )
    {
        _err->setIerr(-1);
        return *_err;
    }

    if (rank < _rank)
    {
        _rank = rank;
        // --------------------------------------------------------------------
        // 1 Initialisation
        if (_rank > n) _rank = n;
        if (_rank > m) _rank = m;

        // _qrAH.zeros(n,n);

        // --------------------------------------------------------------------
        // 1.1 Special case m = n = 1
        if ( (m == 1) && (n == 1) )
        {
            _pivot(1) = 1;
            _diag(1) = _qrA(1,1);
            if (_diag(1) == 0.0)
            {
                _rank = 0;
                _cond = 1.0 / EPMACH;
                // _valid = true;
            }
            else
            {
                _cond = 1.0;
                // _valid = true;
            }
            _valid = true;

            _err->setIerr(0);
            return *_err;
        }

        // --------------------------------------------------------------------
        // 3 Rank-deficient pseudo-inverse
        _valid = computePseudoInverse();

        // --------------------------------------------------------------------
        // 9 Exit
        if (!_valid) return *_err;
        // 9.2 Subcondition of least squares part
        if ( _diag(_rank) != 0.0 )
            _cond = std::fabs( _diag(1) / _diag(_rank) );
        else
            _cond = 0.0;
    }

    _err->setIerr(0);
    return *_err;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
Real
QRDecomp::projectIntoSubspace(Vector const& v) const
{
    Vector  w;
    long    n = _qrA.nc();
//    Real    del, s, s2;
//    long    rk1, n = _qrA.nc();

    // if ( !isValid() ) return 0.0;

    w.zeros(n);
    for (long j = 1; j <= n; ++j) { w(j) = v(_pivot(j)); }

    return _pinv -> projectIntoSubspace(_rank, w);

    /*

    rk1 = _rank + 1;
    del = 0.0;

    for (long j = rk1; j <= n; ++j)
    {
        s2 = 0.0;
        for (long l = 1; l < j; ++l) s2 += _qrAH(l,j)*w(l);

        s    = ( w(j) - s2 )/ _diag(j);
        del += s*s;
        w(j) = s;
    }

    return del;

    */
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
bool
QRDecomp::computeHouseholderReflections(long m, long n)
{
	const Real  reduce = 0.05;
	bool        jd, data, qloop;
	int         stage;
	long        rankh, rk1, k1 = 0, k = 1;
	Real        s, s1, t, h, hmax, dd = 0.0;

	// --------------------------------------------------------------------
	// 1 Initialisation
	if (_rank > n) _rank = n;
	if (_rank > m) _rank = m;

    //_qrAH.zeros(n,n);
    t = 0.0;

	// --------------------------------------------------------------------
	// 1.1 Special case m = n = 1
	if ( (m == 1) && (n == 1) )
	{
		_pivot(1) = 1;
		_diag(1) = _qrA(1,1);
		if (_diag(1) == 0.0)
		{
			_rank = 0;
			_cond = 1.0 / EPMACH;
			// _valid = false;
		}
		else
		{
			_cond = 1.0;
			// _valid = true;
		}
		_valid = true;

        _err->setIerr(0);
		return _valid;
	}

	// --------------------------------------------------------------------
	// 1.2 Initialise pivot-vector
	for (long j = 1; j <= n; ++j)  _pivot(j) = j;

	// --------------------------------------------------------------------
	// 2 Householder triangulation
	jd = true;
	stage = 0;
    rankh = _rank;
    // mh = m;
    data = true;
	rk1 = _rank;
	for (/*long*/ k = 1; k <= rk1; ++k) // for each column
	{
		qloop = true;
		while (qloop)
		{
		qloop = false;
		stage = 1;
		if (k != n)
		{
			long jj;

			k1 = k + 1;
			jd = true;
			while (jd)
			{
				// --------------------------------------------
				// 2.0 Column squared norms
				for (long j = k; j <= n; ++j)
				{
					s = 0.0;
					for (long l = k; l <= m; ++l)
						s += ( _qrA(l,j) * _qrA(l,j) );
					_diag(j) = s;
				}

				// --------------------------------------------
				// 2.1 Column pivoting
				s1 = _diag(k);
				jj = k;
				for (long l = k; l <= n; ++l)
				{
					if (_diag(l) > s1)
					{
						s1 = _diag(l);
						jj = l;
					}
				}
				h = _diag(jj);
				if ( jd ) hmax = h / std::max( 10.0, _cond*reduce );
				jd = false;
				if (h < hmax) jd = true;                 // ???
			} // while jd

			// ----------------------------------------------------
			// 2.2 Column swap
			if (jj != k)
			{
				// SWAP( _pivot(jj), _pivot(k) );
				std::swap( _pivot(jj), _pivot(k) );
				_diag(jj) = _diag(k);                    // ???
				for (long l = 1; l <= m; ++l)
					// SWAP( _A(l,jj), _A(l,k) );
					std::swap( _qrA(l,jj), _qrA(l,k) );
			}
		} // if (k != n)

		h = 0.0;
		for (long l = k; l <= m; ++l)
		{
			h += ( _qrA(l,k) * _qrA(l,k) );
		}
		t = std::sqrt(h);

		// ------------------------------------------------------------
		// 2.3.0 A priori test on pseudo-rank
		if ( k == 1 )  dd = t / _cond;
		if ( t <= dd || k > rankh )
		{
			// ----------------------------------------------------
			// 2.3.1 Rank reduction
			rankh = k - 1;
			if ( data )
			{
				_rank = rankh;
				stage = 3;
			}
		}
		} // while qloop

		// ------------------------------------------------------------
		// 2.4 Householder step
		if (stage == 1)
		{
			s = _qrA(k,k);
			// t = -SIGN( t, s );
			t = (s==0.0) ? -std::fabs(t)
			             : -std::fabs(t)*((s < 0.0) ? -1.0 : 1.0);
			_diag(k) = t;
			_qrA(k,k) = s - t;
			if (k != n)
			{
				t = 1.0/(h - s*t);
				for (long j = k1; j <= n; ++j)
				{
					s = 0.0;
					for (long l = k; l <= m; ++l)
					{
						s += _qrA(l,k) * _qrA(l,j);
					}
					s *= t;
					s1 = -s;
					if (s != 0.0)
					{
						// Update the subcolumns
						for (long l = k; l <= m; ++l)
							_qrA(l,j) += _qrA(l,k)*s1;
					}
					// Update the norm of subcolumns
					_diag(j) -= ( _qrA(k,j) * _qrA(k,j) );
				}
				if ( k == rk1 ) stage = 3;
			}
			else
			{
				stage = 4;
			}
		} // if stage == 1, i.e. Householder step

		// Exit for-loop if ...
		if (stage > 1) break;

	} // for each column k

    _err->setIerr(0);
    _valid = true;
	// --------------------------------------------------------------------
	// 3 Rank-deficient pseudo-inverse
	if (stage == 3)
	{
		_valid = computePseudoInverse();
	}

	// --------------------------------------------------------------------
	// 9 Exit
	if (!_valid) return false;
	// 9.2 Subcondition of least squares part
	if (k == _rank)
	{
		t = _diag(_rank);
	}
	if ( t != 0.0 )
	{
		_cond = std::fabs( _diag(1) / t );
	}
	else
	{
		_cond = 0.0;
	}

	return _valid;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
bool QRDecomp::computePseudoInverse()
{
    return _pinv -> prepare( _qrA, _diag, _rank );
}
//-----------------------------------------------------------------------------
