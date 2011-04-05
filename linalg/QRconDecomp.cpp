// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-11-05 td
// last changed:
//
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

#include "QRconDecomp.h"

using namespace PARKIN;

//-----------------------------------------------------------------------------
void
QRconDecomp::getFirstFactor(Matrix& Q) const
{
    if ( !isValid() ) return;
    return;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void
QRconDecomp::getSecondFactor(Matrix& R) const
{
    if ( !isValid() ) return;
    return;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
ErrType const&
QRconDecomp::decompose(Matrix const& A)
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
QRconDecomp::solve(Matrix const& A, Vector const& b, Vector& x)
{
    Vector c = b;
    decompose(A);

    if ( isValid() ) solve(c,x);

    return *_err;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
ErrType const&
QRconDecomp::solve(Vector& b, Vector& x) const
{
    Vector  v;
    Real    s;  //, sh;
    long    mh; //, rk1, j1;
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
	/*
	if ( (_rank <= _rankc) && (_rank != n) )
	{
	    // long ranc1 = _rankc + 1;
	    // v.set_row(ranc1,n) = zeros(n-ranc1+1,1);
	}
	*/

	// --------------------------------------------------------------------
	// 2 Constrained Householder transformation of right-hand side
    if ( (m != 1) || (n != 1) )
    {
        mh = ( _rankc == 0 ) ? m : _mcon;
        for (long j = 1; j <= _rank; ++j)
        {
            s = 0.0;
            for (long l = j; l <= mh; ++l)
                s += _qrA(l,j) * b(l);
            s /= _diag(j)*_qrA(j,j);

            for (long l = j; l <= m; ++l)
                b(l) += _qrA(l,j)*s;

            if ( j == _rankc ) mh = m;
        }
    }

	// --------------------------------------------------------------------
	// 3 Solution for upper triangular system

	_pinv -> solveR( _qrA, _diag, _rank, b, v);

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
	// 3.1 Computation of best constrained least-squares solution
	if ( (_rank != n) && (_rank != _rankc) )
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
QRconDecomp::solveR(Vector const& b, Vector& x) const
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
	/*
	if ( (_rank <= _rankc) && (_rank != n) )
	{
	    // long ranc1 = _rankc + 1;
	    // v.set_row(ranc1,n) = zeros(n-ranc1+1,1);
	}
	*/

	// --------------------------------------------------------------------
	// 2 Constrained Householder transformation of right-hand side
	/*
    if ( (kred >= 0) && (m != 1) || (n != 1) )
    {
        mh = ( _rankc == 0 ) ? m : _mcon;
        for (long j = 1; j <= _rank; ++j)
        {
            s = 0.0;
            for (long l = j; l <= mh; ++l)
                s += _qrA(l,j) * b(l);
            s /= _diag(j)*_qrA(j,j);

            for (long l = j; j <= m; ++l)
                b(l) += _qrA(l,j)*s;

            if ( j == _rankc ) mh = m;
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
	// 3.1 Computation of best constrained least-squares solution
	if ( (_rank != n) && (_rank != _rankc) )
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
QRconDecomp::setNewRank(long rank)
{
    long m = _qrA.nr();
    long n = _qrA.nc();

    if ( !isValid() )
    {
        _err->setIerr(-1);
        return *_err;
    }

    if ( rank < _rank )
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

        _err->setIerr(0);
        // --------------------------------------------------------------------
        // 3 Rank-deficient pseudo-inverse
        _valid = computePseudoInverse();

        // --------------------------------------------------------------------
        // 9 Exit
        if (!_valid) return *_err;
        // 9.1 Subcondition of constrained part
        if ( _rankc != 0 )
        {
            _condc = _diag(_rankc);
            if ( _condc != 0.0 ) _condc = std::fabs( _diag(1)/_condc );
        }
        else
        {
            _condc = 0.0;
        }
        // 9.2 Subcondition of least squares part
        /*
        if ( -1 == _rank ) t = _diag(_rank);
        if ( (_rankc+1 <= _rank) && ( t != 0.0 )
            _cond = std::fabs( _diag(_rankc+1)/t );
        else
        */
        if ( _diag(_rank) != 0.0)
            _cond = std::fabs( _diag(_rankc+1) / _diag(_rank) );
        else
            _cond = 0.0;
    }

    return *_err;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
bool QRconDecomp::computeHouseholderReflections(long m, long n)
{
	const Real  reduce = 0.05;
	bool        jd, data, qloop;
	int         stage;
	long        mh, rankh, ranc1, rk1, k1 = 0, k = 1;
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
	ranc1 = _rankc + 1;
	mh = _mcon;
	rankh = _rankc;
	data = false;
	if ( mh == 0 )
	{
	    rankh = _rank;
	    mh = m;
	    data = true;
	}
	stage = 0;
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
					for (long l = k; l <= mh; ++l)
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
		for (long l = k; l <= mh; ++l)
		{
			h += ( _qrA(l,k) * _qrA(l,k) );
		}
		t = std::sqrt(h);

		// ------------------------------------------------------------
		// 2.3.0 A priori test on pseudo-rank
		if ( (k == 1) || (k == ranc1) )  dd = t / _cond;
		if ( (t <= dd) || (k > rankh) )
		{
			// ----------------------------------------------------
			// 2.3.1 Rank reduction
			rankh = k - 1;
			if ( (mh != _mcon) || data )
			{
				_rank = rankh;
				stage = ( _rankc == _rank ) ? 4 : 3;
			}
			else
			{
			    _rankc = rankh;
			    if ( _rankc != _mcon )
			    {
			        mh = m;
			        rankh = _rank;
			        jd = data = qloop = true;
			    }
			    else
			    {
			        _err->setIerr(-999);
			        return false;
			    }
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
					for (long l = k; l <= mh; ++l)
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
				if ( k == _rankc ) { mh = m; rankh = _rank; jd = true; }
				if ( k == rk1    ) { stage = 3; }
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
	// 9.1 Subcondition of constrained part
	if ( _rankc != 0 )
	{
	    _condc = _diag(_rankc);
	    if ( _condc != 0.0 )
            _condc = std::fabs( _diag(1) / _condc );
	}
	else
	{
	    _condc = 0.0;
	}
	// 9.2 Subcondition of least squares part
	if (k == _rank)
	{
		t = _diag(_rank);
	}
	if ( (_rankc+1 <= _rank) && (t != 0.0) )
	{
		_cond = std::fabs( _diag(_rankc+1) / t );
	}
	else
	{
		_cond = 0.0;
	}

	return _valid;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
bool QRconDecomp::computePseudoInverse()
{
    return _pinv -> prepare( _qrA, _diag, _rank );
}
//-----------------------------------------------------------------------------
