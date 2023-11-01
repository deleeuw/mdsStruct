
#include "library.h"
#include "mds.h"

double mds( const size_t n, double** delta, const size_t p, double** z, int** fz, double** d, const size_t MAXITER, const double FCRIT, size_t* lastiter, double* lastdif )
{
  const double EPS = DBL_EPSILON;                                          // 2.2204460492503131e-16
  const double TOL = sqrt( EPS );                                          // 1.4901161193847656e-08
  const double CRIT = sqrt( TOL );                                         // 0.00012207031250000000
  const double TINY = pow( 10.0, ( log10( EPS ) + log10( TOL ) ) / 2.0 );  // 1.8189894035458617e-12

  // allocate memory
  double* imb = getvector( n, 0.0 );
  double* s = getvector( p, 0.0 );

  // initialization
  double scale = 0.0;
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= n; j++ ) {
      const double work = delta[i][j];
      scale += work * work;
    }
  }
  int nfz = 0;
  for ( size_t i = 1; i <= n; i++ ) for ( size_t k = 1; k <= p; k++ ) nfz += fz[i][k];
  double h = ( double )( n - 1 );

  // update distances and calculate normalized stress
  euclidean( n, p, z, d );
  double fold = 0.0;
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= n; j++ ) {
      const double work = delta[i][j] - d[i][j];
      fold += work * work;
    }
  }
  fold /= scale;
  double fnew = 0.0;

  // start main loop
  size_t iter = 0;
  for ( iter = 1; iter <= MAXITER; iter++ ) {

    // compute update for z
    for ( size_t k = 1; k <= p; k++ ) {
      double work = 0.0;
      for ( size_t i = 1; i <= n; i++ ) work += z[i][k];
      s[k] = work;
    }
    for ( size_t i = 1; i <= n; i++ ) {
      double work = 0.0;
      for ( size_t j = 1; j <= n; j++ ) work += imb[j] = ( d[i][j] < TINY ? 0.0 : -1.0 * delta[i][j] / d[i][j] );
      imb[i] = -1.0 * work;
      for ( size_t k = 1; k <= p; k++ ) if ( fz[i][k] == 0 ) {  // update free coordinate
        double t = 0.0;
        for ( size_t j = 1; j <= n; j++ ) t += imb[j] * z[j][k];
        z[i][k] = ( t + s[k] - z[i][k] ) / h;
      }
    }

    // update distances and calculate normalized stress
    euclidean( n, p, z, d );
    fnew = 0.0;
    for ( size_t i = 1; i <= n; i++ ) {
      for ( size_t j = 1; j <= n; j++ ) {
        const double work = delta[i][j] - d[i][j];
        fnew += work * work;
      }
    }
    fnew /= scale;

    // check divergence and convergence
    ( *lastdif ) = fold - fnew;
    if ( ( *lastdif ) <= -1.0 * CRIT ) break;
    const double fdif = 2.0 * ( *lastdif ) / ( fold + fnew );
    if ( fdif <= FCRIT ) break;
    fold = fnew;
  }
  ( *lastiter ) = iter;

  // rotate to principal axes of x when no fixed coordinates are in play
  if ( nfz == 0 ) rotate( n, p, z );
 
  // de-allocate memory
  freevector( imb );
  freevector( s );

  return( fnew );
} // mds

void Cmds( int* rn, double* rdelta, int* rp, double* rz, int* rfz, double* rd, int* rmaxiter, double* rfdif, double* rfvalue )
// Function Cmds() performs multidimensional scaling.
// Copyright (C) 2020 Frank M.T.A. Busing (e-mail: busing at fsw dot leidenuniv dot nl)
// This function is free software:
// you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// You should have received a copy of the GNU General Public License along with this function.
// If not, see <https://www.gnu.org/licenses/>.
{
  // transfer to C
  size_t n = *rn;
  size_t p = *rp;
  size_t MAXITER = *rmaxiter;
  double** delta = getmatrix( n, n, 0.0 );
  for ( size_t j = 1, k = 0; j <= n; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) delta[i][j] = rdelta[k];
  double** z = getmatrix( n, p, 0.0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) z[i][j] = rz[k];
  int** fz = getimatrix( n, p, 0 );
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) fz[i][j] = rfz[k];
  double** d = getmatrix( n, n, 0.0 );
  double FCRIT = *rfdif;

  // run function
  size_t lastiter = 0;
  double lastdif = 0.0;
  double fvalue = mds( n, delta, p, z, fz, d, MAXITER, FCRIT, &lastiter, &lastdif );

  // transfer to R
  for ( size_t j = 1, k = 0; j <= p; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rz[k] = z[i][j];
  for ( size_t j = 1, k = 0; j <= n; j++ ) for ( size_t i = 1; i <= n; i++, k++ ) rd[k] = d[i][j];
  ( *rmaxiter ) = ( int ) ( lastiter );
  ( *rfdif ) = ( int ) ( lastdif );
  ( *rfvalue ) = fvalue;

  // de-allocate memory
  freematrix( delta );
  freematrix( z );
  freeimatrix( fz );
  freematrix( d );

} // Cmds
