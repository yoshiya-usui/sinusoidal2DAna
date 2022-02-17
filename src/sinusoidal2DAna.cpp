//-------------------------------------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) 2021 Yoshiya Usui
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//-------------------------------------------------------------------------------------------------------
#define _USE_MATH_DEFINES 
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <iomanip>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <complex>
#include <cmath>
#include "mkl_lapacke.h"

const static int m_numRead = 100;
const static int m_num = 100;
const static int m_nDof = 2 * m_num;
static double m_T = 0.0;

std::vector< std::complex<double> > m_Itheta1;
std::vector< std::complex<double> > m_Itheta2;
std::vector<double> m_obsLocX; 

std::complex<double> getItheta1( const int k, const int n );
std::complex<double> getItheta2( const int k, const int n );
int getPositionInMatrix( const int row, const int col );
void readModfiedBesselFunction( const std::string& inFileName );
void readObsLocX( const std::string& inFileName );
void calcBetaGammaTM();
void calcBetaGammaTE();
void calcMapOfByEx();

int main( int argc, char* argv[] ){

	if( argc < 3 ){
		std::cerr << "You must specify input file names and period(s)" << std::endl;
		exit(1);
	}
	const std::string fileNameBessel = argv[1];
	readModfiedBesselFunction(fileNameBessel);

	const std::string fileNameObsLoc = argv[2];
	readObsLocX(fileNameObsLoc); 

	m_T = atof(argv[3]);
	std::cout << "Period(s) : " << m_T << std::endl;

	calcBetaGammaTM();

	calcBetaGammaTE();

#ifdef _TEST
	calcMapOfByEx();
#endif

	return 0;

}

void readModfiedBesselFunction( const std::string& inFileName ){

	std::ifstream inFile( inFileName.c_str(), std::ios::in );
	if( inFile.fail() ){
		std::cerr << "File open error !!" << std::endl;
		exit(1);
	}
	std::string str;
	inFile >> str >> str >> str >> str;
	for( int i = 0; i < m_numRead; ++i ){
		for( int j = 0; j < m_numRead; ++j ){
			int ibuf(0);
			inFile >> ibuf;
			assert( ibuf == i );
			inFile >> ibuf;
			assert( ibuf == j );
			double dbuf1(0.0);
			double dbuf2(0.0);
			inFile >> dbuf1;
			inFile >> dbuf2;
#ifdef _FLAT
			if( i == 0 ){
				dbuf1 = 1.0;
				dbuf2 = 0.0;
			}else{
				dbuf1 = 0.0;
				dbuf2 = 0.0;
			}
#endif
			m_Itheta1.push_back( std::complex<double>(dbuf1,dbuf2) );
		}
	}

	inFile >> str >> str >> str >> str;
	for( int i = 0; i < m_numRead; ++i ){
		for( int j = 0; j < m_numRead; ++j ){
			int ibuf(0);
			inFile >> ibuf;
			assert( ibuf == i );
			inFile >> ibuf;
			assert( ibuf == j );
			double dbuf1(0.0);
			double dbuf2(0.0);
			inFile >> dbuf1;
			inFile >> dbuf2;
#ifdef _FLAT
			if( i == 0 ){
				dbuf1 = 1.0;
				dbuf2 = 0.0;
			}else{
				dbuf1 = 0.0;
				dbuf2 = 0.0;
			}
#endif
			m_Itheta2.push_back( std::complex<double>(dbuf1,dbuf2) );
		}
	}
	inFile.close();

}

void readObsLocX( const std::string& inFileName ){

	std::ifstream inFile( inFileName.c_str(), std::ios::in );
	if( inFile.fail() ){
		std::cerr << "File open error !!" << std::endl;
		exit(1);
	}

	int nLocX(0);
	inFile >> nLocX;
	for( int i = 0; i < nLocX; ++i ){
		double dbuf(0.0);
		inFile >> dbuf;
		m_obsLocX.push_back(dbuf);
	}

	inFile.close();
}

std::complex<double> getItheta1( const int k, const int n ){
	const int ipos = k * m_numRead + n;
	return m_Itheta1[ipos];
}

std::complex<double> getItheta2( const int k, const int n ){
	const int ipos = k * m_numRead + n;
	return m_Itheta2[ipos];
}

int getPositionInMatrix( const int row, const int col ){
	const int ipos = row + col * ( 2 * m_num );
	return ipos;
}

void calcBetaGammaTE(){

	const double sigma1 = 3.0;
	const double sigma2 = 0.01;
#ifdef _FLAT
	const double delta = 0.0;
#else
	const double delta = 100.0;
#endif
	const double lambda = 1000.0;
	const double f = 1.0 / m_T;
	const double omega = 2.0 * M_PI * f;
	const double nu = 2.0 * M_PI / lambda;
	const double mu0 = 4.0 * M_PI * 1.0e-7;
	const std::complex<double> E0 = std::complex<double>(1.0, 0.0);
	const std::complex<double> imagUnit = std::complex<double>(0.0, 1.0);

	std::complex<double>* matrix = new std::complex<double>[m_nDof*m_nDof];
	std::complex<double>* rhsVector = new std::complex<double>[m_nDof];
	
	// Initialize
	for( int i = 0; i < m_nDof; ++i ){
		rhsVector[i] = std::complex<double>(0.0, 0.0);
	}
	for( int j = 0; j < m_nDof; ++j ){
		const int offset = j * m_num;
		for( int i = 0; i < m_nDof; ++i ){
			matrix[offset+i] = std::complex<double>(0.0, 0.0);
		}
	}

	// Construct matrix and right-hand-side vector
	for( int k = 0; k < m_num; ++k ){
		{// 1st term of Ey equation
			double factor = 2.0;
			if( k == 0 ){
				factor = 1.0;
			}
			rhsVector[ 2 * k ] -= E0 * pow(-1.0, k) * factor * getItheta1(k,0);
		}
		// 2nd term of Ey equation
		for( int n = 0; n < m_num; ++n ){
			double factor = 1.0;
			if( k == 0 ){
				factor = 0.5;
			}
			if( 2 * abs(n + k) < m_nDof ){
				const std::complex<double> val = factor * getItheta1(k, n); 
				const int ipos = getPositionInMatrix( 2 * abs(n + k), 2 * n );
				matrix[ipos] += val; 
				if( 2 * abs(n - k) < m_nDof ){
					const int ipos = getPositionInMatrix( 2 * abs(n - k), 2 * n );
					matrix[ipos] += val;
				}
			}
		}
		// 3rd term of Ey equation
		for( int n = 0; n < m_num; ++n ){
			double factor = 1.0;
			if( k == 0 ){
				factor = 0.5;
			}
			if( 2 * abs(n + k) < m_nDof ){
				const std::complex<double> val = pow(-1.0, k) * factor * getItheta2(k, n); 
				const int ipos = getPositionInMatrix( 2 * abs(n + k), 2 * n + 1 );
				matrix[ipos] -= val; 
				if( 2 * abs(n - k) < m_nDof ){
					const int ipos = getPositionInMatrix( 2 * abs(n - k), 2 * n + 1 );
					matrix[ipos] -= val; 
				}
			}
		}
		{// 1st term of Bx equation
			const std::complex<double> theta = sqrt(std::complex<double>( 0.0, omega * mu0 * sigma1 ));
			double factor = 2.0;
			if( k == 0 ){
				factor = 1.0;
			}
			rhsVector[2 * k + 1] += theta / std::complex<double>(0.0, omega) * E0 * pow(-1.0, k) * factor * getItheta1(k,0);
		}
		// 2nd term of Bx equation
		for( int n = 0; n < m_num; ++n ){
			const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
			double factor = 1.0;
			if( k == 0 ){
				factor = 0.5;
			}
			if( 2 * abs(n + k) + 1 < m_nDof ){
				const std::complex<double> val = theta / std::complex<double>(0.0, omega) * factor * getItheta1(k, n);
				const int ipos = getPositionInMatrix( 2 * abs(n + k) + 1, 2 * n );
				matrix[ipos] += val; 
				if( 2 * abs(n - k) + 1 < m_nDof ){
					const int ipos = getPositionInMatrix( 2 * abs(n - k) + 1, 2 * n );
					matrix[ipos] += val; 
				}
			}
		}
		// 3rd term of Bx equation
		for( int n = 0; n < m_num; ++n ){
			const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
			double factor = 1.0;
			if( k == 0 ){
				factor = 0.5;
			}
			if( 2 * abs(n + k) + 1 < m_nDof ){
				const std::complex<double> val = theta / std::complex<double>(0.0, omega) * pow(-1.0, k) * factor * getItheta2(k, n);
				const int ipos = getPositionInMatrix( 2 * abs(n + k) + 1, 2 * n + 1 );
				matrix[ipos] += val; 
				if( 2 * abs(n - k) + 1 < m_nDof ){
					const int ipos = getPositionInMatrix( 2 * abs(n - k) + 1, 2 * n + 1 );
					matrix[ipos] += val;
				}
			}
		}
	}

	lapack_complex_double* matrixLapack = new lapack_complex_double[m_nDof*m_nDof];
	for( int i = 0; i < m_nDof * m_nDof; ++i ){
		matrixLapack[i].real = matrix[i].real();
		matrixLapack[i].imag = matrix[i].imag();
	}
	lapack_complex_double* resultVectorLapack = new lapack_complex_double[m_nDof];
	for( int i = 0; i < m_nDof; ++i ){
		resultVectorLapack[i].real = rhsVector[i].real();
		resultVectorLapack[i].imag = rhsVector[i].imag();
	}

	int ierr(0);
	int* ipiv = new int[m_nDof];
	int lda = m_nDof;
	int ldb = m_nDof;

	// Factorize
	ierr = LAPACKE_zgetrf( LAPACK_COL_MAJOR, m_nDof, m_nDof, matrixLapack, lda, ipiv );
	if( ierr > 0 ) {
		std::cerr << "Error : Matrix is singular. ierr = " << ierr << std::endl;
		exit(1);
	}else if( ierr < 0 ){
		std::cerr << "Error : " << -ierr << "-th parameter has illegal value." << std::endl;
		exit(1);
	}

	// Solve
	ierr = LAPACKE_zgetrs( LAPACK_COL_MAJOR, 'N', m_nDof, 1, matrixLapack, lda, ipiv, resultVectorLapack, ldb );
	if( ierr < 0 ){
		std::cerr << "Error : " << -ierr << "-th parameter has illegal value." << std::endl;
		exit(1);
	}

	std::complex<double>* beta = new std::complex<double>[m_num];
	std::complex<double>* gamma = new std::complex<double>[m_num];
	for( int i = 0; i < m_num; ++i ){
		beta[i]  = std::complex<double>( resultVectorLapack[ 2 * i     ].real, resultVectorLapack[ 2 * i     ].imag );
		gamma[i] = std::complex<double>( resultVectorLapack[ 2 * i + 1 ].real, resultVectorLapack[ 2 * i + 1 ].imag );
	}

	const int numObsLoc = static_cast<int>( m_obsLocX.size() );
	std::complex<double>* Ey1 = new std::complex<double>[numObsLoc];
	std::complex<double>* Bx1 = new std::complex<double>[numObsLoc];
	std::complex<double>* Bz1 = new std::complex<double>[numObsLoc];
	std::complex<double>* Ey2 = new std::complex<double>[numObsLoc];
	std::complex<double>* Bx2 = new std::complex<double>[numObsLoc];
	std::complex<double>* Bz2 = new std::complex<double>[numObsLoc];
	for( int iObs = 0; iObs <numObsLoc; ++iObs ){
		const double x = m_obsLocX[iObs];
		const double z = delta * cos(nu * x);
		{// Upper region
			const std::complex<double> theta = sqrt(std::complex<double>( 0.0, omega * mu0 * sigma1 ));
			// Ey1
			Ey1[iObs] = E0 * exp( - theta * z );
			for( int n = 0; n < m_num; ++n ){
				const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
				Ey1[iObs] += beta[n] * exp( theta * z ) * cos( n * nu * x );
			}
			// Bx1
			Bx1[iObs] = - E0 * theta / std::complex<double>(0.0, omega) * exp( - theta * z );
			for( int n = 0; n < m_num; ++n ){
				const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
				Bx1[iObs] += theta / std::complex<double>(0.0, omega) * beta[n] * exp( theta * z ) * cos( n * nu * x ); 
			}
			// Bz1
			Bz1[iObs] = std::complex<double>(0.0, 0.0);
			for( int n = 0; n < m_num; ++n ){
				const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
				Bz1[iObs] += n * nu / std::complex<double>(0.0, omega) * beta[n] * exp( theta * z ) * sin( n * nu * x );
			}
		}
		{// Lower region
			const std::complex<double> theta = sqrt(std::complex<double>( 0.0, omega * mu0 * sigma2 ));
			// Ey2
			Ey2[iObs] = std::complex<double>(0.0,0.0);
			for( int n = 0; n < m_num; ++n ){
				const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
				Ey2[iObs] += gamma[n] * exp( - theta * z ) * cos( n * nu * x );
			}
			// Bx2
			Bx2[iObs] = std::complex<double>(0.0,0.0);
			for( int n = 0; n < m_num; ++n ){
				const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
				Bx2[iObs] -= theta / std::complex<double>(0.0, omega) * gamma[n] * exp( - theta * z ) * cos( n * nu * x );
			}
			// Bz2
			Bz2[iObs] = std::complex<double>(0.0,0.0);
			for( int n = 0; n < m_num; ++n ){
				const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
				Bz2[iObs] += n * nu / std::complex<double>(0.0, omega) * gamma[n] * exp( - theta * z ) * sin( n * nu * x );
			}
		}
	}

	const double rad2deg = 180.0 / M_PI;
	double* app1 = new double[numObsLoc];
	double* phs1 = new double[numObsLoc];
	std::complex<double>* vtf1 = new std::complex<double>[numObsLoc];
	double* app2 = new double[numObsLoc];
	double* phs2 = new double[numObsLoc];
	std::complex<double>* vtf2 = new std::complex<double>[numObsLoc];
	for( int iObs = 0; iObs <numObsLoc; ++iObs ){
		const std::complex<double> z1 = Ey1[iObs] / Bx1[iObs];
		app1[iObs] = std::norm(z1)* mu0 / omega;
		phs1[iObs] = atan2( z1.imag(), z1.real() ) * rad2deg + 180.0;
		vtf1[iObs] =  Bz1[iObs] / Bx1[iObs];
		const std::complex<double> z2 = Ey2[iObs] / Bx2[iObs];
		app2[iObs] = std::norm(z2) * mu0 / omega;
		phs2[iObs] = atan2( z2.imag(), z2.real() ) * rad2deg + 180.0;
		vtf2[iObs] =  Bz2[iObs] / Bx2[iObs];
	}

	std::ofstream oFile( "TE.txt");
	if( oFile.fail() ){
		std::cerr << "File open error !!" << std::endl;
		exit(1);
	}
	oFile.precision(6);
	for( int iObs = 0; iObs <numObsLoc; ++iObs ){
		const double x = m_obsLocX[iObs];
		const double z = delta * cos(nu * x);
		oFile << std::setw(15) << std::scientific << x;
		oFile << std::setw(15) << std::scientific << z;
		oFile << std::setw(15) << std::scientific << app1[iObs];
		oFile << std::setw(15) << std::scientific << phs1[iObs];
		oFile << std::setw(15) << std::scientific << app2[iObs];
		oFile << std::setw(15) << std::scientific << phs2[iObs];
		oFile << std::setw(15) << std::scientific << vtf1[iObs].real();
		oFile << std::setw(15) << std::scientific << vtf1[iObs].imag();
		oFile << std::setw(15) << std::scientific << vtf2[iObs].real();
		oFile << std::setw(15) << std::scientific << vtf2[iObs].imag();
		oFile << std::endl;
	}

	oFile.close();

	std::cout << "Finish" << std::endl;

}

void calcBetaGammaTM(){

	const double sigma1 = 3.0;
	const double sigma2 = 0.01;
#ifdef _FLAT
	const double delta = 0.0;
#else
	const double delta = 100.0;
#endif
	const double lambda = 1000.0;
	const double f = 1.0 / m_T;
	const double omega = 2.0 * M_PI * f;
	const double nu = 2.0 * M_PI / lambda;
	const double mu0 = 4.0 * M_PI * 1.0e-7;
	const std::complex<double> B0 = std::complex<double>(1.0, 0.0);

	std::complex<double>* matrix = new std::complex<double>[m_nDof*m_nDof];
	std::complex<double>* rhsVector = new std::complex<double>[m_nDof];
	
	// Initialize
	for( int i = 0; i < m_nDof; ++i ){
		rhsVector[i] = std::complex<double>(0.0, 0.0);
	}
	for( int j = 0; j < m_nDof; ++j ){
		const int offset = j * m_num;
		for( int i = 0; i < m_nDof; ++i ){
			matrix[offset+i] = std::complex<double>(0.0, 0.0);
		}
	}

	// Construct matrix and right-hand-side vector
	for( int k = 0; k < m_num; ++k ){
		{// 1st term of By equation
			double factor = 2.0;
			if( k == 0 ){
				factor = 1.0;
			}
			rhsVector[ 2 * k ] -= B0 * pow(-1.0, k) * factor * getItheta1(k,0);
		}
		// 2nd term of By equation
		for( int n = 0; n < m_num; ++n ){
			double factor = 1.0;
			if( k == 0 ){
				factor = 0.5;
			}
			if( 2 * abs(n + k) < m_nDof ){
				const int ipos = getPositionInMatrix( 2 * abs(n + k), 2 * n );
				matrix[ipos] += factor * getItheta1(k, n); 
				if( 2 * abs(n - k) < m_nDof ){
					const int ipos = getPositionInMatrix( 2 * abs(n - k), 2 * n );
					matrix[ipos] += factor * getItheta1(k, n);
				}
			}
		}
		// 3rd term of By equation
		for( int n = 0; n < m_num; ++n ){
			double factor = 1.0;
			if( k == 0 ){
				factor = 0.5;
			}
			if( 2 * abs(n + k) < m_nDof ){
				const int ipos = getPositionInMatrix( 2 * abs(n + k), 2 * n + 1 );
				matrix[ipos] -= pow(-1.0, k) * factor * getItheta2(k, n); 
				if( 2 * abs(n - k) < m_nDof ){
					const int ipos = getPositionInMatrix( 2 * abs(n - k), 2 * n + 1 );
					matrix[ipos] -= pow(-1.0, k) * factor * getItheta2(k, n); 
				}
			}
		}
		{// 1st term of Ex equation
			const std::complex<double> theta = sqrt(std::complex<double>( 0.0, omega * mu0 * sigma1 ));
			double factor = 2.0;
			if( k == 0 ){
				factor = 1.0;
			}
			rhsVector[2 * k + 1] -= theta / ( mu0 * sigma1 ) * B0 * pow(-1.0, k) * factor * getItheta1(k,0);
		}
		// 2nd term of Ex equation
		for( int n = 0; n < m_num; ++n ){
			const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
			double factor = 1.0;
			if( k == 0 ){
				factor = 0.5;
			}
			if( 2 * abs(n + k) + 1 < m_nDof ){
				const int ipos = getPositionInMatrix( 2 * abs(n + k) + 1, 2 * n );
				matrix[ipos] -= theta / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
				if( 2 * abs(n - k) + 1 < m_nDof ){
					const int ipos = getPositionInMatrix( 2 * abs(n - k) + 1, 2 * n );
					matrix[ipos] -= theta / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
				}
			}
		}
		// 3rd term of Ex equation
		for( int n = 0; n < m_num; ++n ){
			const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
			double factor = 1.0;
			if( k == 0 ){
				factor = 0.5;
			}
			if( 2 * abs(n + k) + 1 < m_nDof ){
				const int ipos = getPositionInMatrix( 2 * abs(n + k) + 1, 2 * n + 1 );
				matrix[ipos] -= theta / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
				if( 2 * abs(n - k) + 1 < m_nDof ){
					const int ipos = getPositionInMatrix( 2 * abs(n - k) + 1, 2 * n + 1 );
					matrix[ipos] -= theta / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
				}
			}
		}
		//{// 1st term of Ez equation
		//	const std::complex<double> theta = sqrt(std::complex<double>( 0.0, omega * mu0 * sigma1 ));
		//	double factor = 2.0;
		//	if( k == 0 ){
		//		factor = 1.0;
		//	}
		//	if( 2 * abs(k + 2) + 1 < m_nDof && 2 * abs(k - 2) + 1 < m_nDof ){
		//		rhsVector[2 * k + 1         ] +=  0.5 * theta / ( mu0 * sigma1 ) * pow(delta * nu, 2) * B0 * pow(-1.0, k) * factor * getItheta1(k,0);
		//		rhsVector[2 * abs(k + 2) + 1] -= 0.25 * theta / ( mu0 * sigma1 ) * pow(delta * nu, 2) * B0 * pow(-1.0, k) * factor * getItheta1(k,0);
		//		rhsVector[2 * abs(k - 2) + 1] -= 0.25 * theta / ( mu0 * sigma1 ) * pow(delta * nu, 2) * B0 * pow(-1.0, k) * factor * getItheta1(k,0);
		//	}
		//}
		//// 2nd term of Ez equation
		//for( int n = 0; n < m_num; ++n ){
		//	const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
		//	double factor = 1.0;
		//	if( k == 0 ){
		//		factor = 0.5;
		//	}
		//	if( 2 * abs(n + k    ) + 1 < m_nDof && 2 * abs(n - k    ) + 1 < m_nDof && 
		//		2 * abs(n + k + 2) + 1 < m_nDof && 2 * abs(n + k - 2) + 1 < m_nDof &&
		//		2 * abs(n - k + 2) + 1 < m_nDof && 2 * abs(n - k - 2) + 1 < m_nDof ){
		//		int ipos(-1);
		//		ipos = getPositionInMatrix( 2 * abs(n + k    ) + 1, 2 * n );
		//		matrix[ipos] +=  0.5 * pow(delta * nu, 2) * theta / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
		//		ipos = getPositionInMatrix( 2 * abs(n - k    ) + 1, 2 * n );
		//		matrix[ipos] +=  0.5 * pow(delta * nu, 2) * theta / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
		//		ipos = getPositionInMatrix( 2 * abs(n + k + 2) + 1, 2 * n );
		//		matrix[ipos] -= 0.25 * pow(delta * nu, 2) * theta / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
		//		ipos = getPositionInMatrix( 2 * abs(n + k - 2) + 1, 2 * n );
		//		matrix[ipos] -= 0.25 * pow(delta * nu, 2) * theta / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
		//		ipos = getPositionInMatrix( 2 * abs(n - k + 2) + 1, 2 * n );
		//		matrix[ipos] -= 0.25 * pow(delta * nu, 2) * theta / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
		//		ipos = getPositionInMatrix( 2 * abs(n - k - 2) + 1, 2 * n );
		//		matrix[ipos] -= 0.25 * pow(delta * nu, 2) * theta / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
		//	}
		//}
		// 3rd term of Ez equation
		for( int n = 0; n < m_num; ++n ){
			const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
			double factor = 2.0 * static_cast<double>(n);
			if( k == 0 ){
				factor = 1.0 * static_cast<double>(n);
			}
			//if( 2 * abs(n + k + 1) + 1 < m_nDof && 2 * abs(n + k - 1) + 1 < m_nDof && 
			//	2 * abs(n - k + 1) + 1 < m_nDof && 2 * abs(n - k - 1) + 1 < m_nDof ){
			//	int ipos(-1);
			//	ipos = getPositionInMatrix( 2 * abs(n + k + 1) + 1, 2 * n );
			//	matrix[ipos] -= 0.25 * delta * pow(nu, 2) / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
			//	ipos = getPositionInMatrix( 2 * abs(n + k - 1) + 1, 2 * n );
			//	matrix[ipos] += 0.25 * delta * pow(nu, 2) / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
			//	ipos = getPositionInMatrix( 2 * abs(n - k + 1) + 1, 2 * n );
			//	matrix[ipos] -= 0.25 * delta * pow(nu, 2) / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
			//	ipos = getPositionInMatrix( 2 * abs(n - k - 1) + 1, 2 * n );
			//	matrix[ipos] += 0.25 * delta * pow(nu, 2) / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
			//}
			if( 2 * abs(n + k + 1) + 1 < m_nDof && 2 * abs(n + k - 1) + 1 < m_nDof ){
				int ipos(-1);
				ipos = getPositionInMatrix( 2 * abs(n + k + 1) + 1, 2 * n );
				matrix[ipos] -= 0.25 * delta * pow(nu, 2) / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
				ipos = getPositionInMatrix( 2 * abs(n + k - 1) + 1, 2 * n );
				matrix[ipos] += 0.25 * delta * pow(nu, 2) / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
			}
			if( 2 * abs(n - k + 1) + 1 < m_nDof && 2 * abs(n - k - 1) + 1 < m_nDof ){
				int ipos(-1);
				ipos = getPositionInMatrix( 2 * abs(n - k + 1) + 1, 2 * n );
				matrix[ipos] -= 0.25 * delta * pow(nu, 2) / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
				ipos = getPositionInMatrix( 2 * abs(n - k - 1) + 1, 2 * n );
				matrix[ipos] += 0.25 * delta * pow(nu, 2) / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
			}
		}
		//// 4th term of Ez equation
		//for( int n = 0; n < m_num; ++n ){
		//	const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
		//	double factor = 1.0;
		//	if( k == 0 ){
		//		factor = 0.5;
		//	}
		//	if( 2 * abs(n + k    ) + 1 < m_nDof && 2 * abs(n - k    ) + 1 < m_nDof && 
		//		2 * abs(n + k + 2) + 1 < m_nDof && 2 * abs(n + k - 2) + 1 < m_nDof && 
		//		2 * abs(n - k + 2) + 1 < m_nDof && 2 * abs(n - k - 2) + 1 < m_nDof ){
		//		int ipos(-1);
		//		ipos = getPositionInMatrix( 2 * abs(n + k    ) + 1, 2 * n + 1 );
		//		matrix[ipos] +=  0.5 * pow(delta * nu, 2) * theta / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
		//		ipos = getPositionInMatrix( 2 * abs(n - k    ) + 1, 2 * n + 1 );
		//		matrix[ipos] +=  0.5 * pow(delta * nu, 2) * theta / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
		//		ipos = getPositionInMatrix( 2 * abs(n + k + 2) + 1, 2 * n + 1 );
		//		matrix[ipos] -= 0.25 * pow(delta * nu, 2) * theta / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
		//		ipos = getPositionInMatrix( 2 * abs(n + k - 2) + 1, 2 * n + 1 );
		//		matrix[ipos] -= 0.25 * pow(delta * nu, 2) * theta / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
		//		ipos = getPositionInMatrix( 2 * abs(n - k + 2) + 1, 2 * n + 1 );
		//		matrix[ipos] -= 0.25 * pow(delta * nu, 2) * theta / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
		//		ipos = getPositionInMatrix( 2 * abs(n - k - 2) + 1, 2 * n + 1 );
		//		matrix[ipos] -= 0.25 * pow(delta * nu, 2) * theta / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
		//	}
		//}
		// 5th term of Ez equation
		for( int n = 0; n < m_num; ++n ){
			const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
			double factor = 2.0 * static_cast<double>(n);
			if( k == 0 ){
				factor = 1.0 * static_cast<double>(n);
			}
			//if( 2 * abs(n + k + 1) + 1 < m_nDof && 2 * abs(n + k - 1) + 1 < m_nDof && 
			//	2 * abs(n - k + 1) + 1 < m_nDof && 2 * abs(n - k - 1) + 1 < m_nDof ){
			//	int ipos(-1);
			//	ipos = getPositionInMatrix( 2 * abs(n + k + 1) + 1, 2 * n + 1 );
			//	matrix[ipos] += 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
			//	ipos = getPositionInMatrix( 2 * abs(n + k - 1) + 1, 2 * n + 1 );
			//	matrix[ipos] -= 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
			//	ipos = getPositionInMatrix( 2 * abs(n - k + 1) + 1, 2 * n + 1 );
			//	matrix[ipos] += 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
			//	ipos = getPositionInMatrix( 2 * abs(n - k - 1) + 1, 2 * n + 1 );
			//	matrix[ipos] -= 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
			//}
			if( 2 * abs(n + k + 1) + 1 < m_nDof && 2 * abs(n + k - 1) + 1 < m_nDof ){
				int ipos(-1);
				ipos = getPositionInMatrix( 2 * abs(n + k + 1) + 1, 2 * n + 1 );
				matrix[ipos] += 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
				ipos = getPositionInMatrix( 2 * abs(n + k - 1) + 1, 2 * n + 1 );
				matrix[ipos] -= 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
			}
			if( 2 * abs(n - k + 1) + 1 < m_nDof && 2 * abs(n - k - 1) + 1 < m_nDof ){
				int ipos(-1);
				ipos = getPositionInMatrix( 2 * abs(n - k + 1) + 1, 2 * n + 1 );
				matrix[ipos] += 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
				ipos = getPositionInMatrix( 2 * abs(n - k - 1) + 1, 2 * n + 1 );
				matrix[ipos] -= 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
			}
		}
	}

	lapack_complex_double* matrixLapack = new lapack_complex_double[m_nDof*m_nDof];
	for( int i = 0; i < m_nDof * m_nDof; ++i ){
		matrixLapack[i].real = matrix[i].real();
		matrixLapack[i].imag = matrix[i].imag();
	}
	lapack_complex_double* resultVectorLapack = new lapack_complex_double[m_nDof];
	for( int i = 0; i < m_nDof; ++i ){
		resultVectorLapack[i].real = rhsVector[i].real();
		resultVectorLapack[i].imag = rhsVector[i].imag();
	}

	int ierr(0);
	int* ipiv = new int[m_nDof];
	int lda = m_nDof;
	int ldb = m_nDof;

	// Factorize
	ierr = LAPACKE_zgetrf( LAPACK_COL_MAJOR, m_nDof, m_nDof, matrixLapack, lda, ipiv );
	if( ierr > 0 ) {
		std::cerr << "Error : Matrix is singular. ierr = " << ierr << std::endl;
		exit(1);
	}else if( ierr < 0 ){
		std::cerr << "Error : " << -ierr << "-th parameter has illegal value." << std::endl;
		exit(1);
	}

	// Solve
	ierr = LAPACKE_zgetrs( LAPACK_COL_MAJOR, 'N', m_nDof, 1, matrixLapack, lda, ipiv, resultVectorLapack, ldb );
	if( ierr < 0 ){
		std::cerr << "Error : " << -ierr << "-th parameter has illegal value." << std::endl;
		exit(1);
	}

#if 0
	std::complex<double> resultVector[m_nDof];
	for( int i = 0; i < m_nDof; ++i ){
		resultVector[i] = std::complex<double>( resultVectorLapack[i].real, resultVectorLapack[i].imag );
	}

	std::complex<double> tempVector[m_nDof];
	for( int i = 0; i < m_nDof; ++i ){
		tempVector[i] = std::complex<double>(0.0, 0.0);
	}
	for( int i = 0; i < m_nDof; ++i ){
		for( int j = 0; j < m_nDof; ++j ){
			const int ipos = getPositionInMatrix( i, j );
			tempVector[i] += matrix[ipos] * resultVector[j];
		}
	}
	double diffMax(0.0);
	for( int i = 0; i < m_nDof; ++i ){
		const std::complex<double> diff = tempVector[i] - rhsVector[i];
		if( std::abs(diff) > diffMax ){
			diffMax = std::abs(diff);
		}
	}
	std::cout << "diffMax(TM mode) = " << diffMax << std::endl;
#endif

	std::complex<double>* beta = new std::complex<double>[m_num];
	std::complex<double>* gamma = new std::complex<double>[m_num];
	for( int i = 0; i < m_num; ++i ){
		beta[i]  = std::complex<double>( resultVectorLapack[ 2 * i     ].real, resultVectorLapack[ 2 * i     ].imag );
		gamma[i] = std::complex<double>( resultVectorLapack[ 2 * i + 1 ].real, resultVectorLapack[ 2 * i + 1 ].imag );
	}

	const int numObsLoc = static_cast<int>( m_obsLocX.size() );
	std::complex<double>* By1 = new std::complex<double>[numObsLoc];
	std::complex<double>* Ex1 = new std::complex<double>[numObsLoc];
	std::complex<double>* Ez1 = new std::complex<double>[numObsLoc];
	std::complex<double>* Et1 = new std::complex<double>[numObsLoc];
	std::complex<double>* In1 = new std::complex<double>[numObsLoc];
	std::complex<double>* By2 = new std::complex<double>[numObsLoc];
	std::complex<double>* Ex2 = new std::complex<double>[numObsLoc];
	std::complex<double>* Ez2 = new std::complex<double>[numObsLoc];
	std::complex<double>* Et2 = new std::complex<double>[numObsLoc];
	std::complex<double>* In2 = new std::complex<double>[numObsLoc];

#ifdef _TEST
	const int num = 100;
#endif

	for( int iObs = 0; iObs <numObsLoc; ++iObs ){
		const double x = m_obsLocX[iObs];
		const double z = delta * cos(nu * x);
		{// Upper region
			const std::complex<double> theta = sqrt(std::complex<double>( 0.0, omega * mu0 * sigma1 ));
			// By1
			By1[iObs] = B0 * exp( - theta * z );
#ifdef _TEST
			for( int n = 0; n < num; ++n ){
#else
			for( int n = 0; n < m_num; ++n ){
#endif
				const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
				By1[iObs] += beta[n] * exp( theta * z ) * cos( n * nu * x );
			}
			// Ex1
			Ex1[iObs] = B0 * theta / ( mu0 * sigma1 ) * exp( - theta * z );
#ifdef _TEST
			for( int n = 0; n < num; ++n ){
#else
			for( int n = 0; n < m_num; ++n ){
#endif
				const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
				Ex1[iObs] -= theta / ( mu0 * sigma1 ) * beta[n] * exp( theta * z ) * cos( n * nu * x ); 
			}
			// Ez1
			Ez1[iObs] = std::complex<double>(0.0,0.0);
#ifdef _TEST
			for( int n = 0; n < num; ++n ){
#else
			for( int n = 0; n < m_num; ++n ){
#endif
				const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
				Ez1[iObs] -= n * nu * beta[n] * exp( theta * z ) * sin( n * nu * x ) / ( mu0 * sigma1 );
			}
//			// Ez1
//			Ez1[iObs] = theta * delta * nu / ( mu0 * sigma1 ) * sin( nu * x ) * B0 * exp( -theta * z );
//#ifdef _TEST
//			for( int n = 0; n < num; ++n ){
//#else
//			for( int n = 0; n < m_num; ++n ){
//#endif
//				const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
//				Ez1[iObs] -= theta * beta[n] * delta * nu * sin( nu * x ) * exp( theta * z ) * cos( n * nu * x ) / ( mu0 * sigma1 );
//				Ez1[iObs] -= n * nu * beta[n] * exp( theta * z ) * sin( n * nu * x ) / ( mu0 * sigma1 );
//			}
			// Et1
			const double phi = atan2( -nu * delta * sin(nu * x), 1.0 );
			//Et1[iObs] = Ex1[iObs] * cos(phi) + Ez1[iObs] * sin(phi);@todo Œã‚Å–ß‚·
			Et1[iObs] = Ex1[iObs] + Ez1[iObs] * ( -nu * delta * sin(nu * x) );
			// In1
			In1[iObs] = sigma1 * ( Ex1[iObs] * sin(phi) - Ez1[iObs] * cos(phi) );
		}
		{// Lower region
			const std::complex<double> theta = sqrt(std::complex<double>( 0.0, omega * mu0 * sigma2 ));
			// By2
			By2[iObs] = std::complex<double>(0.0,0.0);
#ifdef _TEST
			for( int n = 0; n < num; ++n ){
#else
			for( int n = 0; n < m_num; ++n ){
#endif
				const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
				By2[iObs] += gamma[n] * exp( - theta * z ) * cos( n * nu * x );
			}
			// Ex2
			Ex2[iObs] = std::complex<double>(0.0,0.0);
#ifdef _TEST
			for( int n = 0; n < num; ++n ){
#else
			for( int n = 0; n < m_num; ++n ){
#endif
				const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
				Ex2[iObs] += theta / ( mu0 * sigma2 ) * gamma[n] * exp( - theta * z ) * cos( n * nu * x );
			}
			// Ez2
			Ez2[iObs] = std::complex<double>(0.0,0.0);
			for( int n = 0; n < m_num; ++n ){
				const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
				Ez2[iObs] -= n * nu * gamma[n] * exp( - theta * z ) * sin( n * nu * x ) / ( mu0 * sigma2 );
			}
			//// Ez2
			//Ez2[iObs] = std::complex<double>(0.0,0.0);
			//for( int n = 0; n < num; ++n ){
			//	const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
			//	Ez2[iObs] += theta * gamma[n] * delta * nu * sin( nu * x ) * exp( - theta * z ) * cos( n * nu * x ) / ( mu0 * sigma2 );
			//	Ez2[iObs] -= n * nu * gamma[n] * exp( - theta * z ) * sin( n * nu * x ) / ( mu0 * sigma2 );
			//}
			// Et2
			const double phi = atan2( -nu * delta * sin(nu * x), 1.0 );
			//Et2[iObs] = Ex2[iObs] * cos(phi) + Ez2[iObs] * sin(phi);@todo Œã‚Å–ß‚·
			Et2[iObs] = Ex2[iObs] + Ez2[iObs] * ( -nu * delta * sin(nu * x) );
			// In2
			In2[iObs] = sigma2 * ( Ex2[iObs] * sin(phi) - Ez2[iObs] * cos(phi) );
		}
	}

#ifdef _TEST
	for( int iObs = 0; iObs <numObsLoc; ++iObs ){
		const double x = m_obsLocX[iObs];
		const double z = delta * cos(nu * x);
		const int num = 2;
		// Ez2
		std::complex<double> Ez21 = std::complex<double>(0.0,0.0);
		const double factor = -nu * delta * sin(nu * x);
		for( int n = 0; n < num; ++n ){
			const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
			Ez21 -= n * nu * gamma[n] * exp( - theta * z ) * sin( n * nu * x ) / ( mu0 * sigma2 );
		}
		Ez21 *= factor;
		std::complex<double> Ez22 = std::complex<double>(0.0,0.0);
		for( int n = 0; n < num; ++n ){
			const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
			for( int k = 0; k < m_num; ++k ){
				double factor = 2.0 * static_cast<double>(n);
				if( k == 0 ){
					factor = 1.0 * static_cast<double>(n);
				}
				if( 2 * abs(n + k + 1) + 1 < m_nDof && 2 * abs(n + k - 1) + 1 < m_nDof ){
					int ipos(-1);
					ipos = getPositionInMatrix( 2 * abs(n + k + 1) + 1, 2 * n + 1 );
					Ez22 -= 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n) * gamma[n] * cos( abs(n + k + 1) * nu * x ); 
					ipos = getPositionInMatrix( 2 * abs(n + k - 1) + 1, 2 * n + 1 );
					Ez22 += 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n) * gamma[n] * cos( abs(n + k - 1) * nu * x ); 
				}
				if( 2 * abs(n - k + 1) + 1 < m_nDof && 2 * abs(n - k - 1) + 1 < m_nDof ){
					int ipos(-1);
					ipos = getPositionInMatrix( 2 * abs(n - k + 1) + 1, 2 * n + 1 );
					Ez22 -= 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n) * gamma[n] * cos( abs(n - k + 1) * nu * x ); 
					ipos = getPositionInMatrix( 2 * abs(n - k - 1) + 1, 2 * n + 1 );
					Ez22 += 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n) * gamma[n] * cos( abs(n - k - 1) * nu * x ); 
				}
			}
		}
		std::cout << iObs << " " << Ez21 << " " << Ez22 << std::endl; 

		const int  n = 1;
		const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
		const std::complex<double> exp21 = exp( - theta * z );
		std::complex<double> exp22 = std::complex<double>(0.0,0.0);
		for( int k = 0; k < m_num; ++k ){
			double factor = 2.0;
			if( k == 0 ){
				factor = 1.0;
			}
			if( 2 * abs(n + k + 1) + 1 < m_nDof && 2 * abs(n + k - 1) + 1 < m_nDof ){
				exp22 += 0.25 * pow(-1.0, k) * factor * getItheta2(k, n) * cos( k * nu * x );
				exp22 += 0.25 * pow(-1.0, k) * factor * getItheta2(k, n) * cos( k * nu * x );
			}
			if( 2 * abs(n - k + 1) + 1 < m_nDof && 2 * abs(n - k - 1) + 1 < m_nDof ){
				exp22 += 0.25 * pow(-1.0, k) * factor * getItheta2(k, n) * cos( k * nu * x );
				exp22 += 0.25 * pow(-1.0, k) * factor * getItheta2(k, n) * cos( k * nu * x );
			}
		}
		std::cout << iObs << " " << exp21 << " " << exp22 << std::endl; 
	}

	std::ofstream testFile( "ByEx.txt");
	if( testFile.fail() ){
		std::cerr << "File open error !!" << std::endl;
		exit(1);
	}
	testFile.precision(6);
	for( int iObs = 0; iObs <numObsLoc; ++iObs ){
		const double x = m_obsLocX[iObs];
		const double z = delta * cos(nu * x);
		testFile << std::setw(15) << std::scientific << x;
		testFile << std::setw(15) << std::scientific << z;
		testFile << std::setw(15) << std::scientific << std::abs(By1[iObs]);
		testFile << std::setw(15) << std::scientific << std::abs(Ex1[iObs]);
		testFile << std::setw(15) << std::scientific << std::abs(By2[iObs]);
		testFile << std::setw(15) << std::scientific << std::abs(Ex2[iObs]);
		testFile << std::setw(15) << std::scientific << std::abs(Ez1[iObs]);
		testFile << std::setw(15) << std::scientific << std::abs(Ez2[iObs]);
		testFile << std::setw(15) << std::scientific << std::abs(Et1[iObs]);
		testFile << std::setw(15) << std::scientific << std::abs(Et2[iObs]);
		testFile << std::setw(15) << std::scientific << std::abs(In1[iObs]);
		testFile << std::setw(15) << std::scientific << std::abs(In2[iObs]);
		testFile << std::setw(15) << std::scientific << std::abs(In1[iObs]/sigma1);
		testFile << std::setw(15) << std::scientific << std::abs(In2[iObs]/sigma2);
		testFile << std::endl;
	}
	testFile.close();

	std::ofstream testFile2( "test.txt");
	if( testFile2.fail() ){
		std::cerr << "File open error !!" << std::endl;
		exit(1);
	}
	testFile2.precision(6);
	double xPre = m_obsLocX[0];
	double zPre = delta * cos(nu * xPre);
	std::complex<double> By1Pre = By1[0];
	std::complex<double> By2Pre = By2[0];
	for( int iObs = 1; iObs <numObsLoc; ++iObs ){
		const double x = m_obsLocX[iObs];
		const double z = delta * cos(nu * x);
		const double distance = sqrt( (xPre - x) * (xPre - x) + (zPre - z) * (zPre - z) );
		const std::complex<double> dBy1dt = ( By1[iObs] - By1Pre ) / distance;
		const std::complex<double> dBy2dt = ( By2[iObs] - By2Pre ) / distance;
		const std::complex<double> In1    = dBy1dt / mu0;
		const std::complex<double> In2    = dBy2dt / mu0;
		testFile2 << std::setw(15) << std::scientific << 0.5 * (x + xPre);
		testFile2 << std::setw(15) << std::scientific << 0.5 * (z + zPre);
		testFile2 << std::setw(15) << std::scientific << distance;
		//testFile2 << std::setw(30) << std::scientific << dBy1dt;
		//testFile2 << std::setw(30) << std::scientific << dBy2dt;
		//testFile2 << std::setw(30) << std::scientific << In1;
		//testFile2 << std::setw(30) << std::scientific << In2;
		//testFile2 << std::setw(30) << std::scientific << In1 - In2;
		testFile2 << std::setw(15) << std::scientific << std::abs(In1);
		testFile2 << std::setw(15) << std::scientific << std::abs(In2);
		const double phi = atan2( -nu * delta * sin(nu * x), 1.0 );
		//testFile2 << std::setw(15) << std::scientific << phi * 180 / 3.14;
		testFile2 << std::setw(15) << std::scientific << sin(phi);
		testFile2 << std::setw(15) << std::scientific << cos(phi);
		testFile2 << std::setw(30) << std::scientific << sigma1 * Ex1[iObs] * sin(phi);
		testFile2 << std::setw(30) << std::scientific << sigma1 * Ez1[iObs] * cos(phi);
		testFile2 << std::setw(30) << std::scientific << sigma2 * Ex2[iObs] * sin(phi);
		testFile2 << std::setw(30) << std::scientific << sigma2 * Ez2[iObs] * cos(phi);
		testFile2 << std::endl;
		xPre = x;
		zPre = z;
		By1Pre = By1[iObs];
		By2Pre = By2[iObs];
	}
	testFile2.close();
#endif

	const double rad2deg = 180.0 / M_PI;
	double* app1 = new double[numObsLoc];
	double* phs1 = new double[numObsLoc];
	double* app2 = new double[numObsLoc];
	double* phs2 = new double[numObsLoc];
	for( int iObs = 0; iObs <numObsLoc; ++iObs ){
		const std::complex<double> z1 = Ex1[iObs] / By1[iObs];
		app1[iObs] = std::norm(z1)* mu0 / omega;
		phs1[iObs] = atan2( z1.imag(), z1.real() ) * rad2deg;
		const std::complex<double> z2 = Ex2[iObs] / By2[iObs];
		app2[iObs] = std::norm(z2) * mu0 / omega;
		phs2[iObs] = atan2( z2.imag(), z2.real() ) * rad2deg;
	}

	std::ofstream oFile( "TM.txt");
	if( oFile.fail() ){
		std::cerr << "File open error !!" << std::endl;
		exit(1);
	}
	oFile.precision(6);
	for( int iObs = 0; iObs <numObsLoc; ++iObs ){
		const double x = m_obsLocX[iObs];
		const double z = delta * cos(nu * x);
		oFile << std::setw(15) << std::scientific << x;
		oFile << std::setw(15) << std::scientific << z;
		oFile << std::setw(15) << std::scientific << app1[iObs];
		oFile << std::setw(15) << std::scientific << phs1[iObs];
		oFile << std::setw(15) << std::scientific << app2[iObs];
		oFile << std::setw(15) << std::scientific << phs2[iObs];
		oFile << std::endl;
	}

	oFile.close();

	std::cout << "Finish" << std::endl;

}

void calcMapOfByEx(){

	const double sigma1 = 3.0;
	const double sigma2 = 0.01;
#ifdef _FLAT
	const double delta = 0.0;
#else
	const double delta = 100.0;
#endif
	const double lambda = 1000.0;
	const double f = 1.0 / m_T;
	const double omega = 2.0 * M_PI * f;
	const double nu = 2.0 * M_PI / lambda;
	const double mu0 = 4.0 * M_PI * 1.0e-7;
	const std::complex<double> B0 = std::complex<double>(1.0, 0.0);

	std::complex<double>* matrix = new std::complex<double>[m_nDof*m_nDof];
	std::complex<double>* rhsVector = new std::complex<double>[m_nDof];
	
	// Initialize
	for( int i = 0; i < m_nDof; ++i ){
		rhsVector[i] = std::complex<double>(0.0, 0.0);
	}
	for( int j = 0; j < m_nDof; ++j ){
		const int offset = j * m_num;
		for( int i = 0; i < m_nDof; ++i ){
			matrix[offset+i] = std::complex<double>(0.0, 0.0);
		}
	}

	// Construct matrix and right-hand-side vector
	for( int k = 0; k < m_num; ++k ){
		{// 1st term of By equation
			double factor = 2.0;
			if( k == 0 ){
				factor = 1.0;
			}
			rhsVector[ 2 * k ] -= B0 * pow(-1.0, k) * factor * getItheta1(k,0);
		}
		// 2nd term of By equation
		for( int n = 0; n < m_num; ++n ){
			double factor = 1.0;
			if( k == 0 ){
				factor = 0.5;
			}
			if( 2 * abs(n + k) < m_nDof ){
				const int ipos = getPositionInMatrix( 2 * abs(n + k), 2 * n );
				matrix[ipos] += factor * getItheta1(k, n); 
				if( 2 * abs(n - k) < m_nDof ){
					const int ipos = getPositionInMatrix( 2 * abs(n - k), 2 * n );
					matrix[ipos] += factor * getItheta1(k, n);
				}
			}
		}
		// 3rd term of By equation
		for( int n = 0; n < m_num; ++n ){
			double factor = 1.0;
			if( k == 0 ){
				factor = 0.5;
			}
			if( 2 * abs(n + k) < m_nDof ){
				const int ipos = getPositionInMatrix( 2 * abs(n + k), 2 * n + 1 );
				matrix[ipos] -= pow(-1.0, k) * factor * getItheta2(k, n); 
				if( 2 * abs(n - k) < m_nDof ){
					const int ipos = getPositionInMatrix( 2 * abs(n - k), 2 * n + 1 );
					matrix[ipos] -= pow(-1.0, k) * factor * getItheta2(k, n); 
				}
			}
		}
		{// 1st term of Ex equation
			const std::complex<double> theta = sqrt(std::complex<double>( 0.0, omega * mu0 * sigma1 ));
			double factor = 2.0;
			if( k == 0 ){
				factor = 1.0;
			}
			rhsVector[2 * k + 1] -= theta / ( mu0 * sigma1 ) * B0 * pow(-1.0, k) * factor * getItheta1(k,0);
		}
		// 2nd term of Ex equation
		for( int n = 0; n < m_num; ++n ){
			const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
			double factor = 1.0;
			if( k == 0 ){
				factor = 0.5;
			}
			if( 2 * abs(n + k) + 1 < m_nDof ){
				const int ipos = getPositionInMatrix( 2 * abs(n + k) + 1, 2 * n );
				matrix[ipos] -= theta / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
				if( 2 * abs(n - k) + 1 < m_nDof ){
					const int ipos = getPositionInMatrix( 2 * abs(n - k) + 1, 2 * n );
					matrix[ipos] -= theta / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
				}
			}
		}
		// 3rd term of Ex equation
		for( int n = 0; n < m_num; ++n ){
			const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
			double factor = 1.0;
			if( k == 0 ){
				factor = 0.5;
			}
			if( 2 * abs(n + k) + 1 < m_nDof ){
				const int ipos = getPositionInMatrix( 2 * abs(n + k) + 1, 2 * n + 1 );
				matrix[ipos] -= theta / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
				if( 2 * abs(n - k) + 1 < m_nDof ){
					const int ipos = getPositionInMatrix( 2 * abs(n - k) + 1, 2 * n + 1 );
					matrix[ipos] -= theta / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
				}
			}
		}
		// 3rd term of Ez equation
		for( int n = 0; n < m_num; ++n ){
			const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
			double factor = 2.0 * static_cast<double>(n);
			if( k == 0 ){
				factor = 1.0 * static_cast<double>(n);
			}
			if( 2 * abs(n + k + 1) + 1 < m_nDof && 2 * abs(n + k - 1) + 1 < m_nDof ){
				int ipos(-1);
				ipos = getPositionInMatrix( 2 * abs(n + k + 1) + 1, 2 * n );
				matrix[ipos] -= 0.25 * delta * pow(nu, 2) / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
				ipos = getPositionInMatrix( 2 * abs(n + k - 1) + 1, 2 * n );
				matrix[ipos] += 0.25 * delta * pow(nu, 2) / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
			}
			if( 2 * abs(n - k + 1) + 1 < m_nDof && 2 * abs(n - k - 1) + 1 < m_nDof ){
				int ipos(-1);
				ipos = getPositionInMatrix( 2 * abs(n - k + 1) + 1, 2 * n );
				matrix[ipos] -= 0.25 * delta * pow(nu, 2) / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
				ipos = getPositionInMatrix( 2 * abs(n - k - 1) + 1, 2 * n );
				matrix[ipos] += 0.25 * delta * pow(nu, 2) / ( mu0 * sigma1 ) * factor * getItheta1(k, n); 
			}
		}
		// 5th term of Ez equation
		for( int n = 0; n < m_num; ++n ){
			const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
			double factor = 2.0 * static_cast<double>(n);
			if( k == 0 ){
				factor = 1.0 * static_cast<double>(n);
			}
			if( 2 * abs(n + k + 1) + 1 < m_nDof && 2 * abs(n + k - 1) + 1 < m_nDof ){
				int ipos(-1);
				ipos = getPositionInMatrix( 2 * abs(n + k + 1) + 1, 2 * n + 1 );
				matrix[ipos] += 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
				ipos = getPositionInMatrix( 2 * abs(n + k - 1) + 1, 2 * n + 1 );
				matrix[ipos] -= 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
			}
			if( 2 * abs(n - k + 1) + 1 < m_nDof && 2 * abs(n - k - 1) + 1 < m_nDof ){
				int ipos(-1);
				ipos = getPositionInMatrix( 2 * abs(n - k + 1) + 1, 2 * n + 1 );
				matrix[ipos] += 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
				ipos = getPositionInMatrix( 2 * abs(n - k - 1) + 1, 2 * n + 1 );
				matrix[ipos] -= 0.25 * delta * pow(nu, 2) / ( mu0 * sigma2 ) * pow(-1.0, k) * factor * getItheta2(k, n); 
			}
		}
	}

	lapack_complex_double* matrixLapack = new lapack_complex_double[m_nDof*m_nDof];
	for( int i = 0; i < m_nDof * m_nDof; ++i ){
		matrixLapack[i].real = matrix[i].real();
		matrixLapack[i].imag = matrix[i].imag();
	}
	lapack_complex_double* resultVectorLapack = new lapack_complex_double[m_nDof];
	for( int i = 0; i < m_nDof; ++i ){
		resultVectorLapack[i].real = rhsVector[i].real();
		resultVectorLapack[i].imag = rhsVector[i].imag();
	}

	int ierr(0);
	int* ipiv = new int[m_nDof];
	int lda = m_nDof;
	int ldb = m_nDof;

	// Factorize
	ierr = LAPACKE_zgetrf( LAPACK_COL_MAJOR, m_nDof, m_nDof, matrixLapack, lda, ipiv );
	if( ierr > 0 ) {
		std::cerr << "Error : Matrix is singular. ierr = " << ierr << std::endl;
		exit(1);
	}else if( ierr < 0 ){
		std::cerr << "Error : " << -ierr << "-th parameter has illegal value." << std::endl;
		exit(1);
	}

	// Solve
	ierr = LAPACKE_zgetrs( LAPACK_COL_MAJOR, 'N', m_nDof, 1, matrixLapack, lda, ipiv, resultVectorLapack, ldb );
	if( ierr < 0 ){
		std::cerr << "Error : " << -ierr << "-th parameter has illegal value." << std::endl;
		exit(1);
	}

	std::complex<double>* beta = new std::complex<double>[m_num];
	std::complex<double>* gamma = new std::complex<double>[m_num];
	for( int i = 0; i < m_num; ++i ){
		beta[i]  = std::complex<double>( resultVectorLapack[ 2 * i     ].real, resultVectorLapack[ 2 * i     ].imag );
		gamma[i] = std::complex<double>( resultVectorLapack[ 2 * i + 1 ].real, resultVectorLapack[ 2 * i + 1 ].imag );
	}

	const int num = 100;

	const int nx = 100;
	const double startX = -1000.0;
	const double endX = 1000.0;
	const int nz = 100;
	const double startZ = -500.0;
	const double endZ = 500.0;

	const double dx = ( endX - startX ) / static_cast<double>(nx);
	const double dz = ( endZ - startZ ) / static_cast<double>(nz);

	std::complex<double>* By1 = new std::complex<double>[ (nx + 1) * (nz + 1) ];
	std::complex<double>* Ex1 = new std::complex<double>[ (nx + 1) * (nz + 1) ];
	std::complex<double>* Ez1 = new std::complex<double>[ (nx + 1) * (nz + 1) ];
	std::complex<double>* Et1 = new std::complex<double>[ (nx + 1) * (nz + 1) ];
	std::complex<double>* In1 = new std::complex<double>[ (nx + 1) * (nz + 1) ];
	std::complex<double>* By2 = new std::complex<double>[ (nx + 1) * (nz + 1) ];
	std::complex<double>* Ex2 = new std::complex<double>[ (nx + 1) * (nz + 1) ];
	std::complex<double>* Ez2 = new std::complex<double>[ (nx + 1) * (nz + 1) ];
	std::complex<double>* Et2 = new std::complex<double>[ (nx + 1) * (nz + 1) ];
	std::complex<double>* In2 = new std::complex<double>[ (nx + 1) * (nz + 1) ];

	for( int ix = 0; ix < nx + 1; ++ix ){
		const double x = startX + dx * static_cast<double>(ix);
		for( int iz = 0; iz < nz + 1; ++iz ){
			const double z = startZ + dz * static_cast<double>(iz);
			const int index = (nz + 1) * ix + iz;
			{// Upper region
				const std::complex<double> theta = sqrt(std::complex<double>( 0.0, omega * mu0 * sigma1 ));
				// By1
				By1[index] = B0 * exp( - theta * z );
				for( int n = 0; n < num; ++n ){
					const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
					By1[index] += beta[n] * exp( theta * z ) * cos( n * nu * x );
				}
				// Ex1
				Ex1[index] = B0 * theta / ( mu0 * sigma1 ) * exp( - theta * z );
				for( int n = 0; n < num; ++n ){
					const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
					Ex1[index] -= theta / ( mu0 * sigma1 ) * beta[n] * exp( theta * z ) * cos( n * nu * x ); 
				}
				// Ez1
				Ez1[index] = std::complex<double>(0.0,0.0);
				for( int n = 0; n < num; ++n ){
					const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
					Ez1[index] -= n * nu * beta[n] * exp( theta * z ) * sin( n * nu * x ) / ( mu0 * sigma1 );
				}
				//// Ez1
				//Ez1[index] = theta * delta * nu / ( mu0 * sigma1 ) * sin( nu * x ) * B0 * exp( -theta * z );
				//for( int n = 0; n < num; ++n ){
				//	const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma1 ));
				//	Ez1[index] -= theta * beta[n] * delta * nu * sin( nu * x ) * exp( theta * z ) * cos( n * nu * x ) / ( mu0 * sigma1 );
				//	Ez1[index] -= n * nu * beta[n] * exp( theta * z ) * sin( n * nu * x ) / ( mu0 * sigma1 );
				//}
				// Et1
				const double phi = atan2( -nu * delta * sin(nu * x), 1.0 );
				Et1[index] = Ex1[index] * cos(phi) + Ez1[index] * sin(phi);
				// In1
				In1[index] = sigma1 * ( Ex1[index] * sin(phi) - Ez1[index] * cos(phi) );
			}
			{// Lower region
				const std::complex<double> theta = sqrt(std::complex<double>( 0.0, omega * mu0 * sigma2 ));
				// By2
				By2[index] = std::complex<double>(0.0,0.0);
				for( int n = 0; n < num; ++n ){
					const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
					By2[index] += gamma[n] * exp( - theta * z ) * cos( n * nu * x );
				}
				// Ex2
				Ex2[index] = std::complex<double>(0.0,0.0);
				for( int n = 0; n < num; ++n ){
					const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
					Ex2[index] += theta / ( mu0 * sigma2 ) * gamma[n] * exp( - theta * z ) * cos( n * nu * x );
				}
				// Ez2
				Ez2[index] = std::complex<double>(0.0,0.0);
				for( int n = 0; n < num; ++n ){
					const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
					Ez2[index] -= n * nu * gamma[n] * exp( - theta * z ) * sin( n * nu * x ) / ( mu0 * sigma2 );
				}
				//// Ez2
				//Ez2[index] = std::complex<double>(0.0,0.0);
				//for( int n = 0; n < num; ++n ){
				//	const std::complex<double> theta = sqrt(std::complex<double>( pow( n * nu, 2 ), omega * mu0 * sigma2 ));
				//	Ez2[index] += theta * gamma[n] * delta * nu * sin( nu * x ) * exp( - theta * z ) * cos( n * nu * x ) / ( mu0 * sigma2 );
				//	Ez2[index] -= n * nu * gamma[n] * exp( - theta * z ) * sin( n * nu * x ) / ( mu0 * sigma2 );
				//}
				// Et2
				const double phi = atan2( -nu * delta * sin(nu * x), 1.0 );
				Et2[index] = Ex2[index] * cos(phi) + Ez2[index] * sin(phi);
				// In2
				In2[index] = sigma2 * ( Ex2[index] * sin(phi) - Ez2[index] * cos(phi) );
			}
		}
	}
	
	std::ofstream oFile( "MapOfByEx.txt");
	if( oFile.fail() ){
		std::cerr << "File open error !!" << std::endl;
		exit(1);
	}
	oFile.precision(6);
	for( int ix = 0; ix < nx + 1; ++ix ){
		const double x = startX + dx * static_cast<double>(ix);
		for( int iz = 0; iz < nz + 1; ++iz ){
			const double z = startZ + dz * static_cast<double>(iz);
			const int index = (nz + 1) * ix + iz;
			oFile << std::setw(15) << std::scientific << x;
			oFile << std::setw(15) << std::scientific << z;
			oFile << std::setw(15) << std::scientific << std::abs(By1[index]);
			oFile << std::setw(15) << std::scientific << std::abs(Ex1[index]);
			oFile << std::setw(15) << std::scientific << std::abs(By2[index]);
			oFile << std::setw(15) << std::scientific << std::abs(Ex2[index]);
			oFile << std::setw(15) << std::scientific << std::abs(Ez1[index]);
			oFile << std::setw(15) << std::scientific << std::abs(Ez2[index]);
			oFile << std::setw(15) << std::scientific << std::abs(Et1[index]);
			oFile << std::setw(15) << std::scientific << std::abs(Et2[index]);
			oFile << std::setw(15) << std::scientific << std::abs(In1[index]);
			oFile << std::setw(15) << std::scientific << std::abs(In2[index]);
			oFile << std::setw(15) << std::scientific << std::abs(In1[index]/sigma1);
			oFile << std::setw(15) << std::scientific << std::abs(In2[index]/sigma2);
			oFile << std::endl;
		}
	}
	oFile.close();

#ifdef _TEST
	std::ofstream testFile( "testEx.txt");
	if( testFile.fail() ){
		std::cerr << "File open error !!" << std::endl;
		exit(1);
	}
	testFile.precision(6);
	for( int ix = 0; ix < nx + 1; ++ix ){
		const double x = startX + dx * static_cast<double>(ix);
		const int iz = 0;
		double zPre = startZ + dz * static_cast<double>(iz);
		const int index = (nz + 1) * ix + iz;
		std::complex<double> By1Pre = By1[index];
		std::complex<double> By2Pre = By2[index];
		for( int iz = 1; iz < nz + 1; ++iz ){
			const double z = startZ + dz * static_cast<double>(iz);
			const int index = (nz + 1) * ix + iz;
			const double distance = z - zPre;
			const std::complex<double> dBy1dz = ( By1[index] - By1Pre ) / distance;
			const std::complex<double> dBy2dz = ( By2[index] - By2Pre ) / distance;
			const std::complex<double> Ex1    = - dBy1dz / ( sigma1 * mu0 );
			const std::complex<double> Ex2    = - dBy2dz / ( sigma2 * mu0 );
			testFile << std::setw(15) << std::scientific << x;
			testFile << std::setw(15) << std::scientific << 0.5 * (z + zPre);
			testFile << std::setw(15) << std::scientific << distance;
			testFile << std::setw(15) << std::scientific << std::abs(Ex1);
			testFile << std::setw(15) << std::scientific << std::abs(Ex2);
			testFile << std::endl;
			zPre = z;
			By1Pre = By1[index];
			By2Pre = By2[index];
		}
	}
	testFile.close();

	std::ofstream testFile2( "testEz.txt");
	if( testFile2.fail() ){
		std::cerr << "File open error !!" << std::endl;
		exit(1);
	}
	testFile2.precision(6);
	for( int iz = 0; iz < nz + 1; ++iz ){
		const double z = startZ + dz * static_cast<double>(iz);
		const int ix = 0;
		double xPre = startX + dx * static_cast<double>(ix);
		const int index = (nz + 1) * ix + iz;
		std::complex<double> By1Pre = By1[index];
		std::complex<double> By2Pre = By2[index];
		for( int ix = 1; ix < nx + 1; ++ix ){
			const double x = startX + dx * static_cast<double>(ix);
			const int index = (nz + 1) * ix + iz;
			const double distance = x - xPre;
			const std::complex<double> dBy1dx = ( By1[index] - By1Pre ) / distance;
			const std::complex<double> dBy2dx = ( By2[index] - By2Pre ) / distance;
			const std::complex<double> Ez1    = dBy1dx / ( sigma1 * mu0 );
			const std::complex<double> Ez2    = dBy2dx / ( sigma2 * mu0 );
			testFile2 << std::setw(15) << std::scientific << 0.5 * (x + xPre);;
			testFile2 << std::setw(15) << std::scientific << z;
			testFile2 << std::setw(15) << std::scientific << distance;
			testFile2 << std::setw(15) << std::scientific << std::abs(Ez1);
			testFile2 << std::setw(15) << std::scientific << std::abs(Ez2);
			testFile2 << std::endl;
			xPre = x;
			By1Pre = By1[index];
			By2Pre = By2[index];
		}
	}
	testFile2.close();
#endif
	//const double rad2deg = 180.0 / M_PI;
	//double* app1 = new double[numObsLoc];
	//double* phs1 = new double[numObsLoc];
	//double* app2 = new double[numObsLoc];
	//double* phs2 = new double[numObsLoc];
	//for( int iObs = 0; iObs <numObsLoc; ++iObs ){
	//	const std::complex<double> z1 = Ex1[iObs] / By1[iObs];
	//	app1[iObs] = std::norm(z1)* mu0 / omega;
	//	phs1[iObs] = atan2( z1.imag(), z1.real() ) * rad2deg;
	//	const std::complex<double> z2 = Ex2[iObs] / By2[iObs];
	//	app2[iObs] = std::norm(z2) * mu0 / omega;
	//	phs2[iObs] = atan2( z2.imag(), z2.real() ) * rad2deg;
	//}

	//std::ofstream oFile( "TM.txt");
	//if( oFile.fail() ){
	//	std::cerr << "File open error !!" << std::endl;
	//	exit(1);
	//}
	//oFile.precision(6);
	//for( int iObs = 0; iObs <numObsLoc; ++iObs ){
	//	const double x = m_obsLocX[iObs];
	//	const double z = delta * cos(nu * x);
	//	oFile << std::setw(15) << std::scientific << x;
	//	oFile << std::setw(15) << std::scientific << z;
	//	oFile << std::setw(15) << std::scientific << app1[iObs];
	//	oFile << std::setw(15) << std::scientific << phs1[iObs];
	//	oFile << std::setw(15) << std::scientific << app2[iObs];
	//	oFile << std::setw(15) << std::scientific << phs2[iObs];
	//	oFile << std::endl;
	//}

	//oFile.close();

	//std::cout << "Finish" << std::endl;

}