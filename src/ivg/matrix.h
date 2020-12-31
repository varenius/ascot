#ifndef MATRIX_H
#define MATRIX_H

/**
*
* @brief Matrix class
* @author bakkari developer team
* @date 2015-03-24
* @version 0.1
*/

// ===========================================================================
// Matrix-Klasse
// -------------
// 2009-11-30 - 	erste Schritte: Uebung 4 NiC++
// 2009-12-01 - 	Kopierkonstruktor, Zuweisungsoperator, (), +, -, * Operator
// 2009-12-14 - 	Operator += (Skalar) += und -= (Matrix)
// 2009-12-21 - 	Trennung in Header und Quelltext + Aufraeumen
// 2009-12-22 - 	Const-Korrektheit
// 2010-01-11 - 	Umstellung auf BLAS Routinen
//
// 2010-12-03 JP - 	neue Methoden fuer Index des min- & max-Elements
// 2010-12-07 JP - 	neue Methode: numerisch sortieren nach beliebiger Spalte
// 2010-12-23 JP & TA - neue Methode: Singulaerwertzerlegung einer allg.
//                      Matrix ( benoetigt "lapack_wrapper.h" )
// 2010-12-23 JP - 	Ueberladung der Methode get_sub(): Aufruf mit 4 Integern
//                 	(Koord. linke obere Ecke, Koord. rechte untere Ecke)
// 2011-01-05 JP - 	neue Methode: save_matrix() erzeugt eine Ascii-Datei mit
//                 	den Daten im Matrix-Format (nuetzlich z.B. fuer image-plot)
// 2011-01-17 TA - 	Ueberladung der Methode append_rows(): double-Wert an Vektor
// 2011-05-10 TA - 	Operator ^ und pow fuer elementweises "potenzieren"
// 2011-05-17 TA - 	inv repariert und Test in chol_inv, ob Matrix symmetrisch
// 2011-05-31 JP - 	svd: Option die Zeilenanzahl der Matrix anzugeben
// 2011-05-31 TA & JP - rem_r: Option mehrere Zeilen gleichzeitig zu loeschen
// 2012-05-18 ME - 	spezielle Matrix Multiplikationen ergaenzt
// 2013-06-18 SH - 	new methods: rand_norm, norm
// 2013-06-18 SH - 	comment source code / clean-up matrix class
// 2013-06-18 SH -	new reference documentation (as html-file; using doxygen)
// 2013-08-XX SH - 	new method: estimate_cpwlf
// 2013-09-09 SH - 	new methods: min( int &idx ), max( int &idx )
// 2014-12-08 SH -  new namespace ivg for bakkari software
// ===========================================================================

#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <limits>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <numeric>
#include <cmath>
#include <cstdarg>
#include <boost/bind.hpp>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <boost/random/normal_distribution.hpp>
#include "auxfunc.h"
#include <matio.h>
#include "gsl/gsl_errno.h"
#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"

extern "C"
{
#include "cblas.h"
}
#include "lapack_wrapper.h"
#include "funother.h"

using namespace std;

// read sparse matrix from MATLAB binaries
extern "C"
{
    void sparse2array( matvar_t *matvar, double *out );
};


namespace ivg
{

// ===========================================================================
class Matrix
    // ===========================================================================
{
    public:
        // ==============================================
        // =============== Constructors: ================
        // ==============================================
        /**
         *  \b Description: \n
         *        Default constructor
         *  \param [in] no input parameters needed
         *  \return An instance of the class 'Matrix'
         */
        Matrix( );

        /**
         *  \b Description: \n
         *        Constructor using three parameters: number of rows and colomns respectively, and a value
         *  \param [in] [int] number of rows; number of columns
         *               [double] value
         *  \return An instance of the class 'Matrix'
         */
        Matrix( int r, int c, double v = numeric_limits<double>::quiet_NaN() );

        /**
         *  \b Description: \n
         *        Constructor using the number of rows and columns respectively, and two iterators
         *  \param [in] [int] number of rows; number of columns
         *               [vector<double>::iterator]
         *  \return An instance of the class 'Matrix'
         */
        Matrix( vector<double>::iterator begin, vector<double>::iterator end,
                int r, int c );

        /**
         *  \b Description: \n
         *        Copy constructor
         *  \param [in] [Matrix] other matrix
         *  \return An instance of the class 'Matrix'
         */
        Matrix( const Matrix &other );

        /**
         *  \b Description: \n
         *        Constructor using std::vectors for (row,column)-index and values; the dimension of the matrix is [r, c]
         *  \param [in] [vector<int>] vectors with row and column index
         *               [vector<double>] vector with values
         *               [int] r: number of rows; c: number of columns
         *               [double] v: value; if not selected, v = 0
         *  \return An instance of the class 'Matrix'
         */
        Matrix(const vector<int> &rows,const vector<int> &cols,
               const vector<double> &values,
               int r, int c, double v=0);

        /**
         *  \b Description: \n
         *        Constructor for a [n x 1] matrix with regulary increments
         *  \param [in] [double] start: start value; end: end value
         *               [double] schrittweite: increment
         *               [int] cols: number of columns
         *  \return An instance of the class 'Matrix'
         */
        Matrix(double start, double schrittweite, double ende, int cols);

        /**
         *  \b Description: \n
         *        Constructor for a [n x 1] matrix with double-values
         *  \param [in] [vector<double>] vector that needed to convert to a matrix
         *  \return An instance of the class 'Matrix'
         */
        Matrix(const vector<double> &vecvec);

        /**
         *  \b Description: \n
         *        Constructor for a [n x 1] matrix with integer-values
         *  \param [in] [vector<int>] vector that needed to convert to a matrix
         *  \return An instance of the class 'Matrix'
         */
        Matrix(const vector<int> &vecvec);

        Matrix(const vector<Matrix> &vMat) ;

        // Aus string Matrix machen: "[1,2;3,4]"
        //Matrix(const string &strIn);

        //Matrix(int amount, ...);


        // ==============================================
        // =============== Destruktoren: ================
        // ==============================================
        /**
         *  \b Description: \n
         *        Default deconstructor
         *  \param [in] no input parameters needed
         */
        ~Matrix( );



        // ==============================================
        // ================ Operators: ==================
        // ==============================================

        // ***** assignment operators (Zuweisungsoperator) *****
        /**
         *  \b Description: \n
         *        assignment operator '=' using another matrix
         *  \param [in] [Matrix] &other: other matrix
         */
        Matrix & operator=( const Matrix &other ); 		// ... using a matrix

        /**
         *  \b Description: \n
         *        assignment operator '=' using a std::vector
         *  \param [in] [vector<int>] &vec: std::vector
         */
        Matrix & operator=(const vector<int> &vec);	// ...using a std::vector



        // ***** operators for multiplication *****
        /**
         *  \b Description: \n
         *        operator for multiplying a scalar value '*'
         *  \param [in] [double] v: value
         */
        Matrix   operator*(const double &v ) const ;	// multiply with scalar value

        /**
         *  \b Description: \n
         *        operator for multiplying a scalar value + assignment to this-matrix '*='
         *  \param [in] [double] v: value
         */
        Matrix & operator*=( double
                             v );           	// multiply with scalar value + assignment

        /**
         *  \b Description: \n
         *        operator for a matrix multiplication '*'
         *  \param [in] [Matrix] m2: other matrix
         */
        Matrix   operator*( const Matrix &m2 ) const;  // matrix multiplication

        /**
         *  \b Description: \n
         *        operator for a matrix multiplication + assignment to this-matrix '*='
         *  \param [in] [Matrix] m2: other matrix
         */
        Matrix & operator*=( const Matrix & m2 ); // matrix mmultiplication + assignment



        // ***** operators for division *****
        /**
         *  \b Description: \n
         *        operator for dividing by a scalar value '/'
         *  \param [in] [double] v: value
         */
        Matrix operator/( double v );            		// Division by scalar value

        /**
         *  \b Description: \n
         *        operator for an element by element devision in a matrix '/'
         *  \param [in] [Matrix] m2: other matrix
         */
        Matrix operator/(const Matrix &m2 );			// element-by-element division



        // ***** operators for addition *****
        /**
         *  \b Description: \n
         *        operator for adding a scalar value to each element in a matrix '+'
         *  \param [in] [double] v: value
         */
        Matrix operator+( double v )
        const;      		// add scalar value to each element in a matrix

        /**
         *  \b Description: \n
         *        operator for adding a scalar value to each element in a matrix + assignment to this-matrix '+='
         *  \param [in] [double] v: value
         */
        Matrix operator+=( double
                           v );					// add scalar value to each element in a matrix + assignment

        /**
         *  \b Description: \n
         *        operator for adding a matrix '+'
         *  \param [in] [Matrix] m2: other matrix
         */
        Matrix   operator+( const Matrix &m2 ) const;  // add matrix

        /**
         *  \b Description: \n
         *        operator for adding a matrix + assignment to this-matrix '+='
         *  \param [in] [Matrix] m2: other matrix
         */
        Matrix & operator+=( const Matrix &m2 );   		// add matrix + assignment



        // ***** operators for subtraction *****
        /**
         *  \b Description: \n
         *        operator to subtract a scalar value from each element in a matrix '-'
         *  \param [in] [double] v: value
         */
        Matrix operator-( double v )
        const;      		// subtract a scalar value from each element of a matrix

        /**
         *  \b Description: \n
         *        operator to subtract a scalar value from each element in a matrix + assignment to this-matrix '-='
         *  \param [in] [double] v: value
         */
        Matrix operator-=( double
                           v );					// subtract a scalar value from each element of a matrix + assignment

        /**
         *  \b Description: \n
         *        operator to subtract a matrix '-'
         *  \param [in] [Matrix] m2: other matrix
         */
        Matrix operator-( const Matrix &m2 ) const;    // subtract matrix

        /**
         *  \b Description: \n
         *        operator to subtract a matrix + assignment to this-matrix '-='
         *  \param [in] [Matrix] m2: other matrix
         */
        Matrix & operator-=( const Matrix &m2 );   		// subtract matrix + assignment



        // ***** operators for exponentiation *****
        /**
         *  \b Description: \n
         *        operator for an element-by-element exponentiation in a matrix '^'
         *  \param [in] [double] v: value
         */
        Matrix operator^( double v ) const;      // element-by-element exponentiation

        /**
         *  \b Description: \n
         *        operator for an element-by-element exponentiation in a matrix '^'
         *  \param [in] [int] v: value
         */
        Matrix operator^( int v ) const;         // element-by-element exponentiation



        // ***** operators using brackets ******
        /**
         *  \b Description: \n
         *        operator for the use of brackets '()'
         *        read and write access
         *  \param [in] [int] i: row; j: column
         *  \return [double] element of the (i,j)th position
         */
        double & operator()( int i,
                             int j );      // read and write access using row/column

        /**
         *  \b Description: \n
         *        operator for the use of brackets '()'
         *        read and write access; linear adress
         *  \param [in] [int] i: row
         *  \return [double] element of the (i,i)th position
         */
        double & operator()( int
                             i );              // read and write access using a linear adress

        /**
         *  \b Description: \n
         *        operator for the use of brackets '()'
         *        read-only access
         *  \param [in] [int] i: row; j: column
         *  \return [double] element of the (i,j)th position
         */
        double operator()( int i, int j ) const;  // read-only access using row/column

        /**
         *  \b Description: \n
         *        operator for the use of brackets '()'
         *        read-only access; linear adress
         *  \param [in] [int] i: row
         *  \return [double] element of the (i,i)th position
         */
        double operator()( int i )
        const;        	// read-only access using a linear adress

        /**
         *  \b Description: \n
         *        operator for the use of brackets '()'
         *        read-only access; return of the ith column of a matrix
         *  \param [in] [char] c: use of ':' for the whole column
         *               [int] i: row
         *  \return [Matrix] ith column of a matrix
         */
        Matrix operator()( char c, int i ) const;			// selection of a column

        
        /**
         *  \b Description: \n
         *        operator for the use of brackets '()'
         *        read-only access; return of several columns of a matrix
         *  \param [in] [char] c: use of ':' for the whole column
         *              [std:vector<int>] idx: rows
         *  \return [Matrix] columns of a matrix
         */
        Matrix operator()( char c, std::vector<int> idx ) const;			// selection of a several columns

      
        /**
         *  \b Description: \n
         *        operator for the use of brackets '()'
         *        read-only access; return of several columns of a matrix
         *  \param [in] [char] c: use of ":" for the whole column
         *              [std:vector<int>] idx: rows
         *  \return [Matrix] columns of a matrix
         */
        Matrix operator()( std::string c, std::vector<int> idx ) const;			// selection of a several columns
        
        /**
         *  \b Description: \n
         *        operator for the use of brackets '()'
         *        read-only access; return of the ith column of a matrix
         *  \param [in] [string] c: use of ":" for the whole column
         *               [int] i: row
         *  \return [Matrix] ith column of a matrix
         */
        Matrix operator()(const string c, int i ) const; 	// selection of a column

        /**
         *  \b Description: \n
         *        operator for the use of brackets '()':
         *        get sub-matrix using two std::vectors
         *  \param [in] [vector<int>] &idx_r: index vector for rows
         *               [vector<int>] &idx_c: index vector for columns
         *  \return [Matrix] sub-matrix
         */
        Matrix operator()(const vector<int> &idx_r,
                          const vector<int> &idx_c) const;		// get sub-matrix using index vectors

        /**
         *  \b Description: \n
         *        operator for the use of brackets '()':
         *        get sub-matrix using start and end tags for the row and the column of the matrix
         *  \param [in] [int] startR: start tag (row); endR: end tag (row)
         *               [int] startC: start tag (column); endC: end tag (column)
         *  \return [Matrix] sub-matrix
         */
        Matrix operator()(int startR, int endR, int startC ,
                          int endC ) const;			// get sub-matrix using start and end tags


        Matrix operator()(const vector<int> &idx) const;
        Matrix operator()(const Matrix &idx) const;
        Matrix operator()(const Matrix &idx_r, const string c) const;
        Matrix operator()(const vector<int> &idx_r, const string c) const;
        //operator double() const;
        // schreibender (":",i) operator ???
        //Matrix & operator()( char c, int i );




        // ==============================================
        // =========== MEMBER-functions: ================
        // ==============================================

        // ***** resize matrix dimension *****
        /**
         *  \b Description: \n
         *        Method to change the matrix dimension
         *  \param [in] [int] &i: number of rows; &j: number of columns
         *               [double] v: value
         *  \return [double] element of the (i,i)th position
         */
        void resize(const int &i,const int &j,
                    const double &v =numeric_limits<double>::quiet_NaN());

        // ***** get dimension  *****
        /**
         *  \b Description: \n
         *        Method to get the number of rows
         *  \param [in] no input parameters needed
         *  \return [int] number of rows
         */
        int rows() const          // get dimension (#rows)
        {
            return _rows;
        };

        /**
         *  \b Description: \n
         *        Method to get the number of columns
         *  \param [in] no input parameters needed
         *  \return [int] number of columns
         */
        int cols() const          // get dimension (#columns)
        {
            return _cols;
        };

        /**
         *  \b Description: \n
         *        Method to get the maximal dimension (rows or columns)
         *  \param [in] no input parameters needed
         *  \return [int] max. dimension
         */
        int length() const;      // get the maximal dimension (rows or columns)

        /**
         *  \b Description: \n
         *        Method to get the size (rows or columns) of a matrix
         *  \param [in] [int] dim: dim = 1 >>> rows; dim = 2 >>> columns
         *  \return [int] size of a matrix
         */
        int size(int dim) const;

        /**
         *  \b Description: \n
         *        Method to get the size of a matrix
         *  \param [in] no input parameters needed
         *  \return [Matrix] [2 x 1] Matrix with the number of rows and columns respectively
         */
        Matrix size() const;

        /**
         *  \b Description: \n
         *        Method to get the number of elements of a matrix
         *  \param [in] no input parameters needed
         *  \return [int] number of elements
         */
        int numel() const  		// number of elements
        {
            return _rows*_cols;
        }



        // ***** mathematical fundamentals *****
        /**
         *  \b Description: \n
         *        Method to get the minimal element
         *  \param [in] no input parameters needed
         *  \return [double] min. element in matrix
         */
        double min() const;       // get the minimum of the matrix

        /**
         *  \b Description: \n
         *        Method to get the maximal element
         *  \param [in] no input parameters needed
         *  \return [double] max. element in matrix
         */
        double max() const;       // get the maximum of the matrix

        /**
         *  \b Description: \n
         *        Method to get the minimal element and its index
         *  \param [in] no input parameters needed
         *  \return [double] min. element in matrix
         */
        double min( int & idx );       // get the minimum of the matrix

        /**
         *  \b Description: \n
         *        Method to get the maximal element and its index
         *  \param [in] no input parameters needed
         *  \return [double] max. element in matrix
         */
        double max( int & idx );       // get the maximum of the matrix

        /**
         *  \b Description: \n
         *        Method to get the index of the minimal element
         *  \param [in] no input parameters needed
         *  \return [double] index of the min. element in matrix
         */
        int minIdx() const;		// get index of the minimum of the matrix

        /**
         *  \b Description: \n
         *        Method to get the index of the maximal element
         *  \param [in] no input parameters needed
         *  \return [double] index of the max. element in matrix
         */
        int maxIdx() const;		// get index of the maximum of the matrix

        // element-by-element trigonometrical functions
        /**
         *  \b Description: \n
         *        Method to calculate the sinus of each element
         *  \param [in] no input parameters needed
         *  \return [Matrix] new matrix concerning the sinus of each element of the this-matrix
         */
        Matrix sin() const;					// element-by-element sinus

        /**
         *  \b Description: \n
         *        Method to calculate the sinus of each element (as reference)
         *  \param [in] no input parameters needed
         */
        void sin(Matrix &m2) const;		// element-by-element sinus (reference)

        /**
         *  \b Description: \n
         *        Method to calculate the cosine of each element
         *  \param [in] no input parameters needed
         *  \return [Matrix] new matrix concerning the cosine of each element of the this-matrix
         */
        Matrix cos() const;					// element-by-element cosine

        /**
         *  \b Description: \n
         *        Method to calculate the cosine of each element (as reference)
         *  \param [in] no input parameters needed
         */
        void cos(Matrix &m2) const;		// element-by-element cosine (reference)

        /**
         *  \b Description: \n
         *        Method to calculate the tangent of each element
         *  \param [in] no input parameters needed
         *  \return [Matrix] new matrix concerning the tangent of each element of the this-matrix
         */
        Matrix tan() const;					// element-by-element tangent

        /**
         *  \b Description: \n
         *        Method to calculate the tangent of each element (as reference)
         *  \param [in] no input parameters needed
         */
        void tan(Matrix &m2) const;		// element-by-element tangent (reference)


        /**
         *  \b Description: \n
         *        Method to calculate the arc tangent of each element (as reference)
         *        this = atan2( m1,m2 )
         *  \param [in] no input parameters needed
         */
        void atan2(const Matrix &m1, const Matrix &m2);  // arc tangent (reference)

        /**
         *  \b Description: \n
         *        Method to calculate the arc cosine of each element
         *  \param [in] no input parameters needed
         *  \return [Matrix] new matrix concerning the arc cosine of each element of the this-matrix
         */
        Matrix acos() const;					// element-by-element arc cosine

        /**
         *  \b Description: \n
         *        Method to calculate the arc cosine of each element (as reference)
         *  \param [in] no input parameters needed
         */
        void acos(Matrix &m2) const;		// element-by-element arc cosine (reference)
        
        /**
         *  \b Description: \n
         *        Method to calculate the square root of each element
         *  \param [in] no input parameters needed
         *  \return [Matrix] new matrix concerning the square root of each element of the this-matrix
         */
        Matrix sqrt() const;				// element-by-element square root

        /**
         *  \b Description: \n
         *        Method to calculate the square root of each element (as reference)
         *  \param [in] no input parameters needed
         */
        void sqrt(Matrix &m2) const;		// element-by-element square root (reference)

        /**
         *  \b Description: \n
         *        Method to calculate the exponential function of each element
         *  \param [in] no input parameters needed
         *  \return [Matrix] new matrix concerning the exp of each element of the this-matrix
         */
        Matrix exp() const;					// element-by-element exp

        /**
         *  \b Description: \n
         *        Method to calculate the exponential function of each element (as reference)
         *  \param [in] no input parameters needed
         */
        void exp( Matrix &m2) const;		// element-by-element exp (reference)

        /**
         *  \b Description: \n
         *        Method to calculate the natural logarithm of each element
         *  \param [in] no input parameters needed
         *  \return [Matrix] new matrix concerning the exp of each element of the this-matrix
         */
        Matrix log() const;

        /**
         *  \b Description: \n
         *        Method to calculate the natural logarithm of each element (as reference)
         *  \param [in] no input parameters needed
         */
        void log( Matrix &m2) const;       
        
        
        /**
         *  \b Description: \n
         *        Method to get the modulo after division (MATLAB version, as a reference)
         *  \param [in] [double] &v: value
         */
        void mod(Matrix &m2, const double v)
        const;	// modulo after division (MATLAB version, reference)

        /**
         *  \b Description: \n
         *        Method to get the modulo after division (MATLAB version)
         *  \param [in] [double] &v: value
         *  \return [Matrix] new matrix concerning the modulo of each element of the this-matrix after division by value
         */
        Matrix mod(const double &v) const;			// modulo after division (MATLAB version)

        /**
         *  \b Description: \n
         *        Method to get the modulo after division (as a reference)
         *  \param [in] [double] &v: value
         */
        void fmod(Matrix &m2, const double v)
        const;	// modulo after division (reference)

        /**
         *  \b Description: \n
         *        Method to get the modulo after division
         *  \param [in] [double] &v: value
         *  \return [Matrix] new matrix concerning the modulo of each element of the this-matrix after division by value
         */
        Matrix fmod(const double &v) const;			// modulo after division
        
        
         /**
         *  \b Description: \n
         *        Method to get fractional and integral part for every element.
         *  \param [in] [Matrix] &fractpart: matrix containing the fractional part
         *              [Matrix] &fintpart : matrix containing the integral part 
         */
        void modf(Matrix &fractpart, Matrix &intpart) const;
        
         /**
         *  \b Description: \n
         *        Method to get fractional and integral part for every element.
         *  \param [in] [Matrix] &fintpart : matrix containing the integral part 
         * \return [Matrix] &fractpart: matrix containing the fractional part
         */
        Matrix modf(  Matrix & intpart) const;
        

        /**
         *  \b Description: \n
         *        Method to round toward negative infinity (according to MATLAB, as reference)
         *  \param [in] no input parameters needed
         *  \return
         */
        void floor(Matrix &m2)
        const;					// round toward negative infinity (reference)

        /**
         *  \b Description: \n
         *        Method to round toward negative infinity (according to MATLAB)
         *  \param [in] no input parameters needed
         *  \return
         */
        Matrix floor() const;							// round toward negative infinity

        /**
         *  \b Description: \n
         *        Method to round nearest integer (according to MATLAB, as reference)
         *  \param [in] no input parameters needed
         *  \return
         */
        void round(Matrix &m2) const;					// round nearest integer (reference)

        /**
         *  \b Description: \n
         *        Method to round to nearest integer (according to MATLAB)
         *  \param [in] no input parameters needed
         *  \return
         */
        Matrix round() const;							// round nearest integer

        /**
         *  \b Description: \n
         *        Method to get the remainder after division (according to MATLAB, as reference)
         *  \param [in] [double] mean value
         *  \return
         */
        void rem(Matrix &m2,const double v)
        const;	// get the remainder after division (reference)

        /**
         *  \b Description: \n
         *        Method to get the remainder after division (according to MATLAB)
         *  \param [in] [double] mean value
         *  \return
         */
        Matrix rem(const double v) const;				// get the remainder after divvision



        // ***** statistical functions *****
        /**
         *  \b Description: \n
         *        Method to calculate the mean value
         *  \param [in] no input parameters needed
         *  \return [double] mean value
         */
        double meanD() const;     // mean value

        /**
         *  \b Description: \n
         *        Method to calculate the median value
         *  \param [in] no input parameters needed
         *  \return [double] median value
         */
        double median() const;    // median

        /**
         *  \b Description: \n
         *        Method to calculate the standard deviation
         *  \param [in] no input parameters needed
         *  \return [double] standard deviation
         */
        double stdD() const;     // standard deviation



        // ***** matrix operations *****
        /**
         *  \b Description: \n
         *        Method to get the transpose of a matrix
         *  \param [in] no input parameters needed
         *  \return [Matrix] transpose of a matrix
         */
        Matrix transpose() const; 			// returns the transpose of a matrix (as a copy)

        /**
         *  \b Description: \n
         *        Method to get the transpose of a matrix
         *  \param [in] [Matrix] arbitrary (empty) matrix
         *  \return [Matrix] transpose of a matrix >>> &mOut
         */
        void transpose(Matrix &mOut) const;

        /**
         *  \b Description: \n
         *        Method to extract the diagonal of a matrix or create a diagonal matrix
         *  \param [in] no input parameters needed
         *  \return [Matrix] diagonal of a matrix or diagonal matrix
         */
        Matrix diag( )
        const;     			// extract the diagonal of a matrix or create a diagonal matrix
        /**
         *  \b Description: \n
         *        Method to extract the diagonal of a matrix or create a diagonal matrix
         *  \param [in] [Matrix] arbitrary (empty) matrix
         *  \return [Matrix] diagonal of a matrix or diagonal matrix >>> &mOut
         */
        void diag(Matrix &mOut)
        const;     // extract the diagonal of a matrix or create a diagonal matrix

         /**
         *  \b Description: \n
         *        In comparison to diag() optimized method to generate a diagonal matrix based on a column- or row-vector.
         *        *this will be transformed. 
         */
        void to_diag_matrix();
        
                 /**
         *  \b Description: \n
         *        In comparison to diag() optimized method to generate a column-vector based on a diagonal matrix.
         *        *this will be transformed. Furthermore it can be definied if the diagonal elements will be the reciprocals or not.
         */
        void from_diag_matrix(bool reciprocal = false);
        
        /**
         *  \b Description: \n
         *        Method to create an identity matrix of size (dim x dim)
         *  \param [in] [int] dim: dimension
         *  \return [Matrix] identity matrix (dim x dim)
         */
        void eye( int dim);		// erstellt EinheitsMatrix

        /**
         *  \b Description: \n
         *        Method to create a matrix with (pseudo-) random numbers
         *  \param [in] no input parameters needed; matrix has to be initialized
         *  \return [Matrix] matrix with (pseudo-) random numbers
         */
        void random( );           // bereits initialisierte Matrix mit Zufalls

        /**
         *  \b Description: \n
         *  Method to create an 3D rotation matrix with angle psi around x-axis
         *  \param [in] [double] psi: rotation angle
         *  \return [Matrix] rotation matrix (dim x dim)
         */
        void rot3D_x( double psi );

        /**
         *  \b Description: \n
         *        Method to create an 3D rotation matrix with angle psi around y-axis
         *  \param [in] [double] psi: rotation angle
         *  \return [Matrix] rotation matrix (dim x dim)
         */
        void rot3D_y( double psi );

        /**
         *  \b Description: \n
         *   Method to create an 3D rotation matrix with angle psi around z-axis
         *  \param [in] [double] psi: rotation angle
         *  \return [Matrix] rotation matrix (dim x dim)
         */
        void rot3D_z( double psi );
        
        /**
        *  \b Description: \n
        *        This method creates a [size x 1] vector
        *        of zero mean Gaussian random numbers with variance one
        *  \param [in] no input parameters needed; matrix has to be initialized
        */
        void rand_norm( double mean = 0.0, double std = 1.0 );

        /**
         *  \b Description: \n
         *        Method for a element-by-element exponentiation of a matrix
         *  \param [in] [double] power
         *  \return [Matrix] new matrix
         */
        Matrix pow(double v) const;				// Element by element exponentiation

        /**
         *  \b Description: \n
         *        Method for a element-by-element division
         *  \param [in] [Matrix] &m2: other Matrix
         *  \return [Matrix] new matrix
         */
        Matrix div_elem( const Matrix &m2 ) const;  // Element by element division

        /**
         *  \b Description: \n
         *        Method for a element-by-element multiplication
         *  \param [in] [Matrix] &m2: other Matrix
         *  \return [Matrix] new matrix
         */
        Matrix mult_elem( const Matrix &m2 )
        const; // Element by element multiplication

        /**
         *  \b Description: \n
         *        Method to calculate column totals
         *  \param [in] no input parameters needed
         *  \return [Matrix] column totals
         */
        Matrix sum_col() const;    					// column totals (Spaltensumme)

        /**
         *  \b Description: \n
         *        Method to calculate column product
         *  \param [in] no input parameters needed
         *  \return [Matrix] column product
         */
        Matrix prod_col() ;    					// column totals (Spaltensumme)

        // +++ new: SH
        // norm and rand_norm class
        /**
        *  \b Description: \n
        *        This method calculates the norm of a [n x 1] matrix
        *  \param [in] no input parameters needed
        *  \return [matrix] [1 x 1] norm
        */
        Matrix norm( ) ;
        // --- SH


        void mult(double v);						// not implemented?


        /**
         *  \b Description: \n
         *        Method for a cholesky decomposition;
         *        using 'clapack_dpotrf' (Double GEneral TRiangular Factorize)
         *  \param [in] [Matrix] &A: other Matrix;
         *               calculate A = R'*R and save R in this-matrix;
         *               matrix A has to be positive definit
         */
        void chol( const Matrix &A , const char type = 'U');				// cholesky decomposition

        /**
         *  \b Description: \n
         *        Method to calculate the inverse of a matrix using a cholesky decomposition;
         *  \param [in] no input parameters needed
         *  \return [Matrix] inverse
         */
        Matrix chol_inv( )
        const;					// calculation of the inverse of a matrix using a cholesky decomposition

        /**
         *  \b Description: \n
         *        Method to solve a positive definite system of equations;
         *        using 'clapack_dposv' (Double symmetric POsitive definite)
         *  \param [in] [Matrix] &N
         */
        void solve_neq( const Matrix
                        &N );			// solve a positive definite system of equations


        /**
         *  \b Description: \n
         *        Method to solve a general system of equations;
         *        using 'clapack_dgesv' (Double GEneral SolVe)
         *  \param [in] [Matrix] &A
         */
        void solve_geeqs( const Matrix &A );		// solve a general system of equations

        void solve_qr( const Matrix &A, Matrix &Sxx );

        void solve_lu( const Matrix &A );

        void solve_svd( const Matrix &A );

        /**
         *  \b Description: \n
         *        Method to solve an Ordinary Least Squares system;
         *        using 'clapack_dposv' (Double GEneral SolVe)
         *  \param [in] [Matrix] &N
         */
        void ols( const Matrix &A );				// solve an Ordinary Least Squares system

        /**
         *  \b Description: \n
         *        Method to solve an Ordinary Least Squares system;
         *        using 'clapack_dgels' (Double GEneral Least Squares) 
         *        which makes use of an QR decomposition.
         *        If the system is underdetermined, the minimum norm solution is performed
         *  \param [in] [Matrix] &N
         */
        // void ols_qr( const Matrix &A );				// solve an Ordinary Least Squares system


        /**
         *  \b Description: \n
         *        Method to calculate the inserve a matrix;
         *        using 'clapack_dgetrf' (Double GEneral TRiangular Factorize)
         *        and 'clapack_dgetri' (Double GEneral TRiangular Invert)
         *  \param [in] no input parameters needed
         */
        void inv( );								// inverse of a matrix

        /**
         *  \b Description: \n
         *        Method to calculate the pseudo-inserve of a matrix;
         *        using SVD
         *  \param [in] no input parameters needed
         */
        Matrix pinv( );								// pseudo-inverse of a matrix

        
        /**
         *  \b Description: \n
         *        Method to calculate the inserve a matrix;
         *        using 'clapack_dgetrf' (Double GEneral TRiangular Factorize)
         *        and 'clapack_dgetri' (Double GEneral TRiangular Invert)
         *  \param [in] no input parameters needed
         */
        void inv_scal( );							// scaled inverse of a matrix

        /**
         *  \b Description: \n
         *        Method to do a singular value decomposition (SVD);
         *        A = USV'
         *        using 'dgesvd_' (LAPACK)
         *  \param [in] [Matrix] (empty) matrices U, S, VT
         */
        void svd( Matrix &U, Matrix &S, Matrix &VT,
                  int n = -999 ) const;	// singular value decomposition

        /**
         *  \b Description: \n
         *        Method to do a QR decomposition (QR);
         *        A = Q*R
         *        using 'dgeqrf_' (LAPACK)
         *  \param [in] [Matrix] (empty) matrices Q, R
         */
        void qr( Matrix &Q, Matrix &R ) const; 	// QR decomposition
        
        ivg::Matrix calc_impact_factors();

        void eig( Matrix &d, Matrix &U ) const;	// not implemented?

        /**
         *  \b Description: \n
         *        Method to calculate condition based on singular value decomposition
         *  \param [in] no input parameters needed
         */
        double cond()
        const;						// calculate condition based on singular value decomposition

        /**
         *  \b Description: \n
         *        Method to calculate rank based on singular value decomposition
         *  \param [in] no input parameters needed
         */
        int rank( double tol =
                      0.0 );		// calculate rank based on singular value decomposition

        /**
         *  \b Description: \n
         *        Method to calculate the trace of a matrix
         *  \param [in] no input parameters needed
         */
        double trace( ) const;				// calculate the trace of a matrix

        void abs(Matrix &m2) const;		// gets the absolute value

        Matrix abs() const;					// gets the absolute value

        Matrix absD() const;

        Matrix diff() const;

        Matrix sign() const;


        void solve_tri( Matrix &rhs ) const;		// TRiangular matrix Solve Matrix

        /**
         *  \b Description: \n
         *        Method to calculate the symmetric product;
         *        this = A^T * A
         *        using 'cblas dsyrk' (SYmmetric RanK update)
         *  \param [in] [Matrix] &A: matrix
         */
        void is_AtA_of( const Matrix &A ); 					// this = A^T * A

        /**
         *  \b Description: \n
         *        Method to calculate the product A^T * v
         *        using 'cblas dgemv' (GEnerel Matrix Vector product)
         *  \param [in] [Matrix] &A: matrix
         *               [Matrix] &v: vector (of type) matrix
         */
        void is_Atv( const Matrix &A, const Matrix &v );	// this = A^T * v

        /**
         *  \b Description: \n
         *        Method to calculate the product A^T * v
         *        using 'cblas dgemv' (GEnerel Matrix Vector product)
         *  \param [in] [Matrix] &A: matrix
         *               [Matrix] &v: vector (of type) matrix
         */
        void is_Atv2( const Matrix &A, const Matrix &v );	// this = A^T * v

        /**
         *  \b Description: \n
         *        Method to do a matrix multiplication;
         *        this = A(')*B(')
         *        using 'cblas dgemm' (GEnerel Matrix Matrix product)
         *  \param [in] [Matrix] &A: matrix
         *               [Matrix] &B: matrix
         */
        void is_product_of(const Matrix & A,
                           const Matrix & B, 								// this= A(')*B(')
                           CBLAS_TRANSPOSE TransA=CblasNoTrans, CBLAS_TRANSPOSE TransB=CblasNoTrans);

        /**
         *  \b Description: \n
         *        Method to do a matrix multiplication and add to this-matrix;
         *        this += A(')*B(')
         *        using 'cblas dgemm' (GEnerel Matrix Matrix product)
         *  \param [in] [Matrix] &A: matrix
         *               [Matrix] &B: matrix
         */
        void plus_product_of(const Matrix & A,
                             const Matrix & B,								// this+ = A(')*B(')
                             CBLAS_TRANSPOSE TransA=CblasNoTrans, CBLAS_TRANSPOSE TransB=CblasNoTrans);



        // ***** extended matrix operations ******
        /**
         *  \b Description: \n
         *        Method to estimate Continuous PiceWise Linear Functions (CPWLF) \n
         *        Call method with: this = vector (type Matrix) of observations
         *  \param [in] [Matrix] x: epochs of observations
         *               [Matrix] x0: vector (type Matrix) of parameters
         *               [double] w: weights (for constraints), [optional]
         */
        Matrix estimate_cpwlf( Matrix x, Matrix x0, double w = 1.0e-04 ) ;


        /**
         *  \b Description: \n
         *        Calculates the gamma function of each element of a 1xn matrix
         *  \param [in] no input parameters needed
         *  \return [Matrix] Matrix with the gamma function of each element
         */
        Matrix gamma_fct( ) ;


        /**
         *  \b Description: \n
         *        Method to interpolate data
         *  \param [in] [Matrix] data_epochs: all data epochs
         *              [Matrix] epoch:       epoch
         *              [std::string] interpolation_type: interpolation_type \n
         *                              "linear"      linear interpolation
         *                              "polynomial"  polynomial interpolation
         *                              "cspline"     cubic spline interpolation
         */
        Matrix interpolate( Matrix data_epochs, double epoch,
                            std::string interpolation_type ) const;



        // ***** input and output *****
        /**
         *  \b Description: \n
         *        Formatted output on the console
         *  \param [in] no input parameters needed
         */
        void show(const int &nk = 8) const;        // formatted output

        /**
         *  \b Description: \n
         *        Formatted output on the console
         *  \param [in] no input parameters needed
         */
        void disp() const;							 // formatted output

        /**
         *  \b Description: \n
         *        Formatted output of the dimension on the console with cout
         *  \param [in] no input parameters needed
         */
        void cout_size(const string = "Matrix" )
        const;  // formatted output of the dimension (r,c)
        
        /**
         *  \b Description: \n
         *        Formatted output of the dimension on the console with cerr
         *  \param [in] no input parameters needed
         */
        void cerr_size(const string = "Matrix" )
        const;  // formatted output of the dimension (r,c)

        /**
         *  \b Description: \n
         *        Get matrix from binary file
         *  \param [in] [string] &file: binary file name
         */
        void load_bin(const string &file );        // get Matrix from binary file

        /**
         *  \b Description: \n
         *        Write matrix to binary file
         *  \param [in] [string] &file: binary file name
         */
        void save_bin(const string &file ) const; // write Matrix to binary file

        /**
         *  \b Description: \n
         *        Get matrix from ascii file
         *  \param [in] [string] &file: ascii file name
         */
        void load_ascii(const string &file);     	// get Matrix from ascii file

        /**
         *  \b Description: \n
         *        Get matrix from ascii file
         *  \param [in] [string] &file: ascii file name
         */
        void load_vec( const string &file );   	// get Matrix from ascii file


        /**
         *  \b Description: \n
         *        Get one single matrix from MATLAB binary file. Several full or sparse matrices might
        *        be stored in the mat-file, however, only double values are allowed.
        *        matlab> save -mat 'matfile.mat' A b c
         *  \param [in] [string] &varname: matrixname used in MATLAB
         *              [string] &matfile: binary file name
         */
        void  load_mat_bin( string varname,
                            string matfile );  // get matrix from MATLAB binary

        /**
         *  \b Description: \n
         *        Write matrix to ascii file
         *  \param [in] [string] &file: ascii file name
         */
        void save_ascii(const string &file) const; 	// save Matrix to ascii file

        /**
         *  \b Description: \n
         *        Get matrix from ascii file with #r #c value in each row
         *  \param [in] [string] &file: ascii file name
         */
        void save_ascii_idx(const string &file)
        const;	// save Matrix to ascii file with #r #c value in each row

        /**
         *  \b Description: \n
         *        Get matrix from ascii file (in matrix format)
         *  \param [in] [string] &file: ascii file name
         */
        void save_matrix(const string &file,
                         const int &nk=10) const;	// save Matrix to ascii file in Matrix-Format



        // ***** set-functions / matrix modifications *****
        /**
         *  \b Description: \n
         *        This method creates a large matrix (>> this-matrix) consisting of an i-by-j tiling of copies of this matrix
         *  \param [in] [int] i: rows; j, columns
         */
        void repmat(const Matrix & m2,  int i, int j );		// replicate and tile matrix

        //Matrix repmat(const Matrix & m2,  int i, int j ) const;
        
         /**
         *  \b Description: \n
         *        This method reshape a vector to  matrix with the same number of elements as the vector
         *  \param [in] no input parameters needed
         */
        void vec2Mat(int c);

        
        /**
         *  \b Description: \n
         *        This method creates a large vector of type matrix with the dimension [number-of-elements x 1]
         *  \param [in] no input parameters needed
         */
        void vec();								// large vector with dimension [number-of-elements x 1]

        //Matrix vec() const;

        /**
         *  \b Description: \n
         *        This method removes a selected row in a matrix
         *  \param [in] [int] idx: index of the row to be removed
         */
        void rem_r( int idx, int n=1 );      // remove row

        /**
         *  \b Description: \n
         *        This method removes a selected column in a matrix
         *  \param [in] [int] idx: index of the column to be removed
         */
        void rem_c( int idx, int n=1 );      // remove column
        
        
        /**
         *  \b Description: \n
         *        This method removes several selected columns defined by a int-vector
         *  \param [in] [vector<int>] idx: vector of indexes of the colomns to be removed
         */
        void rem_c( vector<int> idxs );
        
        /**
         *  \b Description: \n
         *        This method removes several selected rows defined by a int-vector
         *  \param [in] [vector<int>] idx: vector of indexes of the rows to be removed
         */
        void rem_r( vector<int> idxs );

        /**
         *  \b Description: \n
         *        This method removes a selected row and column in a matrix
         *  \param [in] [int] idx: index of the row/column to be removed
         */
        void rem_rc( int idx );     		// remove row and column

        /**
         *  \b Description: \n
         *        This method adds a sub-matrix to this-matrix
         *  \param [in] [vector<int>] &rows: index vector (row); &col: index vector (column)
         *               [Matrix] &m2: sub-matrix to be added
         */
        void add_sub(const vector<int> &rows, const vector<int> &cols,
                     const Matrix &m2 );		// add sub-matrix

        /**
         *  \b Description: \n
         *        This method sets a sub-matrix to this-matrix
         *  \param [in] [vector<int>] &rows: index vector (row); &col: index vector (column)
         *               [Matrix] &m2: sub-matrix
         */
        void set_sub(const vector<int> &rows, const vector<int> &cols,
                     const Matrix &m2 );		// set sub-matrix using index vectors

        /**
         *  \b Description: \n
         *        This method sets a sub-matrix to this-matrix
         *  \param [in] [int] i: row; j: column (index for the top-left element in the matrix)
         *               [Matrix] &other: sub-matrix
         */
        void set_sub(int i, int j,
                     const Matrix &other );										// set sub-matrix using indices of the top-left element

        /**
         *  \b Description: \n
         *        This method appends a sub-matrix to the end of this-matrix (rows)
         *  \param [in] [Matrix] &m2: sub-matrix to be appended
         */
        void append_rows( const Matrix & m2 ); 		// append new rows on a matrix

        /**
         *  \b Description: \n
         *        This method appends a value to the end of this-matrix (only if matrix is [n x 1])
         *  \param [in] [double] v: value
         *
         */
        void append_rows( double
                          v );          		// append a value to the end of this-matrix

        /**
         *  \b Description: \n
         *        This method appends a sub-matrix to the end of this-matrix (columns)
         *  \param [in] [Matrix] &m2: sub-matrix to be appended
         */
        void append_cols( const Matrix & m2 ); 		// append new columns on a matrix


        /**
         *  \b Description: \n
         *        This method changes the position of a row in a matrix
         *  \param [in] [int] start: initial position; dest: destination
         */
        void change_row( int start,
                         int dest );		// change the position of a row in a matrix

        /**
         *  \b Description: \n
         *        This method changes the position of a column in a matrix
         *  \param [in] [int] start: initial position; dest: destination
         */
        void change_col( int start,
                         int dest );		// change the position of a column in a matrix


        /**
         *  \b Description: \n
         *        This method performs a numerical sort of the columns in a matrix
         *  \param [in] [int] col: column to sort by
         */
        void sort_cols( int c );     				// numerical sort of the columns in a matrix


        /**
         *  \b Description: \n
         *        This method numerically sorts a column in a matrix
         *  \param [in] [int] dim: dimension of the matrix
         */
        void sort_col( int col );     				// numerical sort of a column in a matrix


        /**
         *  \b Description: \n
         *        This method numerically sorts all columns in a matrix
         *  \param [in] [int] dim: dimension of the matrix
         */
        void sort_cols( );     						// numerical sort of all columns in a matrix


        /**
         *  \b Description: \n
         *        This method sets a column of a matrix
         *  \param [in] [int] col: column to be modified
         *               [Matrix] &m2: matrix with new column
         */
        void set_col(int i, const Matrix
                     &m2);		// modify a column of a matrix;  this(:,i)= m2;

        /**
         *  \b Description: \n
         *        This method eliminates the row/column given by an index vector
         *        this(idx) = [];
         *  \param [in] [vector<int>] index vector
         */
        void elim_idx(vector<int>
                      idx);				// eliminate the row/column given by an index vector

        /**
         *  \b Description: \n
         *        This method modifies a sub-matrix of this-matrix (using index std::vectors)
         *  \param [in] [vector<int>] &idx: index vector
         *               [Matrix] &m2: new sub-matrix
         */
        void setIdx(const vector<int> &idx,const Matrix &m2);

        /**
         *  \b Description: \n
         *        This method modifies a sub-matrix of this-matrix (using index vectors of type Matrix)
         *  \param [in] [Matrix] &idx: index vector
         *               [Matrix] &m2: new sub-matrix
         */
        void setIdx(const Matrix &idx,const Matrix &m2);

        /**
         *  \b Description: \n
         *        This method sets any number of values of this-matrix to a given value (using index std::vector)
         *  \param [in] [vector<int>] &idx: index vector
         *               [double] value
         */
        void setIdx(const vector<int> &idx,const double &v);

        /**
         *  \b Description: \n
         *        This method sets any number of values of this-matrix to a given value (using index vector of type matrix)
         *  \param [in] [Matrix] &idx: index vector
         *               [double] value
         */
        void setIdx(const Matrix &idx,const double &v);

        void set(const vector<Matrix> &vMat);		// join any number of matrices

        //void set(const string &strIn);

        /**
         *  \b Description: \n
         *        This method sets all elements to zero
         *  \param [in] [in] no input parameters needed
         */
        void zero( );
        
        /**
         *  \b Description: \n
         *        creates block toepliz Matrix
         *  \param [in] [const std::vector<ivg::Matrix>&] vector of quadratic Matrices
         *              [const char] symmetric = 'S' , not symmetric = 'N'
         */
        void blockToepliz( const std::vector<ivg::Matrix>& mats, const char type = 'S');


        // ***** index search / find-methods *****
        /**
         *  \b Description: \n
         *        This method gets the indices of the elements which are equal to a selected value
         *  \param [in] [double] v: value
         */
        vector<int> find_idx( double v ) const;

        // find: Idex finden mit (>,<,>=,<=,(==))
        vector<int> find_idx(bool (*pfunc)(double i, double j),
                             const double val) const;

        // mit zwei Bedingungen &&
        vector<int> find_idx(bool (*pfunc1)(double i, double j),const double val1,
                             bool (*pfunc2)(double i, double j),const double val2 ) const;

        // find(this_nXm >,<,>=,<=,(==) vals_nXm ): Idex finden mit n x m Werten
        // Idex-vektor besitzt immer die dim: [found,1]
        vector<int> find_idx(bool (*pfunc)(double i, double j),
                             const Matrix &vals) const;

        vector<int> find_idx(bool (*pfunc1)(double i, double j),const Matrix &vals1,
                             bool (*pfunc2)(double i, double j),const Matrix &vals2 ) const;


        /**
         *  \b Description: \n
         *        Iterator of the first element which is equal to a selected value
         *  \param [in] [vector<double>::iterator] start: start index
         *               [vector<double>::iterator] end: end index
         *               [double] v: value
         */
        vector<double>::iterator find_elem( vector<double>::iterator start,
                                            vector<double>::iterator end,
                                            double v ) const;

        // Elemente finden mit Pointer auf Bool funktion (>,<,>=,<=,(==))
        Matrix find_elem(bool (*pfunc)(double i, double j),const double val) const;

        Matrix find_elem(bool (*pfunc)(double i, double j),const double val,
                         vector<int> &idx) const;

        void find_elem(bool (*pfunc)(double i, double j),const double val,
                       vector<int> &idx,Matrix &m2 ) const;

        // find_elem mit zwei Bedingungen
        void find_elem(bool (*pfunc1)(double i, double j),const double val1,
                       bool (*pfunc2)(double i, double j),const double val2,
                       vector<int> &idx,Matrix &m2 ) const;

        Matrix find_elem(bool (*pfunc1)(double i, double j),const double val1,
                         bool (*pfunc2)(double i, double j),const double val2,
                         vector<int> &idx) const;

        Matrix find_elem(bool (*pfunc1)(double i, double j),const double val1,
                         bool (*pfunc2)(double i, double j),const double val2) const;

        // Indexvektor erzeugen
        void indexvec(double start, double schrittweite,double ende);



        // ***** get-functions *****
        /**
         *  \b Description: \n
         *        Method to get a vector with all matrix elements (>> private: _data)
         *  \param [in] no input parameters needed
         *  \return [vector<double>] vector with all matrix elements
         */
        vector<double> get_vec()
        const; 				// get double vector with all matrix elements

        /**
         *  \b Description: \n
         *        Method to get a vector with all matrix elements (>> private: _data)
         *  \param [in] no input parameters needed
         *  \return [vector<double>] vector with all matrix elements
         */
        vector<double>
        get_data_vec( );					// get double vector with all matrix elements

        /**
         *  \b Description: \n
         *        Method to get a vector with all matrix elements
         *  \param [in] [string] art: type; set art = "int"
         *  \return [vector<double>] vector with all matrix elements
         */
        vector<int> get_vec(const string art)
        const;	// get integer vector with all matrix elements

        /**
         *  \b Description: \n
         *        Method to get a vector with the elements given by an index vector (as reference)
         *  \param [in] [vector<int>] &idx: index vector
         *  \return [Matrix] matrix with all elements given by the index vector (as a reference)
         */
        void get_vec(const vector<int> &idx,
                     Matrix &m2) const; 	// get vector (type matrix) with elements of a given index vector (reference)

        /**
         *  \b Description: \n
         *        Method to get a vector with the elements given by an index vector
         *  \param [in] [vector<int>] &idx: index vector
         *  \return [Matrix] matrix vector with all elements given by the index vector
         */
        Matrix get_vec(const vector<int> &idx)
        const;				// get vector (type matrix) with elements of a given index vector

        /**
         *  \b Description: \n
         *        Method to get a vector with the elements given by an index vector (as reference)
         *  \param [in] [matrix] &idx: matrix with indices
         *  \return [Matrix] matrix with all elements given by the index vector (as a reference)
         */
        void get_vec(const Matrix &idx,
                     Matrix &m2) const;			// get vector (type matrix) with elements of a given index matrix (reference)

        /**
         *  \b Description: \n
         *        Method to get a vector with the elements given by an index vector
         *  \param [in] [matrix] &idx: matrix with indices
         *  \return [Matrix] matrix vector with all elements given by the index vector
         */
        Matrix get_vec(const Matrix &idx)
        const;					// get vector (type matrix) with elements of a given index matrix

        /**
         *  \b Description: \n
         *        Method to get a row of a matrix
         *  \param [in] [int] &row: index of a row
         *  \return [Matrix] matrix with only one row of the this-matrix
         */
        Matrix get_row(const int &row) const;						// get a row of a matrix

        /**
         *  \b Description: \n
         *        Method to get any number of rows of a matrix
         *  \param [in] [vector<int>] &row: index vector
         *  \return [Matrix] matrix with any number of rows of the this-matrix
         */
        Matrix get_rows(const vector<int> &rows )
        const;			// get any number of rows of a matrix

        /**
         *  \b Description: \n
         *        Method to get any number of (following) rows of a matrix
         *  \param [in] [int] start: first index; end: last index
         *  \return [Matrix] matrix with any number of (following) rows of the this-matrix
         */
        Matrix get_rows(int start,
                        int end ) const;					// get any number of (following) rows of a matrix

        /**
         *  \b Description: \n
         *        Method to get a column of a matrix
         *  \param [in] [int] &col: index of a row
         *  \return [Matrix] matrix with only one column of the this-matrix
         */
        Matrix get_col(const int &col) const; 						// get a column of a matrix

        /**
         *  \b Description: \n
         *        Method to get any number of column of a matrix
         *  \param [in] [vector<int>] &col: index vector
         *  \return [Matrix] matrix with any number of column of the this-matrix
         */
        Matrix get_cols(const vector<int> &rows )
        const;			// get any number of column of a matrix

        /**
         *  \b Description: \n
         *        Method to get any number of (following) column of a matrix
         *  \param [in] [int] start: first index; end: last index
         *  \return [Matrix] matrix with any number of (following) column of the this-matrix
         */
        Matrix get_cols(int start,
                        int end ) const;					// get any number of (following) column of a matrix

        /**
         *  \b Description: \n
         *        Method to get a sub-matrix of this-matrix
         *  \param [in] [vector<int>] &row: index vector for rows; &col: index vector for columns
         *  \return [Matrix] sub-matrix of this-matrix
         */
        Matrix get_sub(const vector<int> &rows,
                       const  vector<int> &cols ) const; 	// get sub-matrix using index vectors

        /**
         *  \b Description: \n
         *        Method to get a sub-matrix of this-matrix
         *  \param [in] [int] startR: first row index;    endR: last row index
         *               [int] startC: first column index; endC: last column index
         *  \return [Matrix] sub-matrix of this-matrix
         */
        Matrix get_sub(int startR, int startC, int endR,
                       int endC ) const;			// get sub-matrix using start and end tags

        /**
         *  \b Description: \n
         *        interator to the first element of _data
         *  \param [in] no input parameters needed
         *  \return [vector<double>::iterator] iterator to the first element in _data
         */
        vector<double>::iterator begin()
        {
            return _data.begin();
        }

        /**
         *  \b Description: \n
         *        interator to the  past-the-end element of _data
         *  \param [in] no input parameters needed
         *  \return [vector<double>::iterator] iterator to the  past-the-end element in _data
         */
        vector<double>::iterator end()
        {
            return _data.end();
        }
        
        /**
         *  \b Description: \n
         *       refernce to the last element of _data
         *  \param [in] no input parameters needed
         *  \return [vector<double>::iterator] iterator to the  past-the-end element in _data
         */
        double back()
        {
            return _data.back();
        }
        
         /**
         *  \b Description: \n
         *       refernce to the last element of _data
         *  \param [in] no input parameters needed
         *  \return [vector<double>::iterator] iterator to the  past-the-end element in _data
         */
        const double back() const
        {
            return  _data.back();
        }


        /**
         *  \b Description: \n
         *        interator to the first element of _data
         *  \param [in] no input parameters needed
         *  \return [vector<double>::iterator] iterator to the first element in _data
         */
        vector<double>::const_iterator begin( ) const
        {
            return _data.begin();
        }

        /**
         *  \b Description: \n
         *        interator to the  past-the-end element of _data
         *  \param [in] no input parameters needed
         *  \return [vector<double>::iterator] iterator to the  past-the-end element in _data
         */
        vector<double>::const_iterator end( ) const
        {
            return _data.end();
        }

        /**
         *  \b Description: \n
         *        pointer for _data: read and write access
         *  \param [in] no input parameters needed
         *  \return [double*] pointer to the first element in _data
         */
        double * data_ptr()
        {
            return &_data.at(0);    // pointer for _data: read and write access
        }
        
        vector<double> * data_vec_ptr()
        {
            return &_data;    // pointer for _data: read and write access
        }

        /**
         *  \b Description: \n
         *        pointer for _data: read-only access
         *  \param [in] no input parameters needed
         *  \return [double*] pointer to the first element in _data
         */
        const double * data_ptr() const
        {
            return &_data.at(0);    // pointer for _data: read-only access
        }
        
        bool is_vector(){if(rows()==1||cols()==1) return true;else return false;};

        void cross(const ivg::Matrix &a,const ivg::Matrix &b);
    private:
        

        // ==============================================
        // ======== private MEMBER-Functions: ===========
        // ==============================================

        /**
         *  \b Description: \n
         *        Method to get the (i,j)th element of the matrix
         *  \param [in] [int] i: row; j: column
         *  \return [double] element of the (i,j)th position
         */
        double _get(const int &i,
                    const int &j) const; 			// returns the (i,j)th element of a matrix

        /**
         *  \b Description: \n
         *        Method to set the (i,j)th element of the matrix to a given value
         *  \param [in] [int] i: row; j: column
         *               [double] v: value
         */
        void _set(const int &i,const int &j,
                  const double &v); 	// sets the (i,j)th element of a matrix to v

        /**
         *  \b Description: \n
         *        Method to set the (i,i)th element of the matrix to a given value
         *  \param [in] [int] i: row
         *               [double] v: value
         */
        void _set(const int &i, const double
                  &v );				// sets the (i,i)th element of a matrix to v

        /**
         *  \b Description: \n
         *        Method to multiply each element of the matrix with a scalar value
         *  \param [in] [double] v: value
         */
        void _mult(const double
                   &v);    // multiplies each element with a scalar value


        // ==============================================
        // ======== class variables / attributes: =======
        // ==============================================

        size_t _rows;         // number of rows
        size_t _cols;         // number of columns

        vector<double> _data; // data structure for matrix elements: std::vector
        
        static unsigned int _rng_calls;
};

}

#endif  // MATRIX_H


