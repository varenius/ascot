#include "matrix.h"

namespace ivg
{

unsigned Matrix::_rng_calls = 0;    
    
// ===========================================================================
//                      Konstruktoren und Destruktor
// ===========================================================================

Matrix::Matrix ()
{
#ifdef Matrix_DEBUG
    cerr << "+++ Konstruktor der Klasse Matrix" << endl;
#endif

    _rows  = 0;
    _cols  = 0;
}

Matrix::Matrix (int r, int c, double v )
{
#ifdef Matrix_DEBUG
    cerr << "+++ Konstruktor der Klasse Matrix" << endl;
#endif

    _rows  = r;
    _cols  = c;

    // Datenstruktur anlegen
    _data.resize( r * c, v );

    if( r*c == 0 )
    {
        _rows = 0;
        _cols = 0;
    }
}
/*
// Konstruktor mit variablem Matrix input zum zusammensetzen
Matrix::Matrix( int amount, ...)
{

  int i;
  Matrix val;
  va_list vl;
  va_start(vl,amount);
  // mit erster Matrix fuellen
  val=va_arg(vl,Matrix);
  *this=val;
  /*
  void* vp;
  while ((vp = va_arg(vl, void*))!= NULL)
   {
	}  * /
  // Schleife beginnt mit 2. Matrix
  for (i=1;i<amount;i++)
  {
    val=va_arg(vl,Matrix);
    append_cols(val);
  }
  va_end(vl);

} */

//mehrere matrizen als vector<Matrix> zusammenfassen
Matrix::Matrix(const vector<Matrix> &vMat)
{
    int num = vMat.size();

    *this= vMat.at( 0 );

    for ( int x = 1; x < num; x++ )        // Loop until all matrices are added
        append_cols( vMat.at( x ) );
    

}

Matrix::Matrix( vector<double>::iterator begin, vector<double>::iterator end,
                int r, int c )
{
    if( end - begin != r*c )
    {
        stringstream errormessage;
        errormessage << "Matrix::Matrix( vector<double>::iterator begin, "
                     << "vector<double>::iterator end, int r, int c ): "
                     << "ERROR: invalid number of elements " << " Exiting";
        throw logic_error( errormessage.str() );
    }

    _data.resize( r * c );
    copy( begin, end, _data.begin() );
    _rows = r;
    _cols = c;
    
}


Matrix::Matrix(const vector<double> &vec)
{
    _data.resize( vec.size() );
    copy( vec.begin(), vec.end(), _data.begin() );
    _rows = vec.size();
    _cols = 1;
    
}

Matrix::Matrix(const vector<int> &vec)
{
    _data.resize( vec.size() );
    copy( vec.begin(), vec.end(), _data.begin() );
    _rows = vec.size();
    _cols = 1;
    
}

// Matrix erzeugen mit Indexvektoren und Werten, Dim: [r,c] mit werten v sonst (optional)
Matrix::Matrix(const vector<int> &rows, const vector<int> &cols,
               const  vector<double> &values,
               int r, int c, double v)
{
    if( rows.size() == cols.size()  && rows.size()==values.size())
    {

        resize( r , c, v );

        for( int i = 0; i < values.size() ; ++i )
        {
            //cout << "groesse; " << _rows << ","<< _cols << endl;
            //cout << "[" << rows.at(i) << "," << cols.at(i) << "," << values.at(i) << "]" << endl;
            _set(rows.at(i), cols.at(i), values.at(i));
        }
        

    }
    else
    {
        stringstream errormessage;
        errormessage << "Matrix::Matrix(vector<int> rows, vector<int> cols,"
                     << "vector<double> values, int r, int c, double v = 0.0 )"
                     << "ERROR: invalid number of elements " << " Exiting";
        throw logic_error( errormessage.str() );
    }
}

Matrix::Matrix(double start, double schrittweite, double ende, int cols)
{
    indexvec(start, schrittweite, ende);
    Matrix m=*this;
    repmat(m,1,cols);
    
    
    return;
}

// this= other
Matrix::Matrix( const Matrix &other)
{
#ifdef Matrix_DEBUG
    cerr << "+++ Kopierkonstruktor der Klasse Matrix" << endl;
#endif

    _rows = other._rows;
    _cols = other._cols;
    _data.resize( _rows * _cols );
    copy ( other._data.begin(), other._data.end(), _data.begin() );
}

/* klappt nicht <= sagt Andi
// Aus string Matrix machen: Matrix("1,2;3,4")
Matrix::Matrix(const string &strIn)
{
    set(strIn);
}
*/

Matrix::~Matrix ()
{
#ifdef Matrix_DEBUG
    cerr << "--- Destruktor der Klasse Matrix" << endl;
#endif
}

// ===========================================================================
//                      Operatoren der Klasse Matrix
// ===========================================================================

// ...........................................................................
// Zuweisungsoperator =				this=other;
Matrix & Matrix::operator=( const Matrix &other )
{
#ifdef Matrix_DEBUG
    cerr << "+++ Zuweisungsoperator der Klasse Matrix" << endl;
#endif

    _rows = other._rows;
    _cols = other._cols;
    _data.resize( _rows * _cols );
    copy ( other._data.begin(), other._data.end(), _data.begin() );

    return *this;
}

Matrix & Matrix::operator=(const vector<int> &vec)
{
    _data.resize( vec.size() );
    copy( vec.begin(), vec.end(), _data.begin() );
    _rows = vec.size();
    _cols = 1;
}

// ...........................................................................
// ()-Operator fuer 2 Indizes - schreibend
double & Matrix::operator()( int i, int j )
{
    //if (i > _rows && j > _cols)
    //{
    //stringstream errormessage;
    //errormessage << "double & Matrix::operator()( int i, int j ): "
    //"ERROR: index  exeeds dimensions (" << i << ","<< j<<  ") Exiting";
    //throw logic_error( errormessage.str() );
    //}
    int idx = j * _rows + i;
    return( _data.at( idx ) );
}

// ()-Operator fuer 2 Indizes - lesend
double Matrix::operator()( int i, int j ) const
{
    int idx = j * _rows + i;
    return( _data.at( idx ) );
}

// ()-Operator fuer 1 Index -  schreibend
double & Matrix::operator()( int i )
{
    //if (i > numel())
    //{
    //stringstream errormessage;
    //errormessage << "double & Matrix::operator()( int i ): "
    //"ERROR: index  exeeds dimensions" << i <<  " Exiting";
    //throw logic_error( errormessage.str() );
    //}
    return( _data.at( i ) );
}

// ()-Operator fuer 1 Index - lesend
double Matrix::operator()( int i ) const
{
    return( _data.at( i ) );
}


// ()-Operator, um i-te Spalte zurueck zugeben (lesend)
Matrix Matrix::operator()( char c, int i ) const
{
    if( c == ':' )
    {
        Matrix vOut( _rows, 1, 0.0 );
        copy ( _data.begin()+i*_rows, _data.begin()+(i+1)*_rows, vOut._data.begin() );
        return vOut;
    }
    else
    {
        stringstream errormessage;
        errormessage << "Matrix Matrix::operator()( char c, int i ) const: "
                     << "ERROR: invalid signs " << c << " Exiting";
        throw runtime_error( errormessage.str() );
    }
}

// ()-Operator, um i-te Spalte zurueck zugeben (lesend)
Matrix Matrix::operator()( char c, std::vector<int> idx ) const
{
    if( c == ':')
    {
        Matrix vOut( _rows, idx.size(), 0.0 );
        for( int i=0;i<idx.size();++i )
            vOut.set_sub( 0,i,this->operator()( ':',idx.at(i) ) );
        return vOut;
    }
    else
    {
        stringstream errormessage;
        errormessage << "Matrix Matrix::operator()( char c, vector<int> idx )"
                     << "const: ERROR: invalid signs " << c << " Exiting";
        throw runtime_error( errormessage.str() );
    }
}

// ()-Operator, um i-te Spalte zurueck zugeben (lesend)
Matrix Matrix::operator()(const string c, int i ) const
{
    return this->operator()( c.c_str()[0],i );
}

Matrix Matrix::operator()( string c, std::vector<int> idx ) const
{
    return this->operator()( c.c_str()[0],idx );
}

//  m2=this(idx);
Matrix Matrix::operator()(const vector<int> &idx) const
{
    Matrix mOut;
    get_vec(idx,mOut);
    return mOut;
}

//  m2=this(idx);
Matrix Matrix::operator()(const Matrix &idx) const
{
    Matrix mOut;
    get_vec(idx,mOut);
    return mOut;
}

//  m2=this(rows,cols);
Matrix Matrix::operator()(const vector<int> &idx_r,
                          const vector<int> &idx_c) const
{
    return get_sub(idx_r, idx_c);
}

//  m2=this(rows,cols);
Matrix Matrix::operator()(const Matrix &idx_r, const string c) const
{
    if (c!=":")
    {
        stringstream errormessage;
        errormessage <<
                     "Matrix Matrix::operator()(const Matrix &idx_r, const string c) const: "
                     << "please use \":\" for columindex. Exiting!";
        throw runtime_error( errormessage.str() );
    }
    if (idx_r.max() > _rows)
    {
        stringstream errormessage;
        errormessage <<
                     "Matrix Matrix::operator()(const Matrix &idx_r, const string c) const: "
                     << "rowindex exceeds dimensions. Exiting!";
        throw runtime_error( errormessage.str() );
    }
    // Spaltenvektor fuer alle anlegen
    vector<int> idx_c(_cols);
    for (int i=0; i < _cols; i++)
        idx_c.at(i)=i;
    //cout << maxVec(idx_r.vec("int")) << ","<< maxVec(idx_c)<< endl;
    return get_sub(idx_r.get_vec("int"), idx_c);
}

//  m2=this(rows,cols);
Matrix Matrix::operator()(const vector<int> &idx_r, const string c) const
{
    if (c!=":")
    {
        stringstream errormessage;
        errormessage <<
                     "Matrix Matrix::operator()(const vector<int> &idx_r, const string c) const: "
                     << "please use \":\" for columindex. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    if (maxVec(idx_r) > _rows)
    {
        stringstream errormessage;
        errormessage <<
                     "Matrix Matrix::operator()(const vector<int> &idx_r, const string c) const: "
                     << "rowindex exceeds dimensions. Exiting!";
        throw runtime_error( errormessage.str() );
    }
    // Spaltenvektor fuer alle anlegen
    vector<int> idx_c(_cols);
    for (int i=0; i < _cols; i++)
        idx_c.at(i)=i;
    //cout << maxVec(idx_r.vec("int")) << ","<< maxVec(idx_c)<< endl;
    return get_sub(idx_r, idx_c);
}

// Operator mit 4 indizes: andere Def wie get_sub: vgl MATLAB: this(startR: endR, startC:endC)
Matrix Matrix::operator()(int startR, int endR, int startC , int endC ) const
{
    return get_sub(startR, startC, endR, endC );
}

/*
Matrix::operator double() const
{
	if (numel()!=1)
	{
      stringstream errormessage;
      errormessage << "double Matrix::operator double() const: "
          "cannot convert Matrix 2 double, numel()=" << numel();
      throw logic_error( errormessage.str() );
   }
	return static_cast<double> (_data.at(0));
}
*/


// ...........................................................................
//                                MULTIPLIKATION
// ...........................................................................
//
// Multiplikation mit einem Skalar
Matrix Matrix::operator*(const double &v  ) const
{
    Matrix mOut(rows(),cols());
    transform (_data.begin() ,_data.end() ,mOut._data.begin() ,
               bind2nd( multiplies < double >() ,v));

    return mOut;
}

Matrix & Matrix::operator*=( double v )
{
    cblas_dscal( _rows*_cols, v, data_ptr(), 1 );
    return *this;
}

Matrix Matrix::operator/( double v )
{
    Matrix mOut = *this;
    v = 1.0 / v;
    cblas_dscal( _rows*_cols, v, mOut.data_ptr(), 1 );

    return mOut;
}

Matrix Matrix::operator/(const Matrix &m2 )
{
    if( _rows!= m2._rows || _cols!= m2._cols )
    {
        stringstream errormessage;
        errormessage << "Matrix Matrix::operator/(const Matrix &m2 ): "
                     << "Matrix dimensions do not agree. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    return div_elem(m2);
}

// ...........................................................................
// Matrizenmultiplikation
Matrix Matrix::operator*( const Matrix &m2 ) const
{
    Matrix  mOut;
    // elementweise multiplikation, wenn Vektor
    if (_cols==m2._cols && _rows==m2._rows  && _cols==1)
        mOut= mult_elem(m2);

    else if( _cols != m2._rows )
    {
        stringstream errormessage;
        errormessage << "Matrix Matrix::operator*( const Matrix &m2 ) const: "
                     << "Matrix dimensions do not agree. Exiting!";
        throw runtime_error( errormessage.str() );
    }
    if( _cols == m2._rows )
    {
        mOut.resize( _rows, m2._cols, 0.0 );

        // level 3 BLAS - Matrix * Matrix
        cblas_dgemm ( CblasColMajor ,CblasNoTrans ,CblasNoTrans ,_rows, m2._cols,
                      _cols, 1.0, data_ptr(), _rows, m2.data_ptr(), m2._rows, 0.0,
                      mOut.data_ptr(), _rows );
    }
    return mOut;
}


Matrix & Matrix::operator*=( const Matrix & m2 )
{
    // elementweise multiplikation, wenn Vektor
    if (_cols==m2._cols && _rows==m2._rows  && _cols==1)
        *this = mult_elem(m2);

    if( _cols != m2._rows )
    {
        stringstream errormessage;
        errormessage << "Matrix Matrix::operator*( const Matrix &m2 ) const: "
                     << "Matrix dimensions do not agree. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    if( _cols == m2._rows )
    {
        // level 3 BLAS - Matrix * Matrix
        Matrix tmp = *this;
        resize( _rows, m2._cols );

        cblas_dgemm ( CblasColMajor ,CblasNoTrans ,CblasNoTrans ,_rows, m2._cols,
                      _cols, 1.0, tmp.data_ptr(), _rows, m2.data_ptr(), m2._rows, 0.0,
                      data_ptr(), _rows );
    }
    return *this;
}

void Matrix::cross(const ivg::Matrix &a,const ivg::Matrix &b)
{
    if(a.numel()!=3 || b.numel()!=3)
    {
        stringstream errormessage;
        errormessage << "Matrix::cross(const ivg::Matrix &a,const ivg::Matrix &b) "
                     << "wrong dimensions: only 3x1/1x3 is possible. Exiting!";
        throw runtime_error( errormessage.str() );
    }
    
    ivg::Matrix tmp(3,1,0.0);
    tmp(0) = a(1)*b(2)-a(2)*b(1);
    tmp(1) = a(2)*b(0)-a(0)*b(2);
    tmp(2) = a(0)*b(1)-a(1)*b(0);
    
    *this = tmp;
}

// ...........................................................................
//                           ADDITION / SUBTRAKTION
// ...........................................................................
//
// Addition eines Skalars
Matrix Matrix::operator+( double v ) const
{
    Matrix mOut(rows(),cols());
    transform (_data.begin() ,_data.end() ,mOut._data.begin() ,
               bind2nd( plus < double >() ,v));

    return mOut;
}

// Addition eines Skalars
// this += Matrix
Matrix Matrix::operator+=( double v )
{
    transform (_data.begin() ,_data.end() , _data.begin() ,
               bind2nd( plus < double >() ,v));
    return *this;
}

// ...........................................................................
// Matrix + Matrix
Matrix Matrix::operator+( const Matrix &m2 ) const
{
    Matrix mOut;

    // pruefen ob beide Matrizen die selbe Dimension haben
    if ( ( _rows == m2._rows ) && ( _cols == m2._cols ) )
    {
        // Matrix zuweisen und zweite addieren
        mOut = *this;
        mOut += m2;
    }
    else
    {
        stringstream errormessage;
        errormessage << " Matrix Matrix::operator+( const Matrix &m2 ) const: "
                     "Matrix Dimensions do not agree: M1(" << _rows << "," << _cols
                     << ") <-> M2(" << m2._rows << "," << m2._cols << ")" << " Exiting";
        throw runtime_error( errormessage.str() );
    }
    return mOut;
}

// this += Matrix
Matrix & Matrix::operator+=( const Matrix &m2 )
{
    // pruefen ob beide Matrizen die selbe Dimension haben
    if ( ( _rows == m2._rows ) && ( _cols == m2._cols ) )
    {
        // level 1 BLAS - Addition zweier Vektoren (fuer Matrix moeglich, da
        // Elemente linear gespeichert)
        cblas_daxpy ( _rows*_cols , 1.0 , m2.data_ptr(), 1, data_ptr(), 1 );
    }
    else
    {
        stringstream errormessage;
        errormessage << " Matrix Matrix::operator+=( const Matrix &m2 ) const: "
                     "Matrix Dimensions do not agree: M1(" << _rows << "," << _cols
                     << ") <-> M2(" << m2._rows << "," << m2._cols << ")" << " Exiting";
        throw runtime_error( errormessage.str() );
    }
    return *this;
}

// ...........................................................................
// Subtraktion eines Skalars
Matrix Matrix::operator-( double v ) const
{
    Matrix mOut(rows(),cols());
    transform (_data.begin() ,_data.end() ,mOut._data.begin() ,
               bind2nd( minus < double >() ,v));

    return mOut;
}

// ...........................................................................
// Matrix - Matrix
Matrix Matrix::operator-( const Matrix &m2 ) const
{
    Matrix mOut;

    // pruefen ob beide Matrizen die selbe Dimension haben
    if ( ( _rows == m2._rows ) && ( _cols == m2._cols ) )
    {
        // Matrix zuweisen und zweite subrathieren
        mOut = *this;
        mOut -= m2;
    }
    else
    {
        stringstream errormessage;
        errormessage << " Matrix Matrix::operator-( const Matrix &m2 ) const: "
                     "Matrix Dimensions do not agree: M1(" << _rows << "," << _cols
                     << ") <-> M2(" << m2._rows << "," << m2._cols << ")" << " Exiting";
        throw runtime_error( errormessage.str() );
    }
    return mOut;
}

// this -= Matrix
Matrix & Matrix::operator-=( const Matrix &m2 )
{
    // pruefen ob beide Matrizen die selbe Dimension haben
    if ( ( _rows == m2._rows ) && ( _cols == m2._cols ) )
    {
        // level 1 BLAS - Addition zweier Vektoren (fuer Matrix moeglich, da
        // Elemente linear gespeichert) mit Multiplikator -1 => Subtraktion
        cblas_daxpy ( _rows*_cols , -1.0 , m2.data_ptr(), 1, data_ptr(), 1 );
    }
    else
    {
        stringstream errormessage;
        errormessage << " Matrix Matrix::operator-=( const Matrix &m2 ) const: "
                     "Matrix Dimensions do not agree: M1(" << _rows << "," << _cols
                     << ") <-> M2(" << m2._rows << "," << m2._cols << ")" << " Exiting";
        throw runtime_error( errormessage.str() );
    }

    return *this;
}

// Addition eines Skalars
// this -= Matrix
Matrix Matrix::operator-=( double v )
{
    transform (_data.begin() ,_data.end() , _data.begin() ,
               bind2nd( minus < double >() ,v));
    return *this;
}

// ...........................................................................
//                                POWER
// ...........................................................................
Matrix Matrix::operator^( double v ) const
{
    return pow( v );
    /*   Matrix mOut = *this;
       mOut._data.erase( mOut.begin(), mOut.end() );
       transform( _data.begin(), _data.end(), back_inserter( mOut._data ),
                  bind2nd( ptr_fun((double (*)(double, double))std::pow),v ) );

       return mOut;
       */
}
Matrix Matrix::operator^( int v ) const
{
    v = double( v );

    return pow( v );
}

Matrix Matrix::pow( double v ) const
{
    Matrix mOut = *this;
    mOut._data.erase( mOut.begin(), mOut.end() );
    transform( _data.begin(), _data.end(), back_inserter( mOut._data ),
               bind2nd( ptr_fun((double (*)(double, double))std::pow),v ) );

    return mOut;
}


// ===========================================================================
//                                 METHODEN
// ===========================================================================

// ...........................................................................
// ein Element n der Matrix ablegen
void Matrix::_set(const int &i,const  int &j,const double &v )
{
    // Test, ob gueltiger Bereich angefordert wurde
    if ( (i < _rows ) && ( j < _cols ) )
    {
        // Ja => Wert setzen
        _data.at(i+j*_rows) = v;

#ifdef Matrix_DEBUG
        cerr << "M(" << i << ", " << j << ") = " << v  << endl;
#endif
    }
    else
    {
        // Ungueltig => Fehlermeldung ausgeben
        stringstream errormessage;
        errormessage << " void Matrix::_set(int i, int j, double v ): "
                     << "ERROR: Index (" << i << ", " << j
                     << ") exceeds Matrix dimensions (" << _rows << ", " << _cols
                     << ") Exiting!";
        throw runtime_error( errormessage.str() );
    }
}

// ...........................................................................
// ein Element n der Matrix ablegen
void Matrix::_set(const int &i, const double &v )
{
    // Test, ob gueltiger Bereich angefordert wurde
    if ( (i < _rows * _cols ) )
    {
        // Ja => Wert setzen
        _data.at(i) = v;

#ifdef Matrix_DEBUG
        cerr << "M(" << i << ", " << j << ") = " << v  << endl;
#endif
    }
    else
    {
        // Ungueltig => Fehlermeldung ausgeben
        stringstream errormessage;
        errormessage << " void Matrix::_set(int i, int j, double v ): "
                     << "ERROR: Index (" << i << ", "
                     << ") exceeds Matrix dimensions (" << _rows << ", " << _cols
                     << ") Exiting!";
        throw runtime_error( errormessage.str() );
    }
}

// ...........................................................................
// Dimension aendern und Klassenvariablen anpassen
void Matrix::resize(const int &i,const int &j,const double &v )
{
    _rows = i;
    _cols = j;
    _data.resize( _rows * _cols, v );
}


// ...........................................................................
// ein Element der Matrix zurueckgeben
double Matrix::_get(const int &i,const int &j) const
{
    // Test, ob gueltiger Bereich angefordert wurde
    if ( ( i < _rows ) && ( j < _cols ) )
    {
        // Ja => Wert zurueckgeben
        return _data.at(i+j*_rows);
    }
    else
    {
        // Ungueltig => Fehlermeldung ausgeben
        stringstream errormessage;
        errormessage << " void Matrix::_get(int i , int j): "
                     << "ERROR: Index (" << i << ", " << j
                     << ") exceeds Matrix dimensions (" << _rows << ", " << _cols
                     << ") Exiting!";
        throw runtime_error( errormessage.str() );
    }
}

// ...........................................................................
// Matrix formatiert ausgeben
void Matrix::show(const int &nk ) const
{

    // Iterator, um ueber Zeilenanfaenge zu iterieren
    vector<double>::const_iterator cIter;
    for (cIter = _data.begin(); cIter < _data.begin()+_rows;
            cIter++)
    {
        // Iterator ueber Elemente einer Zeile
        vector<double>::const_iterator rIter;
        if( cIter == _data.begin() )
            cout << endl << "   [ " << setprecision(nk);
        else
            cout << endl << "     " << setprecision(nk);
        for (rIter = cIter; rIter < _data.end()-_rows;
                rIter += _rows)
        {
            cout << *rIter << ", ";
        }
        cout << *rIter << ";";
    }
    // Ausgabe der Dimension
    cout << "];_[" << _rows << ", " << _cols << "]" << endl;
}

// ...........................................................................
// groesstes und kleinstes Element der Matrix zurueckgeben
double Matrix::max() const
{
    return *max_element(_data.begin(), _data.end());
}


double Matrix::min() const
{
    return *min_element(_data.begin(), _data.end());
}


double Matrix::max( int & idx )
{
    vector<double>::const_iterator idxIt = max_element( _data.begin(),
                                           _data.end() );
    idx = idxIt - _data.begin();
    return *max_element(_data.begin(), _data.end());
}


double Matrix::min( int & idx )
{
    vector<double>::const_iterator idxIt = min_element( _data.begin(),
                                           _data.end() );
    idx = idxIt - _data.begin();
    return *min_element(_data.begin(), _data.end());
}


double Matrix::meanD() const
{
    double sum = accumulate( _data.begin(), _data.end(), 0.0 );
    double average = ( sum / _data.size() );

    return average;
}


double Matrix::median() const
{
    vector<double> dat( _data.size() );
    copy ( _data.begin(), _data.end(), dat.begin() );
    size_t n = dat.size() / 2;
    nth_element( dat.begin(), dat.begin()+n, dat.end() );

    double med = dat.at( n );
    if( n%2 != 0 )
    {
        --n;
        nth_element( dat.begin(), dat.begin()+n, dat.end() );
        med = .5*( med+dat.at( n ) );
    }

    return med;
}


int Matrix::minIdx() const
{
    vector<double>::const_iterator idxIt = min_element( _data.begin(),
                                           _data.end() );
    int idx = idxIt - _data.begin();
    return idx;
}


int Matrix::maxIdx() const
{
    vector<double>::const_iterator idxIt = max_element( _data.begin(),
                                           _data.end() );
    int idx = idxIt - _data.begin();
    return idx;
}


double Matrix::stdD() const
{
    Matrix tmp = *this;

    double ave = tmp.meanD();
    tmp = tmp - ave;

    double x = std::inner_product( tmp.begin(), tmp.end(), tmp.begin(), 0.0 );
    double sigma = std::sqrt( x / ( _data.size() - 1 ) );

    return sigma;
}

Matrix Matrix::sum_col() const
{
    Matrix sum( 1, _cols, 0.0 );
    vector<double>::const_iterator iter = _data.begin();

    int i = 0;
    while( iter < _data.end()-1 )
    {
        sum( i ) = accumulate( iter, iter+_rows, 0.0 );

        iter += _rows;
        i++;
    }

    return sum;
}


// +++ edit SH: 2013-11-20
Matrix Matrix::prod_col()
{
    Matrix tmp = ( *this ) ;
    Matrix prod( 1, tmp.size(2), 0.0 );
    double pro = 1.0 ;

    for( int i = 0; i <= tmp.size(2)-1; i++ )
    {
        for( int j = 0; j <= tmp.size(1)-1; j++ )
        {
            pro *= tmp( j,i );
        }

        prod( 0,i ) = pro ;
        pro = 1.0 ;
    }

    return prod;
}
// --- edit SH: 2013-11-20



// ...........................................................................
// Matrix mit Skalar multiplizieren
void Matrix::_mult(const double &v)
{
    transform(_data.begin(), _data.end(), _data.begin(),
              bind2nd(multiplies<double>(), v));
}

// ...........................................................................
// Element by element division
Matrix Matrix::div_elem( const Matrix & m2 ) const
{
    if( _rows != m2._rows || _cols != m2._cols )
    {
        stringstream errormessage;
        errormessage << "Matrix Matrix::div_elem( const Matrix & m2 ) const: "
                     << "dimension missmatch. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    Matrix out( _rows, _cols, 0.0 );
    for( int i=0; i<_rows*_cols; i++ )
        out( i ) = _data.at( i ) / m2( i );

    return out;
}

// ...........................................................................
// Element by element multiplication
Matrix Matrix::mult_elem( const Matrix & m2 ) const
{
    if( _rows != m2._rows || _cols != m2._cols )
    {
        stringstream errormessage;
        errormessage << "Matrix Matrix::mult_elem( const Matrix & m2 ) const: "
                     << "dimension missmatch. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    Matrix out( _rows, _cols, 0.0 );
    for( int i=0; i<_rows*_cols; i++ )
        out( i ) = _data.at( i ) * m2( i );

    return out;
}

void Matrix::transpose(Matrix &mOut) const
{
    
    if(_cols == 1 || _rows == 1){
        mOut._data = _data;
        mOut._cols = _rows;
        mOut._rows = _cols;
    } else {
        
        mOut.resize (_cols, _rows);
        
        // schleife ueber alle Matrixelemente, auslesen des entsprechenden Wertes
        // und mit verteuschten Zeilen- und Spaltennummer an mOut zuweisen
        for (int r = 0; r < _rows; r++)
        {
            for (int c = 0; c < _cols; c++)
            {
                mOut._set(c,r, _get(r,c));
            }
        }
    }
    return;
}

// ...........................................................................
// Matrix transponieren
Matrix Matrix::transpose() const
{
    Matrix mOut;
    transpose(mOut);
    return mOut;
}

// ...........................................................................
//                        EXTERNE Matrix EINLESEN UND SCHREIBEN
//
// Matrix aus binaerer Datei einlesen
void Matrix::load_bin( const string &file )
{
    ifstream inStream ( file.c_str(), ios::in | ios::binary);

    if( !inStream.is_open() )
    {
        stringstream errormessage;
        errormessage << "void Matrix::load_bin( string file ) : "
                     << "Failed to open file: " << file << " Exiting";
        throw runtime_error( errormessage.str() );
    }

    inStream.read((char*)&_rows, sizeof(int));
    inStream.read((char*)&_cols, sizeof(int));

#ifdef Matrix_DEBUG
    cerr << "M(" << _rows << ", " << _cols << ")" << endl;
#endif

    _data.resize( _rows * _cols, numeric_limits<double>::quiet_NaN() );
    inStream.read( (char *)(&_data.at(0)), _rows*_cols*sizeof(_data.at(0)) );

    inStream.close();
}


// Matrix binaere schreiben
void Matrix::save_bin( const string &file ) const
{
    ofstream outStream ( file.c_str(), ios::out | ios::binary);

    if( !outStream.is_open() )
    {
        stringstream errormessage;
        errormessage << "void Matrix::save_bin( string file ) : "
                     << "Failed to open file: " << file << " Exiting";
        throw runtime_error( errormessage.str() );
    }

    outStream.write((char*)&_rows, sizeof(int));
    outStream.write((char*)&_cols, sizeof(int));

    outStream.write((char*) &_data.at(0), _rows*_cols*sizeof(_data.at(0)));

    outStream.close();
}



// Matrix aus ASCII Datei einlesen
void Matrix::load_ascii( const string &file )
{
    ifstream inStream ( file.c_str(), ios::in);

    if( !inStream.is_open() )
    {
        stringstream errormessage;
        errormessage << "void Matrix::load_ascii( string file ) : "
                     << "Failed to open file: " << file << " Exiting";
        throw runtime_error( errormessage.str() );
    }

    // read Dimension
    inStream >> _rows >> _cols;

#ifdef Matrix_DEBUG
    cerr << "M(" << _rows << ", " << _cols << ")" << endl;
#endif

    // read Matrixelements
    int i = 0;
    _data.resize( _rows * _cols, numeric_limits<double>::quiet_NaN() );

    while( !inStream.eof() && ( i < _rows*_cols ) )
    {
        inStream >> _data.at(i);
        i++;
    }

    inStream.close();
}

// Matrix in ASCII Datei schreiben
void Matrix::save_ascii( const string &file ) const
{
    ofstream outStream ( file.c_str(), ios::out);

    if( !outStream.is_open() )
    {
        stringstream errormessage;
        errormessage << "void Matrix::save_ascii( string file ) : "
                     << "Failed to open file: " << file << " Exiting";
        throw runtime_error( errormessage.str() );
    }

    // zunaechst Anzahl der Zeilen und Spalten schreiben
    outStream << _rows << endl;
    outStream << _cols << endl;

    // nun die Daten (jeweils in neuer Zeile)
    outStream << setprecision(30);
    copy( _data.begin(),_data.end(),ostream_iterator<double>( outStream, "\n") );

    outStream.close();
}

// Matrix in ASCII Datei schreiben
void Matrix::save_ascii_idx(const string &file) const
{
    ofstream outStream ( file.c_str(), ios::out);

    if( !outStream.is_open() )
    {
        stringstream errormessage;
        errormessage << "void Matrix::save_ascii( string file ) : "
                     << "Failed to open file: " << file << " Exiting";
        throw runtime_error( errormessage.str() );
    }

    outStream << setprecision(16);
    for( int c=0; c<_cols; c++ )
    {
        outStream << endl;
        for( int r=0; r<_rows; r++ )
        {
            outStream << r << " " << c << " " << _get( r,c ) << endl;
        }
    }

    outStream.close();
}

// Matrix in ASCII Datei schreiben ( in Matrix-Format )
void Matrix::save_matrix( const string &file,const int &nk) const
{
    ofstream outStream ( file.c_str(), ios::out);

    if( !outStream.is_open() )
    {
        stringstream errormessage;
        errormessage << "void Matrix::save_ascii( string file ) : "
                     << "Failed to open file: " << file << " Exiting";
        throw runtime_error( errormessage.str() );
    }

    Matrix neu = transpose();

    outStream << setprecision(nk);
    for( int i=0; i<neu._rows*neu._cols; i++ )
    {
        if( i % neu._rows == 0 && i != 0 )
            outStream << endl;
        outStream << neu._data.at(i) << " ";
    }
    outStream.close();
}

// Matrix aus ASCII Datei einlesen
void Matrix::load_vec( const string &file )
{
    ifstream inStream ( file.c_str(), ios::in);

    if( !inStream.is_open() )
    {
        stringstream errormessage;
        errormessage << "void Matrix::load_vec( const string &file ) : "
                     << "Failed to open file: " << file << " Exiting";
        throw runtime_error( errormessage.str() );
    }

    int rows=0;
    string a;
    // einmal laenge durchlaufen
    while( !inStream.eof() )
    {
        getline(inStream,a) ;
        rows++;
    }
    inStream.close();
    // read Matrixelements
    resize( rows, 1);

    inStream.open(file.c_str(), ios::in);
    for (int i=0; i<rows; i++ )
    {
        inStream >> _data.at(i);
    }

    inStream.close();
}

// ...........................................................................
int Matrix::length( ) const
{
    // maximale Dimension
    int dim = _rows < _cols ? _cols : _rows;

    return dim;
}

// ...........................................................................
// Hauptdiagonale einer Matrix extrahieren oder DiagonalMatrix aus Vektor
// bauen
Matrix Matrix::diag() const
{
    Matrix mOut;
    diag(mOut);
    return mOut;
}
// ...........................................................................
//
void Matrix::from_diag_matrix(bool reciprocal)
{
    
    if( _cols == _rows)
    {
        vector<double> new_data;
        new_data.resize(_rows,0.0);
        
        if(reciprocal)
            for(int i=0; i<_rows; i++)
                new_data.at(i) = 1.0/_data.at(i * _rows + i);
        else
            for(int i=0; i<_rows; i++)
                new_data.at(i) = _data.at(i * _rows + i);
    
        _data = new_data;
        _cols = 1;        
    }
    else
        throw runtime_error("void Matrix::from_diag_matrix(): Unexpected dimensions. " );
        
}
// ...........................................................................
//
void Matrix::to_diag_matrix()
{
    if( _cols == 1 )
    {
        vector<double> new_data;
        new_data.resize(_rows*_rows,0.0);
        
        for(int i=0; i<_rows; i++)
            new_data.at(i * _rows + i) = _data.at(i);
    
        _data = new_data;
        _cols = _rows;
    }
    else if( _rows == 1 )
    {
        vector<double> new_data;
        new_data.resize(_cols*_cols,0.0);
        _rows = _cols;
        
        for(int i=0; i<_cols; i++)
            new_data.at(i * _rows + i) = _data.at(i);
    
        _data = new_data;
    }
    else
        throw runtime_error("void Matrix::to_diag_matrix(): Unexpected dimensions. " );
}

// ...........................................................................
// Hauptdiagonale einer Matrix extrahieren oder DiagonalMatrix aus Vektor
// bauen
void Matrix::diag(Matrix &mOut ) const
{

    // vector -> Matrix
    if (_cols == 1 || _rows == 1)
    {
        // maximale Dimension
        int dim = _rows < _cols ? _cols : _rows;
        mOut.resize( dim, dim,  0.0 );

        for( int i = 0; i < dim; i++ )
            mOut( i, i ) = operator()( i );
        return;
    }

    // Matrix -> vektor
    else
    {
        int dim = _rows < _cols ? _rows : _cols;

        mOut.resize( dim, 1,  0.0 );

        for( int i = 0; i < dim; i++ )
            mOut( i ) = operator()( i, i );

        return;
    }
}

// ...........................................................................
// EinheitsMatrix erstellen
void Matrix::eye(int dim)
{
    Matrix tmp(dim,1,1.0);
    *this = tmp.diag();
}

// ...........................................................................
// Matrix mit (Pseudo-)Zufallszahlen fuellen
void Matrix::random( )
{
    // initialisieren rand mit der aktuellen Zeit
    srand( (unsigned)time(0) );

    // Schleife ueber den Datenbereich und mit Zufallszahlen zwischen 0 und 1
    // fuellen
    for( int i = 0; i < _data.size(); i++)
    {
        _data.at(i) = (rand() / ((double)RAND_MAX));
    }
}


// ...........................................................................
// Gleichungssystem (in Dreiecksgestalt = this) mit mehreren rechten Seiten
// (rhs) lösen  => rhs = this^-1 * rhs
void Matrix::solve_tri( Matrix & rhs ) const
{
    if( _rows != rhs.rows() )
    {
        stringstream errormessage;
        errormessage << "Matrix Matrix::solve_tri( const Matrix &rhs ) const: "
                     << "Matrix dimensions do not agree. Exiting!";
        throw runtime_error( errormessage.str() );
    }
    cblas_dtrsm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans,
                 CblasNonUnit, _rows, rhs._cols, 1.0, data_ptr(), _rows,
                 rhs.data_ptr(), rhs._rows );
}

// ...........................................................................
// das symmetrische Produkt this = A^T * A bilden
void Matrix::is_AtA_of( const Matrix & A )
{
    resize( A._cols, A._cols, 0.0 );
    cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, A._cols,
                 A._rows, 1.0, A.data_ptr(), A._rows, 0.0, data_ptr(), _rows );

    // unteres Dreieck fuellen
    for( int j = 0; j < _cols; j++ )
        for( int i = j+1; i < _rows; i++ )
            _set( i, j, _get( j, i ) );
}

void Matrix::is_Atv( const Matrix & A, const Matrix & v )
{
    // test, ob v Vektor und Dimensionen passen
    if( A._rows != v.rows() || v.cols() != 1)
    {
        stringstream errormessage;
        errormessage << "void Matrix::is_Atv( const Matrix & A, const Matrix & v ): "
                     << "dimension missmatch. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    resize( A.cols(), 1, 0.0 );

    cblas_dgemv( CblasColMajor, CblasTrans, A._rows, A._cols, 1.0, A.data_ptr(),
                 A._rows, v.data_ptr(), 1, 0.0, data_ptr(), 1 );

}

void Matrix::is_Atv2( const Matrix & A, const Matrix & v )
{
    // test, ob v Vektor und Dimensionen passen
    if( A._rows != v.rows() || v.cols() != 1)
    {
        stringstream errormessage;
        errormessage << "void Matrix::is_Atv2( const Matrix & A, const Matrix & v ): "
                     << "dimension missmatch. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    resize( v.rows(), 1, 0.0 );

    cblas_dgemv( CblasColMajor, CblasTrans, A._cols, A._rows, 0.0, A.data_ptr(),
                 A._rows, v.data_ptr(), 1, 0.0, data_ptr(), 1 );
}

// ...........................................................................
// Cholesky Zerlegung A = R' * R rechnen und R in this speichern
// A muss positiv definit sein
void Matrix::chol( const Matrix & A, const char type )
{
    *this = A;
    int stat;
    dpotrf_( type, _rows, data_ptr(), _rows, stat );

    switch(type){
        case 'U':{
            // unteres Dreieck mit Nullen fuellen
            for( int j = 0; j < _cols; j++ )
                for( int i = j+1; i < _rows; i++ )
                    _set( i, j, 0.0 );
            break;
        }
        case 'L':{
            // oberes Dreieck mit Nullen fuellen
            for( int j = 1; j < _cols; j++ )
                for( int i = 0; i < j; i++ )
                    _set( i, j, 0.0 );
            break;
        }
    }
    

    if( stat != 0 )
    {
        stringstream errormessage;
        errormessage << "Matrix::chol: LAPACK error: " << stat << endl;
        throw runtime_error( errormessage.str() );
    }
}

// ...........................................................................
// loese ein positiv definites Gleichungssystem
// this beinhaltet die rechte(n) Seite(n) und wird mit Loesung(en)
// ueberschrieben
void Matrix::solve_neq( const Matrix & N )
{
    if( _rows != N.rows() )
    {
        stringstream errormessage;
        errormessage << "Matrix Matrix::solve_neq( const Matrix & N ): "
                     << "Matrix dimensions do not agree. Exiting!";
        throw runtime_error( errormessage.str() );
    }
    Matrix tmp = N;

    int stat;
    dposv_( 'U', tmp._rows, _cols, tmp.data_ptr(), tmp._rows, data_ptr(), _rows, stat );

    if( stat != 0 )
    {
        stringstream errormessage;
        errormessage << "Matrix::solve_neq: LAPACK error: " << stat;
        throw runtime_error( errormessage.str() );
    }

}


// ...........................................................................
// solve GENERAL equation system
// this beinhaltet die rechte(n) Seite(n) und wird mit Loesung(en)
// ueberschrieben
void Matrix::solve_geeqs( const Matrix & A )
{
    if( _rows != A.rows() )
    {
        stringstream errormessage;
        errormessage << "Matrix Matrix::solve_geeqs( const Matrix & A ): "
                     << "Matrix dimensions do not agree. Exiting!";
        throw runtime_error( errormessage.str() );
    }
    
    Matrix tmp = A;
    int *ipiv = new int[ _rows ];
    int stat;
    
    dgesv_( _rows, _cols,  tmp.data_ptr(), tmp._rows, ipiv, data_ptr(), _rows, stat);

    if( stat != 0 )
    {
        stringstream errormessage;
        errormessage << "Matrix::solve_geeqs: LAPACK error: " << stat << endl;
        throw runtime_error( errormessage.str() );
    }

}



// ...........................................................................
void Matrix::ols( const Matrix &A )
// ...........................................................................
{
    if( _rows != A.rows() )
    {
        stringstream errormessage;
        errormessage << "Matrix Matrix::ols( const Matrix & A ): "
                     << "Matrix dimensions do not agree. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    Matrix n = A.transpose() * (*this) ;
    Matrix N = A.transpose() * A ;

    n.solve_neq( N ) ;
    (*this) = n ;
}

//// ...........................................................................
//void Matrix::ols_qr( const Matrix &A )
//// ...........................................................................
//{
//    if( _rows != A.rows() )
//    {
//        stringstream errormessage;
//        errormessage << "Matrix Matrix::ols_qr( const Matrix & A ): "
//                     << "Matrix dimensions do not agree. Exiting!";
//        throw runtime_error( errormessage.str() );
//    }
//
//    ivg::Matrix tmp = A;
//
//    int stat;
//    
////    = clapack_dgels( CblasColMajor, CblasNoTrans, tmp._rows, tmp._cols, _cols,
////                              tmp.data_ptr(), tmp._rows, data_ptr(), _rows );
//    
//    dgels_('N', const int& M, const int& N, const int& NRHS, double* A,
//                 const int&LDA, double* B, const int& LDB, double* WORK, const int& LWORK, stat);
//
//    resize( A.cols(), _cols );
//}

// ...........................................................................
void Matrix::solve_lu( const Matrix &A )
// ...........................................................................
{
   ivg::Matrix X = A;
//   int *ipiv = new int[ _rows ];
  
//   // LU decomposition 
//   clapack_dgetrf( CblasColMajor, X._rows, X._cols, X.data_ptr(), X._rows, ipiv );
//   // solve equation system
//   clapack_dgetrs( CblasColMajor, CblasNoTrans, X._rows, _cols, X.data_ptr(),
//                   X._rows, ipiv, data_ptr(), _rows );

//   char trans = 'N';
//   long info = 0;
//   long r = X._rows;
//   long c = X._cols;
//   long rhs = _cols;
//   dgetrs_( &trans, &r, &rhs, X.data_ptr(), &r, ipiv, data_ptr(), &r, &info );   



   int n = _rows;

   int info = 0;
   int ione = 1;
   double *af = new double[ _rows*_rows ];
   double *r  = new double[ _rows ];
   double *c  = new double[ _rows ];
   double *work = new double[ 4*_rows ];
   int *ipiv = new int[ _rows ];
   int *iwork = new int[ _rows ];

   char fact = 'E';
   char trans = 'N';
   char equed;
   double rcond, ferr, berr;

   ivg::Matrix b = *this;

   dgesvx_( &fact, &trans, &n, &ione, X.data_ptr(), &n, af, &n, ipiv,
            &equed, r, c, b.data_ptr(), &n, data_ptr(), &n, &rcond, &ferr, &berr,
            work, iwork, &info );
   cerr << "rcond = " << rcond << endl;
}

// ...........................................................................
void Matrix::solve_qr( const Matrix &A, Matrix &Sxx )
// ...........................................................................
{
    if( _rows != A.rows() )
    {
        stringstream errormessage;
        errormessage << "void Matrix::solve_qr( const Matrix &A, Matrix &Sxx ): "
                     << "Matrix dimensions do not agree. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    ivg::Matrix R = A;


    // solve Q*R=A
    // ===========
    long r = A._rows;
    long c = A._cols;
    long rhs = _cols;
    int  *ldvt = new int[ A._cols ];
    long lwork = -1;
    double iwork;
    long info = 0;
    Matrix tau( _rows,1,0.0 );

    // Explore the optimal allocation of the WORK array by performing a dry run with LWORK=-1.
    dgeqrf_( &r, &c, R.data_ptr(), &r, tau.data_ptr(), &iwork, &lwork, &info );
//    dgeqp3_( &r, &c, R.data_ptr(), &r, ldvt, tau.data_ptr(), &iwork, &lwork, &info );
    R = A;
    
    // perform the final decomposition
    lwork = (long)iwork;
    double* work = new double[lwork];
    dgeqrf_( &r, &c, R.data_ptr(), &r, tau.data_ptr(), work, &lwork, &info );
//    dgeqp3_( &r, &c, R.data_ptr(), &r, ldvt, tau.data_ptr(), work, &lwork, &info );

    bool scale = false;
    ivg::Matrix S;
    if( scale )
    { 
       S = (R^2).sum_col();
       S = S^(-0.5);
       S = S.diag();
    
       R = A*S;
       dgeqrf_( &r, &c, R.data_ptr(), &r, tau.data_ptr(), work, &lwork, &info );
    }

    // Multiply by Q^T to get Q^T*b as in R*x=Q^T*b
    // ============================================
    info = 0;
    lwork = -1;
    char side = 'L';
    char trans = 'T';
    dormqr_( &side, &trans, &r, &rhs, &c, R.data_ptr(), &r, tau.data_ptr(), 
             data_ptr(), &r, &iwork, &lwork, &info );
    lwork = (long)iwork;
    work = new double[lwork];
    dormqr_( &side, &trans, &r, &rhs, &c, R.data_ptr(), &r, tau.data_ptr(), 
             data_ptr(), &r, work, &lwork, &info );

    // ----------------------------------------------------------------------------------------
/*    // classical Gram-Schmidt orthogonalization   (FUNKTIONIERT NICHT!!!)
    ivg::Matrix Q( A._rows,A._rows,0.0 );
    ivg::Matrix T( A._cols,A._cols,0.0 );
    R = T;
    for( int i=0;i<A._cols;++i )
    {
       Q.set_sub( 0,i,A( ":",i ) ); 
       for( int j=0;j<i;++j )
       {
          //R(j,i) = (Q(":",j).transpose()*A(":",i))(0);
          R(j,i) = (Q(":",j).transpose()*Q(":",i))(0);
          Q.set_sub( 0,i, (Q(":",i)-Q(":",j)*R(j,i)) );
       }
       R(i,i) = (Q(":",i).norm())(0);
       Q.set_sub( 0,i, (Q(":",i)*1.0/R(i,i)) );
    }
    // calculate Q^T*b
    *this = Q.transpose()* *this;
*/
    // ----------------------------------------------------------------------------------------

    // Use back-substitution using the upper right part of A (which is still stored in R)
    // =====================
    info = 0;
    char uplo = 'U';
    char diag = 'N';
    trans = 'N';
    dtrtrs_( &uplo, &trans, &diag, &c, &rhs, R.data_ptr(), &r, data_ptr(), 
             &r, &info );

    // calculate VCM via R^-1*R^-T via back-substitution
    Sxx.eye( c );
    dtrtrs_( &uplo, &trans, &diag, &c, &c, R.data_ptr(), &c, Sxx.data_ptr(), 
             &c, &info );



    resize( A.cols(), _cols );
}

// ...........................................................................
void Matrix::solve_svd( const Matrix &A )
// ...........................................................................
{
    if( _rows != A.rows() )
    {
        stringstream errormessage;
        errormessage << "Matrix Matrix::solve_svd( const Matrix & A ): "
                     << "Matrix dimensions do not agree. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    Matrix S;
    Matrix U;
    Matrix VT;

    A.svd( U,S,VT );
    ivg::Matrix Sinv( VT._rows, U._rows,0.0 );
    for( int i=0;i<S._rows;++i )
       Sinv( i,i ) = 1.0/S(i);

    *this = VT.transpose()*Sinv*U.transpose()* *this;
}


// ...........................................................................
Matrix Matrix::chol_inv( ) const
// inversion ueber cholesky-zerlegung
{
    double testMin, testMax;
    {
        Matrix test = *this - transpose();
        testMin = test.min();
        testMax = test.max();
        testMax = std::abs( testMin ) <std:: abs( testMax ) ? std::abs(
                      testMin ) : std::abs( testMax );
    }
    if( testMax > 2e-16 )
    {
        stringstream errormessage;
        errormessage << "Matrix Matrix::chol_inv( ): "
                     << "Matrix has to be symmetric. Max diff: " << testMax << ". Exiting!";
        throw runtime_error( errormessage.str() );
    }

    // bilde EinheitsMatrix
    Matrix Qxx( _rows, 1, 1.0 );
    Qxx = Qxx.diag();

    // loese Gleichungssystem mit EinheitsMatrix als rechter Seite
    // => Inverse
    Qxx.solve_neq( *this );

    return Qxx;
}

// ...........................................................................
void Matrix::inv( )
{
    int *ipiv = new int[ _rows ];
    int info;
    dgetrf_( _rows, _rows, data_ptr(), _rows, ipiv, info);
    
    double* work = new double[ _rows ];
    
    dgetri_( _rows, data_ptr(), _rows, ipiv, work, _rows, info );
    
    delete[] ipiv;
    delete[] work;
    
//    clapack_dgetrf( CblasColMajor, _rows, _rows, data_ptr(), _rows, ipiv );
//    clapack_dgetri( CblasColMajor, _rows, data_ptr(), _rows, ipiv );
}

// ...........................................................................
void Matrix::inv_scal( )
{
    Matrix S = diag();
    for( int i=0; i<S.length(); i++ )
        S( i ) = 1.0 / std::sqrt( S( i ) );

    S = S.diag();
    *this = S **this * S;
    
    int *ipiv = new int[ _rows ];
    int info;
    dgetrf_( _rows, _rows, data_ptr(), _rows, ipiv, info);
    
    double* work = new double[ _rows ];
    
    dgetri_( _rows, data_ptr(), _rows, ipiv, work, _rows, info );
    
    *this = S **this * S;
    
    delete[] ipiv;
    delete[] work;

}
// ...........................................................................
Matrix Matrix::pinv( )
// ...........................................................................
{
    Matrix S;
    Matrix U;
    Matrix VT;

    svd( U,S,VT );
    
    double tol = std::max (_rows, _cols) * S.max() * 1e-16;
    int idx = -1;
    while( ++idx<S._rows && S(idx)>tol )
       S(idx) = 1.0/S(idx);
    ivg::Matrix P = VT.transpose()*S.diag()*U.transpose();
    
    return P;
}
// ...........................................................................
// Singulaerwertzerlegung A = USV'
// A:  m x n, U:  m x m, S:  m x n, V': n x n
void Matrix::svd( Matrix & U, Matrix & S, Matrix & VT, int n ) const
{
    Matrix A;
    if( n == -999 )
    {
        n = _rows;
        A = *this;
    }
    else
        A = get_sub( 0,0,n-1,_cols-1 );
    char job  = 'A';
    long r = long(n);
    long c = _cols;
    // LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
    long tmp = 1;
    long lwork = std::max( tmp, 3*std::min(r,c) + std::max(r,c));
    lwork = std::max(lwork, 5*std::min(r,c)) + 100;
    double work[lwork];
    long info = 0;

    size_t nr = n;
    S.resize(std::min(nr,_cols), 1, 0.0);
    U.resize(n, n);
    VT.resize(_cols, _cols);

    dgesvd_( &job, &job, &r, &c, A.data_ptr(), &r, S.data_ptr(), U.data_ptr(), &r,
             VT.data_ptr(), &c, work, &lwork, &info);

    if( info != 0 )
    {
        stringstream errormessage;
        errormessage << "Matrix::svd: LAPACK error: " << info << endl;
        throw runtime_error( errormessage.str() );
    }
}
// ...........................................................................
ivg::Matrix Matrix::calc_impact_factors()
// ...........................................................................
{
   // Singular Value Decomposition
   ivg::Matrix  U,S,VT;
   (*this).svd( U, S, VT );
   // Data resolution matrix
   ivg::Matrix Ur = U.get_sub( 0, 0, U.rows()-1, S.rows()-1 );
   ivg::Matrix DRM = Ur * Ur.transpose();
   // Impact factors
   DRM.from_diag_matrix();
   
   return DRM;
}
// ...........................................................................
// QR-Zerlegung A = Q*R
void Matrix::qr( Matrix & Q, Matrix & R ) const
{
    R = *this;

    long r = _rows;
    long c = _cols;
    // LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
    long tmp = 1;
    long lwork = std::max( tmp, 3*std::min(r,c) + std::max(r,c));
    lwork = std::max(lwork, 5*std::min(r,c)) + 100;
    double work[lwork];
    long info = 0;
    Matrix tau( _rows,1,0.0 );

    dgeqrf_( &r, &c, R.data_ptr(), &r, tau.data_ptr(), work, &lwork, &info );

    // final (decomposed) Matrix R is stored in the upper triangular and the
    // diagonal entries of R, the Matrix Q has to be determined from the
    // lower triangular and the entries of tau
    long k = std::min(r,c);
    Q.resize( std::max(_rows,_cols),std::max(_rows,_cols),0.0 );
    copy( R.begin(), R.end(), Q.begin() );
    dorgqr_( &r, &r, &k, Q.data_ptr(), &r, tau.data_ptr(), work, &lwork, &info );
    Q.resize( _rows,_rows );

    for( int i=0; i<_cols; i++ )
    {
        for( int j=i+1; j<_rows; j++ )
            R(j,i) = 0.0;
    }

    if( info != 0 )
    {
        stringstream errormessage;
        errormessage << "Matrix::qr: LAPACK error: " << info << endl;
        throw runtime_error( errormessage.str() );
    }
}



// ...........................................................................
// Elemente loeschen
// (1) Spalte
void Matrix::rem_c( int idx, int n )
{
    _data.erase( _data.begin() + ( idx * _rows ),
                 _data.begin() + ( idx * _rows + n*_rows )
               );
    _cols-=n;
}


// (2) Zeile
void Matrix::rem_r( int idx, int n )
{
    // iterator on first element in last column that should be deleted
    // (i.e., last element in first row that should be deleted)
    vector<double>::iterator iterBeg = _data.end() - ( _rows - idx );
    // iterator pointing behind last element that should be deleted
    vector<double>::iterator iterEnd = iterBeg + n;

    for( iterBeg; iterBeg >= _data.begin(); iterBeg -= _rows, iterEnd -= _rows )
        _data.erase( iterBeg, iterEnd );

    _rows-=n;
}

// (3) Zeile und Spalte
void Matrix::rem_rc( int idx )
{
    rem_r( idx );
    rem_c( idx );
}

// ...........................................................................
// ...........................................................................
// Removes Elements
// ...........................................................................
void Matrix::rem_c( vector<int> idxs )
// ...........................................................................
{
    vector<double> new_data;
    int elements = (_cols * _rows) - (idxs.size()*_rows);
    new_data.reserve(elements);
    
    for(int col=0; col<_cols; col++){
        if(find(idxs.begin(), idxs.end(), col) == idxs.end())
            new_data.insert(new_data.end(), _data.begin() + ( col * _rows ), _data.begin() + ( col * _rows + _rows) );
    }
    
    _data = new_data;
    _cols = _cols - idxs.size();
}
// ...........................................................................
void Matrix::rem_r( vector<int> idxs )
// ...........................................................................
{
    // really really slow! who will fix it?!
    // TODO seriously, who did that?  
    *this=transpose();
    rem_c(idxs);
    *this=transpose();
}

// ...........................................................................
Matrix Matrix::get_sub(const vector<int> &rows,
                       const  vector<int> &cols ) const
{
    // dimensionen der indizes
    int nRows = rows.end() - rows.begin();
    int nCols = cols.end() - cols.begin();

    // dimension check
    int maxidxR=*max_element(rows.begin(), rows.end());
    int maxidxC=*max_element(cols.begin(), cols.end());

    //cout << "max idx=" << maxidxR << ","<< maxidxC <<endl;
    //cout << "max mat idx=" << _rows - 1 << ","<< _cols-1 <<endl;

    if( maxidxR > (_rows-1) || maxidxC > (_cols-1))
    {
        stringstream errormessage;
        errormessage <<
                     "Matrix Matrix::get_sub(const vector<int> &rows,const  vector<int> &cols ) const: "
                     << "max index exceeds Matrix dimensions (" << maxidxR << ", "
                     << maxidxC << ") ." << "max possible indices: "<< _rows- 1<< "," <<_cols-1
                     <<". Exiting!";
        throw runtime_error( errormessage.str() );
    }

    Matrix out( nRows, nCols, 0*_get(0,0) );
    
    for( int c=0; c < nCols; c++ )
      for( int r=0; r < nRows; r++ )
	out( r,c ) = _get( rows.at( r ), cols.at( c ) );

    return out;
}

// ...........................................................................
vector<double> Matrix::get_data_vec( )
{
    vector<double> out( _data.size() );
    copy( _data.begin(), _data.end(), out.begin() );

    return out;
}

// ...........................................................................
Matrix Matrix::get_rows(const vector<int> &rows ) const
{
    // dimensionen
    int nRows = rows.end() - rows.begin();
    int nCols = _cols;

    Matrix out( nRows, nCols, 0 );

    for( int r=0; r < nRows; r++ )
        for( int c=0; c < nCols; c++ )
            out( r,c ) = _get( rows.at( r ), c );

    return out;
}

// ...........................................................................
Matrix Matrix::get_rows(int start, int end ) const
{
    // (startR, startC, endR, endC )
    return get_sub(start, 0, end, _cols-1 );
}

// ...........................................................................
Matrix Matrix::get_cols(const vector<int> &cols ) const
{
    // dimensionen
    int nCols = cols.end() - cols.begin();
    int nRows = _rows;

    Matrix out( nRows, nCols, 0 );

    for( int r=0; r < nRows; r++ )
        for( int c=0; c < nCols; c++ )
            out( r,c ) = _get( r, cols.at( c ) );

    return out;
}

// ...........................................................................
Matrix Matrix::get_cols(int start, int end ) const
{
    // (startR, startC, endR, endC )
    return get_sub(0, start, _rows-1, end );
}

// ...........................................................................
Matrix Matrix::get_sub(int startR, int startC, int endR, int endC ) const
{
    // dimensionen
    int nRows = endR - startR + 1;
    int nCols = endC - startC + 1;

    vector<int> Rows( nRows, 0 );
    vector<int> Cols( nCols, 0 );

    for( int i=0; i<nRows; i++ )
        Rows.at( i ) = startR + i;
    for( int j=0; j<nCols; j++ )
        Cols.at( j ) = startC + j;

    Matrix Out = get_sub( Rows, Cols );

    return Out;
}

// ...........................................................................
// A([1,..],[1,3,5]) = ones(3);
void Matrix::set_sub(const vector<int> &rows, const vector<int> &cols,
                     const Matrix & m2 )
{
    int nRows = rows.end() - rows.begin();
    int nCols = cols.end() - cols.begin();

    if( nRows != m2._rows || nCols != m2._cols || nRows > _rows ||
            nCols > _cols )
    {
        stringstream errormessage;
        errormessage <<
                     "void Matrix::set_sub( vector<int> rows, vector<int> cols, Matrix & m2 ): "
                     << "dimension missmatch. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    for( int r=0; r < nRows; r++ )
        for( int c=0; c < nCols; c++ )
            _set( rows.at( r ), cols.at( c ), m2( r, c ) );
}



// ...........................................................................
void Matrix::add_sub( const vector<int> &rows, const vector<int> &cols,
                      const Matrix & m2  )
{
    int nRows = rows.end() - rows.begin();
    int nCols = cols.end() - cols.begin();

    if( nRows != m2._rows || nCols != m2._cols || nRows > _rows ||
            nCols > _cols )
    {
        stringstream errormessage;
        errormessage <<
                     "void Matrix::add_sub( vector<int> rows, vector<int> cols, Matrix & m2 ): "
                     << "dimension missmatch. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    for( int r=0; r < nRows; r++ )
        for( int c=0; c < nCols; c++ )
            _set( rows.at( r ), cols.at( c ), m2( r, c )+_get( rows.at( r ),
                    cols.at( c ) ) );
}


// ...........................................................................
void Matrix::append_cols( const Matrix & m2 )
// m2 als neue Spalten am Ende von this einfuegen
{
    //check if this is already initialized
    if(_rows == 0 && _cols == 0)
    {
        _rows = m2._rows;
        
    }// werfe exception, wenn Anzahl der Zeilen ungleich
    else if( _rows != m2.rows() )
    {
        stringstream errormessage;
        errormessage << "void Matrix::append_cols( Matrix & m2 ): "
                     << "Matrix dimensions do not agree (" << _rows << ", "
                     << m2.rows() << "). Exiting!";
        throw runtime_error( errormessage.str() );
    }
    // this vergroeszern
    resize( _rows, _cols + m2._cols );

    // datenbereich von m2 ans Ende kopieren
    copy( m2._data.begin(), m2._data.end(), _data.end() - m2._rows*m2._cols );
}

// ...........................................................................
void Matrix::append_rows( const Matrix & m2 )
// m2 als neue Zeilen am Ende von this einfuegen
{
    // werfe exception, wenn Anzahl der Spalten ungleich
    if( _rows != 0 && _cols != 0 && _cols != m2.cols() )
    {
        stringstream errormessage;
        errormessage << "void Matrix::append_rows( Matrix & m2 ): "
                     << "Matrix dimensions do not agree. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    // neue Elemente spaltenweise einfuegen (Schleife ueber Spalten von m2)
    vector<double>::iterator iter1;
    int n = 1;
    for( vector<double>::const_iterator iter2 = m2._data.begin(); iter2 < m2._data.end(); iter2 += m2._rows, n++ )
    {
        // setze Iterator auf den naechsten Spaltenanfang
        iter1 = _data.begin() + n * _rows + (n-1) * m2._rows;

        _data.insert( iter1, iter2, iter2 + m2._rows );
    }

    if(_cols == 0)
        _cols = m2._cols;
    
    // Anzahl der Zeilen anpassen
    _rows = _rows + m2._rows;
}


// ...........................................................................
void Matrix::change_row( int start, int dest )
{
    // werfe exception, wenn Zieladresse groeszer als Dimension
    if( _rows < dest )
    {
        stringstream errormessage;
        errormessage << "void Matrix::change_row( int start, int dest ): "
                     << "destination row not within dimension. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    // geht das schneller, wenn man mit eine EinheitsMatrix erstellt, deren
    // Spalte start mit dest vertauscht und anschließend E * this
    // mit lapack rechnet? Entsprechend dann fuer changeCol bzw. changeRowAndCol
    double tmp;
    for( int i=0; i<_cols; i++ )
    {
        tmp = _get( start, i );
        _set( start, i, _get( dest, i ) );
        _set( dest, i, tmp );
    }
}

// ...........................................................................
// double-Wert an Vektor hinzufuegen
void Matrix::append_rows( double v )
// v am Ende von this einfuegen
{
    if( _cols == 0 && _rows == 0 )
    {
        _cols = 1;
    } // werfe exception, wenn kein Vektor
    else if( !( _cols == 1 || _rows == 1 ) )
    {
        stringstream errormessage;
        errormessage << "void Matrix::append_rows( double v ): "
                     << "Matrix dimensions do not agree. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    _data.push_back( v );

    if( _cols == 1 )
        _rows++;
    else
        _cols++;
}

// ...........................................................................
void Matrix::change_col( int start, int dest )
{
    // werfe exception, wenn Zieladresse groeszer als Dimension
    if( _cols < dest )
    {
        stringstream errormessage;
        errormessage << "void Matrix::change_row( int start, int dest ): "
                     << "destination row not within dimension. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    double tmp;
    for( int i=0; i<_rows; i++ )
    {
        tmp = _get( i, start );
        _set( i, start, _get( i, dest ) );
        _set( i, dest, tmp );
    }
}

// ...........................................................................
vector<double>::iterator Matrix::find_elem( vector<double>::iterator start,
        vector<double>::iterator end,
        double v  ) const
{
    vector<double>::iterator iter = find( start, end, v );

    return iter;
}

// Element finden mit Bool funktion
Matrix Matrix::find_elem(bool (*pfunc)(double i, double j),
                         const double val) const
{
    vector<int> idx = find_idx(pfunc,val);
    Matrix mOut=operator()(idx);
    return mOut;
}

// Element finden mit Bool funktion inkl. idx als Referenz
Matrix Matrix::find_elem(bool (*pfunc)(double i, double j),const double val,
                         vector<int> &idx ) const
{
    Matrix mOut;
    find_elem(pfunc,val, idx, mOut);
    return mOut;
}

// Element finden mit Bool funktion inkl. m2, idx als Referenz
void Matrix::find_elem(bool (*pfunc)(double i, double j),const double val,
                       vector<int> &idx,Matrix &m2 ) const
{
    idx = find_idx(pfunc,val);
    m2=operator()(idx);
    return;
}
// find_elem mit zwei Bedingungen
// Element finden mit Bool funktion inkl. m2, idx als Referenz
void Matrix::find_elem(bool (*pfunc1)(double i, double j),const double val1,
                       bool (*pfunc2)(double i, double j),const double val2,
                       vector<int> &idx, Matrix &m2 ) const
{
    find_elem(pfunc1,val1,idx,m2);
    Matrix idxM(idx);
    m2.find_elem(pfunc2,val2,idx,m2);
    // Index aktualisieren
    idxM=idxM(idx);
    idx=idxM.get_vec("int");
    return;
}


// find_elem mit zwei Bedingungen
// Element finden mit Bool funktion inkl.  idx als Referenz
Matrix Matrix::find_elem(bool (*pfunc1)(double i, double j),const double val1,
                         bool (*pfunc2)(double i, double j),const double val2,
                         vector<int> &idx ) const
{
    Matrix mOut;
    find_elem(pfunc1,val1,pfunc2,val2,idx,mOut);
    return mOut;
}

// find_elem mit zwei Bedingungen &&
// Element finden mit Bool funktion inkl.  idx als Referenz
Matrix Matrix::find_elem(bool (*pfunc1)(double i, double j),const double val1,
                         bool (*pfunc2)(double i, double j),const double val2) const
{
    Matrix mOut;
    vector<int> idx;
    find_elem(pfunc1,val1,pfunc2,val2,idx,mOut);
    return mOut;
}


// ...........................................................................
vector<int> Matrix::find_idx(double v ) const// ==
{
    vector<int> out;
    vector<double>::const_iterator iter = _data.begin();
    iter = find( iter, _data.end(), v );

    while( iter < _data.end() )
    {
        out.push_back( iter - _data.begin() );
        iter = find( iter+1, _data.end(), v );
    }

    return out;
}

// ...........................................................................
// ab i0, j0 wird sub Matrix geschrieben
void Matrix::set_sub( int i0, int j0, const Matrix & other )
{
    // teste Datenbereich
    if( i0 + other._rows > _rows || j0 + other._cols > _cols )
    {
        stringstream errormessage;
        errormessage <<
                     "void Matrix::set_sub( int i0, int j0, const Matrix & other ): "
                     << "dimension conflict. Exiting!";
        throw runtime_error( errormessage.str() );
    }

    for( int j=0; j<other.cols(); j++ )
    {
        copy( other._data.begin()+j*other.rows(),
              other._data.begin()+j*other.rows()+other._rows,
              _data.begin() + i0 + _rows * (j+j0) );
    }
}


// ...........................................................................
double Matrix::trace( ) const
{
    Matrix d = diag();
    double sum = 0;

    for( int i=0; i<d.length(); i++ )
        sum += d( i );

    return sum;
}

// Dimension der Matrix
int Matrix::size(int dim) const
{
    if (dim ==1)
    {
        return (_rows);
    }
    if (dim ==2)
    {
        return (_cols);
    }
    else
    {
        cerr << "Matrix existiert nicht size(A) = 0,0" << endl;
        int null = 0;
        return (null);
    }
}

Matrix Matrix::size() const
{
    Matrix S (2,1,0.0);
    S(0) = this->size(1);
    S(1) = this->size(2);
    return (S);
}

// ...........................................................................
// sortiere numerisch nach Spalte c
void Matrix::sort_cols( int c )
{
    int i, j, t;
    Matrix hilf( 1, _cols );

    for( i=_rows-1; i>0; i-- )
    {
        for( j=0; j<i; j++ )
        {
            if( _data.at( (c*_rows)+j ) > _data.at( (c*_rows)+j+1 ) )
            {
                // tauschen erforderlich
                for( t=0; t<_cols; t++)
                {
                    hilf( t ) = _data.at( (t*_rows)+j );
                    _data.at( (t*_rows)+j ) = _data.at( (t*_rows)+j+1 );
                    _data.at( (t*_rows)+j+1 ) = hilf( t );
                }
            }
        }
    }
}


// ...........................................................................
// sort a colimn in a matrix
void Matrix::sort_col( int col )
{

    Matrix A = *this;
    Matrix tmp = A.get_sub( col, 0, col, A.size(2)-1 );

    std::sort ( tmp._data.begin(), tmp._data.end() );

    A.set_sub ( col, 0, tmp );

    (*this) = A ;

}


// ...........................................................................
// sort all columns in a matrix
void Matrix::sort_cols()
{
    Matrix A = *this ;
    Matrix tmp ;

    for ( int i = 0; i < A.size(1)-1; i++ )
    {
        tmp = A.get_sub( i, 0, i, A.size(2)-1 );
        std::sort ( tmp._data.begin(), tmp._data.end() );
        A.set_sub ( i, 0, tmp );
    }

    (*this) = A ;

}


// ...........................................................................
void Matrix::cout_size( string name ) const
{
    cout << name << " (" << _rows << ", " << _cols << ")" << endl;
}

// ...........................................................................
void Matrix::cerr_size( string name ) const
{
    cerr << name << " (" << _rows << ", " << _cols << ")" << endl;
}

// ...........................................................................
// calculate condition based on singular value decomposition
double Matrix::cond() const
{
    Matrix S;
    Matrix U;
    Matrix VT;

    svd( U,S,VT );
    double cond = S.max() / S.min();

    return cond;
}

// ...........................................................................
// calculate rank based on singular value decomposition
int Matrix::rank( double tol )
{
    Matrix S;
    Matrix U;
    Matrix VT;

    svd( U,S,VT );

    if( tol == 0.0 )
        tol = std::max (_rows, _cols) * S.max() * 1e-16;

    cerr << tol << endl;

    int rnk = S.length();
    while( rnk > 0 && S(rnk-1) < tol )
        rnk--;

    return rnk;
}


//----------------------- spezielle Matrix Multiplikationen mit cblas-------------------
// this= A(') * B(')
void Matrix::is_product_of(const Matrix& A,const Matrix& B,
                           CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB)
{
    // Fallunterscheidungen von TransA und TransB
    // => sehr unpraktisch, weil M, N, K bei
    //    "CblasTrans" jeweils geaendert werden muss!
    int r_neu,c_neu, M, N, K;
    if (TransA== CblasNoTrans && TransB==CblasNoTrans)
    {
        r_neu=A.rows();
        c_neu=B.cols();
        M=A.rows();
        N=B.cols();
        K=A.cols();
    }
    else if (TransA== CblasTrans   && TransB==CblasTrans)
    {
        r_neu=A.cols();
        c_neu=B.rows();
        M=A.rows();
        N=B.cols();
        K=A.cols();
    }
    else if (TransA== CblasNoTrans && TransB==CblasTrans)
    {
        r_neu=A.rows();
        c_neu=B.rows();
        M=A.rows();
        N=B.rows();
        K=A.cols();
    }
    else if (TransA== CblasTrans   && TransB==CblasNoTrans)
    {
        r_neu=A.cols();
        c_neu=B.cols();
        K=A.rows();
        N=B.cols();
        M=A.cols();
    }
    resize(r_neu,c_neu);
    cblas_dgemm( CblasColMajor, TransA, TransB, M, N, K, 1.0,
                 A.data_ptr(), A.rows(), B.data_ptr(), B.rows(), 0.0, data_ptr() , _rows );
    return;
}

// this+ = A(') * B(') (wie isroductOf, nur BLAS aufruf 0.0=> 1.0)
void Matrix::plus_product_of(const Matrix & A,const Matrix & B,
                             CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB)
{
    // Fallunterscheidungen von TransA und TransB
    // => sehr unpraktisch, weil M, N, K bei
    //    "CblasTrans" jeweils geaendert werden muss!
    int r_neu,c_neu, M, N, K;
    if (TransA== CblasNoTrans && TransB==CblasNoTrans)
    {
        r_neu=A.rows();
        c_neu=B.cols();
        M=A.rows();
        N=B.cols();
        K=A.cols();
    }
    else if (TransA== CblasTrans   && TransB==CblasTrans)
    {
        r_neu=A.cols();
        c_neu=B.rows();
        M=A.rows();
        N=B.cols();
        K=A.cols();
    }
    else if (TransA== CblasNoTrans && TransB==CblasTrans)
    {
        r_neu=A.rows();
        c_neu=B.rows();
        M=A.rows();
        N=B.rows();
        K=A.cols();
    }
    else if (TransA== CblasTrans   && TransB==CblasNoTrans)
    {
        r_neu=A.cols();
        c_neu=B.cols();
        K=A.rows();
        N=B.cols();
        M=A.cols();
    }
    resize(r_neu,c_neu);
    cblas_dgemm( CblasColMajor, TransA, TransB, M, N, K, 1.0,
                 A.data_ptr(), A.rows(), B.data_ptr(), B.rows(), 1.0, data_ptr() , _rows );
    return;
}

// Matrix anzeigen auf andere Art
void Matrix::disp( ) const
{
    cout << setprecision(4);
    int idx_ij;
    if (_cols==0 && _rows==0)
    {
        cout << "leere Matrix!  " << endl;
    }
    else
    {
        for( int i = 0; i < _rows; ++i )
        {
            for( int j = 0; j < _cols; ++j )
            {
                // Ausgabe des Matrixelements A(i,j)
                idx_ij= i * 1 + j * _rows;
                cout << _data.at(idx_ij) << "\t";
            }
            cout << endl;
        }
    }
# ifdef DEBMODE
    cout << "Dimension: [" << _rows << "," << _cols << "]" << endl;
#endif
    return;
}

// Kopie einer Matrix i,j mal anhaengen, this= repm(m2,i,j), laeuft
void Matrix::repmat(const Matrix & m2,  int i, int j )
{
    *this=m2;
    for( int k = 1; k < j ; ++k )
    {
        append_cols(m2 );
    }
    Matrix m;
    m=*this;
    for( int l = 1; l <  i; ++l )
    {
        append_rows(m);
    }
    return;
}

 void Matrix::vec2Mat(int c){
     
    double rows = (double)_rows/(double)c;
    double fractpart, intpart;

    fractpart = std::modf (rows , &intpart);
    if(fractpart != 0.0){
        throw logic_error( "void Matrix::vec2Mat number of element must bot change" );
    } else {
        _cols = c;
        _rows = rows;
    }
 }


// gibt die Matrix spaltenweise als [numel,1] aus, , laeuft
void Matrix::vec()
{
    _rows  = _rows *_cols;
    _cols =1;
}
// als Vektorklasse
vector<double> Matrix::get_vec() const
{
    return _data;
}

// als Vektorklasse
vector<int> Matrix::get_vec(const string art) const
{
    if (art!="int")
    {
        stringstream errormessage;
        errormessage << "vector<int> Matrix::vec(string art) const: "
                     "ERROR: please set art= \"int\" "<<  " Exiting";
        throw logic_error( errormessage.str() );
    }
    vector<int> vOut(_data.begin(), _data.end());
    return vOut;
}

Matrix Matrix::sin() const
{
    Matrix out;
    sin(out);
    return out;
}

Matrix Matrix::cos() const
{
    Matrix out;
    cos(out);
    return out;
}

Matrix Matrix::tan() const
{
    Matrix out;
    tan(out);
    return out;
}

Matrix Matrix::exp() const
{
    Matrix out;
    exp(out);
    return out;
}

Matrix Matrix::log() const
{
    Matrix out;
    log(out);
    return out;
}

// Wurzel
Matrix Matrix::sqrt() const
{
    Matrix out= operator ^(0.5);
    return out;
}

//sin als referenz
void Matrix::sin(Matrix &m2) const
{
    m2.resize(rows(),cols());
    // sinus elementeweise, laeuft
    //transform (_data.begin() ,_data.end() , m2._data.begin() , std::ptr_fun(sinf)); ungenau, nicht verwenden!!!
    for( int i = 0; i < _data.size() ; ++i )
    {
        m2._set(i, std::sin(_data.at(i)));
    }
    return;
}

void Matrix::cos(Matrix &m2) const
{
    m2.resize(rows(),cols());
    // cosinus elementeweise, laeuft
    //transform (_data.begin() ,_data.end() ,m2._data.begin() , std::ptr_fun(cosf)); ungenau, nicht verwenden!!!
    for( int i = 0; i < _data.size() ; ++i )
    {
        m2._set(i, std::cos(_data.at(i)));
    }
    return;
}

Matrix Matrix::acos() const
{
    Matrix out;
    acos(out);
    return out;
}

void Matrix::acos(Matrix &m2) const
{
    m2.resize(rows(),cols());
    // cosinus elementeweise, laeuft
    //transform (_data.begin() ,_data.end() ,m2._data.begin() , std::ptr_fun(cosf)); ungenau, nicht verwenden!!!
    for( int i = 0; i < _data.size() ; ++i )
        m2._set(i, std::acos(_data.at(i)));
    return;
}

void Matrix::tan(Matrix &m2) const
{
    m2.resize(rows(),cols());
    // tangen elementeweise, laeuft
    //transform (_data.begin() ,_data.end() ,m2._data.begin() , std::ptr_fun(tanf)); ungenau, nicht verwenden!!!
    for( int i = 0; i < _data.size() ; ++i )
    {
        m2._set(i, std::tan(_data.at(i)));
    }
    return;
}

// Wurzel
void Matrix::sqrt(Matrix &m2) const
{
    m2= operator ^(0.5);
    return;
}

void Matrix::exp(Matrix &m2) const
{
    m2.resize(rows(),cols());
    // cosinus elementeweise, laeuft
    //transform (_data.begin() ,_data.end() ,m2._data.begin() , std::ptr_fun(cosf)); ungenau, nicht verwenden!!!
    for( int i = 0; i < _data.size() ; ++i )
    {
        m2._set(i, std::exp(_data.at(i)));
    }
    return;
}


void Matrix::log(Matrix &m2) const
{
    m2.resize(rows(),cols());
    // cosinus elementeweise, laeuft
    //transform (_data.begin() ,_data.end() ,m2._data.begin() , std::ptr_fun(cosf)); ungenau, nicht verwenden!!!
    for( int i = 0; i < _data.size() ; ++i )
    {
        m2._set(i, std::log(_data.at(i)));
    }
    return;
}


// atan2, this wird ueberschrieben
void Matrix::atan2(const Matrix &m1, const Matrix &m2)
{
    if( m1._rows != m1._rows && m1._cols != m1._cols )
    {
        stringstream errormessage;
        errormessage << "Matrix::atan2(Matrix &m1, Matrix &m2): "
                     << "Matrix dimensions do not agree. Exiting!";
        throw runtime_error( errormessage.str() );
    }
    resize(m1._rows,m1._cols);
    // atan2 elementeweise
    //transform (m1._data.begin() , m1._data.end() ,m2._data.begin(), _data.begin() ,std::ptr_fun(atan2f));

    for( int i = 0; i < _data.size() ; ++i )
    {
        _set(i, std::atan2(m1._data.at(i),m2._data.at(i)));
    }

    return;
}


// legt einen Indexvektor am Platz an, vgl. Matlab [start : schrittweite : ende]
void Matrix::indexvec(double start, double schrittweite,double ende)
{
    if( start > ende && schrittweite > 0 )
    {
        stringstream errormessage;
        errormessage <<
                     "Matrix::indexvec(double start, double schrittweite, double ende): "
                     << "start > ende! Exiting!";
        throw runtime_error( errormessage.str() );
    }
    if( start < ende && schrittweite < 0 )
    {
        stringstream errormessage;
        errormessage <<
                     "Matrix::indexvec(double start, double schrittweite, double ende): "
                     << "start < ende! Exiting!";
        throw runtime_error( errormessage.str() );
    }
    if (start==ende)
    {
        resize(1,1,start);
        return;
    }
    
    // es geht abwaerts
    bool up = true;
    if( schrittweite < 0 )
    {
        double tmp = start;
        start = ende;
        ende = tmp;
        schrittweite *= -1.0;
        up = false;
    }
    
    //cout <<  "swe = " << start << "," << schrittweite << ","<< ende <<endl;
    int laenge= int(ceil(fabs((ende - start) /
                              schrittweite) )); // immer aufrunden
    //cout <<" laenge "<<laenge<< endl;
    _rows=laenge;
    _cols=1;
    _data.resize(laenge,1);

    //out <<"[" << _rows <<","<<_cols<< "] " << endl;
    for( int i = 0; i < laenge ; ++i )
    {
        //cout <<"[" << i <<","<<start + i * schrittweite<< "] " << endl;
        _set(i, start + i * schrittweite);
    }

    //cout <<"ende,schrittweite= "<< ende << "," << schrittweite << endl;
    //cout <<" fmod(ende,schrittweite) = "<< fmod(ende,schrittweite) << endl;
    if  ((_get(laenge-1,0) < ende) )
    {
        //cout <<" anhaengen "<<laenge << endl;
        _data.push_back( double(start + laenge * schrittweite) );
        _rows+=1;
        if  ((_get(laenge,0) > ende) )
        {
            //cout <<" wieder entfernen "<<laenge << endl;
            _data.pop_back();
            _rows-=1;
        }
    }
    
    if( !up )        
        std::reverse( _data.begin(),_data.end()); 
    return;
}

// eine Spalte erhalten
Matrix Matrix::get_col(const int &col) const
{
    //vector<int> idx( 1, col );
    //return get_cols(idx);
    if (col > _cols -1 )
    {
        stringstream errormessage;
        errormessage << "Matrix Matrix::get_col(const int &col) const: "
                     "ERROR: invalid col index: " << col <<  " Exiting";
        throw logic_error( errormessage.str() );
    }
    Matrix vOut( _rows, 1, 0.0 );
    copy ( _data.begin()+col*_rows, _data.begin()+(col+1)*_rows,
           vOut._data.begin() );
    return vOut;
}

// eine Zeile erhalten
Matrix Matrix::get_row(const int &row) const
{
    return operator()(':', row );
}

// eine Spalte schreiben  this(:,i)= m2;
void Matrix::set_col(int i,const Matrix &m2)
{
    if( _rows != m2._rows )
    {
        stringstream errormessage;
        errormessage << "Matrix::set_col(int i, Matrix &m2): "
                     << "Matrix dimensions do not agree. Exiting!";
        throw runtime_error( errormessage.str() );
    }
    copy (m2._data.begin(), m2._data.end(), _data.begin()+i*_rows );
    return;
}

// Indizes schreiben  this(idx)= m2;
void Matrix::setIdx(const vector<int> &idx,const Matrix &m2)
{
    if(  idx.size() != m2.numel())
    {
        stringstream errormessage;
        errormessage << "void Matrix::setIdx(vector<int> idx,const Matrix &m2): "
                     << "Index and Matrix must be same size! Exiting!"
                     << "idx.size()=" << idx.size() << ", m2.numel() =" <<m2.numel();
        throw runtime_error( errormessage.str() );
    }
    for (int i = 0; i < idx.size(); i++)
    {
        if(  idx.at(i) > numel()-1)
        {
            stringstream errormessage;
            errormessage << "void Matrix::setIdx(vector<int> idx,const Matrix &m2): "
                         << "Index exeeds dimensions! Exiting!";
            throw runtime_error( errormessage.str() );
        }
        operator()(idx.at(i))= m2(i);
    }
    return;
}

// Indizes schreiben  this(idx)= m2;
void Matrix::setIdx(const vector<int> &idx,const double &v)
{
    for (int i = 0; i < idx.size(); i++)
    {
        if(  idx.at(i) > numel()-1)
        {
            stringstream errormessage;
            errormessage << "setIdx(const vector<int> &idx,const double &v): "
                         << "Index exeeds dimensions! Exiting!";
            throw runtime_error( errormessage.str() );
        }
        operator()(idx.at(i))= v;
    }
    return;
}

// Indizes schreiben  this(idx)= m2;
void Matrix::setIdx(const Matrix &idx,const double &v)
{
    for (int i = 0; i < idx.numel(); i++)
    {
        if(  idx.numel() > numel()-1)
        {
            stringstream errormessage;
            errormessage << "setIdx(const vector<int> &idx,const double &v): "
                         << "Index exeeds dimensions! Exiting!";
            throw runtime_error( errormessage.str() );
        }
        operator()(int(idx(i)))= v;
    }
    return;
}

// Indizes schreiben  this(idx)= m2;
void Matrix::setIdx(const Matrix &idx,const Matrix &m2)
{
    if(  idx.rows() != m2.numel() &&  idx.cols() != 1)
    {
        stringstream errormessage;
        errormessage << "void Matrix::setIdx(const Matrix &idx,const Matrix &m2): "
                     << "Index and Matrix must be same size! Exiting!";
        throw runtime_error( errormessage.str() );
    }
    for (int i = 0; i < idx.rows(); i++)
    {
        if(  idx(i) > numel()-1)
        {
            stringstream errormessage;
            errormessage << "vvoid Matrix::setIdx(const Matrix &idx,const Matrix &m2): "
                         << "Index exeeds dimensions! Exiting!";
            throw runtime_error( errormessage.str() );
        }
        operator()( int(idx(i)) )= m2(i);
    }
    return;
}

//  &m2=this(idx);
void Matrix::get_vec(const vector<int> &idx,Matrix &m2) const
{
    // falls idx leer
    if (idx.size() ==0)
    {
        m2.resize(0,0);
        return;
    }

    // dimension check
    int max_id=*max_element(idx.begin(), idx.end());
    if(  max_id > numel()-1)
    {
        stringstream errormessage;
        errormessage <<
                     "void Matrix::get_vec(const vector<int> idx,Matrix &m2) const: "
                     << "Index " << max_id << " exeeds dimensions! Exiting!";
        throw runtime_error( errormessage.str() );
    }
    // ggf. groesse aendern
    if(  idx.size() != m2.numel())
    {
        m2.resize(idx.size(),1);
    }
    for (int i = 0; i < idx.size(); i++)
    {
        m2(i)=operator()(idx.at(i));
    }
    return;
}

//  m2=this(idx);
Matrix Matrix::get_vec(const vector<int> &idx) const
{
    Matrix mOut;
    get_vec(idx,mOut);
    return mOut;
}

//  &m2=this(idx);
void Matrix::get_vec(const Matrix &idx,Matrix &m2) const
{
    // falls idx leer
    if (idx.numel() ==0)
    {
        m2.resize(0,0);
        return;
    }

    int max_id= int (idx.max());
    if( max_id > numel()-1)
    {
        stringstream errormessage;
        errormessage << "void Matrix::get_vec(const Matrix &idx,Matrix &m2) const: "
                     << "Index " << max_id << " exeeds dimensions! Exiting!";
        throw runtime_error( errormessage.str() );
    }
    // ggf. groesse aendern
    if(  idx.numel() != m2.numel())
    {
        m2.resize(idx.numel(),1);
    }
    for (int i = 0; i < idx.numel(); i++)
    {
        m2(i)=operator()(int(idx(i)));
    }
    return;
}

//  m2=this(idx);
Matrix Matrix::get_vec(const Matrix &idx) const
{
    Matrix mOut;
    get_vec(idx,mOut);
    return mOut;
}

// find: Idex finden mit (>,<,>=,<=,(==))
// Idex-vektor besitzt immer die dim: [found,1]
vector<int> Matrix::find_idx(bool (*pfunc)(double i, double j),
                             const double val ) const
{
    vector<int> idx;
    vector<double>::const_iterator it = _data.begin();

    it = find_if ( it, _data.end(), boost::bind(pfunc, _1, val) );

    while( it < _data.end() )
    {
        idx.push_back( it - _data.begin() );
        it = find_if ( it+1, _data.end(), boost::bind(pfunc, _1, val) );
    }
    return idx;
}

// find: Idex finden mit (>,<,>=,<=,(==)) mit zwei Bedingungen &&
// Idex-vektor besitzt immer die dim: [found,1]
vector<int> Matrix::find_idx(bool (*pfunc1)(double i, double j),
                             const double val1,
                             bool (*pfunc2)(double i, double j),const double val2 ) const
{
    vector<int> idx	;
    Matrix m2=find_elem(pfunc1, val1, idx);
    vector<int> idxOut=m2.find_idx(pfunc2, val2);
    Matrix idxM(idx);
    idxM=idxM(idxOut);
    idxOut=idxM.get_vec("int");
    return idxOut;
}

// find(this_nXm >,<,>=,<=,(==) vals_nXm ): Idex finden mit n x m Werten, elemw. Vergleich von Vektoren
// Idex-vektor besitzt immer die dim: [found,1]
vector<int> Matrix::find_idx(bool (*pfunc)(double i, double j),
                             const Matrix &vals)  const
{
    //cout << "vals:" <<endl;
    //vals.show();

    if(  vals.numel() != numel())
    {
        stringstream errormessage;
        errormessage <<
                     "vector<int> Matrix::find_idx(bool (*pfunc)(double i, double j), const Matrix vals): "
                     << "dimension error, argument vals must be same size as *this! Exiting!";
        throw runtime_error( errormessage.str() );
    }
    vector<int> idx;
    vector<double>::const_iterator it;

    for(int i=0; i< numel(); i++)
    {
        it = find_if ( _data.begin()+i, _data.begin()+i+1, boost::bind(pfunc, _1,
                       vals(i) ) );
        //cout << "_data.at(i)= " << _data.at(i)<<  "  and  vals(i)= " << vals(i);
        //cout << "  it - _data.begin() - i= " << it - _data.begin() << endl;
        if( it ==_data.begin()+i )
        {
            idx.push_back( it - _data.begin() );
        }
    }
    return idx;
}

// find(this_nXm >,<,>=,<=,(==) vals_nXm ): Idex finden mit n x m Werten, elemw. Vergleich von Vektoren
// Idex-vektor besitzt immer die dim: [found,1]
vector<int> Matrix::find_idx(bool (*pfunc1)(double i, double j),
                             const Matrix &vals1,
                             bool (*pfunc2)(double i, double j),const Matrix &vals2 ) const
{
    vector<int> idx=find_idx(pfunc1, vals1);
    Matrix m2= operator ()(idx);
    Matrix vals_neu= vals2(idx);
    vector<int> idxOut=m2.find_idx(pfunc2, vals_neu);
    Matrix idxM(idx);
    idxM=idxM(idxOut);
    idxOut=idxM.get_vec("int");
    return idxOut;
}

//  this(idx)=[];
void Matrix::elim_idx(vector<int> idx)
{
    // ggf. groesse aendern
    if(  _cols== 1 || _rows == 1)
    {
        sort(idx.begin(),
             idx.end()); // Reihenfolge des Outputs wird dadurch nicht veraendert
        vector<int>::iterator it;
        // auf Eindeutigkeit pruefen:
        it = unique (idx.begin(), idx.end());
        if (idx.end()-it != 0)
        {
            stringstream errormessage;
            errormessage << "void Matrix::elim_idx(vector<int> idx): "
                         << "duplicate indices! Exiting!";
            throw runtime_error( errormessage.str() );
        }
        for (int i = 0; i < idx.size(); i++)
        {
            if(  idx.at(i) > numel()-1)
            {
                stringstream errormessage;
                errormessage << "void Matrix::elim_idx(vector<int> idx): "
                             << "Index exeeds dimensions! Exiting!";
                throw runtime_error( errormessage.str() );
            }
            _data.erase (_data.begin()+idx.at(i)-i);
            //   cout << "del idx:" << idx.at(i)-i << endl;
        }
        if (_rows==1)
            _cols= _data.size();
        if (_cols==1)
            _rows= _data.size();
    }
    else
    {
        stringstream errormessage;
        errormessage << "void Matrix::elim_idx(vector<int> idx): "
                     << "data must be a vector (n,1)! Exiting!";
        throw runtime_error( errormessage.str() );
    }
    return;
}

//mehrere matrizen als vector<Matrix> zusammenfassen
void Matrix::set(const vector<Matrix> &vMat)
{
    int num = vMat.size();

    *this= vMat.at( 0 );

    for ( int x = 1; x < num; x++ )        // Loop until all matrices are added
        append_cols( vMat.at( x ) );
    return;
}

// modulo nach MATLAB =! fmod! (fmod erhaelt das Vorzeichen!), referenz
void Matrix::mod(Matrix &m2,const double v) const
{
    m2.resize(rows(),cols());
    // sinus elementeweise, laeuft
    for( int i = 0; i < _data.size() ; ++i )
    {
        // aus funother
        m2._set(i, modulo(_data.at(i),v));
    }
    return;
}

Matrix Matrix::mod(const double &v) const
{
    Matrix out;
    mod(out,v);
    return out;
}

void Matrix::fmod(Matrix &m2,const double v) const
{
    m2.resize(rows(),cols());
    // sinus elementeweise, laeuft
    for( int i = 0; i < _data.size() ; ++i )
    {
        // aus funother
        m2._set(i, std::fmod(_data.at(i),v));
    }
    return;
}


Matrix Matrix::fmod(const double &v) const
{
    Matrix out;
    fmod(out,v);
    return out;
}


void Matrix::modf(Matrix &fractpart, Matrix &intpart) const{
    intpart.resize(rows(), cols());
    fractpart.resize(rows(), cols());
    for( int i = 0; i < _data.size() ; ++i )
    {
            // aus funother
            fractpart._set(i, std::modf(_data.at(i), &intpart(i)) );
        }
    }

Matrix Matrix::modf( Matrix & intpart) const{
    Matrix fractpart;
    modf(fractpart,intpart);
    return fractpart;
}



//abs als referenz
void Matrix::abs(Matrix &m2) const
{
    m2.resize(rows(),cols());
    // sinus elementeweise, laeuft
    for( int i = 0; i < _data.size() ; ++i )
    {
        m2._set(i, std::fabs(_data.at(i)));
    }
    return;
}

Matrix Matrix::abs() const
{
    Matrix out;
    abs(out);
    return out;
}

Matrix Matrix::absD( ) const
{
    Matrix out( _rows, _cols );

    for( int i=0; i<_cols*_rows; i++ )
    {
        out( i ) = fabs( _data.at( i ) );
    }

    return( out );
}


// ...........................................................................
Matrix Matrix::diff( ) const
{
    Matrix out( _rows-1, _cols );

    for( int i=0; i<_cols; i++ )
    {
        for( int j=1; j<_rows; j++ )
            out( j-1,i ) = _get( j,i ) - _get( j-1,i );
    }

    return( out );
}

// ...........................................................................
Matrix Matrix::sign( ) const
{
    Matrix a;
    abs(a);

    Matrix out = *this;
    out = out.div_elem( a );

    return( out );
}

/* klappt nicht sagt Andi
void Matrix::set(const string &strIn)
{
    resize(0,0,0);
    string::size_type pos = 0;
    string::size_type pos2= 0;
    string delimiter = ";";
    string delimiter2 = ",";
    string line ;
    string val;
    double wert;
    int i=0,j=0;
    do
    {
        // Zeilen durchgehen
        line = strIn.substr(pos, strIn.find(delimiter, pos) - pos);
        //cout << line << endl;
        i++;
        pos2 = 0;
        // Zeile auftrennen
        do
        {
            std::istringstream stm;
            stm.str(line.substr(pos2, line.find(delimiter2, pos2) - pos2));
            stm >> wert;
            _data.push_back(wert);
            //cout << wert << endl;
            j++;
        }
        while((pos2 = line.find(delimiter2, pos2))++ != string::npos);
    }
    while((pos = strIn.find(delimiter, pos))++ != string::npos);

    _rows=j/i;
    _cols=i;
    if ( modulo(double(j),double(i))!= 0)
    {
        stringstream errormessage;
        errormessage << "Matrix::Matrix(const string &strIn): "
                     "ERROR: invalid arrangement von elements: Exiting";
        throw logic_error( errormessage.str() );
    }
    *this=transpose();
}
 * */

// abrunden nach MATLAB , referenz
void Matrix::floor(Matrix &m2) const
{
    m2.resize(rows(),cols());
    // sinus elementeweise, laeuft
    for( int i = 0; i < _data.size() ; ++i )
    {
        m2._set(i,std::floor(_data.at(i)));
    }
    return;
}

Matrix Matrix::floor() const
{
    Matrix out;
    floor(out);
    return out;
}

// abrunden nach MATLAB , referenz
void Matrix::round(Matrix &m2) const
{
    m2.resize(rows(),cols());
    // sinus elementeweise, laeuft
    for( int i = 0; i < _data.size() ; ++i )
    {
        m2._set(i,roundD(_data.at(i)));
    }
    return;
}

Matrix Matrix::round() const
{
    Matrix out;
    round(out);
    return out;
}

// abrunden nach MATLAB , referenz
void Matrix::rem(Matrix &m2,const double v) const
{
    m2.resize(rows(),cols());
    // sinus elementeweise, laeuft
    for( int i = 0; i < _data.size() ; ++i )
    {
        m2._set(i,remainder(_data.at(i),v));
    }
    return;
}

Matrix Matrix::rem(const double v) const
{
    Matrix out;
    rem(out,v);
    return out;
}

// .....................................................................
Matrix Matrix::norm( )
// .....................................................................
{
    // get transpose of this-matrix
    Matrix m1 = *this ;
    Matrix trans = m1.transpose();

    // calculate norm as the square root of the transpose of this-matrix times this-matrix
    Matrix m2 = trans * m1 ;
    m2 = m2.sqrt( ) ;

    return m2 ;
}


// .....................................................................
void Matrix::rand_norm( double mean, double std )
// .....................................................................
{
    _rng_calls++;   
    // create vector of zero mean Gaussian random numbers with variance one
    typedef boost::mt19937 RNGType; // mersenne twister generator

    RNGType rng( time(NULL)+_rng_calls ) ;
    
    // < normal distribution with mean of 0.0 and standard deviation of 1.0
    boost::normal_distribution<> rdist( mean,
                                        std ); 

    boost::variate_generator< RNGType, boost::normal_distribution<> > get_rand(
        rng, rdist );

    generate( _data.begin(), _data.end(), get_rand );
    

}



// .....................................................................
void Matrix::zero( )
{
    // .....................................................................

    // set all elements to zero
    Matrix m2( _rows, _cols , 0.0 );

    *this = m2 ;

    return ;
}

void Matrix::blockToepliz( const std::vector<ivg::Matrix>& mats, const char type){
   
    unsigned n = mats[0].rows();
    unsigned nb = mats.size();
    resize(nb*n, nb*n);
    
    switch (type) {
        case 'N': {
            std::vector<unsigned> toe_idx (nb);
            for(int i = 0; i < toe_idx.size(); ++i){
                toe_idx[i] = i;
            }

            unsigned c = 0;
            for(unsigned i = 0; i < nb; ++i){
                for(unsigned j = 0; j < nb; ++j){
                    unsigned idx = toe_idx[c%nb];
                    set_sub( i*n, j*n, mats[idx] );
                    ++c;
                }
                c--;
            }
            break;
        }
        case 'S' : {
            
            for(unsigned i = 0; i < nb; ++i){
                unsigned c = 0;
                for(unsigned j = 0; j < nb; ++j){
                    if(j>=i){
                        if(i==j){
                            set_sub( i*n, j*n, mats[c] );
                        } else {
                            set_sub( i*n, j*n, mats[c] );
                            set_sub( j*n, i*n, mats[c].transpose() );
                        }
                        ++c;
                    }
                }
            }
            
            break;
        }
        default:{
            throw logic_error( "void Matrix::blockToepliz. invalid type " + type);
        }
    }
   
    
}


// .....................................................................
Matrix Matrix::estimate_cpwlf( Matrix x, Matrix x0, double w )
// .....................................................................
{
    Matrix y = *this ;

    int no = x0.size(1);	// no of obs

    // seperate data into segments
    std::vector< std::vector<int> > seg ;
    std::vector<int> found ;

    for( int i = 0; i < no-1; i++ )
    {
        found = x.find_idx( gt, x0(i), le, x0(i + 1) ) ;
        seg.push_back( found ) ;

    }

    // normal equation system
    Matrix N( no, no, 0.0 );
    Matrix n( no, 1, 0.0 );

    // constraints (for stations where no met. data is measured)
    // the corresponding weights are in Matrix W
    Matrix B( no, no, 0.0 );

    Matrix W( no, 1, w ) ;
    W = W.diag( );

    for( int i = 0; i < no; i++ )
    {
        if ( i != 0 )		 // do not for the first time
        {
            if( seg[ i-1 ].begin() != seg[ i-1 ].end() )
            {
                N( i, i-1 ) -= ( (  ( x(seg[i-1]) - x0(i-1) ).mult_elem( x(seg[i-1])-x0(
                                        i) )  ).sum_col() / ( std::pow ( (x0(i) - x0(i-1)), 2.0 ) ) )( 0 );
                N( i, i )   += ( (  ( x(seg[i-1]) - x0(i-1))^2.0  ).sum_col() / ( std::pow( (
                                     x0(i) - x0(i-1)), 2.0 ) ) )( 0 ) ;

                n( i )   	+= ( (  ( x(seg[i-1]).mult_elem( y(seg[i-1]) ) ).sum_col() - ( ( y(
                                        seg[i-1]) ).sum_col() ) * x0(i-1) ) / (x0(i) - x0(i-1)) )( 0 );
            }

        }
        if ( i != no-1 ) 	// do not for the last time
        {
            if( seg[ i ].begin() != seg[ i ].end() )
            {

                N( i, i )   += ( (  ( x(seg[i]) - x0(i+1))^2.0 ).sum_col() / ( std::pow( (x0(
                                     i+1) - x0(i)), 2.0 ) ) )( 0 ) ;
                N( i, i+1 ) -= ( (  ( x(seg[i]) - x0(i) ).mult_elem( x(seg[i])-x0(
                                        i+1) )  ).sum_col() / ( std::pow( (x0(i+1) - x0(i)), 2.0 ) ) )( 0 );

                n( i )      += ( (  ( x(seg[i]).mult_elem( y(seg[i]) ) ).sum_col() *
                                    (-1.0) + ( ( y(seg[i]) ).sum_col() ) * x0(i+1) ) / (x0(i+1) - x0(i)) )( 0 );
            }

            B( i, i ) 	= -1.0 ;
            B( i, i+1 )	=  1.0 ;
        }
    }

    Matrix C = B.transpose() * W * B ;

    // addition of the normal equation system
    // (the addition to the right hand side (n) would be zero)
    Matrix p = ( N + C ).chol_inv() * n;

    return p ;
}




// +++ edit SH: 2013-11-20
// ......................................
Matrix Matrix::gamma_fct( )
// ......................................
{
    Matrix gam = ( *this );
    Matrix M_gamma( gam.size(1), gam.size(2), 0.0 );
    Matrix tmp ;

    if ( _cols == 1 || _rows == 1 )
    {
        if ( _cols == 1 )
        {
            tmp = gam ;
            tmp.transpose();
            gam = tmp ;
        }

        for( int i = 0; i <= gam.size(2)-1; i++ )
        {
            M_gamma(0,i) = tgamma( gam(0,i) ) ;
        }
    }
    else
    {
        stringstream errormessage;
        errormessage << "Matrix:: Matrix::gamma_fct( )"
                     << "ERROR: invalid number of elements " << " Exiting";
        throw logic_error( errormessage.str() );
    }


    return M_gamma ;


}
// --- edit SH: 2013-11-20


// .....................................................................
void Matrix::load_mat_bin( string varname, string matfile )
// .....................................................................
{
    mat_t    *matfp;
    matvar_t *matvar;


    matfp = Mat_Open( matfile.c_str(),MAT_ACC_RDONLY);
    if ( NULL == matfp )
    {
        stringstream errormessage;
        errormessage << "void Matrix::load_mat_bin( string file ) : "
                     << "Failed to open MAT file: " << matfile << ". Exiting";
        throw runtime_error( errormessage.str() );
    }

    matvar = Mat_VarReadNext(matfp);
    while ( NULL != matvar )
    {
        size_t* dim;
        dim = matvar->dims;
        //cout <<matvar->name << " - " << varname << " - " << dim[0] << ", " << dim[1] << endl;

        if( matvar->name == varname )
        {
            resize( (int)dim[0],(int)dim[1],0.0 );

            switch( matvar->class_type )
            {
            case MAT_C_DOUBLE:
            {
                int    start[2]= {0,0};
                int    stride[2]= {1,1};
                int    edge[2];
                edge[0]=dim[0];
                edge[1]=dim[1];

                Mat_VarReadData(matfp,matvar,data_ptr(),start,stride,edge);
                break;
            }
            case MAT_C_SPARSE:
            {
                sparse2array( matvar,data_ptr() );
                break;
            }
            }
        }
        Mat_VarFree(matvar);
        matvar = Mat_VarReadNext(matfp);
    }
    Mat_Close(matfp);
}

void Matrix::rot3D_x( double psi )
{
    resize( 3,3,0.0 );

    _data.at( 0 ) =  1.0;
    _data.at( 4 ) =  std::cos(psi);
    _data.at( 5 ) = -std::sin(psi);
    _data.at( 7 ) =  std::sin(psi);
    _data.at( 8 ) =  std::cos(psi);
}

void Matrix::rot3D_y( double psi )
{
    resize( 3,3,0.0 );

    _data.at( 0 ) =  std::cos(psi);
    _data.at( 2 ) =  std::sin(psi);
    _data.at( 4 ) =  1.0;
    _data.at( 6 ) = -std::sin(psi);
    _data.at( 8 ) =  std::cos(psi);
}

void Matrix::rot3D_z( double psi )
{
    resize( 3,3,0.0 );

    _data.at( 0 ) =  std::cos(psi);
    _data.at( 1 ) = -std::sin(psi);
    _data.at( 3 ) =  std::sin(psi);
    _data.at( 4 ) =  std::cos(psi);
    _data.at( 8 ) =  1.0;
}




// .....................................................................
Matrix Matrix::interpolate( Matrix data_epochs, double epoch,
                            std::string interpolation_type ) const
// .....................................................................
{

    double data_arr[ _data.size() ];
    std::copy( _data.begin(), _data.end(), data_arr );

    // declare accelerator and spline object
    gsl_interp_accel *accel_ptr;
    gsl_spline *spline_ptr;

    int npts = data_epochs.size(1);
    
    if(npts == 1)
        throw runtime_error( "Matrix Matrix::interpolate( Matrix data_epochs, double epoch, std::string interpolation_type ): No interpolation possible. Insufficient number of datapoints.\n" );

    // allocate accelerator and spline object
    accel_ptr = gsl_interp_accel_alloc();
    

    if( interpolation_type.compare( "linear" ) == 0 )
    {
        spline_ptr = gsl_spline_alloc( gsl_interp_linear, npts );
    }
    else if( interpolation_type.compare( "polynomial" ) == 0 )
    {
        spline_ptr = gsl_spline_alloc( gsl_interp_polynomial, npts );
    }
    else if( interpolation_type.compare( "cspline" ) == 0 )
    {
        spline_ptr = gsl_spline_alloc(  gsl_interp_cspline, npts );
    }
    else if( interpolation_type.compare( "neville" ) == 0)
    {
        ivg::Matrix result_row(1,cols(),0.0);
        for(int c=0; c< cols(); c++)
        {
            ivg::Matrix Q(npts,npts,0.0);
            
            if(npts != this->get_col(c).numel())
                throw runtime_error("Matrix Matrix::interpolate(... neville ...): Length of data and data_epochs does not correspond.");
            
            Q.set_col(0,this->get_col(c));
            
            for(int i=2; i<=npts; i++)
            {
                for(int j=2; j<=i; j++)
                {
                    double n1 = (epoch - data_epochs(i -(j-1)-1)) * Q(i-1,j-1-1);
                    double n2 = (epoch - data_epochs(i-1)) * Q(i-1-1,j-1-1);
                    double d = data_epochs(i-1) - data_epochs(i-(j-1)-1);
                    Q(i-1,j-1) = (n1-n2)/d;
                }
            }
            result_row(c) = Q(npts-1,npts-1);
        }    
        
        return result_row;
    }
    else if( interpolation_type.compare( "nearest_neighbor" ) == 0 )
    {
        int min_ix;
        (data_epochs - epoch).abs().min(min_ix);
        vector<int> ix(1);
        ix.at(0) = min_ix;
        return this->operator()(ix,":");
    }
    else
    {
        stringstream errormessage;
        errormessage <<
                     "Matrix Matrix::interpolate( Matrix data_epochs, double epoch, "
                     << "std::string interpolation_type ): "
                     << "ERROR: invalid interpoltion type " << interpolation_type << ". Exiting";
        throw logic_error( errormessage.str() );
    }


    // initialize spline
    ivg::Matrix out( 1,_cols,0.0 );

    int status;
    for( int i = 0; i < _cols; ++i )
    {
        status = gsl_spline_init( spline_ptr, data_epochs.data_ptr(), data_ptr()+i*_rows, npts );
        
        // evaluate spline
        out(i) = gsl_spline_eval( spline_ptr, epoch, accel_ptr );
    }

    // free accelerator and spline object
    gsl_spline_free( spline_ptr );
    gsl_interp_accel_free( accel_ptr );

    return out;

}

} // # namespace ivg


