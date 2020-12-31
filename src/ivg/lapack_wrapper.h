extern"C" {
    void dgtsv_(const long *Np, const long *NRHSp, double *DL,
                double *D, double *DU, double *B, const long *LDBp,
                long *INFOp);
   // void dgesv_(const long *Np, const long *NRHSp, double *A, const long *LDAp,
   //             long *IPIVp, double *B, const long *LDBp, long *INFOp );
    void dgesvd_(const char *JOBUp, const char *JOBVTp, const long *Mp,
                 const long *Np, double *A, const long *LDAp, double *S,
                 double *U, const long *LDUp, double *VT, const long *LDVTp,
                 double *WORK, const long *LWORKp, const long *INFOp );
    void dgeqp3_( const long *Mp, const long *Np, double *A, const long *LDAp, const int *JPVTp, 
                  double *tau, double *WORK, const long *LWORKp, const long *INFOp );
    void dgeqrf_( const long *Mp, const long *Np, double *A, const long *LDAp,
                  double *tau, double *WORK, const long *LWORKp,
                  const long *INFOp );
    void dorgqr_( const long *Mp,  const long *Np, const long *Kp, double *A,
                  const long *LDAp, double *tau, double *WORK,
                  const long *LWORKp, const long *INFOp );
    void dormqr_( const char*  SIDE, const char* TRANS, const long* M, const long* N, const long* K, 
                    double* A, const long* LDA, double* TAU, double* C, const long* LDC,
                    double* WORK, long* LWORK, long* INFO );
    void dtrtrs_( const char* UPLO, const char* TRANS, const char* DIAG, 
                const long* N, const long* NRHS, double* A, const long* LDA, 
                double* B, const long* LDB, long* INFO );
    //void dgetrf_( long *Mp, long *Np, double *A, long *LDAp, int *ipiv, long *INFO );
    void dgetrs_(char *TRANS, long *N, long *NRHS, double *A, 
                      long *LDA, int *IPIV, double *B, long *LDB, long *INFO );    
    void dgesvx_( char *fact, char *trans, int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv,
           char *equed, double *r, double *c, double *b, int *ldb, double *x, int *ldx, double *rcond, double *ferr, double *berr,
            double *work, int *iwork, int *info );
    
     void dpotrf_(const char& UL, const int& N, double* A, const int& lda, int& info);
     
     void dposv_( const char& UPLO, const int& N, const int& NRHS, double* A,
                          const int& LDA, double* B, const int& LDB, int& INFO );
     
     void dgels_(const char& TRANS, const int& M, const int& N, const int& NRHS, double* A,
                 const int& LDA, double* B, const int& LDB, double* WORK, const int& LWORK, int& INFO);
     
     void dgetrf_( const int& M, const int& N, double* A, const int& LDA, int* IPIV, int& INFO);
     
     void dgetri_( const int& N, double* A, const int& LDA, int* IPIV, double* WORK, const int& LWORK, int& INFO );
     
    void dgesv_( const int& N, const int& NRHS, double* A, const int& LDA, int* IPIV, double* B, const int& LDB, int& INFO   );
}
