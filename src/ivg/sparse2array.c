#include <matio.h>

void sparse2array( matvar_t *matvar, double *mat )
{
   int i,j;
   mat_sparse_t *sparse;

   sparse = matvar->data;
   size_t* dim;
   dim = matvar->dims;

   double *data;

   data = sparse->data;
   for ( i = 0; i < sparse->njc-1; i++ ) 
   {
      for (j = sparse->jc[i]; j<sparse->jc[i+1] && j<sparse->ndata;j++ )
         mat[sparse->ir[j]+dim[0]*i] = data[j];
   }
};
