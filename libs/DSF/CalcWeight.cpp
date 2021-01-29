#include "math.h"
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <vector>
#include <iostream>

//using namespace std;

using std::cout;
using std::endl;
using std::vector;

#define matrix(i,j) matrix[i+j*m]

void PlotMatrix (double * matrix, int n, int m) {
    for (int j=0; j<n; j++) {
        for(int i=0; i<m; i++) {
            std::cout << matrix(i,j) << ' ';
        }
        std::cout << '\n';
    }
    return;
}

double CalcNodeWeight(double * matrix,int n,int m, std::vector<int> S, int id)
{  
    
    double w=0.0;
    if (S.size()==1){
        w=1.0;
    }else{        
        //get the position of the i element in the S
        int pos=-1;
        for(int i=0;i<S.size();i++){
            if(S[i]==id){
                pos=i;
                break;
            }
        }

        if (pos==-1) {
            //mexErrMsgTxt("ERROR: Element not found");  
        }else{
            S.erase(S.begin() + pos);
        }
        
        int h=S[0];
        //std::cout << " h=" << h << endl;
        int j=0;
        
        for(int c=0;c<S.size();c++){   
            j=S[c];
            //std::cout << "id = " << id << " h=" << h << " j=" << j << endl;
            //std::cout << " w=" << w << endl;
            //std::cout << "El i,j" << matrix(id,j) << "    El h,j " << matrix(h,j) << endl;
            w=w+(matrix(id,j)-matrix(h,j))*CalcNodeWeight(matrix,n,m,S,j);
            
        }
        /*
        std::cout << "S contains:";
        for (unsigned i=0; i<S.size(); ++i){
            std::cout << ' ' << S[i];
        }
        std::cout << '\n';*/
     }
     return w;   
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    /* Variables */
    const mwSize *dim_array, *dim_G;
    //mwSize* dim_array;
    //mwSize* dim_G;
    int r, c, m, n;
    double *w, *A, *i, *G; 
            
    /* get a pointer to the input matrix */
    A        = mxGetPr (prhs[0]); 
	G        = mxGetPr (prhs[1]); 
	i        = mxGetPr (prhs[2]); 
	
    /* get the number of rows, columns and features */
    dim_array = mxGetDimensions(prhs[0]);
    n         = dim_array[0];
    m         = dim_array[1];
	
	dim_G     = mxGetDimensions(prhs[1]);
    r         = dim_G[0];
    c         = dim_G[1];
        
    /* create the output matrix */
    plhs[0] = mxCreateDoubleScalar(mxREAL);
    /* get a pointer to the real data in the output matrix */
    w = mxGetPr(plhs[0]);
	
    /* create a vector for S */
    std::vector<int> S;    
    for(int j=0;j<c;j++){
        S.push_back(G[j]-1);
    }
    
    int ids=(*i);
    ids=ids-1;
    //std::cout << (*i) << " ->" << ids << '\n';
    /*
    std::cout << "S contains " << S.size() << " elements";
    for (unsigned jj=0; jj<S.size(); ++jj){
        std::cout << ' ' << S[jj];
    }
    std::cout << " ->" << (*i) << '\n';
    */
    
    *w=CalcNodeWeight(A,n,m,S,ids);
    
    //std::cout << (*w) << '\n';
    /* Free Memory */
    //mxFree (P);
    //mxFree (Q);
}

/* copy a vector to a matlab array */
mxArray * getMexArray(const std::vector<double>& v){
    mxArray * mx = mxCreateDoubleMatrix(1,v.size(), mxREAL);
    std::copy(v.begin(), v.end(), mxGetPr(mx));
    return mx;
}

/* const static double eps = 1e-4;

/* equivalent to matlab function sum(sum(a)) 
double SumMatrix (double *matrix, int n, int m, int k) {
    double sum = 0;
    for (int j=0; j<n; j++) {
        for(int i=0; i<m; i++) {
            if (mxIsNaN(matrix(i,j,k))) {
                mexErrMsgTxt("ERROR: Input must not be NaN.");  
            }
            sum += matrix(i,j,k);
        }
    }
    return sum;
}

/* equivalent to matlab function sum(sum(a.*b)) 
double ProdMatrix (double *matrix, int n, int m, int k1, int k2) {
    double prod = 0;
    for (int j=0; j<n; j++) {
        for(int i=0; i<m; i++) {
            prod += matrix(i,j,k1)*matrix(i,j,k2);
        }
    }
    return prod;
}
*/