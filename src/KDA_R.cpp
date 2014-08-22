#include <Rcpp.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <time.h>    
using namespace std;
using namespace Rcpp;

// Function for the quick sort routine
static  int intcompare2(const void *i, const void *j)
{
  int *u,*v;
  u = (int *) i;
  v = (int *) j;
  if (*u > *v)
    return (1);
  if (*u < *v)
    return (-1);
  return (0);
}


long double *K;                                   // K[A]=K(D,A)in Etienne's paper
long double J, SPP;
int *Abund;
int numspecies;



void calcLogKDA() {           //after Jerome Chave

	int i,n,im,s;
   
    //ofstream out(nomfo);    
    //out  <<"S\t J\t Theta\t Std_Theta\t I\t Std_I\t m\t Std_m\t loglike_min\t Theta_Ewens\t loglike_Ewens\t Theta2\t Std_Theta2\t I2\t Std_I2\t m2\t Std_m2\t loglike_min2\n";
    qsort(Abund,numspecies,sizeof(int),intcompare2);
    
    J=0;
    SPP = numspecies;
    for(s=0;s<SPP;s++)
        J += Abund[s];
     //cerr << "Number of individuals: "<< J << endl;
   

    int MaxA = 0;
    MaxA=Abund[(int)SPP-1];

    //cerr <<" "<<endl;
   // cerr << "Maximal abundance: " << MaxA << endl;

   
    // abundance distribution
    int *Phi = new int[MaxA+1];
    for(s=0;s<=MaxA;s++) Phi[s]=0;
    for(s=0;s<SPP;s++) Phi[Abund[s]]++;

    // Number of distinct abundances
    int NDA=0;
    for(s=0;s<=MaxA;s++) if(Phi[s] > 0) {NDA++;}
    
    
    //cerr << "Start computing Stirling numbers ...\n";
    // FIRST STAGE: compute the Stirling numbers
    // I use the relationship S(n,m)=S(n-1,m-1)-(n-1)S(n-1,m)
    // where n >= m >= 1
    // The relation of interest is sum_m S(n,m)/S(n,1)*S(m,1) * x^m
    // defined to be equal to sum_m T(n,m) * x^m
    // The recurrence relation on T(n,m) is 
    // T(n,m)= T(n-1,m) + T(n-1,m-1)*(m-1)/(n-1)
   
    int *f = new int[NDA];
    int *g = new int[NDA];
    i=0;
    for(s=0;s<NDA;s++) {f[s]=0;g[s]=0;}
    for(n=0;n<=MaxA;n++) if(Phi[n] > 0) {             
        f[i] = Phi[n];                                  
        g[i] = n;                                        
        i++;
        }
    long double **T= new long double*[NDA];          // T(n,m) just for the n which are useful
    T[0] = new long double[g[0]+1];
    T[0][0]=0;T[0][1]=1;
    if (g[0]!=1){
        long double *lS2 = new long double[g[0]+1]; 
        lS2[0]=0;lS2[1]=1;
        for (n=2;n<=g[0];n++) {
            long double *lS1 = new long double[n+1];                
            for(im=0;im<=n-1;im++) {
                lS1[im] = lS2[im];
            }
            lS1[n]=0;
            for(im=2;im<=n;im++) {
                lS2[im] = lS1[im]+lS1[im-1]*(im-1)/(n-1); 
            }
            delete[] lS1;            
        }
        for(im=2;im<=g[0];im++) {
            T[0][im]=lS2[im];
        }
        delete[] lS2;
    }
    for (int in=1;in<i;in++){
        T[in]= new long double[g[in]+1];
        T[in][0]=0;T[in][1]=1;
        long double *lS2 = new long double[g[in]+1];         
        for(im=0;im<=g[in-1];im++) {
                lS2[im] = T[in-1][im];
            }
        for (n=g[in-1]+1;n<=g[in];n++) {
            long double *lS1 = new long double[n+1];                
            for(im=0;im<=n-1;im++) {
                lS1[im] = lS2[im];
            }
            lS1[n]=0;
            for(im=2;im<=n;im++) {
                lS2[im] = lS1[im]+lS1[im-1]*(im-1)/(n-1); 
            }
            delete[] lS1;            
        }
        for(im=2;im<=g[in];im++) {
            T[in][im]=lS2[im];
        }
        delete[] lS2;
    }
    // After this stage we have stored in T[i][m] T(g[i],m)
    // with T(n,m) = S(n,m)*S(m,1)/S(n,1) for i>0

    //cerr << "Start computing ln(K(D,A)) ...\n";
    // SECOND STAGE: compute the K(D,A)
    // I follow Etienne's route. Compute the product of polynomials 
    // of length J
    int j,nn,mm;
    K = new long double[int(J)+1];
    long double *poly2 = new long double[int(J)+1];
    for(i=0;i<=J;i++){
        K[i] = poly2[i] = 0.0;
    }
    K[0]=1;
    int degree = 0;
    int spe=0;
    for(i=0;i<NDA;i++) // loop over number of distinct abundances
        for(j=0;j<f[i];j++){ // loop over abundances per class             
            for(nn=0;nn<=degree;nn++)
                for(mm=1;mm<=g[i];mm++){
                    if (K[nn]>0){                   
                       poly2[nn+mm] += T[i][mm]*K[nn];
                    }
                    
                }              
            degree += g[i];
            for(nn=0;nn<=degree;nn++){            
                K[nn] = (poly2[nn]/powl(10,(4500.0/SPP)));
                poly2[nn] = 0.0;
            }
            spe++;
        }
    
    for(i=int(SPP);i<=J;i++){
        K[i] = logl(K[i]);                                    // now K[A]=ln(K(D,A)/10^4500) in Etienne's paper
    }
    for(i=0;i<NDA;i++) delete[] T[i];
	delete[] T;
    delete[] poly2;
    delete[] f;
    delete[] g;

// search of "infinite" values in K[A]
    int borneinf=int(SPP-1);
    int bornesup=int(J+1);
    long double maxlog=11333.2;
    int infinity=0;
    for(i=int(SPP);i<=J;i++){
        if ((K[i]>maxlog)||(K[i]<-maxlog)) {
            infinity=1;
            break;
        }
        borneinf++;
    }  //after that, borneinf=indice next to infinity but before
    for(int i=0;i<=J-SPP;i++){
        if ((K[(int)J-i]>maxlog)||(K[(int)J-i]<-maxlog)) {
            infinity=1;
            break;
        }
        bornesup--;
    }    //after that, bornesup=indice next to infinity but after
  if (infinity==1){
  //cerr << "WARNING : the sample is too large to compute an exact likelihood, the program is thus doing approximations. The following results are to be taken with caution"<<endl;
  //cerr << "Value of A above which K(D,A) is computed approximately ="<<borneinf<<endl;
  //cerr << "Value of A below which K(D,A) is computed approximately ="<<bornesup<<endl;

  //fitting of the infinite values of K[A] by a polynom of degree 3
    //computing of the derivatives at the critic points
    long double Kprimeinf = K[borneinf]-K[borneinf-1];
    long double Kprimesup = K[bornesup+1]-K[bornesup];
    // definition of the parameters of the fitted polynom aX^3+bX^2+cX+d
    long double a,b,c,d;
    //inversion of the linear system of equations (with the Gauss method)
    long double borneinf2=(long double)borneinf*(long double)borneinf;
    long double borneinf3=(long double)borneinf2*(long double)borneinf;
    long double bornesup2=(long double)bornesup*(long double)bornesup;
    long double bornesup3=(long double)bornesup2*(long double)bornesup;
    d=(Kprimesup-3*bornesup2*K[borneinf]/borneinf3+(2*bornesup/(long double)borneinf-3*bornesup2/borneinf2)*(Kprimeinf-3*K[borneinf]/(long double)borneinf)-((1+3*bornesup2/borneinf2-4*bornesup/(long double)borneinf)/(bornesup-2*bornesup2/(long double)borneinf+bornesup3/borneinf2))*(K[bornesup]-bornesup3*K[borneinf]/borneinf3+(bornesup2/(long double)borneinf-bornesup3/borneinf2)*(Kprimeinf-3*K[borneinf]/(long double)borneinf)))/((6*bornesup2/borneinf3)-(6*bornesup/borneinf2)-((1+3*bornesup2/borneinf2-4*bornesup/(long double)borneinf)/(bornesup-2*bornesup2/(long double)borneinf+bornesup3/borneinf2))*(1-3*bornesup2/borneinf2+2*bornesup3/borneinf3));
    c=((K[bornesup]-bornesup3*K[borneinf]/borneinf3+(bornesup2/(long double)borneinf-bornesup3/borneinf2)*(Kprimeinf-3*K[borneinf]/(long double)borneinf))-d*(1-3*bornesup2/borneinf2+2*bornesup3/borneinf3))/(bornesup-2*bornesup2/(long double)borneinf+bornesup3/borneinf2);
    b=(Kprimeinf-3*K[borneinf]/(long double)borneinf+2*c+3*d/(long double)borneinf)/(0.0-(long double)borneinf);
    a=(K[borneinf]-b*borneinf2-c*(long double)borneinf-d)/borneinf3;
    
    //reconstruction of K[A] with the fitted polynom
    for (int i=borneinf+1;i<bornesup;i++) {
     K[i]=(a*i*i*i+b*i*i+c*i+d);
    }
   }
}




// [[Rcpp::export]]
NumericVector calcKDA(NumericVector A)
{
    //convert abundances from A to Species
	numspecies = A.size();
	Abund = new int[numspecies];                          

	int J = 0;
	for(int s = 0; s < numspecies; ++s) {
	   Abund[s] = A[s];
	   J += Abund[s];
	}
	//call calcLogKDA
	calcLogKDA();
	//return K
	   
	int sizeofK =  J + 1; //I hope!
	NumericVector out(sizeofK);
	for(int i = 0; i < sizeofK; ++i) {
	   out[i] = K[i] + 4500.0 * logl(10);
	}


    return out;	
}
