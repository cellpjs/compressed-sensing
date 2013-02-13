/************************************************************************/
/*  suprem algorithm for compressed sensing recovery                    */ 
/*  This is a C source file for creating a mex file                     */
/*  which can be called from MATLAB                                     */
/*  written by Jinsoo Park and Mehmet Akcakaya                          */             
/*  Harvard School of Engineering and Applied Sciences 2009             */
/************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mex.h"

#define vs 0.00000000000000000001

/* comparison function for qsort - descending order */
int compare_doubles (const void *a, const void *b){
   const double *da = (const double *) a;
   const double *db = (const double *) b;
 
   return (*da < *db) - (*da > *db);
}

void decode(double *xhat, double *r, double *E, int K, int M, int N, int dc, int dv, int L, double sd, int n_rounds, int n_iterations, int cstop){
	
	int round, i, j, k, kin, kout, jv, jc;
	int active, count, istop;
	int *inactive;
	double **cin, **cout, **vin, **vout;
	double sum, msum, mu, rdis, rdis_min, delta, threshold;
	double **message, **bel, *beta, *pool, *rhat, *tempo;
	int *bigo, *bigo1, *bigo2, *old;
						
	/* allocate memory */
    
    
    
    
	message=(double **)mxCalloc(K, sizeof(double *));
	for (i=0;i<K;i++){
	  message[i]=(double *)mxCalloc(2, sizeof(double));
	}
	bel=(double **)mxCalloc(M, sizeof(double *));
	for (i=0;i<M;i++){
	  bel[i]=(double *)mxCalloc(2, sizeof(double));
	}
	beta=(double *)mxCalloc(M, sizeof(double));	
	pool=(double *)mxCalloc(M, sizeof(double));
	bigo=(int *)mxCalloc(M, sizeof(int));
	bigo1=(int *)mxCalloc(M, sizeof(int));
	bigo2=(int *)mxCalloc(M, sizeof(int));
	tempo=(double *)mxCalloc(M, sizeof(double));
	old=(int *)mxCalloc(K, sizeof(int));
	rhat=(double *)mxCalloc(N, sizeof(double));
  	inactive=(int *)mxCalloc(dv, sizeof(int));
	cin=(double **)mxCalloc(dc, sizeof(double *));
	for (i=0;i<dc;i++){
	  cin[i]=(double *)mxCalloc(2, sizeof(double));
	}
	cout=(double **)mxCalloc(dc, sizeof(double *));
	for (i=0;i<dc;i++){
	  cout[i]=(double *)mxCalloc(2, sizeof(double));
	}
	vin=(double **)mxCalloc(dv, sizeof(double *));
	for (i=0;i<dv;i++){
	  vin[i]=(double *)mxCalloc(2, sizeof(double));
	}
	vout=(double **)mxCalloc(dv, sizeof(double *));
	for (i=0;i<dv;i++){
	  vout[i]=(double *)mxCalloc(2, sizeof(double));
	}  
    
 	srand((unsigned)time(NULL));
    
	rdis_min=1000000;
	
	for (jc=0;jc<N;jc++){
		rhat[jc]=0;
	}
	
	for (jv=0;jv<M;jv++){
		xhat[jv]=0;
	}
		
	istop=0;
	
	/*--------------------------- begin round --------------------------*/
	for (round=0;round<n_rounds;round++){
		
		count=0;
		
		for (jv=0;jv<M;jv++){
			pool[jv]=0;
		}
		for (j=0;j<K;j++){
			old[j]=0;
		}
	
		/* ----- iteration 0 ----- */
		istop++;
	    // beta and lambda
		for (j=0;j<K;j++){
	        jv=(int)floor((int)E[j]/dv);
	        jc=(int)floor(j/dc);
	        pool[jv]+=(r[jc]-rhat[jc]);
		}
	    for (jv=0;jv<M;jv++){
		    beta[jv]=xhat[jv]*xhat[jv]+pool[jv]*pool[jv]/(double)(dv*dv);
	   		if (beta[jv]<vs){
	        	beta[jv]=vs;
	       	}
	    	for (j=0;j<dv;j++){
	       		message[dv*jv+j][1]=beta[jv];
	           	message[dv*jv+j][0]=0;
	    	}
		}
		
		// bigo1(beta)			
	    for (jv=0;jv<M;jv++){
			bigo1[jv]=0;
			tempo[jv]=fabs(beta[jv]);
		}
		qsort(tempo, M, sizeof(double), compare_doubles);
		threshold=tempo[L];
		k=0;
		for (jv=0;jv<M;jv++){
			if (fabs(beta[jv])>threshold){			
				bigo1[jv]=1;
				k++;
			}
		}
		if (k<L){
			for (jv=0;jv<M;jv++){
				if (fabs(beta[jv])==threshold){			
					bigo1[jv]=1;
					k++;
				}
				if (k==L) break;
			} 
		}
		
		// theta
		for (jc=0;jc<N;jc++){
			for (j=0;j<dc;j++){
				kin=(int)E[dc*jc+j];
				cin[j][0]=message[kin][0];
				cin[j][1]=message[kin][1];
			}
			for (j=0;j<dc;j++){
				cout[j][0]=r[jc];
				cout[j][1]=0;
				for (k=1;k<dc;k++){
					cout[j][0]-=cin[(j+k)%dc][0];
					cout[j][1]+=cin[(j+k)%dc][1];				        
	       		}
	       		cout[j][1]+=sd*sd;
	       		if (cout[j][1]<vs){
	           		cout[j][1]=vs;
	       		}
	   		}
	   		for (j=0;j<dc;j++){
	   			kout=(int)E[dc*jc+j];	
	       		message[kout][0]=cout[j][0];
	       		message[kout][1]=cout[j][1];
	    	}
	   	}    
		
	   	
		/* ------------- iteration --------------- */
		for (i=1;i<n_iterations;i++){
		    istop++;
		    //printf("\n SuPrEM Iteration Number %d", istop);
		    
		    // VARIABLE NODE - beta, message lambda (v to c), belief-xhat
	        for (jv=0;jv<M;jv++){
		    	// read vin<-message    
		        for (j=0;j<dv;j++){
			        vin[j][0]=message[dv*jv+j][0];
			        vin[j][1]=message[dv*jv+j][1];
			        inactive[j]=old[dv*jv+j];
		    	}
		        active=0;
		        sum=0;
		        msum=0;
		        for (j=0;j<dv;j++){
			        if (inactive[j]==0){
				        active++;
			        	sum+=1/vin[j][1];
			        	msum+=vin[j][0]/vin[j][1]; 
		        	}
				}
				if (active>0){
					mu=msum/sum;
					sum+=1/beta[jv];
					if (active==1) {				
						beta[jv]=(pow(mu,2)+1/sum)/3;
					}
					else{
						beta[jv]=(pow(msum/sum,2)+1/sum)/3;
					}
		       		if (beta[jv]<vs){
			        	beta[jv]=vs;
	        	   	}
	       		}
	        	// lambda
	        	if (active>0){
		        	for (j=0;j<dv;j++){
				        sum=0;
				        msum=0;
				        for (k=1;k<dv;k++){
					    	sum+=1/vin[(j+k)%dv][1];
					    	msum+=vin[(j+k)%dv][0]/vin[(j+k)%dv][1];			    	
		           		}
		           		sum+=1/beta[jv];
		           		vout[j][1]=1/sum;
			           	vout[j][0]=msum/sum;
		        	}
	        	}
	        	else{
		        	for (j=0;j<dv;j++){
		           		vout[j][1]=vin[j][1];
			           	vout[j][0]=vin[j][0];
		        	}
	        	}
	        	for(j=0;j<dv;j++){
		        	if (vout[j][1]<vs) {
			        	vout[j][1]=vs;
		        	}
	        	}
	        	/* write vout->message*/
	        	for (j=0;j<dv;j++){
			        message[dv*jv+j][0]=vout[j][0];
			        message[dv*jv+j][1]=vout[j][1];
		    	}
		    	// belief and xhat
		    	sum=0;
		        msum=0;
		        for (j=0;j<dv;j++){
			    	sum+=1/vin[j][1];
			    	msum+=vin[j][0]/vin[j][1];			    	
	       		}
	       		sum+=1/beta[jv];
	       		bel[jv][1]=1/sum;
	           	bel[jv][0]=msum/sum;
	           	xhat[jv]=bel[jv][0];
	    	}
	    			
	    	
		    // bigo2
		    for (jv=0;jv<M;jv++){
				bigo2[jv]=0;
				tempo[jv]=beta[jv];
			}
			qsort(tempo, M, sizeof(double), compare_doubles);
			threshold=tempo[L];
			k=0;
			for (jv=0;jv<M;jv++){
				if (beta[jv]>threshold){			
					bigo2[jv]=1;
					k++;
				}
			}
			// tiebreak
			if (k<L){
				for (jv=0;jv<M;jv++){
					if (beta[jv]==threshold){			
						bigo2[jv]=1;
						k++;
					}
					if (k==L) break;
				} 
			}
			
			
			for (jv=0;jv<M;jv++){
				bigo[jv]=(bigo1[jv]||bigo2[jv]);
			}
					
			// bigo1			
		    for (jv=0;jv<M;jv++){
				bigo1[jv]=0;
				tempo[jv]=0;
			}
			k=0;
			for (jv=0;jv<M;jv++){
				if (bigo[jv]){
					tempo[k]=fabs(xhat[jv]);
					k++;
				}
			}
			qsort(tempo, k, sizeof(double), compare_doubles);
			threshold=tempo[L];
			k=0;
			for (jv=0;jv<M;jv++){
				if (bigo[jv] && fabs(xhat[jv])>threshold){			
					bigo1[jv]=1;
					k++;
				}
			}
			// tiebreak
			if (k<L){
				for (jv=0;jv<M;jv++){
					if (bigo[jv] && fabs(xhat[jv])==threshold){			
						bigo1[jv]=1;
						k++;
					}
					if (k==L) break;
				} 
			}
			
			/* set message-mean and xhat to zero if not in bigo1 */
			for (jv=0;jv<M;jv++){
				if (bigo1[jv]==0){
					for (j=0;j<dv;j++){
			        	message[dv*jv+j][0]=0;
		        	}
		        	xhat[jv]=0;
		        }
			}
			
						
		    // CHECK NODE - message theta (c to v)
		    for (jc=0;jc<N;jc++){
		    	if (rand( )<=(RAND_MAX-1)/2){
					for (j=0;j<dc;j++){
						kin=(int)E[dc*jc+j];
						cin[j][0]=message[kin][0];
						cin[j][1]=message[kin][1];
					}
					for (j=0;j<dc;j++){
						cout[j][0]=r[jc];
						cout[j][1]=0;
						for (k=1;k<dc;k++){
							cout[j][0]-=cin[(j+k)%dc][0];
							cout[j][1]+=cin[(j+k)%dc][1];				        
			       		}
			       		cout[j][1]+=sd*sd;
			       		if (cout[j][1]<vs){
			           		cout[j][1]=vs;
			       		}
		       		}
			    	for (j=0;j<dc;j++){
	   					kout=(int)E[dc*jc+j];
	   					old[kout]=0;	
	       				message[kout][0]=cout[j][0];
	       				message[kout][1]=cout[j][1];
	    			}		       		
	        	}
	        	else{
		        	for (j=0;j<dc;j++){
						kout=(int)E[dc*jc+j];
						old[kout]=1;
		        	}
	        	}
	    	}    	
		    	      	    		
	        // distortion
		    sum=0;
		    for (jc=0;jc<N;jc++){
			    rhat[jc]=0;
			    for (j=0;j<dc;j++){
				    jv=(int)floor((int)E[dc*jc+j]/dv);
				    rhat[jc]+=xhat[jv];
		    	}
		    	delta=fabs(rhat[jc]-r[jc]);
		    	sum+=delta*delta;
			}
			rdis=sqrt(sum);			
			if (rdis<rdis_min){
				rdis_min=rdis;
				count=0;
			}
			else{
				count++;
			}
				    			
			//stopping condition
			if (count==30 && cstop) break;
		}// iteration
		if (count==30 && cstop) break;
	}// round

	for (i=0;i<K;i++){
	  mxFree(message[i]);
	}
	mxFree(message);			
	for (i=0;i<M;i++){
	  mxFree(bel[i]);
	}	
	mxFree(bel);
	mxFree(beta);
	mxFree(pool);
	mxFree(bigo);
	mxFree(bigo1);
	mxFree(bigo2);
	mxFree(tempo);
	mxFree(old);
	mxFree(rhat);
    mxFree(inactive);
  	for (i=0;i<dc;i++){
	  mxFree(cin[i]);
	}	
	mxFree(cin);
  	for (i=0;i<dc;i++){
	  mxFree(cout[i]);
	}	
	mxFree(cout);
  	for (i=0;i<dv;i++){
	  mxFree(vin[i]);
	}	
	mxFree(vin);
  	for (i=0;i<dv;i++){
	  mxFree(vout[i]);
	}	
	mxFree(vout);    
	
	return;
	
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *xhat;
  double *r;
  double *E;
  int K, M, N, dc, dv;
  int L;
  double sd;
  int n_rounds;
  int n_iterations;
  int cstop;
    
  /*  check number of input arguments */
  if (nrhs!=12)
    mexErrMsgTxt("usage: xhat(M,1)=sprecover(r(N,1), E(K,1), K, M, N, dc, dv, L, sd, n_rounds, n_iterations, cstop)");

  /* get input */
  r=mxGetPr(prhs[0]);
  E=mxGetPr(prhs[1]);
  K=(int)mxGetScalar(prhs[2]);
  M=(int)mxGetScalar(prhs[3]);
  N=(int)mxGetScalar(prhs[4]);
  dc=(int)mxGetScalar(prhs[5]);
  dv=(int)mxGetScalar(prhs[6]);
  L=(int)mxGetScalar(prhs[7]);
  sd=mxGetScalar(prhs[8]);
  n_rounds=(int)mxGetScalar(prhs[9]);
  n_iterations=(int)mxGetScalar(prhs[10]);
  cstop=(int)mxGetScalar(prhs[11]);  
  plhs[0]=mxCreateDoubleMatrix(M,1,mxREAL);    
  xhat = mxGetPr(plhs[0]);

  /*  check input argument dimensions */
  if (mxGetM(prhs[0])*mxGetN(prhs[0])!=N)
    mexErrMsgTxt("usage: r should be length N");
  if (mxGetM(prhs[1])*mxGetN(prhs[1])!=K)
    mexErrMsgTxt("usage: E should be length K");
  if (mxGetM(prhs[2])*mxGetN(prhs[2])!=1)  
    mexErrMsgTxt("K should be a scalar");  
  if (mxGetM(prhs[3])*mxGetN(prhs[3])!=1)  
    mexErrMsgTxt("M should be a scalar");
  if (mxGetM(prhs[4])*mxGetN(prhs[4])!=1)  
    mexErrMsgTxt("N should be a scalar");
  if (mxGetM(prhs[5])*mxGetN(prhs[5])!=1)  
    mexErrMsgTxt("dc should be a scalar");
  if (mxGetM(prhs[6])*mxGetN(prhs[6])!=1)  
    mexErrMsgTxt("dv should be a scalar");  
  if (mxGetM(prhs[7])*mxGetN(prhs[7])!=1)  
    mexErrMsgTxt("L should be a scalar");
  if (mxGetM(prhs[8])*mxGetN(prhs[8])!=1)
    mexErrMsgTxt("sd should be a scalar");
  if (mxGetM(prhs[9])*mxGetN(prhs[9])!=1)
    mexErrMsgTxt("n_rounds should be a scalar");
  if (mxGetM(prhs[10])*mxGetN(prhs[10])!=1)
    mexErrMsgTxt("n_iterations should be a scalar");
  if (mxGetM(prhs[11])*mxGetN(prhs[11])!=1)
    mexErrMsgTxt("cstop should be a scalar");

  decode(xhat, r, E, K, M, N, dc, dv, L, sd, n_rounds, n_iterations, cstop);
  
  return;
  
}
