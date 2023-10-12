// ga_bool.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <time.h>                      // define time()
#include "randomc.h"                   // define classes for random number generators
//#include "userintf.cpp"                // define system specific user interface
#include "stdio.h"
#include "stdlib.h"
#include "mersenne.cpp"                // members of class TRandomMersenne
#pragma warning(disable:4996) //crt_warning hatasý engelleme kýsayolu    

#define N 256//67108864//16777216//4194304//65536//262144//1048576//262144//65536//16384//4096//1024//256//67108864//65536//16384//4096//65536//262144//67108864//16777216//4194304//1048576//   //length of truth table
#define n 8//26//24//22//16//18//20//18//16//14//12//10//8//26//16//14//12//16//26//24///22//20		//number of variables
#define NLb  120//33550336//8386560//2096128//32640//130816//523776//130816//32640//8128//2016//496//120//33550336//32640//8128///2016//32640//130816//2096128//523776    // max. nonlinearity
#define nlin 16//8192//4096//2048//256//512//1024//512//256//128//64//32//16//8192//256//128//64//256//512//4096//8192//1024		//2^(n/2)
#define scst 1      // 1: our cost 2: paper's cost
#define SLC	1 // 1: bizim selection 2: sugo'nun selection 
#define NSC 3 // 1: our k-NS 2: sugo's k-NS 3: our k-NSx

#define KNS 2 // number of neighbors
#define psize 1000  //initial population   
#define K  40  // parent size
#define Kb 100  // number of best selections
#define EPR 60  // % olasýlýk (Parent'daki "elite"lerin olasýlýðý)

#define Ka 2000000  // number of runs
#define maxgen 200 // number of iterations (generations)
#define ktour 3 //turnuva büyüklüðü
#define ELT	 EPR*K/100 //elite fonksiyon sayýsý (selection için)
#define KTR  K-ELT //k-turnava elde edilen sayýsý (selection için)       
#define STP 20 // stopping criteria
#define CPR 100 // Crossover probability 
#define nez nlin/2		//2^(n/2-1)
#define nuz N/2+nez //number of zeros of a bent function       

long double mutcost;
//long double *AC=(long double *)malloc(N*sizeof(long double));	
signed char *AT=(signed char *)malloc(N*sizeof(signed char));  
signed char *BT=(signed char *)malloc(N*sizeof(signed char));     
//long double *FW2=(long double *)malloc((N+1)*sizeof(long double));	
//long double *FW4=(long double *)malloc((N+1)*sizeof(long double));   	

int cek=0,BNL=-1,bnl=-1,BAC=1000000,bac=1000000;
int *bix=(int *)malloc(nez*sizeof(int));
int *cix=(int *)malloc(nez*sizeof(int));

int *JB=(int *)malloc(N*sizeof(int));
//int *ID0=(int *)malloc(N*sizeof(int));	
//int *ID1=(int *)malloc(N*sizeof(int));	
signed char *BNT=(signed char *)malloc(N*sizeof(signed char));	
//signed char *C=(signed char *)malloc(N*sizeof(signed char));	
//int *IX=(int *)malloc((nuz+1)*sizeof(int));	
//int *IX=(int *)malloc(N*sizeof(int));	
//int *IY=(int *)malloc((N/2)*sizeof(int));	
//int *IY=(int *)malloc(2*N*sizeof(int));	
//long double *ACt=(long double *)malloc(N*sizeof(long double));

int32 ir,ir1,ir2;                            // random integer number
int32 seed = (unsigned int) time(0);                // random seed

int *BE=(int *)malloc(Kb*nez*sizeof(int));
signed char *H=(signed char *)malloc(N*sizeof(signed char));

FILE *out=fopen("Best_Solutions.txt", "w");
FILE *outs=fopen("stats.txt", "w");

int _tmain(int argc, _TCHAR* argv[])   
{

	if (JB==0 || AT==0 || BT==0 || H==0 || BE==0 || BNT==0 || bix==0 || cix==0)// || FW2==0)// (FW4==0 || ACt==0 || AC==0 || IX==0 || IY==0 || C==0 || || ID0==0 || ID1==0
	{
		printf("Memory Allocation Error1!");
		return 0;
	}

	int weight(signed char *B);     
	void gen_bbool(int *B);   
	int bal_chk(signed char *B);  
	int findmaxwh(long double *tt);
	void gen_bent(signed char *B);
	void make_bal(signed char *B);
	long double sumsse(long double *FW);
	int find_min(long double *COST, int R);
	void tohex(int *TT, unsigned int *tt);
	//void breed(int *P1, int *P2, int *CH);
	//void ACOR(int *FW, int *AC);
	int findmaxac(int *tt);
	void mutation_rnd(signed char *B, long double *FW);
	//long double sumsse_cost3(int *FW);
	void select_ours(int p_size, long double *CSTF, long double *CSTB, int *BF, int *BP);
	//void select_sugo(int p_size, long double *CSTF, long double *CSTB, int *BF, int *BP);
	//void breed_new(int *P1, int *P2, int *CH);

	void fastwh(signed char *T, long double *FW);
	void fastwhld(long double *FW, long double *TTs);
	int acor(long double *FW);
	int anf(signed char *TT);	
	void breed_new(signed char *P1, signed char *P2, signed char *CH);

	time_t rawtime;
	struct tm * timeinfo;

	//signed char x;
	int *BA,*CHA,*BP,*BF;
	int *BB;
	int i,j,k,t,NL,cnt,J,ACR,kt,d;//,ct0,ct1;
	unsigned int cntm=0,idx,cnl=0,KA;
	long double *COST,Cr;//,I;

	int *NLA=(int *)malloc(maxgen*sizeof(int));
	signed char *B=(signed char *)malloc(N*sizeof(signed char));	
	//int *Br=(int *)malloc(N*sizeof(int));	
	long double *CSTB=(long double *)malloc(K*sizeof(long double));	
	long double *CSTC=(long double *)malloc(K*(K-1)/2*sizeof(long double));	
	long double *CSTF=(long double *)malloc((K+K*(K-1)/2)*sizeof(long double));	

	long double *FW=(long double *)malloc(N*sizeof(long double));	
	signed char *P1=(signed char *)malloc(N*sizeof(signed char));	
	signed char *P2=(signed char *)malloc(N*sizeof(signed char));	
	signed char *CH=(signed char *)malloc(N*sizeof(signed char));	

	if (P1==0 || P2==0 || FW==0 || CH==0 || B==0 || NLA==0 || CSTB==0 || CSTC==0 || CSTF==0)
	{
		printf("Memory Allocation Errorx!");
		return 0;
	}

	//int IDX[Kb],chk,p_size;
	//long double cost,CSTE[Kb],Maxi;
	//signed char *BE,IDW[psize];
	//int irk;

	//BE=(signed char *)malloc(Kb*N*sizeof(signed char));
	CHA=(int *)malloc((K+K*(K-1)/2)*nez*sizeof(int));	
	BB=(int *)malloc(maxgen*nez*sizeof(int));	
	BA=(int *)malloc(psize*nez*sizeof(int));	
	BP=(int *)malloc(K*nez*sizeof(int));	
	BF=(int *)malloc((K+K*(K-1)/2)*nez*sizeof(int));	
	COST=(long double *)malloc(psize*sizeof(long double));	

	if (JB==0 || FW==0 || BB==0 || AT==0 || BT==0 || CHA==0 || BA==0 || BP==0 || BF==0 || COST==0 || B==0 || P1==0 || P2==0 || CH==0 || H==0 || BE==0 || BNT==0)// || FW2==0)// (FW4==0 || || ACt==0 AC==0 || IX==0 || IY==0 || C==0 || || ID0==0 || ID1==0
	{
		printf("Memory Allocation Error!");
		return 0;
	}

	//i=0;
	//while (!feof(outh))
	//{
	//	fscanf(outh,"%d ",&x);
	//	*(H+i)=x;
	//	i=i+1;
	//}
	//printf("\ni=%d",i);
	//fclose(outh);


	//hadamard[0][0]=1;
 //   for(int x=1;x<n;x+=x){
 //       for(int i=0;i<x;i++){
 //           for(int j=0;j<x;j++){
 //               hadamard[i+x][j]=hadamard[i][j];
 //               hadamard[i][j+x]=hadamard[i][j];
 //               hadamard[i+x][j+x]=-hadamard[i][j];
 //           }
 //       }
	//
	H[0]=1;
    for(k=1;k<nlin;k+=k){
        for(int i=0;i<k;i++){
            for(int j=0;j<k;j++){
                H[(i+k)*nlin+j]=H[i*nlin+j];
                H[i*nlin+j+k]=H[i*nlin+j];
                H[(i+k)*nlin+j+k]=-H[i*nlin+j];
            }
        }
	}
	for (i=0;i<N;i++)
		H[i]=(1-H[i])/2;


	//I=0;
	//for (i=0;i<N+1;i++){
	//	FW2[i]=(I*I-N)*(I*I-N);
	//	I+=1;
	//}
	//I=0;
	//for (i=0;i<N+1;i++){
	//	FW4[i]=I*I*I*I;
	//	I+=1;
	//}
	for (KA=0;KA<Ka;KA++)
	{
	
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	printf ("\n(Start) Current local time and date: %s", asctime(timeinfo));
	fprintf (out,"\n(Start) Current local time and date: %s", asctime(timeinfo));
	fclose(out);
	out=fopen("Best_Solutions.txt", "a");

	// Initial population (bent or balanced)
	gen_bent(B);
	j=0;
	for (i=0;i<N;i++)
		if (B[i]==0)
			j=j+1;
	printf("\n0's:%d",j);
	j=0;
	for (i=0;i<N;i++)
		if (B[i]==1)
			j=j+1;
	printf("\n1's:%d",j);

	for (i=0;i<N;i++)
		BNT[i]=B[i];
	fastwh(B, FW);

	NL=N/2-findmaxwh(FW)/2;
	if (NL!=NLb)
	{
		printf("\nError..bent %d %d",i,NL);
		return 0;
	}
	//for (i=0;i<nuz;i++)
	//{
	//	printf("\ni=%d",i);
	//	IX[i]=0;
	//}
	//return 0;
	for (i=0;i<psize;i++)
	{
		make_bal(B);
		t=bal_chk(B);
		if (t!=N/2)
		{
			printf("\nError..bal %d; %d",i,weight(B));
			return 0;
		}
		//gen_bbool(B);
		//for (j=0;j<N;j++)
		//	//BA[i][j]=B[j];
		//	*(BA+i*N+j)=B[j];	//Store balanced functions into BA
		
		for (j=0;j<nez;j++)
			*(BA+i*nez+j)=bix[j];	//Store balanced functions into BA
		fastwh(B, FW);

///////////////////////////////////////////////////////

		NL=N/2-findmaxwh(FW)/2;
		if (bnl<NL)
			bnl=NL;

//////////////////////////////////////////////////////

		if (scst==1)
			COST[i]=sumsse(FW);		//Compute and store costs of balanced function into COST
		//else
		//	COST[i]=sumsse_cost3(FW);
		for (j=0;j<N;j++)
			B[j]=BNT[j];
		printf("\ni=%d (psize) => NL=%d cost=%f",i,NL,COST[i]);
	}
	
	//ct0=0;	// number of 0's of bent function
	//ct1=0;	// number of 1's of bent function
	//for (i=0;i<N;i++)
	//{
	//	
	//	if (BNT[i]==0)
	//	{
	//		ID0[ct0]=i;		// Store positions of 0's into ID0
	//		ct0=ct0+1;
	//	}
	//	else
	//	{
	//		ID1[ct1]=i;		// Store positions of 1's into ID1
	//		ct1=ct1+1;
	//	}
	//}

	// Selection
	if (SLC==1)
		select_ours(psize, COST, CSTB, BA, BP);
	//else
	//	select_sugo(psize, COST, CSTB, BA, BP);
	
	// CurrSol'a Parent'ý yerleþtir
	for (i=0;i<K;i++)
	{
		//for (j=0;j<N;j++)
		//	//BF[i][j]=BP[i][j];
		//	*(BF+i*N+j)=*(BP+i*N+j);
		for (j=0;j<nez;j++)
			*(BF+i*nez+j)=*(BP+i*nez+j);

		CSTF[i]=CSTB[i];  
	}
	kt=0;
	// maxgen = number of generations  
	for (J=0;J<maxgen;J++)
	{

		// Child'larý üret (crossover/breed)
		cnt=0;
		for (i=0;i<K-1;i++)
			for (j=i+1;j<K;j++)
			{

				cek=cek+10;
				seed=(unsigned int) time(0)+cek;
				TRandomMersenne rg(seed);            // make instance of random number generator
				ir = rg.IRandom(1,100);

				if (ir>100-CPR)
				{
					for (k=0;k<N;k++)	P1[k]=BNT[k];
					for	(k=0;k<nez;k++)	P1[BP[i*nez+k]]=BNT[BP[i*nez+k]]^1;


					for (k=0;k<N;k++)	P2[k]=BNT[k];
					for	(k=0;k<nez;k++)	P2[BP[j*nez+k]]=BNT[BP[j*nez+k]]^1;


					//for (k=0;k<N;k++)
					//{
					//	//P1[k]=BP[i][k];
					//	P1[k]=*(BP+i*N+k);
					//	//P2[k]=BP[j][k];
					//	P2[k]=*(BP+j*N+k);
					//}
					//breed(P1,P2,CH);
					
					breed_new(P1,P2,CH);
					
					//for (k=0;k<N;k++)
					//	//CHA[cnt][k]=CH[k];
					//	*(CHA+cnt*N+k)=CH[k];	// Store CH into CHA
					for (k=0;k<nez;k++)
						*(CHA+cnt*nez+k)=cix[k];	// Store CH into CHA
					fastwh(CH,FW);

///////////////////////////////////////////////////////

					NL=N/2-findmaxwh(FW)/2;
					if (bnl<NL)
						bnl=NL;

//////////////////////////////////////////////////////

					if (scst==1)
						CSTC[cnt]=sumsse(FW);	// Store cost of CH into CSTC
					//else
					//	CSTC[cnt]=sumsse_cost3(FW);
					cnt=cnt+1;
				}
			}
			
		//Child'lara mutation uygula
		for (i=0;i<cnt;i++)
		{
			for (k=0;k<N;k++)	CH[k]=BNT[k];
			for (k=0;k<nez;k++)	CH[CHA[i*nez+k]]=BNT[CHA[i*nez+k]]^1;
			//for (k=0;k<N;k++)
			//	//CH[k]=CHA[i][k];
			//	CH[k]=*(CHA+i*N+k);
			mutcost=CSTC[i];
			//printf("\nhello");
			mutation_rnd(CH,FW);
			//printf("\nhello");
			//for (k=0;k<N;k++)
			//	//CHA[i][k]=CH[k];
			//	*(CHA+i*N+k)=CH[k];
			for (k=0;k<nez;k++)
				*(CHA+i*nez+k)=cix[k];
			CSTC[i]=mutcost;
		}

		
		// Child'larý CurSol'a aktar (cost dahil)
		for (i=0;i<cnt;i++)
		{
			//for (j=0;j<N;j++)
			//	//BF[i+K][j]=CHA[i][j];
			//	*(BF+(i+K)*N+j)=*(CHA+i*N+j);
			for (j=0;j<nez;j++)
				*(BF+(i+K)*nez+j)=*(CHA+i*nez+j);
			CSTF[i+K]=CSTC[i];
		}

		// Selection
		if (SLC==1)
			select_ours(K+cnt, CSTF, CSTB, BF, BP);
		//else
		//	select_sugo(K+cnt, CSTF, CSTB, BF, BP);

		/////////////////////////////////////////////////////////////

		for (i=0;i<K;i++)
		{
			for (k=0;k<N;k++)	CH[k]=BNT[k];
			for (k=0;k<nez;k++)	CH[BP[i*nez+k]]=BNT[BP[i*nez+k]]^1;
			//for (k=0;k<N;k++)
			//	CH[k]=*(BP+i*N+k);
			fastwh(CH,FW);
			NL=N/2-findmaxwh(FW)/2;
			if (bnl<NL)
				bnl=NL;
		}

		////////////////////////////////////////////


		idx=find_min(CSTB,K);
		if (idx==-1)
		{
			printf("\nError..idx=%d",idx);
			return 0;
		}

		for (i=0;i<N;i++)	B[i]=BNT[i];
		for (i=0;i<nez;i++)	B[BP[idx*nez+i]]=BNT[BP[idx*nez+i]]^1;	// B: best function having minimum cost

		//for (i=0;i<N;i++)
		//	//B[i]=BP[idx][i];
		//	B[i]=*(BP+idx*N+i);		// B: best function having minimum cost

		//öncekilerle ayný mý?
		//for (i=0;i<N;i++)
		//	//BB[J][i]=B[i];
		//	*(BB+J*N+i)=B[i];		// Store B into BB

		for (i=0;i<nez;i++)
			*(BB+J*nez+i)=BP[idx*nez+i];		// Store B into BB

		Cr=CSTB[idx];		//Store min. cost
		//for (i=0;i<N;i++)
		//	Br[i]=B[i];		// Store B as Br
		for (i=0;i<nez;i++)
			cix[i]=BP[idx*nez+i];		// Store B as Br
		if (J>0)
		{
			j=0;
			for (i=0;i<nez;i++)
				if (*(BB+(J-1)*nez+i)==cix[i])
					j=j+1;
			if (j==nez)
				kt=kt+1;		// Number of repetitions
			else
				kt=0;
			if (kt==STP)
				kt=0;
		}

		//if (J>0)
		//{
		//	j=0;
		//	for (i=0;i<N;i++)
		//		//if (BB[J-1][i]==B[i])
		//		if (*(BB+(J-1)*N+i)==B[i])
		//			j=j+1;
		//	if (j==N)
		//		kt=kt+1;		// Number of repetitions
		//	else
		//		kt=0;
		//	if (kt==STP)
		//		kt=0;
		//}

		fastwh(B,FW);
		NL=N/2-findmaxwh(FW)/2;
		NLA[J]=NL;

		if (NL>BNL)
			BNL=NL;

		//ACOR(FW, AC);
		//ACR=findmaxac(AC);
		ACR=acor(FW);
		if (bac>ACR)
			bac=ACR;
		d=anf(B);
		k=0;
		for (i=0;i<K;i++)
			if (CSTB[i]==CSTB[idx])
				k=k+1;					// k: Number of functions with the same best cost

		//if (J%100==0)
		{

			printf("\nRun#:%d, G:%d (Cst=%10.0f; NL=%d AC=%d d=%d sf=%d sc=%d BNL=%d bnl=%d bac=%d): idx=%d",KA,J,CSTB[idx],NL,ACR,d,kt,k,BNL,bnl,bac,idx);
			fprintf(out,"\nRun#:%d, G:%d (Cst=%10.0f; NL=%d AC=%d d=%d sf=%d sc=%d BNL=%d bnl=%d bac=%d): idx=%d\n",KA,J,CSTB[idx],NL,ACR,d,kt,k,BNL,bnl,bac,idx);
		}
			
		if ((n==24 && (NL>=8386522 || ACR<=936)) || (n==18 && (NL>=130800 || ACR<=256)) || (n==20 && (NL>=523756 || ACR<=408)) || (n==22 && (NL>=2096096 || ACR<=600)) || (n==26 && (NL>=33550272 || ACR<=1408)) || (n==16 && (NL>=32628 || ACR<=152)) || (n==14 && (NL>=8120 || ACR<=88)) || (n==12 && (NL>=2010 || ACR<=40)) || (n==10 && (NL>=492 || ACR<=24)) || (n==8 && (NL>116 || ACR<=16))  || (n==6 && (NL>24 || ACR<=16)))
		{
			printf("\nRun#:%d, G:%d (Cst=%10.0f; NL=%d AC=%d d=%d sf=%d sc=%d BNL=%d bnl=%d bac=%d):",KA,J,CSTB[idx],NL,ACR,d,kt,k,BNL,bnl,bac);
			fprintf(out,"\nRun#:%d, G:%d (Cst=%10.0f; NL=%d AC=%d d=%d sf=%d sc=%d BNL=%d bnl=%d bac=%d): \n",KA,J,CSTB[idx],NL,ACR,d,kt,k,BNL,bnl,bac);
			for (i=0;i<N;i++)
				fprintf(out,"%d ",B[i]);
		}

		for (i=0;i<K;i++)    
		{
			//for (j=0;j<N;j++)
			//	//BF[i][j]=BP[i][j];
			//	*(BF+i*N+j)=*(BP+i*N+j); 
			for (j=0;j<nez;j++)
				*(BF+i*nez+j)=*(BP+i*nez+j);
			CSTF[i]=CSTB[i];
		}
		if (k==K || kt==STP)
		{
			for (i=0;i<psize;i++)
			{
				for (j=0;j<N;j++)
					B[j]=BNT[j];
				make_bal(B);
				t=bal_chk(B);
				if (t!=N/2)
				{
					printf("\nError..bal %d; %d",i,weight(B));
					return 0;
				}
				//gen_bbool(B);
				//for (j=0;j<N;j++)
				//	//BA[i][j]=B[j];
				//	*(BA+i*N+j)=B[j];	//Store balanced functions into BA
				for (j=0;j<nez;j++)
					*(BA+i*nez+j)=bix[j];	//Store balanced functions into BA
				fastwh(B, FW);

		///////////////////////////////////////////////////////

				NL=N/2-findmaxwh(FW)/2;
				if (bnl<NL)
					bnl=NL;

		//////////////////////////////////////////////////////

				if (scst==1)
					COST[i]=sumsse(FW);		//Compute and store costs of balanced function into COST
				//else
				//	COST[i]=sumsse_cost3(FW);
				//for (j=0;j<N;j++)
				//	B[j]=BNT[j];
			}
			// Selection
			if (SLC==1)
				select_ours(psize, COST, CSTB, BA, BP);
			//else
			//	select_sugo(psize, COST, CSTB, BA, BP);
			CSTB[K-1]=Cr;
			i=K-1;
			for (j=0;j<nez;j++)
				*(BP+i*nez+j)=cix[j];
			// CurrSol'a Parent'ý yerleþtir
			for (i=0;i<K;i++)
			{
				//for (j=0;j<N;j++)
				//	//BF[i][j]=BP[i][j];
				//	*(BF+i*N+j)=*(BP+i*N+j);
				for (j=0;j<nez;j++)
					*(BF+i*nez+j)=*(BP+i*nez+j);
				CSTF[i]=CSTB[i];  
			}
		}

	}

	j=-1;
	for (i=0;i<maxgen;i++)
		if (NLA[i]>j)
			j=NLA[i];

	if (j==116 && n==8)
		cnl=cnl+1;
	else if (j==492 && n==10)
		cnl=cnl+1;
	else if (j==2010 && n==12)
		cnl=cnl+1;
	else if (j==8120 && n==14)
		cnl=cnl+1;
	else if (j==32628 && n==16)
		cnl=cnl+1;
	else if (j==130800 && n==18)
		cnl=cnl+1;
	else if (j==523756 && n==20)
		cnl=cnl+1;
	else if (j==2096096 && n==22)
		cnl=cnl+1;
	else if (j==8386522 && n==24)
		cnl=cnl+1;
	else if (j==33550272 && n==26)
		cnl=cnl+1;
		

	printf("\nKA=%d Best NL=%d (%d)",KA,j,cnl);

	fprintf(outs,"%d ",j);
	fclose(outs);
	outs=fopen("stats.txt", "a");

	fclose(out);
	out=fopen("Best_Solutions.txt", "a");


	}
	fclose(out);
	fclose(outs);
	return 0;
}

int anf(signed char *TT)
{
	int i,j;
	for (i=0;i<N;i++)	AT[i]=TT[i];
	int t1=1,t2,Zx,l1,l2,lx,u1,u2,ux;
	for (int k=0;k<n;k++)
	{
		t1=t1*2;
		t2=t1/2;
		Zx=N/t1;
		for (int b=0;b<Zx;b++)
		{
			i=2*b;
			l1=t2*i;	u1=t2*(i+1);
			l2=t2*(i+1);u2=t2*(i+2);
			lx=t1*b;	ux=t1*(b+1);
			for (int j=l1;j<u1;j++)
			{
				BT[lx]=AT[j];
				lx=lx+1;
			}
			for (int j=l2;j<u2;j++)
			{
				BT[lx]=AT[j]^AT[l1];
				lx=lx+1;l1=l1+1;
			}
		}
		for (j=0;j<N;j++)
			AT[j]=BT[j];
	}
	int ad=0,d;
	for (i=0;i<N;i++)
		if (AT[i]!=0)
		{
			d=0;
			for (j=0;j<n;j++)
				if ((i&(1<<j))!=0)
					d=d+1;
			if (d>ad)
				ad=d;

		}
	return ad;
}

void fastwh(signed char *T, long double *FW)
{	
	int i,j,i1,i2,i3,k1=2,k2=N/2,k3=1,L1;
	long double temp1,temp2;
	for (i=0;i<N;i++)
		FW[i]=1-2*((long double) T[i]);

	for (i1=0;i1<n;i1++)  
	{
	   L1=1;
	   for (i2=0;i2<k2;i2++)
	   {
		  for (i3=0;i3<k3;i3++)
		  {
			 i=i3+L1-1; j=i+k3; 
		     temp1= FW[i]; temp2 = FW[j]; 
			 FW[i]=temp1+temp2;
		     FW[j]=temp1-temp2;
		  }
	      L1=L1+k1; 
	   }
	   k1=k1*2; k2=k2/2; k3=k3*2;
	}
}

void fastwhld(long double *FW, long double *TTs)
{	
	int i,j,i1,i2,i3,k1=2,k2=N/2,k3=1,L1;
	long double temp1,temp2;
	for (i=0;i<N;i++)
		FW[i]=TTs[i];
	for (i1=0;i1<n;i1++)  
	{
	   L1=1;
	   for (i2=0;i2<k2;i2++)
	   {
		  for (i3=0;i3<k3;i3++)
		  {
			 i=i3+L1-1; j=i+k3; 
		     temp1= FW[i]; temp2 = FW[j]; 
			 FW[i]=temp1+temp2;
		     FW[j]=temp1-temp2;
		  }
	      L1=L1+k1; 
	   }
	   k1=k1*2; k2=k2/2; k3=k3*2;
	}
}

int acor(long double *FW)
{
	int i,j,i1,L1,i2,i3,k1=2,k2=N/2,k3=1;
	long double temp1,temp2;
	for (i=0;i<N;i++)
		FW[i]=FW[i]*FW[i];

	for (i1=0;i1<n;i1++)  
	{
	   L1=1;
	   for (i2=0;i2<k2;i2++)
	   {
		  for (i3=0;i3<k3;i3++)
		  {
			 i=i3+L1-1; j=i+k3; 
		     temp1= FW[i]; temp2 = FW[j]; 
			 FW[i]=temp1+temp2;
		     FW[j]=temp1-temp2;
		  }
	      L1=L1+k1; 
	   }
	   k1=k1*2; k2=k2/2; k3=k3*2;
	}

	long double D,Maxi=-1;
	for (i=1;i<N;i++)
	{
		D=FW[i];
		if (FW[i]<0)
			D=-FW[i];
		if (D>Maxi)
			Maxi=D;
	}
	Maxi=Maxi/N;
	return ((int) Maxi);
}


//int acor(int *FW)
//{
//	int i;
//	for (i=0;i<N;i++)
//		ACt[i]=((long double) FW[i])*((long double) FW[i]);
//	fastwhld(AC,ACt);
//
//	long double D,Maxi=-1;
//	for (i=1;i<N;i++)
//	{
//		D=AC[i];
//		if (AC[i]<0)
//			D=-AC[i];
//		if (D>Maxi)
//			Maxi=D;
//	}
//	Maxi=Maxi/N;
//	return ((int) Maxi);
//}


int find_min(long double *COST, int R)
{
	int i,idxmin=-1;
	long double Maxi=1.6e+50;
	for (i=0;i<R;i++)
		if (COST[i]<Maxi)
		{
			Maxi=COST[i];
			idxmin=i;
		}
	return idxmin;
}

//void select_sugo(int p_size, long double *CSTF, long double *CSTB, signed char *BF, signed char *BP)
//{
//	int i,j,k,IDW[psize],IDX[ELT],chk,irk;
//	long double Maxi;
//	
//	for (i=0;i<p_size;i++)
//		IDW[i]=0;
//
//	// select best ELT functions (elitism - % 60)
//	for (i=0;i<ELT;i++)
//	{
//		IDX[i]=find_min(CSTF,p_size);
//		if (IDX[i]==-1)
//		{
//			printf("\nError..IDX[%d]=%d",i,IDX[i]);
//			return;
//		}
//		CSTB[i]=CSTF[IDX[i]]; // Store costs of elits into CSTB
//		CSTF[IDX[i]]=1.7e+50;
//		//for (j=0;j<N;j++)
//		//	//BP[i][j]=BF[IDX[i]][j];
//		//	*(BP+i*N+j)=*(BF+IDX[i]*N+j);
//		for (j=0;j<nez;j++)
//			*(BP+i*nez+j)=*(BF+IDX[i]*nez+j);
//		IDW[IDX[i]]=1;
//	}
//	// make k-tournament
//	for (i=0;i<KTR;i++)
//	{
//		Maxi=1.6e+50;
//		for (k=0;k<ktour;k++)
//		{
//			chk=0;
//			while (chk==0)
//			{
//				cek=cek+10;
//				seed=(unsigned int) time(0)+cek;
//				TRandomMersenne rg(seed);            // make instance of random number generator
//				ir = rg.IRandom(0,p_size-1);
//				if (IDW[ir]==0)
//					chk=1;
//			}
//			if (CSTF[ir]<Maxi)
//			{
//				Maxi=CSTF[ir];
//				irk=ir;
//			}
//		}
//		//for (j=0;j<N;j++)
//		//	//BP[i+ELT][j]=BF[irk][j];
//		//	*(BP+(i+ELT)*N+j)=*(BF+irk*N+j);
//
//		for (j=0;j<nez;j++)
//			*(BP+(i+ELT)*nez+j)=*(BF+irk*nez+j);
//
//		CSTB[i+ELT]=Maxi;	// Add costs of k-tournament result into CSTB
//		IDW[irk]=1;
//	}
//}

void select_ours(int p_size, long double *CSTF, long double *CSTB, int *BF, int *BP)
{
	int i,j,k,IDX[Kb],chk,irk;

	long double CSTE[Kb],Maxi;
	int IDW[Kb];

	for (i=0;i<Kb;i++)
		IDW[i]=0;
	// select best Kb functions (elitism)
	for (i=0;i<Kb;i++)
	{
		//printf("\nhx i=%d",i);
		IDX[i]=find_min(CSTF,p_size);
		if (IDX[i]==-1)
		{
			printf("\nError..IDX[%d]=%d",i,IDX[i]);
			return;
		}
		CSTE[i]=CSTF[IDX[i]];
		CSTF[IDX[i]]=1.7e+50;
		//for (j=0;j<N;j++)
		//	//BE[i][j]=BF[IDX[i]][j];
		//	//BE[i][j]=*(BF+IDX[i]*N+j);			//MEMORY????
		//	*(BE+i*N+j)=*(BF+IDX[i]*N+j);			//MEMORY????
		for (j=0;j<nez;j++)
			*(BE+i*nez+j)=*(BF+IDX[i]*nez+j);			//MEMORY????
	}
	// make k-tournament
	for (i=0;i<K;i++)
	{
		Maxi=1.6e+50;
		for (k=0;k<ktour;k++)
		{
			chk=0;
			while (chk==0)
			{
				cek=cek+10;
				seed=(unsigned int) time(0)+cek;
				TRandomMersenne rg(seed);            // make instance of random number generator
				ir = rg.IRandom(0,Kb-1);
				if (IDW[ir]==0)
					chk=1;
			}
			if (CSTE[ir]<Maxi)
			{
				Maxi=CSTE[ir];
				irk=ir;
			}
		}
		//for (j=0;j<N;j++)
		//	//BP[i][j]=BE[irk][j];
		//	//*(BP+i*N+j)=BE[irk][j];
		//	*(BP+i*N+j)=*(BE+irk*N+j);
		for (j=0;j<nez;j++)
			*(BP+i*nez+j)=*(BE+irk*nez+j);
		CSTB[i]=Maxi;
		IDW[irk]=1;
	}
}

int bal_chk(signed char *B)
{
	int i,s=0;
	for (i=0;i<N;i++)
		s=s+B[i];
	return s;
}

int weight(signed char *B)
{
	int i,s=0;
	for (i=0;i<N;i++)
		s=s+B[i];
	return s;
}

void breed_new(signed char *P1, signed char *P2, signed char *CH)
{
	int i,j,k,I1[nez],I2[nez],PI[nez],c1=0,c2=0,chk,kt;
	for (i=0;i<N;i++)
	{
		if (P1[i]==P2[i])
			CH[i]=P1[i];
		else
			CH[i]=2;
	}
	for (i=0;i<N;i++)
	{
		if (P1[i]==1 && BNT[i]==0)
		{
			I1[c1]=i;
			//printf("\nc1=%d I1[c1]=%d",c1,I1[c1]);
			c1=c1+1;
		}
		if (P2[i]==1 && BNT[i]==0)
		{
			I2[c2]=i;
			//printf("\nc1=%d I1[c1]=%d",c1,I1[c1]);
			c2=c2+1;
		}
	}
	//printf("\nc1=%d c2=%d",c1,c2);
	kt=0;
	for (i=0;i<nez;i++)
	{
		//printf("\ni=%d",i);
		chk=0;
		for (j=0;j<nez;j++)
			if (I1[i]==I2[j])
				chk=1;
		if (chk==0)
		{
			//printf(" kt=%d; %d",kt,I1[i]);
			//printf("kt=%d; %d %d %d",kt,JB[kt],I1[i],CH[I1[i]]);
			JB[kt]=I1[i];
			kt=kt+1;
			CH[I1[i]]=1;
		}
	}
	
	for (i=0;i<nez;i++)
	{
		chk=0;
		for (j=0;j<nez;j++)
			if (I2[i]==I1[j])
				chk=1;
		if (chk==0)
		{
			JB[kt]=I2[i];
			kt=kt+1;
			CH[I2[i]]=1;
		}
	}
	
	for (i=0;i<kt/2;i++)
	{
		chk=0;
		while (chk==0)
		{
			chk=1;
			cek=cek+10;
			seed=(unsigned int) time(0)+cek;
			TRandomMersenne rg(seed);            // make instance of random number generator
			ir = rg.IRandom(0,kt-1);
			for (k=0;k<i;k++)
				if (PI[k]==ir)
					chk=0;
			if (chk==1)	PI[i]=ir;
		}
	}
	for (i=0;i<kt/2;i++)
		CH[JB[PI[i]]]=0;
	for (i=0;i<N;i++)
		if (CH[i]==2)
			printf("\nError...breed_new");
	j=bal_chk(CH);
	if (j!=N/2)
	{
		printf("\nError..bal...breed_new %d; %d",i,weight(CH));
		return;
	}
	j=0;
	for (i=0;i<N;i++)
		if (CH[i]!=BNT[i])
		{
			//printf("j=%d",j);
			cix[j]=i;
			j=j+1;
		}
	if (j!=nez)
		printf("\nError..cix");

}

//void breed(int *P1, int *P2, int *CH)
//{
//	int hd=0,i,cnt0=0,cnt1=0;
//	for (i=0;i<N;i++)
//		hd=hd+P1[i]^P2[i];
//	//printf("\nhd=%d",hd);
//	if (hd>N/2)
//		for (i=0;i<N;i++)
//			P1[i]=P1[i]^1;
//	for (i=0;i<N;i++)
//	{
//		if (P1[i]==P2[i])
//		{
//			CH[i]=P1[i];
//			if (CH[i]==1)
//				cnt1=cnt1+1;
//			else
//				cnt0=cnt0+1;
//		}
//		else
//			CH[i]=2;
//	}
//	//if (cnt0+cnt1>200)
//	//	printf("\np1==p2 => %d",cnt0+cnt1);
//
//	for (i=0;i<N;i++)
//	{
//		if (CH[i]==2)
//		{
//			if (cnt0==(N/2))
//				CH[i]=1;
//			else if (cnt1==(N/2))
//				CH[i]=0;
//			else
//			{
//				cek=cek+10;
//				seed=(unsigned int) time(0)+cek;
//				TRandomMersenne rg(seed);            // make instance of random number generator
//				ir = rg.IRandom(0,1);
//				CH[i]=ir;
//			}
//			if (CH[i]==1)
//				cnt1=cnt1+1;
//			else
//				cnt0=cnt0+1;
//			//printf("\ni=%d cnt0=%d cnt1=%d",i,cnt0,cnt1);
//		}
//	}
//}

long double sumsse(long double *FW){
	
	int i;
	long double sum=0,tmp;
	
	for (i=0;i<N;i++){
		tmp=FW[i]*FW[i];
		sum=sum+(tmp-N)*(tmp-N);
		//if (FW[i]<0)
		//	sum=sum+FW2[-FW[i]];
		//else
		//	sum=sum+FW2[FW[i]];
	}
	return sum;
	
}

//long double sumsse_cost3(int *FW){
//	
//	int i;
//	long double sum=0;
//	
//	for (i=0;i<N;i++){
//		if (FW[i]<0)
//			sum=sum+FW4[-FW[i]];
//		else
//			sum=sum+FW4[FW[i]];
//	}
//	return sum;
//	
//}

void gen_bbool(int *B)
{	
	int i,k,PI[N/2],chk;

	for (i=0;i<N/2;i++)
	{
		chk=0;
		while (chk==0)
		{
			chk=1;
			cek=cek+10;
			seed=(unsigned int) time(0)+cek;
			TRandomMersenne rg(seed);            // make instance of random number generator
			ir = rg.IRandom(0,N-1);
			if (ir>=N)	chk=0;
			for (k=0;k<i;k++)
				if (PI[k]==ir)
					chk=0;
			if (chk==1)	PI[i]=ir;
		}
	}
	for (i=0;i<N;i++)
		B[i]=0;
	for (i=0;i<N/2;i++)
		B[PI[i]]=1;
}

//void make_bal(signed char *B)
//{
//	int i,k,PI[nez],chk;
//	//printf("\nhello");
//	for (i=0;i<nez;i++)
//	{
//		chk=0;
//		while (chk==0)
//		{
//			chk=1;
//			cek=cek+10;
//			seed=(unsigned int) time(0)+cek;
//			TRandomMersenne rg(seed);            // make instance of random number generator
//			ir = rg.IRandom(0,nuz-1);
//			for (k=0;k<i;k++)
//				if (PI[k]==ir)
//					chk=0;
//			if (chk==1)	PI[i]=ir;
//		}
//	}
//	k=0;
//	for (i=0;i<N;i++)
//		if (B[i]==0)
//		{
//			//printf("\n%d",k);
//			IX[k]=i;
//			k=k+1;
//		}
//	//printf("\nhello");
//	for (i=0;i<nez;i++)
//	{
//		bix[i]=IX[PI[i]];
//		B[IX[PI[i]]]=1;
//	}
//
//}

void make_bal(signed char *B)
{
	int i,chk,cbx=0;//,PI[nez],k;
	//printf("\nhello");
	for (i=0;i<nez;i++)
	{
		chk=0;
		while (chk==0)
		{
			chk=1;
			cek=cek+10;
			seed=(unsigned int) time(0)+cek;
			TRandomMersenne rg(seed);            // make instance of random number generator
			ir = rg.IRandom(0,N-1);
			if (B[ir]==0)
			{
				B[ir]=1;
				bix[cbx]=ir;
				cbx=cbx+1;
			}
			else
				chk=0;
			//for (k=0;k<i;k++)
			//	if (PI[k]==ir)
			//		chk=0;
			//if (chk==1)	PI[i]=ir;
		}
	}
	//k=0;
	//for (i=0;i<N;i++)
	//	if (B[i]==0)
	//	{
	//		//printf("\n%d",k);
	//		IX[k]=i;
	//		k=k+1;
	//	}
	//printf("\nhello");
	//for (i=0;i<nez;i++)
	//{
	//	bix[i]=IX[PI[i]];
	//	B[IX[PI[i]]]=1;
	//}

}

void gen_bent(signed char *B)
{
	int i,j,k,PI[nlin],chk;

	for (i=0;i<nlin;i++)
	{
		chk=0;
		while (chk==0)
		{
			chk=1;
			cek=cek+10;
			seed=(unsigned int) time(0)+cek;
			TRandomMersenne rg(seed);            // make instance of random number generator
			ir = rg.IRandom(0,nlin-1);
			for (k=0;k<i;k++)
				if (PI[k]==ir)
					chk=0;
			if (chk==1)	PI[i]=ir;
		}
	}
	//printf("\n");
	//for (i=0;i<16;i++)
	//	printf("%d ",PI[i]);
	for (i=0;i<nlin;i++)
		for (j=0;j<nlin;j++)
			B[i*nlin+j]=H[PI[i]*nlin+j];
}


int findmaxwh(long double *tt)
{
	int i;
	long double D,Maxi=-1;
	for (i=0;i<N;i++)
	{
		D=tt[i];
		if (tt[i]<0)
			D=-tt[i];
		if (D>Maxi)
			Maxi=D;
	}
	return ((int) Maxi);
}

void mutation_rnd(signed char *B, long double *FW)
{
	int i,j,chk,NL;

	cek=cek+10;
	seed=(unsigned int) time(0)+cek;
	TRandomMersenne rg(seed);            // make instance of random number generator
	ir = rg.IRandom(0,N-1);

	if (ir<2)
	{
		chk=0;
		while (chk==0)
		{
			chk=1;
			cek=cek+10;
			seed=(unsigned int) time(0)+cek;
			TRandomMersenne rgx(seed);            // make instance of random number generator
			ir1 = rgx.IRandom(0,N-1);

			cek=cek+10;
			seed=(unsigned int) time(0)+cek;
			TRandomMersenne rgy(seed);            // make instance of random number generator
			ir2 = rgy.IRandom(0,N-1);

			if ((B[ir1]==0 && B[ir2]==1) || (B[ir1]==1 && B[ir2]==0))
			{
				B[ir1]=B[ir1]^1;
				B[ir2]=B[ir2]^1;
			}
			else
				chk=0;


		}

		fastwh(B,FW);

///////////////////////////////////////////////////////
		NL=N/2-findmaxwh(FW)/2;
		if (bnl<NL)
			bnl=NL;

//////////////////////////////////////////////////////

		if (scst==1)
			mutcost=sumsse(FW);
		//else
		//	mutcost=sumsse_cost3(FW);
	}
	j=0;
	for (i=0;i<N;i++)
		if (B[i]!=BNT[i])
		{
			cix[j]=i;
			j=j+1;
		}

}



//void mutation_rnd(signed char *B, long double *FW)
//{
//	int i,j,cx=0,cy=0,IZ[nez];
//
//	for (i=0;i<N;i++)
//	{
//		if (B[i]==1 && BNT[i]==0)
//		{
//			IZ[cx]=i;
//			cx=cx+1;
//		}
//		else if (B[i]==0 && BNT[i]==0)
//		{
//			IY[cy]=i;
//			cy=cy+1;
//		}
//		//printf("\ncx=%d cy=%d",cx,cy);
//	}
//	for (i=0;i<N;i++)
//		C[i]=B[i];
//
//	cek=cek+10;
//	seed=(unsigned int) time(0)+cek;
//	TRandomMersenne rg(seed);            // make instance of random number generator
//	ir = rg.IRandom(0,N-1);
//
//	if (ir<2)
//	{
//		cek=cek+10;
//		seed=(unsigned int) time(0)+cek;
//		TRandomMersenne rgx(seed);            // make instance of random number generator
//		ir = rgx.IRandom(0,nez-1);
//		B[IZ[ir]]=B[IZ[ir]]^1;
//
//		cek=cek+10;
//		seed=(unsigned int) time(0)+cek;
//		TRandomMersenne rgy(seed);            // make instance of random number generator
//		ir = rgy.IRandom(0,N/2-1);
//		B[IY[ir]]=B[IY[ir]]^1;
//		fastwh(B,FW);
//
/////////////////////////////////////////////////////////
//		int NL;
//		NL=N/2-findmaxwh(FW)/2;
//		if (bnl<NL)
//			bnl=NL;
//
////////////////////////////////////////////////////////
//
//		if (scst==1)
//			mutcost=sumsse(FW);
//		//else
//		//	mutcost=sumsse_cost3(FW);
//	}
//	j=0;
//	for (i=0;i<N;i++)
//		if (B[i]!=BNT[i])
//		{
//			cix[j]=i;
//			j=j+1;
//		}
//
//}


void tohex(int *TT, unsigned int *tt){
	
	int i,j;
	
	for (i=0;i<N/32;i++){
		tt[i]=0;
		for (j=31;j>=0;j--)
			tt[i]=(TT[i*32+j]<<(31-j))^tt[i];
	}
}

