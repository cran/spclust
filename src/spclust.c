#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void distbc(int *n, int *nmrk, double *prmat1, double *prmat2, double *distmat)
{
	double k, tot;
	int a, b, c;
	for (a=0; a < *n; a++)
	for (b=0; b < *n; b++)
	{
		k=0;
		tot=0;
		for (c=0; c < *nmrk; c++) 
		if ((prmat1[a*(*nmrk)+c] >= 0) && (prmat1[b*(*nmrk)+c] >= 0)) { 
			k += prmat1[a*(*nmrk)+c]*prmat1[b*(*nmrk)+c] + prmat2[a*(*nmrk)+c]*prmat2[b*(*nmrk)+c];
			tot++;
		}
	      	distmat[a*(*n)+b] = sqrt(1.0-k/tot);
	}
}

void distf2(int *n, int *nmrk, double *prmat1, double *prmat2, double *prmat3, double *distmat)
{
	double k, tot;
	int a, b, c;
	for (a=0; a < *n; a++)
	for (b=0; b < *n; b++)
	{ 
		k=0;
		tot=0;
		for (c=0; c < *nmrk; c++)
		if ((prmat1[a*(*nmrk)+c] >= 0) && (prmat1[b*(*nmrk)+c] >= 0))
		{
		  k += 2*prmat1[a*(*nmrk)+c]*prmat1[b*(*nmrk)+c] + 2*prmat2[a*(*nmrk)+c]*prmat2[b*(*nmrk)+c] + 2*prmat3[a*(*nmrk)+c]*prmat3[b*(*nmrk)+c]
		  	+ prmat1[a*(*nmrk)+c]*prmat2[b*(*nmrk)+c] + prmat1[b*(*nmrk)+c]*prmat2[a*(*nmrk)+c]
			+ prmat2[a*(*nmrk)+c]*prmat3[b*(*nmrk)+c] + prmat2[b*(*nmrk)+c]*prmat3[a*(*nmrk)+c];
		  tot++;
		}
		distmat[a*(*n)+b] = sqrt(1-k/2.0/tot);
	}
}

void distmagic(int *n, int *nmrk, int *nfounders, double *prmat, double *distmat)
{
	double k, tot;
	int a, b, c;
	for (a=0; a < *n; a++)
	for (b=0; b < *n; b++)
	{ 
		k=0;
		tot=0;
		for (c=0; c < *nmrk; c++)
		if ((prmat[a*(*nmrk)+c] >= 0) && (prmat[b*(*nmrk)+c] >= 0))
		{
			k += prmat[a*(*nmrk)+c]*prmat[b*(*nmrk)+c];
			tot++;
		}
		distmat[a*(*n)+b] = sqrt(1.0-(*nfounders)*k/tot);
	}
}
