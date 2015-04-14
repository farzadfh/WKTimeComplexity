/* Comuptue weighted distance from and correlation with identity for all permutations of a given length */

#include<stdio.h>
#include<stdlib.h>
#include <iostream>
using namespace std;

double dist(const int *,const double *,const int );
int getNextPerm(int *, const int);
void swap(int *, const int, const int);
void printIntArray(const int *, const int);
int find(const int *, const int, const int );

int main(int argc, char *argv[])
{
	// Determine permutation length
	int n;		// length of permutation
	if (argc != 2) {
		n = 4;	// default is 4 if length not given
	} else {
		n = atoi(argv[1]);
	}

	// Setting up variables
	int i,j;
	double *w = new double[n-1]; 	// The weight function. w[0], w[1], ..., w[n-2]
	for ( i = 0; i < n-1; ++i) { 	// Setup weight function: w[i] = 1/(i+1)
		w[i] = 1.0/(i+1);
		// w[i] = (i < 2) ? 1 : 0;
		printf("%f \n", w[i]);
	}
	int *pMax = new int[n];	// Permutation that gives maximum weighted Kendall Distance;
	int *p = new int[n];

	for ( i = 0; i < n ; ++i) {	// Setup pMax and p
		p[i] = i;
		pMax[i] = n-i-1;
	}
	double WKDMax = dist(pMax,w,n);	// Maximum weighted Kendall Distance
	double WKD;
	double WKDCorr;

	double corrTotal = 0;
	int count = 0;
	int done = 0;
	cout << "\nperm, distance, correlation\n";
	while (!done) {
		WKD = dist(p,w,n);
		WKDCorr = 1 - 2 * WKD / WKDMax;
		corrTotal += WKDCorr;
		++count;
		printIntArray(p,n);
		printf(", %f, %+f\n", WKD,WKDCorr);
		done = getNextPerm(p,n);	// Get next permutation
	}
	cout << corrTotal/count << "\n";
	delete p;
	delete pMax;
	delete w;
	return 0;
}

double dist(const int *p, const double *w, const int n)	
// Compute the distance between p (of length n) and id according to w
{	
	double d = 0;	// distance
	int *q = new int[n]; 	// copy of p
	for (int i = 0; i<n; ++i){
		q[i] = p[i];
	}

	for (int i = 0; i<n; ++i) {	
	// starting from i=0, each i is pushed to its proper position  
		int pos = find(q, i, n);	// find current position of i
		for (int j = pos; j>i; j--) {	// push it to its correct position
			swap(q, j, j-1);	// swap elements in positions j and j-1
			d = d + w[j-1];		// add the cost to distance
		}
	}
	delete q;	
	return d;
}

int getNextPerm(int *p, const int n)	
// find the permutation that lexigraphically follows p
{
    	int i = n - 1;
    	while ( (p[i-1] >= p[i]) && (i>0) ) 
	  	i = i-1;
	if ( i == 0 )
		return 1;
    	int j = n;
    	while (p[j-1] <= p[i-1]) 
	  	j = j-1;
      	
    	swap(p,i-1, j-1);    
	
    	i++; j = n;
    	while (i < j)
    	{
	  	swap(p,i-1, j-1);
	  	i++;
	  	j--;
    	}
	return 0;
}
				
void swap(int *p, const int i, const int j)
// Swaps elements in positions i and j in array p
{
	int temp = p[i];
	p[i] = p[j];
	p[j] = temp;
}

void printIntArray(const int *a, const int n){
// Print an array a of length n
	for (int j = 0; j < n; ++j){
		printf("%i", a[j]);
	}
}

int find(const int *a, const int v, const int n){
// Find element v in an array a of length n
	int i = 0;
	for ( ; i<n; i++) {
		if (a[i] == v)
			break;
	}
	return i;
}
