/* Compute weighted distance from and correlation with identity for all permutations of a given length */

#include<iostream>
using namespace std;
#include <algorithm>	// std::random_shuffle
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <stdint.h>

#ifdef WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

double dist(const int *, const int *, const double *,const int );
double dist2(const int *, const int *, const double *,const int );
int getNextPerm(int *, const int);
void printIntArray(const int *, const int);
int find(const int *, const int, const int );
int *getInv(const int *, const int);
int *permMultiply(const int *, const int *, const int);
void mergeAndCount(int *, int *, const int, const int);
void mergeSortAndCount(const int *p, int *invCount, const int n);
uint64_t GetTimeMs64();

int myrandom (int i) { return rand()%i;}

int main(int argc, char *argv[])
{
	int n = 5;
	if (argc>1) 
		n = atoi(argv[1]);

	// Setting up variables
	srand ( unsigned ( time(0) ) );
	int i,j;
	int *p = new int[n];
	int *q = new int[n];
	int *id = new int[n];
	for (i = 0; i<n; ++i) {id[i] = p[i] = q[i] = i;}
	random_shuffle (p, p+n, myrandom);
	random_shuffle (q, q+n, myrandom);
	int *reverse = new int[n];
	for ( i = 0; i<n; ++i)
		reverse[i] = id[n-1-i];
	double w[n-1]; 	// The weight function. w[0], w[1], ..., w[n-2]
	for ( i = 0; i < n-1; ++i){  	// Setup weight function: w[i] = 1/(i+1)
		w[i] = 1.0/(i+1);
	}
	
	double WKDMax;	// Maximum weighted Kendall Distance
	double WKD;
	double WKDCorr;
	
	uint64_t T0,T1,T2;
	T0 = GetTimeMs64();
	//WKDMax = dist(id,reverse,w,n);
	//WKD = dist(p,q,w,n);
	//WKDCorr = 1 - 2 * WKD / WKDMax;
	//printIntArray(p,n);
	//printf(", ");
	//printIntArray(q,n);
	//printf(": Distance = %f, Max Distance = %f, Correlation = %f\n", WKD,WKDMax,WKDCorr);
	
	T1 = GetTimeMs64();
	WKDMax = dist2(id,reverse,w,n);
	WKD = dist2(p,q,w,n);
	WKDCorr = 1 - 2 * WKD / WKDMax;
	printIntArray(p,n);
	printf(", ");
	printIntArray(q,n);
	printf(": Distance = %f, Max Distance = %f, Correlation = %f\n", WKD,WKDMax,WKDCorr);
	T2 = GetTimeMs64();

	//cout << "Simple: " << T1-T0 << "\n";
	cout << "Mergesort Time: " << T2-T1 << "\n";

	delete reverse;
	delete p;
	delete q;
	delete id;
	return 0;
}

double dist(const int *p1, const int *p2, const double *w, const int n)	
// Compute the distance between p1 and p2 (of length n) according to w: O(n^2)
{	
	double d = 0;	// distance
	int q[n];	// copy of p2
	for (int i = 0; i<n; ++i){
		q[i] = p2[i];
	}

	for (int i = 0; i<n; ++i) {	
	// starting from i=0, each p1[i] is pushed to its proper position  
		int pos = find(q, p1[i], n);	// find current position of p1[i] in q
		for (int j = pos; j>i; j--) {	// push it to its correct position
			swap(q[j], q[j-1]);	// swap elements in positions j and j-1
			d = d + w[j-1];		// add the cost to distance
		}
	}	
	return d;
}


double dist2(const int *p, const int *q, const double *w, const int n)	
// Compute the distance between p1 and p2 (of length n) according to w: O(n lg n)
{	double W[n]; 	// The sum weight function. w[1], w[2], ..., w[n-1]
	W[0] = 0;
	int i;
	for ( i = 1; i < n; ++i){ 
		W[i] = W[i-1] + w[i-1];
	}

	double d = 0;	// distance
	int L_newq_i;
	int *invp = getInv(p,n);
	int *newq = permMultiply(invp,q,n);
	//cout<<"p: ";printIntArray(p,n);cout<<"\n";
	//cout<<"q: ";printIntArray(q,n);cout<<"\n";
	//cout<<"inv p: ";printIntArray(invp,n);cout<<"\n";
	//cout<<"new q: ";printIntArray(q,n);cout<<"\n";
	int *invCount = new int[n];
	for (i = 0; i<n; ++i) invCount[i] = 0;
	
	mergeSortAndCount(newq, invCount, n);
	//cout<<"invCount: ";printIntArray(invCount,n);cout<<"\n";
	for ( i = 0; i<n; ++i){
		L_newq_i = (newq[i] + i + invCount[newq[i]]) / 2;
		d += ( 2.0 * W[L_newq_i] - W[i] - W[newq[i]]) / 2.0;
		//printf("i = %i, L_newq_i = %i, d=%f\n",i,L_newq_i,d);
	} 

	delete invp;
	delete newq;
	delete invCount;	
	return d;
}

				

void printIntArray(const int *a, const int n)
// Print an array a of length n
{
	for (int j = 0; j < n; ++j){
		printf("%i", a[j]);
	}
}

int find(const int *a, const int v, const int n)
// Find element v in an array a of length n
{
	int i = 0;
	for ( ; i<n; i++) {
		if (a[i] == v)
			break;
	}
	return i;
}

int *getInv(const int *p, const int n)
// compute the inverse of permutation p of length n
{
	int *inv = new int[n];
	for ( int i = 0; i<n; i++) 
		inv[p[i]] = i;
	return inv;
}


int *permMultiply(const int *p, const int *s, const int n)
// multiply p and s: i --> p[s[i]]
{
	int *q = new int[n];
	for (int i = 0; i < n; i++){
		q[i] = p[s[i]];
	}
	return q;
}


void mergeSortAndCount_(int *, int *, const int);
void mergeSortAndCount(const int *p, int *invCount, const int n)
{
	int *p_copy = new int[n];
	for (int i = 0; i<n; ++i) p_copy[i] = p[i];
	mergeSortAndCount_(p_copy, invCount, n);
	delete p_copy;

}

void mergeSortAndCount_(int *p, int *invCount, const int n)
{
	if ( n == 1 ) return;
	mergeSortAndCount_(p, invCount, n/2);
	mergeSortAndCount_(p+n/2, invCount, n-n/2);
	mergeAndCount(p, invCount, n/2,n);
}

void mergeAndCount(int *p, int *invCount, const int n1, const int n2)
// merge two lists pn[0..n1-1] and p[n1..n2-1] and add to count the number of inversions for each element caused by merging
{
	int i1 = 0;
	int i2 = n1;
	int i = 0;
	int j;
	int *tmp = new int[n2];
	while ( (i1 < n1) && (i2 < n2) ){
		if ( p[i1] < p[i2] ) {
			j = p[i1];
			invCount[j] += i2 - n1;   // how many have we had from the second list
			i1++;
		} else {
			j = p[i2];
			invCount[j] += n1 - i1;	 // how many are left in the first list
			i2++;
		}
		tmp[i++] = j;
	}
	while (i1 < n1) {
		j = p[i1];
		invCount[j] += n2 - n1;   // how many have we had from the second list
		tmp[i] = j;
		i++; i1++;
	}
	while (i2 < n2) {
		j = p[i2];
		tmp[i] = j;
		i++; i2++;
	}
	for (int i = 0; i < n2; ++i)
		p[i] = tmp[i];
	delete tmp;
}
		


uint64_t GetTimeMs64()
/* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both windows and linux. */
{
#ifdef WIN32
 /* Windows */
 FILETIME ft;
 LARGE_INTEGER li;

 /* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
  * to a LARGE_INTEGER structure. */
 GetSystemTimeAsFileTime(&ft);
 li.LowPart = ft.dwLowDateTime;
 li.HighPart = ft.dwHighDateTime;

 uint64_t ret = li.QuadPart;
 ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
 //ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */
 ret /= 10; /* From 100 nano seconds (10^-7) to 1 microsecond (10^-6) intervals */

 return ret;
#else
 /* Linux */
 struct timeval tv;

 gettimeofday(&tv, NULL);

 uint64_t ret = tv.tv_usec;
 //ret /= 1000; /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */

 //ret += (tv.tv_sec * 1000); // Adds the seconds (10^0) after converting them to milliseconds (10^-3)
 ret += (tv.tv_sec * 1000000); // Adds the seconds (10^0) after converting them to milliseconds (10^-6)
 return ret;
#endif
}		
