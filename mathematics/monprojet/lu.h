#include <stdio.h>
#include <stdlib.h>
//la méthode de résolution LU
int n, k=1 ; // taille du système
float *A ;  // matrice A
float *b ;    // vecteur b
float *x ;
float *L ;
float *y ;

int numero(int i, int j){

	return (i-1)*n + j-1;
}

void affiche_matrice(int n, float *A){

	int i;
	int j;
   printf("Affichage de la Matrice A : \n");
	for(i = 1; i<= n; i++){

		for(j = 1; j <= n; j++){
           printf("\t%.2f", A[numero(i,j)]);
		   printf("  ");
		}
		printf("\n");
	}

}

void affiche_vecteur(int n, float *b){
	int i;

  printf("Affichage du vecteur b : \n");

	for (i = 1; i<= n; i++){

		printf("\t%.2f",b[i]);
		printf("\n");
	}

}

void affiche_solution(int n, float *x){
	int i;

printf("les solutions de l'équation sont : \n");
	for (i = 1; i<= n; i++){

		printf("\tx%d = %.2f",i, x[i]);
		printf("\n");
	}

}

/* Saisie des éléments de la matrice A et b */
void saisie(int n,float *b,float *A )
{
     int i , j ;
     printf("Saisie de la matrice A ligne par ligne \n");

     for(i = 1 ; i <= n ; i++)
     {
        for( j = 1 ; j <= n ; j++)
        {
           printf("  A[%d][%d] : ",i,j);
           scanf("%f",&A[numero(i,j)]);
        }
     printf("\n");
     }

     printf("Saisie du vecteur b : \n");

     for(i = 1 ; i <= n ; i++)
     {
        printf("  b[%d] : ",i);
        scanf("%f",&b[i]);
        printf("\n");
     }
}

void descente(int n, float *y, float *b, float *L){

	int i;
	int j;
	y[1]=b[1];
	for (i= 2; i <=n; i++){
	y[i]= b[i];
		for (j=1 ; j <= i-1; j++)

			y[i]-=(L[numero(i,j)]*y[j]);

	}
}

void remonte(int n, float *x, float *b, float *A){

	int i;
	int j;
	x[n]=b[n]/A[numero(n,n)];
	for (i= n-1; i >=1; i--){
	x[i]= b[i];
			for (j=i+1 ; j <= n; j++){

				x[i]= x[i]-(A[numero(i,j)]*x[j]);
			}


	       x[i] = x[i]/A[numero(i,i)];

	}
}

void elimination(int k, float *A){

	int i,j;
	float r;
	for(i=k+1; i<=n; i++)
	{
		r= A[numero(i,k)]/A[numero(k,k)];
		for (j = k; j <= n; j++)
		{
			A[numero(i,j)]-=r*A[numero(k,j)];

		}

	}



}

void remplir_L(int k, float *A){

	int i;
	for(i=k+1; i<=n; i++)
		L[numero(i,k)]= A[numero(i,k)]/A[numero(k,k)];
	L[numero(k,k)]=1;

}








	






