#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int n; // taille du système
float *A ;  // matrice A
float *b ;    // vecteur b
float *x ;
float *B ;
float *y ;

int num(int i, int j){

	return (i-1)*n + j-1;
}

void afficher_matrice(int n, float *A){

	int i;
	int j;
	for(i = 1; i<= n; i++){

		for(j = 1; j <= n; j++){
			printf("\t%.2f", A[num(i,j)]);
		   printf("  ");
		}
		printf("\n");
	}
}

void afficher_vecteur(int n, float *b){
	int i;

  printf("\nAffichage du vecteur b : \n");

	for (i = 1; i<= n; i++){
		printf("\t%.2f",b[i]);
		printf("\n");
	}
}

void affiche_sol(int n, float *x){
	int i;

printf("\nles solutions de l'équation sont : \n");
	for (i = 1; i<= n; i++){

		printf("\tx%d = %.2f",i, x[i]);
		printf("\n");
	}
}
void saisir(int n,float *b,float *A )
{
     int i , j ;
     printf("\nsaisir de la matrice A ligne par ligne : \n");

     for(i = 1 ; i <= n ; i++)
     {
        for( j = 1 ; j <= n ; j++)
        {
           printf("  A[%d][%d] : ",i,j);
           scanf("%f",&A[num(i,j)]);
        }
     printf("\n");
     }

     printf("\nsaisir du vecteur b : \n");

     for(i = 1 ; i <= n ; i++)
     {
        printf("  b[%d] : ",i);
        scanf("%f",&b[i]);
        printf("\n");
     }
}

void descente_cholesky(int n, float *y, float *b, float *B){

	int i;
	int j;
	y[1]=b[1]/B[num(1,1)];
	for (i= 2; i <=n; i++)
	{
		y[i]= b[i];
		for (j=1 ; j <= i-1; j++)
			y[i]-=B[num(i,j)]*y[j];

		y[i]=y[i]/B[num(i,i)];
	}
}
void remonte_cholesky(int n, float *x, float *y, float *B){

	int i;
	int j;
	x[n]=y[n]/B[num(n,n)];
	for (i= n-1; i >=1; i--)
	{
		x[i]= y[i];
		for (j=i+1 ; j <= n; j++)
			x[i]-= B[num(j,i)]*x[j];

	    x[i] = x[i]/B[num(i,i)];

	}
}



void facto_cholesky(float *A, float *B){

	int i,k,j;

		for (j=1; j<=n; j++)
		{
			B[num(j,j)]=A[num(j,j)];

			for (k = 1; k <= j-1; k++)
				B[num(j,j)] -= pow(B[num(j,k)],2);

		B[num(j,j)] = sqrt(B[num(j,j)]);
		for (i = j+1; i <=n; i++)
		{
			B[num(i,j)]=A[num(i,j)];
			for (k= 1; k <= j-1; k++)
				B[num(i,j)] -= B[num(i,k)]*B[num(j,k)];

			B[num(i,j)]= B[num(i,j)]/B[num(j,j)];
		}
	}

}

/*void Calcul_de_B(int j, float *A, float *B){

	int i,k;
	float sommeCarree=0, somme=0;
	if (j==1)
	{
		for (i = 1; i <=n; i++)
			B[num(i,j)] = A[num(i,j)] / sqrt(A[num(j,j)]);
	}
	else if (j>1)
	{

		for (k = 1; k <= j-1; k++)
		{
			sommeCarree += pow(B[num(j,k)],2);
		}
		B[num(j,j)] = sqrt(A[num(j,j)] - sommeCarree);
		for (i = 1; i <=n; i++)
		{
			somme=0;
			for (k= 1; k <= j-1; k++)
			{
				somme += B[num(i,k)]*B[num(j,k)];
			}

			B[num(i,j)]= (A[num(i,j)]- somme)/B[num(j,j)];
		}
	}

}*/





	






