#include<stdio.h>
#include<conio.h>
#include<windows.h>
#include"lu.h"
#include"cholesky.h"

#define ENABLE_QUICK_EDIT_MODE  0x0040
#define ENABLE_INSERT_MODE      0x0020

#define NO_MOUSE_BUTTON         0
#define MOUSE_BUTTON_LEFT       1
#define MOUSE_BUTTON_RIGHT      2
#define t 3

#define numeroLigneOption1      5
#define numeroLigneOption2      7
#define numeroLigneOption3      9
#define numeroLigneOption4      11
#define numeroLigneOption5      13


float **matriceId;
float **matrice1;
float **NouvelleMatrice;

void ajouter(float *T)
{
    int i;
    for(i=0;i<t;i++)
    {
        printf("donner les valeurs de l'equation\n");
        scanf("%f",&T[i]);
    }
}
void AjouterMatrice()
{

    int i,j;
    //SAISIE DES ELEMENTS DE LA MATRICE
    printf("Enter les elements de la matrice\n");
    for (i=0; i<t; i++){
        for (j=0; j<t; j++){
        printf("Saisir la valeur de mat[%d][%d] : ", i, j);
        scanf("%f", &matrice1[i][j]);
        printf("\n");
        }
    }
}

void afficherMatrice()
{
    int i,j;
    printf("Affichage de la matrice \n");
    for (i=0;i<t;i++)
    {
        printf("  ");
        for (j=0;j<t;j++)
        {
            printf("%f\t  ",matrice1[i][j]);
        }
        printf(" )");
        printf("\n" );
    }
}

void afficherMatriceIdentite()
{
    int i,j;
    for (i=0;i<t;i++)
    {
        for (j=0;j<t;j++)
        {
            printf("%f  ",matriceId[i][j]);
        }
        printf("\n" );
    }
}


void afficherMatriceInverse()
{
    int i,j;
    float elem;
    for (i=0;i<t;i++)
    {
        printf("  (");
        for (j=t;j<2*t;j++)
        {
            printf("%f\t  ",NouvelleMatrice[i][j]);
        }
        printf(" ) ");
        printf("\n" );
    }
}
void creernouvellemat()
{
    int i,j;
    for (i=0;i<t;i++)
    {
        for (j=0;j<t;j++)
        {
            if (i==j)
            {
                matriceId[i][j] = 1;
            }
            else
            {
                matriceId[i][j] = 0;
            }
        }
    }
}
void definirNouvelleMatrice()
{
    int i,j;
    i=j=0;
    for (i=0;i<t;i++)
    {
        for (j=0;j<2*t;j++)
        {
            if (j<t)
            {
                NouvelleMatrice[i][j] = matrice1[i][j];
            }
            else
            {
                NouvelleMatrice[i][j] = matriceId[i][j-t];
            }
        }
    }
}
 void modifierMatrice()
{
    creernouvellemat();
    definirNouvelleMatrice();

}
int MethodeGauss()
{
    int inversible = 1;
    int k,i,colonne,colonnebis;
    float var,var1;
    k=0;
    while((inversible == 1)&&(k<t))
    {
            if (NouvelleMatrice[k][k] != 0)
            {
                var = NouvelleMatrice[k][k];
                for (colonne=0;colonne<2*t;colonne++)
                {
                    NouvelleMatrice[k][colonne] = NouvelleMatrice[k][colonne]/var;  //Normalisation de la ligne contenant l'élément diagonal
                }
                for (i=0;i<t;i++)
                {
                    if (i != k)
                    {
                        var1=NouvelleMatrice[i][k];
                        for (colonnebis=0;colonnebis<2*t;colonnebis++)
                        {
                            NouvelleMatrice[i][colonnebis] = NouvelleMatrice[i][colonnebis] - NouvelleMatrice[k][colonnebis]*var1;
                        }
                    }
                }
                k++;
            }
            else
            {
                inversible = 0;
            }
    }
    return inversible;
}
void affichage_systeme(float *b ,int n)
{
	int i , j ;
	printf(" Affichage du systeme : \n\n\n");

	for(i = 0 ; i < n ; i++)
	{
		printf("  (");
		for(j = 0 ; j < n ; j++)
		{
			printf("  %.3f  ",matrice1[i][j]);
		}
		printf(" )    (X%d)   =",i+1);
		printf("\t%.3f",b[i]);
		printf("\n\n");
	}
}


void affichage_solution(float *x, int n)
{
    int i ;
	printf("Affichage de la solution : \n\n\n");

	for(i = 0 ; i < n ; i++)
	{
        printf("(X%d)   =",i+1);
		printf("\t%.6f",x[i]);
		printf("\n\n");
	}
}


void gauss(float *b, float *x, int n)
{
 int i, j, k ;
     int iminimal ;
     float p ;
     float somme, valeur_minim, tempon1, tempon2 ;

     for(k = 0 ; k < n-1 ; k++)
     {
        /* Dans un premier temps, on cherche l'élément minimum (non */
        /* nul) en valeur absolue dans la colonne k et d'indice i   */
        /* supérieur ou égal à k.                                   */

        valeur_minim = matrice1[k][k] ; iminimal = k ;
        for(i = k+1 ; i < n ; i++)
        {
           if (valeur_minim != 0)
           {
              if (abs(matrice1[i][k]) < abs(valeur_minim) && matrice1[i][k] != 0)
              {
                 valeur_minim = matrice1[i][k] ;
                 iminimal = i ;
              }
           }
           else
           {
                 valeur_minim = matrice1[i][k] ;
                 iminimal = i ;
           }
        }

        /* Si l'élément minimum est nul, on peut en déduire */
        /* que la matrice est singulière. Le pogramme est   */
        /* alors interrompu.                                */

        if (valeur_minim == 0.)
        {
           printf("\n\n\nAttention! Matrice singuliere!\n\n\n") ;
           exit( EXIT_FAILURE ) ;
        }

        /* Si la matrice n'est pas singulière, on inverse    */
        /* les éléments de la ligne imax avec les éléments   */
        /* de la ligne k. On fait de même avec le vecteur b. */

        for(j = 0 ; j < n ; j++)
        {
           tempon1 = matrice1[iminimal][j] ;
           matrice1[iminimal][j] = matrice1[k][j] ;
           matrice1[k][j] = tempon1 ;
        }

        tempon2 = b[iminimal] ;
        b[iminimal] = b[k] ;
        b[k] = tempon2 ;


        /* On procède à la réduction de la matrice par la */
        /* méthode d'éliminimalation de Gauss. */

        for(i = k+1 ; i < n ; i++)
        {
           p = matrice1[i][k]/matrice1[k][k] ;

           for(j = 0 ; j < n ; j++)
           {
              matrice1[i][j] = matrice1[i][j] - p*matrice1[k][j] ;
           }

           b[i] = b[i] - p*b[k] ;
        }
     }



     if (matrice1[n-1][n-1] == 0)
     {
        printf("\n\n\nAttention! Matrice singuliere!\n\n\n") ;
        exit( EXIT_FAILURE ) ;
     }



     x[n-1] = b[n-1]/matrice1[n-1][n-1] ;

     for(i = n-2 ; i > -1 ; i--)
     {
           somme = 0 ;

           for(j = n-1 ; j > i ; j--)
           {
              somme = somme + matrice1[i][j]*x[j] ;
           }
           x[i] = (b[i] - somme)/matrice1[i][i] ;
     }
}
float determinant (float **mat, int n)
{
	float result;

	float **mat_n_moins_1;
	int i, j, k;
	double l, puiss=0;
    mat_n_moins_1 = (float **)malloc(t*sizeof(float *));
    for (i=0; i<t; i++)
    {
        mat_n_moins_1[i]= (float *)malloc(t*sizeof(float));
    }
	if (n==1)
		return mat [0][0];
	else
	{
		result = 0;
		for (i=0; i<n; i++)
		{
			for (j=0; j<n-1; j++)
			{
				for (k=0; k<n-1; k++)
				{
					mat_n_moins_1 [j][k] = mat [j+1][k+(k>=i)];
				}
			}
			l=i;

			result = result + pow(-1, l) *mat[0][i] * determinant(mat_n_moins_1, n-1);
		}
		return result;
	}
}


typedef struct _DATA_MOUSE{
        unsigned char buttonPressed;
        COORD coordinates;
}DMOUSE, *PDMOUSE;



void mouse(PDMOUSE m){
     HANDLE in = GetStdHandle(STD_INPUT_HANDLE);
     DWORD mode;
     DWORD read;
     INPUT_RECORD go[1];
     int nbre = 1,col = 16;

     BOOL buttonMouseValidate = FALSE;
     BOOL keyValidate = FALSE;

     GetConsoleMode(in, &mode);
     SetConsoleMode(in, (mode | ENABLE_MOUSE_INPUT) & ~ENABLE_QUICK_EDIT_MODE & ~ENABLE_INSERT_MODE);

     while(!buttonMouseValidate){
        ReadConsoleInput(in, go, 1, &read);

        if(go[0].EventType == MOUSE_EVENT){
               if(buttonMouseValidate = go[0].Event.MouseEvent.dwButtonState == FROM_LEFT_1ST_BUTTON_PRESSED)
                    m->buttonPressed = MOUSE_BUTTON_LEFT;
               else if(buttonMouseValidate = go[0].Event.MouseEvent.dwButtonState == RIGHTMOST_BUTTON_PRESSED)
                    m->buttonPressed = MOUSE_BUTTON_RIGHT;

               if(buttonMouseValidate){
                  m->coordinates.X = go[0].Event.MouseEvent.dwMousePosition.X;
                  m->coordinates.Y = go[0].Event.MouseEvent.dwMousePosition.Y;
               }
        }else if(go[0].EventType == KEY_EVENT){
            if(go[0].Event.KeyEvent.wVirtualKeyCode == VK_DOWN){
                    nbre++;
                    switch(nbre){
                        case 1:
                            gotoxy(numeroLigneOption1,col);
                            break;
                        case 2:
                            gotoxy(numeroLigneOption2,col);
                            break;
                        case 3:
                            gotoxy(numeroLigneOption3,col);
                            break;
                        case 4:
                            gotoxy(numeroLigneOption4,col);
                            break;
                        case 5:
                            gotoxy(numeroLigneOption5,col);
                            break;
                        default:nbre--;break;
                    }
            }else if(go[0].Event.KeyEvent.wVirtualKeyCode == VK_UP){
                nbre--;
                    switch(nbre){
                        case 1:
                            gotoxy(numeroLigneOption1,col);
                            break;
                        case 2:
                             gotoxy(numeroLigneOption2,col);
                             break;
                        case 3:
                             gotoxy(numeroLigneOption3,col);
                             break;
                        case 4:
                             gotoxy(numeroLigneOption4,col);
                             break;
                        case 5:
                             gotoxy(numeroLigneOption5,col);
                             break;
                        default:nbre++;break;
                    }
            }
        }
     }

     SetConsoleMode(in, mode);
}
void color(int couleurDuTexte,int couleurDeFond){  // fonction d'affichage de couleurs
        HANDLE H=GetStdHandle(STD_OUTPUT_HANDLE);
        SetConsoleTextAttribute(H,couleurDeFond*15+couleurDuTexte);
}
void gotoxy(int y, int x){
  COORD coord;
  coord.X = x;
  coord.Y = y;
  SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), coord);
}
void menu(){
    printf("\n\n\n");
    color(0,15);
    printf("\t\t\t\t ********************************** MENU PRINCIPALE ************************************ \t\t\t\t\n");
    printf("\t\t\t\t\t\t\t\t\n");
    printf("\t\t   Systeme lineaire             \n");
    printf("\t\t                                                \n");
    printf("\t\t  Inverse dune matrice                          \n");
    printf("\t\t                                                \n");
    printf("\t\t  Determination du determinant                  \n");
    printf("\t\t                                                \n");
    printf("\t\t  La methode LU                          \n");
    printf("\t\t                                                \n");
    printf("\t\t  La methode  cholesky                          \n");
}
void modePleinEcran(int Mode){ // parametre Mode : 1=plein ecran et 2=mode fenetre
    typedef BOOL WINAPI (*SetConsoleDisplayModeT)(HANDLE,DWORD,DWORD*);
    SetConsoleDisplayModeT SetConsoleDisplayMode;

   HINSTANCE hLib=LoadLibrary("KERNEL32.DLL");
   SetConsoleDisplayMode=(SetConsoleDisplayModeT)
      GetProcAddress(hLib,"SetConsoleDisplayMode");

   HANDLE h=CreateFile("CONOUT$",GENERIC_WRITE|GENERIC_READ,FILE_SHARE_READ |
      FILE_SHARE_WRITE,NULL,OPEN_EXISTING,0,0);

   DWORD oldmode;

   SetConsoleDisplayMode(h,Mode,&oldmode);
}
int main(){
    DMOUSE m;
    int i;
     //pour activer le mode plein écran
     float * b ;
    float **tmp;
    float * x ;
    float **tab;
    matriceId = (float **)malloc(t*sizeof(float *));
    for (i=0; i<t; i++)
    {
        matriceId[i]= (float *)malloc(t*sizeof(float));
    }
    tab = (float **)malloc(t*sizeof(float *));
    for (i=0; i<t; i++)
    {
        tab[i]= (float *)malloc(t*sizeof(float));
    }

    tmp = (float **)malloc(t*sizeof(float *));
    for (i=0; i<t; i++)
    {
        tmp[i]= (float *)malloc(t*sizeof(float));
    }

    matrice1 = (float **)malloc(t*sizeof(float *));
    for (i=0; i<t; i++)
    {
        matrice1[i]= (float *)malloc(t*sizeof(float));
    }
    NouvelleMatrice = (float **)malloc(t*sizeof(float *));
    for (i=0; i<t; i++)
    {
        NouvelleMatrice[i]= (float *)malloc((2*t)*sizeof(float));
    }
     b = (float *) malloc (sizeof (float *) * t) ;
    x = (float *) malloc (sizeof (float *) * t) ;

    menu();

    while(1){
        mouse(&m);
        color(10,3);
        if(m.buttonPressed==1 && (m.coordinates.X >=0 && m.coordinates.X<=100) && m.coordinates.Y == numeroLigneOption1){
            system("cls");
            AjouterMatrice();
            afficherMatrice();
            ajouter(b);
            modifierMatrice();

            if (MethodeGauss() == 1)
            {
                printf("\n\n******* Affichage de la Matrice inverse\n" );

            }
            else
            {
                printf("La matrice n'est pas inversible\n" );
            }
            tmp = matrice1;

            affichage_systeme(b,t) ;
            gauss(b,x,t) ;


            affichage_systeme(b,t) ;


            affichage_solution(x,t) ;
            printf("Appuyer sur entrer pour retourner au menu");
            getch();
            system("cls");
            main();
        }
        if(m.buttonPressed==1 && (m.coordinates.X >=0 && m.coordinates.X<=100) && m.coordinates.Y == numeroLigneOption2){
            printf("Vous avez choisi le menu 2 avec la souris\n");
            AjouterMatrice();
            system("cls");
            afficherMatrice();
            modifierMatrice();
            if (MethodeGauss() == 1)
            {
                printf("Matrice inverse\n" );
                afficherMatriceInverse();
            }
            else
            {
                printf("La matrice n'est pas inversible\n" );
            }
            printf("Appuyer sur entrer pour retourner au menu");
            getch();
            system("cls");
            main();
        }
        if(m.buttonPressed==1 && (m.coordinates.X >=0 && m.coordinates.X<=100) && m.coordinates.Y == numeroLigneOption3){
            system("cls");
            AjouterMatrice();
            tab = matrice1;
            afficherMatrice();
            printf ("determinant : %f\n", determinant (tab, 3));
            printf("Appuyer sur entrer pour retourner au menu");
            getch();
            system("cls");
            main();
        }

        if(m.buttonPressed==1 && (m.coordinates.X >=0 && m.coordinates.X<=100) && m.coordinates.Y == numeroLigneOption4){
           /* system("cls");
            exit(0);*/
            printf("Les tailles des matrices et vecteurs  : ");
	       scanf("%d", &n) ;
	       printf("\n") ;
	A = (float *) malloc (n*n* sizeof (float));
	L = (float *) malloc (n*n* sizeof (float));
	b = (float *) malloc ((n+1)* sizeof (float));
	x = (float *) malloc ((n+1)* sizeof (float));
	y = (float *) malloc ((n+1)* sizeof (float));
  	//Saisie de la matrice A et du vecteur b
	saisie(n,b,A);

  	// Affichage de la matrice
	affiche_matrice(n, A);
	// Affichage du vecteur
	affiche_vecteur(n,b);
	// application de la méthode de gauss

	while((A[numero(k,k)]!=0) && (k<=n))
	{

		remplir_L(k,A);
		elimination(k,A);
		affiche_matrice(n,A);
		k++;
	}

	descente(n,y,b,L);
	affiche_matrice(n,L);
	remonte(n,x,y,A);
	affiche_matrice(n,A);
	affiche_solution(n,x);

	free(A);
	free(b);
	free(L);
	free(y);
	free(x);

        }
        if(m.buttonPressed==1 && (m.coordinates.X >=0 && m.coordinates.X<=100) && m.coordinates.Y == numeroLigneOption5){
            //printf("ok");
            // Saisie de la taille du système d'équation
	printf("\nEntrez les dimension du systeme   : ");
	scanf("%d", &n) ;
	printf("\n") ;
	A = (float *) malloc (n*n* sizeof (float));
	B = (float *) malloc (n*n* sizeof (float));
	b = (float *) malloc ((n+1)* sizeof (float));
	x = (float *) malloc ((n+1)* sizeof (float));
	y = (float *) malloc ((n+1)* sizeof (float));
  	//Saisie de la matrice A et du vecteur b
	saisie(n,b,A);

  	// Affichage de la matrice
  	printf("\nAffichage de la Matrice: A \n");
	afficher_matrice(n, A);
	// Affichage du vecteur
	afficher_vecteur(n,b);
	// application de la méthode de cholesky

	facto_cholesky(A,B);
	descente_cholesky(n,y,b,B);
	remonte_cholesky(n,x,y,B);
	affiche_sol(n,x);

	          free(A);
	          free(b);
	          free(B);
	          free(y);
	          free(x);
          }

    }

}
