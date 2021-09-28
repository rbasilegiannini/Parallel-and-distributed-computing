#include <stdlib.h>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

/**************PROTOTIPI************/
int* matxvet(int mat_rig, int mat_col, int* x, int* A, double* tdiff);
/***********************************/

int main(int argc, char *argv[]) {
	
//0. Init
	int mat_col, mat_rig; 	//dimensioni della matrice
	
	int* A;				//matrice
	int* x;				//vettore
	int* res;				//vettore risultato
	
	int i,j;

	double tdiff;			//differenza tempo finale - tempo iniziale

	//lettura di mat_rig e mat_col
	mat_rig = atoi(argv[1]);
	mat_col = atoi(argv[2]);

//1. Controllo input
	
	//controllo validità mat_rig
	if (mat_rig < 1) {
		printf ("[E0] Numero %d di righe non valido. Valore default: 20.\n", mat_rig);
		mat_rig = 20;
	}
	
	//controllo validità mat_col
	if (mat_col < 1) {
		printf ("[E1] Numero %d di colonne non valido. Valore default: 20.\n", mat_col);
		mat_col = 20;
	}
		
	//Allocazioni
	A = 	(int*) malloc (mat_rig*mat_col*sizeof(int));		//Allocazione della matrice
	x = 	(int*) malloc (mat_col*sizeof(int));			//Allocazione del vettore

//2. Inizializzazione di A e x
	printf("Generazione matrice %dx%d...\n", mat_rig, mat_col);
	
	srand(time(NULL));
	for (i=0; i<mat_rig*mat_col; i++) {
		A[i] = 1+rand()%9;
	}
	
	for (i=0; i<mat_rig; i++){
		for (j=0; j<mat_col; j++) {
			printf("%d ", A[i*mat_col+j]);
		}
		printf("\n");
	}
	
	printf("Generazione vettore...\n");
	
	for (i=0; i<mat_col; i++) {
		x[i] = 1;
	}
	
	for (i=0; i<mat_col; i++) {
		printf("%d \n", x[i]);
	}
	
//3. Calcolo Ax = b
	printf("Calcolo Ax = b...\n");
	res = matxvet(mat_rig, mat_col, x, A, &tdiff);

//4. Stampa risultato
	printf("Risultato:\n");
	for (i=0; i<mat_rig; i++) {
		printf("%d \n", res[i]);
	}

//5. Stampa dei tempi
	printf("Tempo di esecuzione con %e s threads: %ld \n", omp_get_max_threads(), tdiff);
 
//6. Finalizzazione
	//deallocazioni
	free(A);
	free(x);
	free(res);
	
	return 0;
}

/**********IMPLEMENTAZIONI*************/

int* matxvet(int mat_rig, int mat_col, int* x, int* A, double* tdiff) {

	int i,j;
	int *res;

	//variabili per le misure
	struct timeval time;		//tempo 
	double ti;					//tempo iniziale
	double tf;					//tempo finale
	
	res = (int*) calloc (mat_rig, sizeof(int));			//Inizializzazione del vettore risultato a 0
	
	//inizio regione parallela
	printf("Inizio regione parallela...\n");	  
	
	gettimeofday(&time, NULL);
	ti = time.tv_sec + (time.tv_usec/1000000.0);
	#pragma omp parallel for default(none) shared(mat_rig, mat_col, A, x, res) private(i,j)
	
		for (i=0; i<mat_rig; i++) {
			for (j=0; j<mat_col; j++)
				res[i] += A[i*mat_col+j] * x[j];
		}
	gettimeofday(&time, NULL);
	tf = time.tv_sec + (time.tv_usec/1000000.0);

	//fine regione parallela
	printf("Fine regione parallela.\n");

	*tdiff = (tf - ti);	

	return res;

}

/********************************/

