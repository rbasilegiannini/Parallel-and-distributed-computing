#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "mpi.h"

//Roberto Basile Giannini
//N97000340

/**************PROTOTIPI FUNZIONI**************/

void 	crea_griglia(	MPI_Comm *griglia, MPI_Comm *grigliar, MPI_Comm *grigliac, 
						int menum, int rig, int col, int *coordinate);	

void distribuzione_mat(	int menum, int nproc, int* mat, int* blocco, int mat_rig, int mat_col, 
					  	int blocco_col, int blocco_rig, int p, int q, MPI_Comm comm);

/**************FINE PROTOTIPI***************/

/****************MAIN*****************/
int main (int argc, char *argv[]) 
{
//0. Init
	
	//variabili
	int i, j, k, h, nproc, menum;
	int *mat_A, *mat_B, *blocco_A, *blocco_B,	//matrici e sottoblocchi
		*blocco_A_bcast, *blocco_C;
	int mat_col_A, mat_rig_A;					//dimensioni di A
	int mat_col_B, mat_rig_B;					//dimensioni di B
	int m_input;								//dimensione in input
	int root_bcast;								//sorgente Bcast in BMR
	int tag_send, tag_recv;						//tag comunicazione
	int dest, mitt;								//destinatario e mittente comunicazione
	int comm_grid_c_id;							//id nel communicator delle colonne della griglia

	double t0, t1, tdiff;						//variabili utilizzate per la misura dei tempi
	double timetot;
	
	MPI_Request request;
	MPI_Status status;
			
	int blocco_rig_A, blocco_col_A;			//numero di righe e colonne del blocco A
	int blocco_rig_B, blocco_col_B;			//numero di righe e colonne del blocco B
	int blocco_rig_C, blocco_col_C;			//numero di righe e colonne del blocco C

	int p;									//numero di righe della griglia
	int q;									//numero di colonne della griglia
	int *coordinate;						//coordinate della griglia
	MPI_Comm comm_grid, comm_grid_r, 		//communicator della griglia
		comm_grid_c;
	//
	
	//inizializzazione MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &menum);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	//
//

//1. Lettura input ed inizializzazioni
	
	if (menum == 0) {
		p = atoi(argv[1]);
		
		//controllo validità di p
		if (p < 1) {
			printf ("[E0] Dimensione griglia %d non valida. Valore default: 2 \n", p);
			p = 2;
		}
	}
	
	MPI_Bcast (&p, 1, MPI_INT, 0, MPI_COMM_WORLD); 			//comunicazione in BCast di p
	q = nproc/p;
	
	//controllo griglia quadrata
	if (p != q) {
		if (menum == 0)
			printf("[E1] p e q sono diversi: griglia non quadrata. Terminazione... \n");
		
		MPI_Finalize();
		return -1;
	} 
	
	else {
		if (menum == 0)
			printf("[P0] dimensioni griglia: %dx%d. \n", p, q);
		
	}
    
	//Lettura dimensione matrici: basta solo una dimensione	
	if (menum == 0) {
		m_input = atoi(argv[2]);
		
		//controllo validità di m_input
		if (m_input < 1) {
			printf ("[E2] Dimensione input %d non valido. Valore default: 16. \n", m_input);
			m_input = 16;
		}

		printf("[P0] Dimensione in input: %d.\n", m_input);
		
		//controllo m_input multiplo di p
		if ( (m_input%p) != 0) {
			printf ("[E3] Dimensione input %d non e' multiplo di p = %d. Valore di default: p*10. \n", m_input, p);
			m_input = p*10;
		}

		//costruzione matrici A e B	
		
		//dimensioni matrici quadrate
		mat_rig_A = m_input;
		mat_col_A = m_input;
		mat_rig_B = m_input;
		mat_col_B = m_input;
		
		//allocazione matrici A e B
		mat_A = (int*) malloc (mat_rig_A * mat_col_A * sizeof(int));
		mat_B = (int*) malloc (mat_rig_B * mat_col_B * sizeof(int));

		srand(time(NULL));		
		//inizializzazione matrice A
		for (i=0; i<mat_rig_A*mat_col_A; i++) {
			mat_A[i] = 1+rand()%9;
		}
		
		//inizializzazione matrice B
		for (i=0; i<mat_rig_B*mat_col_B; i++) {
			mat_B[i] = 1+rand()%9;
		}
		
		//stampa matrice A
		printf ("[P0] Matrice A (%dx%d):\n", mat_rig_A, mat_col_A);
		
		for (i=0; i<mat_rig_A; i++){
			for (j=0; j<mat_col_A; j++) {
				printf("%d ", mat_A[i*mat_col_A+j]);
			}	
			printf("\n");
		}
		
		//stampa matrice B
		printf ("[P0] Matrice B (%dx%d):\n", mat_rig_B, mat_col_B);
		
		for (i=0; i<mat_rig_B; i++){
			for (j=0; j<mat_col_B; j++) {
				printf("%d ", mat_B[i*mat_col_B+j]);
			}	
			printf("\n");
		}
		
		printf ("\n");
			
		//dimensioni blocchi quadrati
		blocco_rig_A = mat_rig_A / p;
		blocco_col_A = mat_col_A / q;
		blocco_rig_B = mat_rig_B / p;
		blocco_col_B = mat_col_B / q;
	}
	
	MPI_Bcast (&blocco_rig_A, 1, MPI_INT, 0, MPI_COMM_WORLD); 			//comunicazione in BCast di blocco_rig_A
	MPI_Bcast (&blocco_col_A, 1, MPI_INT, 0, MPI_COMM_WORLD); 			//comunicazione in BCast di blocco_col_A
	MPI_Bcast (&blocco_rig_B, 1, MPI_INT, 0, MPI_COMM_WORLD); 			//comunicazione in BCast di blocco_rig_B
	MPI_Bcast (&blocco_col_B, 1, MPI_INT, 0, MPI_COMM_WORLD); 			//comunicazione in BCast di blocco_col_B
	
	//dimensioni blocco_C
	blocco_rig_C = blocco_rig_A;
	blocco_col_C = blocco_col_B;
	
	//allocazione blocchi
	blocco_A 		= (int*) malloc (blocco_rig_A * blocco_col_A * sizeof(int));
	blocco_A_bcast 	= (int*) malloc (blocco_rig_A * blocco_col_A * sizeof(int));
	blocco_B 		= (int*) malloc (blocco_rig_B * blocco_col_B * sizeof(int));
	blocco_C 		= (int*) calloc (blocco_rig_C * blocco_col_C , sizeof(int));

//

//2. Creazione griglia	

	//Init coordinate
	coordinate = (int*) calloc (2, sizeof(int)); 			//coordinate bidimensionali, ogni proc ha le sue coordinate
	MPI_Barrier(MPI_COMM_WORLD);
	
	//routine di creazione griglia...
	crea_griglia(&comm_grid, &comm_grid_r, &comm_grid_c, menum, p, q, coordinate);
	MPI_Comm_rank(comm_grid_c, &comm_grid_c_id);
	
	MPI_Barrier(comm_grid);	
//
	
//3. Distribuzione matrici

	//distribuzione matrice A
	distribuzione_mat(	menum, nproc, mat_A, blocco_A, mat_rig_A, mat_col_A, blocco_col_A, 
						blocco_rig_A, p, q, comm_grid);

	MPI_Barrier (comm_grid);
	
	//distribuzione matrice B
	distribuzione_mat(	menum, nproc, mat_B, blocco_B, mat_rig_B, mat_col_B, blocco_col_B, 
						blocco_rig_B, p, q, comm_grid);

	MPI_Barrier (comm_grid);
	//A questo punto, i processori P_ij hanno a propria disposizione i blocchi A_ij e B_ij
	
	//stampe dei blocchi A e B
    for (k = 0; k<nproc; k++) {
		
		if (menum == k) {
		
			printf("[P%d] Il mio primo blocco A:\n", menum);
			for (i=0; i<blocco_rig_A; i++) {
		        
				for (j=0; j<blocco_col_A; j++) {
		    	    printf("%d ",blocco_A[i*blocco_col_A+j]);
		        }
		    
			    printf("\n");
		    }
		   
			printf("[P%d] Il mio primo blocco B:\n", menum);
			for (i=0; i<blocco_rig_B; i++) {
		        
				for (j=0; j<blocco_col_B; j++) {
		    	    printf("%d ",blocco_B[i*blocco_col_B+j]);
		        }
		    
			    printf("\n");
		    }
			
			printf("\n");
		}
					
		MPI_Barrier(comm_grid);
	}
//

	MPI_Barrier(comm_grid);
	t0 = MPI_Wtime();				//Viene segnato t0 da tutti i processori

//4.BMR

	//bisogna ripetere BMR p volte, dove p è il numero di diagonali
	for (k=0; k<p; k++) {
		
	//Broadcast			
			
		//se si tratta di un processore sulla k-esima diagonale
		if (coordinate[1] == (coordinate[0] + k)%p) {
					
			//predisponi il blocco da spedire
			for (j=0; j<blocco_rig_A*blocco_col_A;j++)
				blocco_A_bcast[j] = blocco_A[j];
			
		}
		
		root_bcast = (coordinate[0] + k)%p;			//Sorgente del Bcast: è l'ID del comm di riga relativo all'elemento diagonale della rispettiva riga

		//Bcast solo su riga di root_bcast
		MPI_Bcast(blocco_A_bcast, blocco_rig_A*blocco_col_A, MPI_INT, root_bcast, comm_grid_r);
	//

	//Multiply
		//prodotto blocco_A_bcast * blocco_B		
		for (i=0; i<blocco_rig_C; i++){
			for (j=0; j<blocco_col_C; j++) {
				for (h=0; h<blocco_rig_B; h++) {
				
					blocco_C[i*blocco_col_C+j] += blocco_A_bcast[i*blocco_col_A+h] * blocco_B[h*blocco_col_B+j];
				
				}
			}
		}
	//
	
		MPI_Barrier(comm_grid);
				
	//Rolling				
				
		dest = (comm_grid_c_id - 1 + p)%p;
		mitt = (comm_grid_c_id + 1)%p;
		tag_send = 50+dest;
		tag_recv = 50+comm_grid_c_id;
				
		MPI_Isend(	blocco_B, blocco_rig_B*blocco_col_B, MPI_INT, dest, tag_send,
					comm_grid_c, &request);
					
		MPI_Barrier(comm_grid);

		MPI_Recv(	blocco_B, blocco_rig_B*blocco_col_B, MPI_INT, mitt, tag_recv,
					comm_grid_c, &status);
	//
		MPI_Barrier(comm_grid);
	}
//

	t1 = MPI_Wtime();	//viene segnato t1 da tutti i processori
	tdiff = t1 - t0;		//tempo impiegato da ogni processo

	//P0 ora conosce il max tra i tempi impiegati	
	MPI_Reduce(&tdiff, &timetot, 1, MPI_DOUBLE, MPI_MAX, 0, comm_grid);

//5. Stampa risultato e tempi
	for (h = 0; h<nproc; h++) {
			
		if (menum == h) {
			printf("[P%d] Il mio blocco risultato C:\n", menum);
			for (i=0; i<blocco_rig_C; i++) {
		        
				for (j=0; j<blocco_col_C; j++) {
		    	    printf("%d ",blocco_C[i*blocco_col_C+j]);
		        }
		    
			    printf("\n");
		    }
			    
		}
						
		MPI_Barrier(comm_grid);
		
	}

	//segnalazione del timetot max
	if (menum == 0)
		printf("[TIMEMAX][P0] Tempo massimo impiegato: %e s\n", timetot);

//6. Finalizzazione
	MPI_Finalize();

	free(blocco_A);
	free(blocco_A_bcast);
	free(blocco_B);
	free(blocco_C);
	free(coordinate);
	if (menum == 0) {
		free(mat_A);
		free(mat_B);
	}

	return 0;
}
/*******************FINE MAIN************************/

/******************IMPLEMENTAZIONE FUNZIONI******************/

void crea_griglia(MPI_Comm *griglia, MPI_Comm *grigliar, MPI_Comm *grigliac, 
				  int menum, int grid_rig, int grid_col, int *coordinate) {

//Output : griglia, grigliar, grigliac, coordinate

	int dim, reorder;
	int *ndim, *period, *vc;
	
	dim		= 2; //griglia bidimensionale
	reorder = 0; //no riordinamento degli identificativi
	
	ndim	= (int*) calloc (dim, sizeof(int));
	period	= (int*) calloc (dim, sizeof(int));
	vc		= (int*) calloc (dim, sizeof(int));
	
	ndim[0] = grid_rig; //numero di righe
	ndim[1] = grid_col; //numero di colonne
	
	period[0] = 1; //periodicità righe
	period[1] = 1; //periodicità colonne
	
	MPI_Cart_create (MPI_COMM_WORLD, dim, ndim, period, reorder, griglia);	//creazione della griglia
	MPI_Cart_coords (*griglia, menum, dim, coordinate);						// assegnazione coordinate

	vc[0] = 0; //per cancellare il collegamento tra le righe
	vc[1] = 1; //per mantenere il collegamento tra le colonne
	
	MPI_Cart_sub (*griglia, vc, grigliar); //divisione in righe del comm
	
	vc[0] = 1; //per mantenere il collegamento tra le righe
	vc[1] = 0; //per cancellare il collegamento tra le colonne
	
	MPI_Cart_sub (*griglia, vc, grigliac); //divisione in colonne del comm
	
	//deallocazioni
	free(ndim);
	free(vc);
	free(period);
	
	return;
}


void distribuzione_mat(	int menum, int nproc, int* mat, int* blocco, int mat_rig, int mat_col, 
					  	int blocco_col, int blocco_rig, int p, int q, MPI_Comm comm) {

//Output : blocco
		  		
    MPI_Datatype type;
    MPI_Datatype type2;	
	
	int *disps, *counts;			
	int i, j;
	
	//allocazioni vettori
	disps 	= 	(int*) malloc (nproc*sizeof(int)); 					//displacement
	counts 	= 	(int*) malloc (nproc*sizeof(int));					//sizecounts
		
	//costruzione del blocco (ho provato MPI_Type_vector)
    MPI_Type_vector(blocco_rig, blocco_col, mat_col, MPI_INT, &type2);	//blocco_rig righe di blocco_col elementi con stride mat_col
   	MPI_Type_create_resized( type2, 0, sizeof(int), &type);				//reset del lower bound e dell'upper bound
   	MPI_Type_commit(&type);												//commit del datatype
	
	//inizializzazione dei displacements e del sendcounts. nproc = p * q
   	for (i=0; i<p; i++) {
       	for (j=0; j<q; j++) {
       		disps[i*q+j] = i*mat_col*blocco_rig+j*blocco_col;
       		counts [i*q+j] = 1;
        }
    }
	
	//ad ogni processore viene spedito un blocco tramite scatterv
    MPI_Scatterv(mat, counts, disps, type, blocco, blocco_rig*blocco_col, MPI_INT, 0, comm);

	if (menum == 0) {
		free(disps);
		free(counts);
	}

	return;	  		
}

/*************************************************************/


