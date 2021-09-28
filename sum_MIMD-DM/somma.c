#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h> 
#include "mpi.h"

//Roberto Basile Giannini
//N97000340

/**************PROTOTIPI FUNZIONI**************/

//distribuzione dati
void distribuzione_dati(int me, float *data, float *dataloc, int nl, int np, int res, MPI_Status *stat);	

//strategie
void I_strategia(int me, int ID_s, float *sum_p, float *sum, int np, MPI_Status *stat);
void II_strategia(int me, float *sum_p, float *sum_tot, int idx, int *p_2, int *p_2_plus_1, MPI_Status *stat);
void III_strategia(int me, float *sum_p, float *sum_tot, int idx, int *p_2, int *p_2_plus_1, MPI_Status *stat);

//stampa risultato
void stampa_ris(int me, float sum_tot, int ID_s, int N_str, MPI_Status *stat);

/**************FINE PROTOTIPI***************/

/****************MAIN*****************/
int main (int argc, char *argv[]) 
{
//0. Init
	//variabili
	int menum, nproc;
	int n, nloc, tag, 
		rest, i, j, 
		start, N_strat_temp, 
		N_strat, ID_stamp, ID_stamp_time,
		index;
	
	int *pow_2, *pow_2_plus_1;
	
	float sum, sum_parz;
	float *x, *xloc;
	
	double t0, t1, time;
	double timetot;

	MPI_Status status;
	//
	
	//inizializzazione MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &menum);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	//
//

//1. lettura e distribuzione dei dati
	if (menum == 0) {	//P0 demandato al compito di lettura dei dati
		
		//lettura e controllo di N	
		n = atoi(argv[1]);
		
		//ERR: addendi non sufficienti
		if ((n > argc-4) && (n <= 20)) {//Numero addenti non sufficienti
			printf ("[ERROR E1] Addendi non sufficienti. Devono essere almeno %d. Verranno considerati solo quelli presenti.\n", n);
			n = argc-4;
			//In caso il numero di addendi sia superiore ad N, considero solo i primi N
		}
		
		//ERR: N non valido
		if (n < 1) {//Se N non è valido
			printf ("[ERROR E2] Numero %d di addendi non valido. Valore default: 21\n", n);
			n = 21;
		}
		
		//ERR: N non sufficiente
		if ((n < argc-4) && (n <= 20)) {
			printf ("[ERROR E3] N non sufficientemente grande. Verranno considerati solo i primi %d elementi.\n", n);
		}
		
		//allocazione di x con n
		x = (float*) malloc(n*sizeof(float));
		
		//lettura input
		j=0;
		for (i=4; i<n+4; i++) {
			if (n <= 20) {//N <= 20 e Numero addendi sufficiente

				x[j] = atof(argv[i]);  //prelievo input da linea di comando
				j++;
				
			}
			else {//N > 20
				
				x[j] = 1+rand()%100;  //assegnazione randomica
				j++;
				
			}
		}
	}
	
	//Broadcast tra i processi: P0 informa i restanti processi quanti addendi ci sono
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	nloc = n/nproc;
	rest = n%nproc;
	
	//gestione del resto
	if (menum < rest)
		nloc = nloc+1; //necessario per una corretta allocazione di xloc qualora ci fosse del resto
	
	if(menum != 0) { //non necessario per P0
		xloc = (float*) malloc(nloc*sizeof(float));
	}	
	
	//distribuzione dei dati
	distribuzione_dati(menum, x, xloc, nloc, nproc, rest, &status);	
	
//

//2. Controllo degli input di controllo
	ID_stamp = atoi(argv[3]); //prelievo ID del processore di stampa
	
	//ERR: controllo ID_stamp
	if (ID_stamp < -1 || ID_stamp >= nproc) {//In caso di ID_stamp non valido
		if (menum == 0) {//P0 demandato alla segnalazione di errori
			printf ("[ERROR E4] ID_stamp %d non valido, deve essere compreso tra [-1, %d]. Valore default: 0\n", ID_stamp, nproc-1);
		}	
		ID_stamp = 0;
	}
	
	N_strat = 1; //valore default
	N_strat_temp = atoi(argv[2]);
	
	//ERR: numero di strategia non valido
	if ((N_strat_temp < 1) || (N_strat_temp > 3)) {
		if (menum == 0) {//P0 demandato alla segnalazione di errori
			printf("[ERROR E5] Strategia %d non valida, possibili solo: 1, 2 o 3. Valore default: 1\n", N_strat_temp);
		}
	}
	
	else if (N_strat_temp != 1){//Numero strategia valido
	
			//Verifica applicabilità strategia II o III
			if((n != 0) && ((n &(n - 1)) == 0)){ //Se la strategia scelta è la II o III, verifica se N è potenza di 2.
					N_strat = N_strat_temp;		
			}
			
			else {//ERR: non applicabilità strategia 2 o 3
				if(menum == 0) //P0 demandato alla segnalazione di errori
					printf("[ERROR E6] Strategia %d applicabile solo se N e' potenza di 2. Strategia di default: 1\n", N_strat_temp);
			}	
			
		}	//Altrimenti N_strat = N_strat_temp che è uguale a 1, ma 1 è già il valore di default di N_strat  
	

	
	//ERR: se la strategia non è la III e viene richiesta la stampa di tutti
	if ((N_strat != 3) && (ID_stamp == -1)) { 
		if(menum == 0) {//P0 demandato alla segnalazione degli errori
			printf("[ERROR E7] Strategia %d incompatibile con la stampa di tutti. Strategia imposta: 3\n", N_strat);
		}
		N_strat = 3;
	}	
//

//3. Calcolo variabili di supporto

	if (N_strat == 2 || N_strat == 3) {
	
		index = log(nproc)/log(2); //per le strategie II e III
		pow_2 = (int*) malloc(index*sizeof(int));
		pow_2_plus_1 = (int*) malloc(index*sizeof(int));
		
		for (i=0; i<index; i++) {
			pow_2[i] = (int)pow(2,i);
			pow_2_plus_1[i] = (int)pow(2,i+1);
		}
	
	}
	
	MPI_Barrier (MPI_COMM_WORLD); //Tutti i processi del communicator specificato si fermano all'istruzione Barrier per poi riprendere tutti insieme
	t0 = MPI_Wtime(); //Viene segnato t0 da tutti i processi

//4. Calcolo somma locale
	sum = 0;
	
	if (menum == 0) {//P0
		for (i=0; i<nloc; i++) {
			sum += x[i];
		}
	}
	else {//P1,...,Pm-1
		for (i=0; i<nloc; i++) {
			sum += xloc[i];
		}
	}
//


//5. Calcolo somma totale
	switch (N_strat) {
		
		case 1: 
			I_strategia(menum, ID_stamp, &sum_parz, &sum, nproc, &status);	
			break;
			
		case 2:
			II_strategia(menum, &sum_parz, &sum, index, pow_2, pow_2_plus_1, &status);	
			break;
			
		case 3:
			III_strategia(menum, &sum_parz, &sum, index, pow_2, pow_2_plus_1, &status);	
			break;
		
		//ERR: Numero strategia non valido già gestito
		default: 
			break;
	}
	
	t1 = MPI_Wtime(); //viene segnato t1 da tutti i processi
	time = t1 - t0; //tempo impiegato da ogni processo
	
	//Comunicazione collettiva tra i processi
	if (ID_stamp == -1) //nel caso si voglia che tutti i processi stampino i risultati
		ID_stamp_time = 0; //solo P0 stamperà il timetot massimo
	else
		ID_stamp_time = ID_stamp;
	
	//PID_stamp_time ora conosce il max tra i tempi impiegati	
	MPI_Reduce(&time, &timetot, 1, MPI_DOUBLE, MPI_MAX, ID_stamp_time, MPI_COMM_WORLD);
		
//

//6. Stampa risultato
	stampa_ris(menum, sum, ID_stamp, N_strat, &status);
//
	
//7. Stampa tempi dopo i risultati
	
	//segnalazione del tempo di ogni processo
	printf("[TIME][P%d] Tempo impiegato: %e s.\n", menum, time);
		
	MPI_Barrier (MPI_COMM_WORLD);
	//segnalazione del timetot max
	if (menum == ID_stamp_time) 
		printf("[TIMEMAX][P%d] Tempo massimo impiegato: %e s\n", menum, timetot);


//8. Finalizzazione
	if (menum == 0)
		free(x);
	else 
		free(xloc);
		
	if (N_strat == 2 || N_strat == 3){
	
		free (pow_2);
		free (pow_2_plus_1);
	
	}
	
	MPI_Finalize();
//

	return 0;
}
/*******************FINE MAIN************************/

/******************IMPLEMENTAZIONE FUNZIONI******************/


void distribuzione_dati(int me, float *data, float *dataloc, int nl, int np, int res, MPI_Status *stat) {	
	//distribuzione dei dati
		
	int menum;
	float *x;
	float *xloc;
	int nloc;
	int nproc;
	int i;
	int tag;
	int tmp;
	int start;
	int rest;
	MPI_Status *status;
	
	menum = me;
	x = data;
	xloc = dataloc;
	nloc = nl;
	nproc = np;
	rest = res;
	status = stat;
	
	if (menum == 0) { //P0 demandato al compito di distribuzione dei dati
		
		tmp = nloc;	//gestione del resto
		start = 0;
		
		for (i=1; i<nproc; i++) { //send ai processi, da P0
			start = start+tmp;
			tag = 22+i;
			if (i == rest)
				tmp = tmp-1; //gestione del resto
			MPI_Send(&x[start], tmp, MPI_FLOAT, i, tag, MPI_COMM_WORLD);		
		}
	}
	
	else {
		//recv ai processi, da P0
		tag = 22+menum;
		MPI_Recv(xloc, nloc, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, status);	
		
	}
}


void I_strategia(int me, int ID_s, float *sum_p, float *sum_tot, int np, MPI_Status *stat) {
	//I strategia

	int i;
	int sender;
	int tag;
	int menum;
	int ID_stamp;
	float *sum_parz;
	float *sum;
	int nproc;
	MPI_Status *status;
	
	menum = me;
	ID_stamp = ID_s;
	sum_parz = sum_p;
	sum = sum_tot;
	nproc = np;
	status = stat;	
	
	if (menum == ID_stamp) {//P_IDs riceve da P_Senders
		*sum_parz = 0;
		for (i=menum+1; i<nproc+menum; i++) {
		
			tag = 80+ (i%nproc);
			sender = i % nproc;
			MPI_Recv(sum_parz, 1, MPI_FLOAT, sender, tag, MPI_COMM_WORLD, status);
			*sum += *sum_parz;
				
		}
	} 
	else {//P_Senders spediscono a P_IDs
		tag = 80+menum;
		MPI_Send(sum, 1, MPI_FLOAT, ID_stamp, tag, MPI_COMM_WORLD);
	}
}


void II_strategia(int me, float *sum_p, float *sum_tot, int idx, int *p_2, int *p_2_plus_1, MPI_Status *stat) {
	//II strategia

	int i;
	int tag;
	int index;
	int menum;
	float *sum_parz;
	float *sum;
	int *pow_2;
	int *pow_2_plus_1;
	MPI_Status *status;
	
	menum = me;
	sum_parz = sum_p;
	sum = sum_tot;
	status = stat;	
		
	pow_2 = p_2;
	pow_2_plus_1 = p_2_plus_1;
	

	*sum_parz = 0;
	index = idx;
	for (i=0; i<index; i++) {//Numero passi di comunicazione
		
		if ((menum % pow_2[i]) == 0) {//partecipanti al passo i
			
			if((menum % pow_2_plus_1[i]) == 0) {//ricevi somma parziale
		
				tag = 80 + menum + pow(2,i);
				MPI_Recv(sum_parz, 1, MPI_FLOAT, menum + pow_2[i], tag, MPI_COMM_WORLD, status);
				*sum += *sum_parz;
		
			}
			else {//spedisci somma parziale
				
				tag = 80 + menum;
				MPI_Send(sum, 1, MPI_FLOAT, menum - pow_2[i], tag, MPI_COMM_WORLD);

			}
		}
	}	
}


void III_strategia(int me, float *sum_p, float *sum_tot, int idx, int *p_2, int *p_2_plus_1, MPI_Status *stat) {
	//III strategia

	int i;
	int tag;
	int index;
	int menum;
	float *sum_parz;
	float *sum;
	int *pow_2;
	int *pow_2_plus_1;
	
	MPI_Status *status;
	
	menum = me;
	sum_parz = sum_p;
	sum = sum_tot;
	status = stat;
	
	pow_2 = p_2;
	pow_2_plus_1 = p_2_plus_1;
	
	*sum_parz = 0;
	index = idx;
	for (i=0; i<index; i++) {//# passi di comunicazione
		if((menum % pow_2_plus_1[i]) < (pow_2[i])) {
			
			//Send precede Recv per evitare deadlock
			tag = 80 + menum + pow(2,i); //tag spedizione
			MPI_Send(sum, 1, MPI_FLOAT, menum + pow_2[i], tag, MPI_COMM_WORLD);
	
			tag = 80 + menum; //tag ricezione
			MPI_Recv(sum_parz, 1, MPI_FLOAT, menum + pow_2[i], tag, MPI_COMM_WORLD, status);
					
			*sum += *sum_parz;
							
		}
		else {
					
			tag = 80 + menum - pow_2[i]; //tag spedizione
			MPI_Send(sum, 1, MPI_FLOAT, menum - pow_2[i], 80 + menum - pow_2[i], MPI_COMM_WORLD);
			
			tag = 80+menum; //tag ricezione
			MPI_Recv(sum_parz, 1, MPI_FLOAT, menum - pow_2[i], 80 + menum, MPI_COMM_WORLD, status);
	
			*sum += *sum_parz;			
					
		}
	}
}


void stampa_ris(int me, float sum_tot, int ID_s, int N_str, MPI_Status *stat) {
	//stampa del risultato
	
	int i;
	int tag;
	int ID_stamp;
	int menum;
	float sum;
	int N_strat;
	MPI_Status *status;
	
	menum = me;
	sum = sum_tot;
	N_strat = N_str;
	status = stat;
	ID_stamp = ID_s;

	if (ID_stamp == -1) { //tutti devono stampare
		
		if (N_strat == 3)
			printf("[P%d] Somma: %f, strategia: %d.\n", menum, sum, N_strat);
		
		else {
			if (menum == 0) 
				printf("[P0] Non è possibile effettuare la stampa. \n");
		} 
					
	} 	
	
	else if (ID_stamp > -1) { //ID_stamp deve stampare
		
		if ((N_strat != 3) && (N_strat != 1)) { //Se non III e non è I, sono nella II strategia
			
			if (ID_stamp == 0) { //P0 ha già il riusltato,
				if (menum == ID_stamp) //può già stampare
					printf ("[P%d] Somma: %f, strategia: %d.\n", menum, sum, N_strat);
			}
			
			
			//altrimenti, se ID_stamp != P0
			else if (menum == 0) {	//P0 deve spedire sum a ID_stamp
					tag = 90 + ID_stamp;
					MPI_Send(&sum, 1, MPI_FLOAT, ID_stamp, tag, MPI_COMM_WORLD);
			
			}
			
				else if (menum == ID_stamp) { //Solo ID_stamp deve ricevere
						tag = 90 + menum;
						MPI_Recv(&sum, 1, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, status);	
						printf ("[P%d] Somma: %f, strategia: %d.\n", menum, sum, N_strat);
					
				}
		}
		
		else if (menum == ID_stamp)
				printf ("[P%d] Somma: %f, strategia: %d.\n", menum, sum, N_strat);

	}	
	
}

/*************************************************************/

