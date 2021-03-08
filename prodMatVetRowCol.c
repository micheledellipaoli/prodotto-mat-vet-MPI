#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h> 

#define limiteStampaM 20
#define limiteStampaN 20

//Funzione che computa il prodotto matrice-vettore
void prod_mat_vett(double *a, double *v, int ROWS, int COLS, double *w){
    int i, j;
    
    for(i=0;i<ROWS;i++){
        w[i]=0;
        for(j=0;j<COLS;j++){ 
            w[i] += a[i*COLS+j]* v[j]; 
        } 
    }    
}


int main(int argc, char **argv) {

	int nproc;						// Numero di processi totale            	
	int me, me_grid, me_row, me_col;      			// id del processore attuale nel contesto: globale, griglia, riga, colonna

	int dim, *ndim, reorder, *period;			//Variabili utili alle operazioni relative alla Griglia di processori
	int *coord_grid, coord_row, coord_col;
	int *belongs;       	

	int p, q;						// Indici della griglia di processori

	double *A;						// Matrice Originale m*n
	int m, n;  						// Dimensione della matrice A originale	

	double *v;						// Vettore da moltiplicare di dimensione n
	double *w;						// Vettore risultato di dimensione m

	double *local_A, *local_AT, *local_w;			// Vettori locali ai processori appartenenti alla PRIMA COLONNA. 
							// local_A e` formato da (m/p) RIGHE e n COLONNE, local_AT e` formato da n RIGHE e (m/p) COLONNE, local_w e` di dimensione m/p

	double *local_v, *local_result; 			// Vettori locali a ciascun processore della griglia

	double *local_blockT;					// Sottomatrice di local_AT di (m/p)*(n/q) elementi.
							// local_block T e` locale ad ogni processore della griglia, e viene ricevuta dal processore con indice 0 della propria riga.

	double *local_block;					// Trasposta di local_blockT

	int local_m, local_n;					// Dimensione dei dati locali per la matrice local_A       
	int local_mT, local_nT;  				// Dimensione dei dati locali per la matrice local_AT    
	
	int i, j;                				// Iteratori vari     	
	
	double T_inizio, T_fine, T_max; 			//Variabili per la misurazione del tempo di esecuzione


	/*Attiva MPI*/
	MPI_Init(&argc, &argv);

	/*Trova il numero totale dei processi*/
	MPI_Comm_size (MPI_COMM_WORLD, &nproc);
	/*Da ad ogni processo il proprio numero identificativo*/
	MPI_Comm_rank (MPI_COMM_WORLD, &me);

	// Se sono il processore root
	if(me == 0){

		printf("Inserire numero di righe m della Matrice: \n"); 
    		scanf("%d",&m); 
    
    		printf("Inserire numero di colonne n della Matrice: \n"); 
    		scanf("%d",&n);

    		printf("Considerando che nProc=%d, inserire il numero di righe p della griglia della Topologia: \n", nproc); 
    		scanf("%d",&p);

    		// Calcolo di q, ovvero indice di righe della Topologia di processori.
    		q = nproc/p;

    		// Numero di righe della matrice locale local_A da inviare a ciascun processore della PRIMA COLONNA 
    		local_m = m/p;

    		// Numero di elementi del vettore locale local_v da inviare a ciascun processore della PRIMA RIGA 
    		local_n = n/q;
    
    		// Allocazione spazio di memoria
    		A = (double*) malloc(m*n*sizeof(double));
    		v = (double*) malloc(n*sizeof(double));
    		w = (double*) malloc(m*sizeof(double)); 
    
    		//Inizializzazione e stampa della matrice A con m RIGHE e n COLONNE
    		printf("\nA = \n"); 
    		for (i=0;i<m;i++){
       			for(j=0;j<n;j++){
            			if (j==0){
                			A[i*n+j]= 1.0/(i+1)-1;
            			}else{
                			A[i*n+j]= 1.0/(i+1)-pow(1.0/2.0,j); 
            			}
				printf("%f ", A[i*n+j]);
        		}
        		printf("\n");
    		}

    		//Inizializzazione e stampa del vettore v
    		printf("\nv = \n");
    		for (j=0;j<n;j++){
    			v[j]=j; 
    			printf("%f ", v[j]);
    		}
    		printf("\n");

	} // fine me==0

	// Spedisco m, n, local_m, local_n
	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);  
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&local_m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&local_n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	p = m/local_m;
	q = n/local_n;

	// Allocazione spazio e inizializzazione del vettore contenente le dimensioni della griglia della Topologia
	dim = 2;
    	ndim = (int*) malloc(dim*sizeof(int));
    	ndim[0] = p;
    	ndim[1] = q;
	
	coord_grid = (int*) malloc(dim*sizeof(int));	
	
	belongs = (int*) malloc(dim*sizeof(int));

    	period=(int*)calloc(dim,sizeof(int));
    	period[0]=0;
    	period[1]=0;
    	reorder=0;

    	// Definizione della Griglia bidimensionale di processori.
    	// Si ottiene in output il contesto comm_grid
	MPI_Comm comm_grid;
    	MPI_Cart_create(MPI_COMM_WORLD, dim, ndim, period, reorder, &comm_grid);
	
    	// Ottenimento dell' id del processore attuale nel Contesto Griglia appena creato.
    	MPI_Comm_rank(comm_grid, &me_grid);
	
    	// Ottenimento delle coordinate del processore attuale me_grid nella Griglia il cui contesto e` comm_grid
    	MPI_Cart_coords(comm_grid, me_grid, dim, coord_grid);

	// Creazione Contesto dei processori appartenenti alla STESSA RIGA della Griglia (comm_row)
	MPI_Comm comm_row;
	belongs[0] = 0;
    	belongs[1] = 1;
	MPI_Cart_sub(comm_grid, belongs, &comm_row);
	MPI_Comm_rank(comm_row, &me_row);
	MPI_Cart_coords(comm_row, me_row, 1, &coord_row);

	// Creazione Contesto dei processori appartenenti alla STESSA COLONNA della Griglia (comm_col)
        MPI_Comm comm_col;
	belongs[0] = 1;
    	belongs[1] = 0;
	MPI_Cart_sub(comm_grid, belongs, &comm_col);
	MPI_Comm_rank(comm_col, &me_col);
	MPI_Cart_coords(comm_col, me_col, 1, &coord_col);

	// La Barrier assicura che tutti conoscano le proprie coordinate prima di continuare
	MPI_Barrier(MPI_COMM_WORLD); 

	// Allocazione del vettore local_v di dimensione local_n=(n/q)
	local_v = (double*) malloc(local_n*sizeof(double));

	MPI_Barrier(MPI_COMM_WORLD); 

	//Codice eseguito solo dai processori appartenenti alla PRIMA RIGA della Griglia
	if(coord_grid[0] == 0){

		// Distribuzione del vettore v da P00 ai processori della PRIMA RIGA che ricevono la propria porzione in local_v
		MPI_Scatter(v, local_n, MPI_DOUBLE, local_v, local_n, MPI_DOUBLE, 0, comm_row);

		// Stampa del vettore local_v ricevuto da ciascun processore della PRIMA RIGA
		if(m<=limiteStampaM && n<=limiteStampaN){
			printf("\nlocal_v [%d] = \n", me);
			for(i=0; i<local_n; i++){
				printf("%f\t", local_v[i]);
			}
			printf("\n");
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD); 

	//Codice eseguito solo dai processori appartenenti alla PRIMA COLONNA della Griglia
	if(coord_grid[1] == 0){
	
		// Allocazione del vettore local_A formato da local_m=(m/p) RIGHE e n COLONNE.
		local_A = (double*) malloc(local_m*n*sizeof(double));

		// Il processore P00 invia un blocco (di local_m=(m/p) RIGHE e n COLONNE, ovvero m/p*n ELEMENTI) ai processori della PRIMA COLONNA
		MPI_Scatter(A, local_m*n, MPI_DOUBLE, local_A, local_m*n, MPI_DOUBLE, 0, comm_col);

		// Stampa di local_A
		if(m<=limiteStampaM && n<=limiteStampaN){
			printf("\nlocal_A [%d] = \n", me); 
			for(i=0; i<local_m; i++){
	    			for(j=0; j<n; j++){
	        			printf("%f\t", local_A[i*n+j]);
	    			}
	    			printf("\n");
			}
			printf("\n");
		}

		// Allocazione del vettore local_AT di n RIGHE e local_m=(m/p) COLONNE.
		local_AT = (double*) malloc(n*local_m*sizeof(double));
	
		int local_mT = n;
		int local_nT = local_m;

		// Calcolo della matrice trasposta local_AT
		for(i=0; i<local_mT; i++){
	    		for(j=0; j<local_nT; j++){
				local_AT[i*local_nT+j] = local_A[j*n+i];
	    		}
		}
		
		// Stampa di local_AT
		if(m<=limiteStampaM && n<=limiteStampaN){
			printf("\nlocal_AT [%d] = \n", me); 
			for(i=0; i<local_mT; i++){
	    			for(j=0; j<local_nT; j++){
	        			printf("%f\t", local_AT[i*local_nT+j]);
	    			}
	    			printf("\n");
			}	
			printf("\n");
		}

	
	}
	
	MPI_Barrier(MPI_COMM_WORLD); 

	// Invio di local_v da ciascun processore P0j ai processori P0j, P1j,..., Pp-1,j per (0<=j<q).
	// In altre parole, ogni processore della PRIMA RIGA della Griglia invia una copia di local_v ai processori appartenenti alla sua STESSA COLONNA.
	MPI_Bcast(local_v, local_n, MPI_DOUBLE, 0, comm_col);	
	
	// Stampa del vettore local_v ricevuto da ciascun processore della Griglia.
	if(m<=limiteStampaM && n<=limiteStampaN){
		printf("\nlocal_v [%d] = \n", me);
		for(i=0; i<local_n; i++){
			printf("%f\t", local_v[i]);
		}
		printf("\n");
	}
	
	MPI_Barrier(MPI_COMM_WORLD); 	
	
	// Allocazione del vettore local_blockT che rappresenta una sottomatrice di local_AT, di local_n=(n/q) RIGHE e local_m=(m/p) COLONNE.
	local_blockT = (double*) malloc(local_n*local_m*sizeof(double));
 
	//Distribuzione della matrice trasposta local_AT dai processori della prima COLONNA (Pi0 con 0<= i < m) ai processori appartenenti alla STESSA RIGA (Pij con 0<= j < n).
	//Ciascun processore riceve un blocco di local_n=(n/q) RIGHE e local_m=(m/p) COLONNE che memorizza in local_blockT.
	MPI_Scatter(local_AT, local_n*local_m, MPI_DOUBLE, local_blockT, local_n*local_m, MPI_DOUBLE, 0, comm_row);

	MPI_Barrier(MPI_COMM_WORLD);
	
	// Stampa di local_blockT.
	if(m<=limiteStampaM && n<=limiteStampaN){
		printf("\nlocal_blockT [%d] = \n", me); 
		for (i=0; i<local_n; i++){
   			for(j=0; j<local_m; j++){
        			printf("%f\t", local_blockT[i*local_m+j]);
    			}
    			printf("\n");
		}
		printf("\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Allocazione del vettore local_block di local_m=(m/p) RIGHE e local_n=(n/q) COLONNE.
	local_block = (double*) malloc(local_m*local_n*sizeof(double));

	// Calcolo di local_block, ottenuto dalla trasposizione di local_blockT.
	for (i=0; i<local_m; i++){
   		for(j=0; j<local_n; j++){
        		local_block[i*local_n+j] = local_blockT[j*local_m+i];
    		}
	}

	// Stampa di local_block
	if(m<=limiteStampaM && n<=limiteStampaN){
		printf("\nlocal_block [%d] = \n", me); 
		for (i=0; i<local_m; i++){
   			for(j=0; j<local_n; j++){
        			printf("%f\t", local_block[i*local_n+j]);
    			}
    			printf("\n");
		}
		printf("\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	// Allocazione del vettore local_result di local_m ELEMENTI.
	local_result = (double*) malloc(local_m*sizeof(double));

	// Codice eseguito dai processori appartenenti alla PRIMA COLONNA della Griglia.
	if(coord_grid[1] == 0){

		// Allocazione del vettore local_w di local_m=(m/p) ELEMENTI in cui raggruppare, tramite Reduce, i vettori local_w computati da ciascun processore.
		local_w = (double*) malloc(local_m*sizeof(double));
	}

	MPI_Barrier(MPI_COMM_WORLD);
	T_inizio=MPI_Wtime();
	
	// Computazione del prodotto local_block*local_v, memorizzato nel vettore local_result.
	prod_mat_vett(local_block, local_v, local_m, local_n, local_result);

	// Ogni processore ha computato local_result. 
	// Bisogna sommare i vettori local_result di tutti i processori appartenenti alla STESSA RIGA, ed inviare il risultato al PRIMO processore della RIGA.
	// Tutti i processori della PRIMA COLONNA raccolgono tale somma nel vettore local_w.
	MPI_Reduce(local_result, local_w, local_m, MPI_DOUBLE, MPI_SUM, 0, comm_row);
	
	MPI_Barrier(MPI_COMM_WORLD);
	T_fine=MPI_Wtime()-T_inizio;

	// Stampa del vettore local_result di dimensione local_m=(m/p) computato da ciascun processore effettuando il prodotto tra local_block e local_v.
	if(m<=limiteStampaM && n<=limiteStampaN){
		printf("\nlocal_result [%d] = \n", me);
		for(i=0; i<local_m; i++){
			printf("%f\t", local_result[i]);
		}
		printf("\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Codice eseguito dai processori appartenenti alla PRIMA COLONNA della Griglia.
	if(coord_grid[1] == 0){

		// Stampa del vettore local_w
		if(m<=limiteStampaM && n<=limiteStampaN){
			printf("\nlocal_w [%d]: \n", me);
			for(i=0; i<local_m; i++){
				printf("%f\t", local_w[i]);
			}
			printf("\n");
		}
	}

	printf("\n");

	MPI_Barrier(MPI_COMM_WORLD);
	

	/* Raccoglie il tempo di esecuzione Massimo e lo invia al processore 0*/
	MPI_Reduce(&T_fine, &T_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	
	// Codice eseguito solo dai processori appartenenti alla PRIMA COLONNA della Griglia
	if(coord_grid[1] == 0){
		// Bisogna raggruppare i vettori local_w, di dimensione local_m=(m/p), dei processori della PRIMA COLONNA nel vettore w, di dimensione m, del processore radice P00. 
		MPI_Gather(local_w, local_m, MPI_DOUBLE, w, local_m, MPI_DOUBLE, 0, comm_col);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Il processore 0 stampa la soluzione
	if(me==0){ 
    		printf("\nw = \n"); 
    		for(i = 0; i < m; i++){
        		printf("%f ", w[i]);
    		}
		printf("\n");
   
    		printf("\nTempo calcolo locale: %lf\n", T_fine);
    		printf("\nMPI_Reduce max time: %f\n",T_max);

    	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	
	free(local_blockT);
	free(local_block);
	free(local_result);
	free(local_v);
	free(ndim);
	free(belongs);
	free(period);
	if(coord_grid[1] == 0){
		free(local_A);
		free(local_AT);
		free(local_w);
	}
	if(me==0){
		free(A);
    		free(v);
    		free(w);
	}

	MPI_Finalize (); /* Disattiva MPI */
	return 0;  
}
