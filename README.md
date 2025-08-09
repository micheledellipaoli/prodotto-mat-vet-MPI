# Prodotto matrice-vettore utilizzando MPI
Questo repository contiene degli algoritmi paralleli sviluppati utilizzando MPI (Message Passing Interface).

Nello specifico, gli algoritmi implementati effettuano il prodotto tra:
* una Matrice A di dimensione m*n; 
* un Vettore v di dimensione n.

Il risultato ottenenuto è un Vettore w di dimensione m.

<br>

Il prodotto matrice-vettore viene effettuato utilizzando diverse strategie, implementata dai seguenti algoritmi:
* algoritmo sequenziale (il cui file è chiamato prodMatVetSeq.c);
* algoritmo parallelo che implementa la strategia a BLOCCHI di COLONNE (prodMatVetCol.c);
* algoritmo parallelo che implementa la strategia a BLOCCHI di RIGHE e COLONNE (prodMatVetRowCol.c).

<br><br>

Per la descrizione e l'analisi dei primi due algoritmi, consultare il file "Esercitazione2_Michele_Delli_Paoli.pdf".<br>
Per la descrizione e l'analisi dell'algoritmo che implementa l'ultima strategia, consultare il file "Esercitazione3_Michele_Delli_Paoli.pdf".
