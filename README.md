# Prodotto Matrice Vettore utilizzando MPI

<br>

Gli algoritmi implementati effettuano il prodotto tra:
* una Matrice A di dimensione m*n; 
* un Vettore v di dimensione n.

Il risultato ottenenuto è un Vettore w di dimensione m.

<br>

Il prodotto Matrice-Vettore viene effettuato utilizzando diverse strategie, implementata dai seguenti algoritmi:
* algoritmo sequenziale (il cui file è chiamato prodMatVetSeq.c);
* algoritmo parallelo che implementa la strategia a BLOCCHI di COLONNE (prodMatVetCol.c);
* algoritmo parallelo che implementa la strategia a BLOCCHI di RIGHE e COLONNE (prodMatVetRowCol.c).

Per la descrizione e l'analisi dei primi due algoritmi, consultare il file "Esercitazione2_Michele_Delli_Paoli.pdf".
Per la descrizione e l'analisi dell'algoritmo che implementa l'ultima strategia, consultare il file "Esercitazione3_Michele_Delli_Paoli.pdf".
