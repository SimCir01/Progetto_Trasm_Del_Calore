## CORSO DI TRASMISSIONE DEL CALORE IN APPLICAZIONI BIOMEDICALI (A.A. 2023-24)
## Progetto finale: 
## Analisi termica del transitorio 1D relativo ad una massa tumorale nel tessuto mammario a seguito di una termoablazione a microonde

## Autori
* **Benedetta Masone** 
* **Simone Cirnelli** 

## Descrizione
Questo progetto MATLAB è stato sviluppato per analizzare il comportamento termico del tessuto mammario
contenente una massa tumorale durante una procedura di termoablazione a microonde. 


Il progetto è composto da un file principale (main.m) e diversi file funzione che supportano i calcoli e la generazione dei risultati.

**************************************** 'main' ****************************************
Tale file rappresenta il punto di ingresso principale per l'esecuzione delle simulazioni.
Esegue i seguenti passaggi:
    - Definizione dello studio: Permette all'utente di scegliere l'obiettivo dello studio (campo di temperatura, stima della potenza, verifica della stabilità).
    - Raccolta parametri: Richiede all'utente di inserire i parametri necessari per la simulazione (numero di nodi, durata della simulazione, passo temporale).
    - Risoluzione numerica: Esegue il calcolo del campo di temperatura e genera grafici e video dei risultati.

Il programma permette di selezionare tra tre diversi obiettivi:

    [0] Campo di temperatura
    [1] Stima della massima potenza fissata una percentuale di danno
    [2] Verifica condizioni di stabilità per il metodo esplicito

Durante l'esecuzione, il programma richiederà all'utente di inserire vari parametri,
come il numero di nodi, la durata della simulazione, e il passo temporale.
In base all'obiettivo selezionato, verranno richiesti ulteriori parametri specifici.

Il programma genera vari grafici e video che mostrano l'andamento della temperatura nel tempo e nello spazio, oltre alla percentuale di danno nei vari tessuti.
Grafici e Video Generati

    - Progressione nel tempo della distribuzione spaziale di temperatura
    - Progressione nello spazio della distribuzione temporale di temperatura
    - Plot 3D della temperatura nello spazio e nel tempo
    - Rappresentazione grafica strati di tessuti

Le funzioni utilizzate nel progetto includono:

- La funzione "Arr" per il calcolo della percentuale di danno all'interno di un tessuto generico;
- La funzione "cond" consente, conoscendo la posizione nel dominio  spaziale (quindi i e dx), di definire in quali tessuti vengono a trovarsi
  i nodi i-esimo, (i-1)-esimo e (i+1)-esimo;
- La funzione "h_aria" consente il calcolo del coefficiente convettivo nell'interazione termica con l'aria (convezione naturale);
- La funzione "h_sangue" consente il calcolo del coefficiente convettivo % nell'interazione termica con il sangue (convezione forzata);
- La funzione "matrix" serve per il calcolo delle matrice A  (matrice coefficienti) e b (termini noti) e 
  il campo di temperatura T(x,t);
- La funzione "SAR" consente di ottenere, scelta in input la posizione spaziale, il corrispettivo valore SAR 
  in base all'andamento di riferimento;
- La funzione "Sigma" consente, conoscendo la posizione nel dominio  spaziale (quindi i e dx), di definire gli eventuali spessori
  dei nodi (i-1)-esimo e (i+1)-esimo, in modo alla fine % da poter ricavare il parametro Sigma necessario nei bilanci;
- La funzione "T_sup_tumore" serve per il calcolo delle temperatura lungo la superficie esterna del tumore e per l'impostazione della 
  condizione temporale sui suoi bordi;
- La funzione "tissue" consente di plottare l'immamgine dei vari tessuti che compongono il dominio spaziale 
  con la discretizzazione impostata dei nodi