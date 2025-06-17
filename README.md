# WSL-git

Benvenuto in **WSL-git**, il repository che raccoglie tutti gli esercizi del corso di Laboratorio di Simulazione Numerica (LSN). Ogni esercitazione è organizzata in una cartella a sé stante, contenente il codice sorgente C++, eventuali notebook Jupyter per l’analisi e i file di supporto (Primes, seed.in). Il flusso tipico è: compilare con `make`, eseguire il programma C++ e infine lanciare il notebook per visualizzare grafici e risultati.

---

## 📁 Struttura del repository

```plaintext
WSL-git/
├── Esercitazione1/
├── Esercitazione2/
├── Esercitazione3/
├── Esercitazione4/
├── Esercitazione6/
├── Esercitazione7/
├── Esercitazione8/
├── Esercitazione9/
├── Esercitazione10/
├── Esercitazione11/
├── Esercitazione12/
└── README.md
```

* **EsercitazioneX/**: contiene tutto il necessario per svolgere l’esercizio X.
* Il file **README.md** spiega struttura, comandi di compilazione ed esecuzione.

---

## 🚀 Guida rapida all’esecuzione

Per tutte le esercitazioni (eccetto 11 e 12), il workflow è:

1. **Compilazione**: entra nella cartella e lancia `make`.
2. **Esecuzione**: esegui l’eseguibile generato (o `mpirun` per MPI).
3. **Analisi**: apri il notebook Jupyter per visualizzare i risultati.

Segue il dettaglio per ciascuna.

### Esercitazione 1, 2, 3, 8

Queste esercitazioni seguono lo stesso schema:

```bash
cd EsercitazioneX
make               # compila il programma
./EsX.Y.exe        # esegui ogni sotto-esercizio Y
jupyter notebook  # apri i notebook EsX.Y.ipynb
```

*Sostituisci X con 1, 2, 3 o 8 e Y con il numero del sotto-esercizio.*

### Esercitazione 4, 6, 7

Queste simulazioni richiedono solo la compilazione del simulatore e il lancio:

```bash
cd EsercitazioneX/SOURCE
make               # compila il simulatore
./simulator.exe    # esegui la simulazione
jupyter notebook   # apri il notebook Es0X.ipynb sotto Python/
```

*Dove X è 4, 6 o 7.*

### Esercitazione 9

```bash
cd Esercitazione9
make               # compila l’eseguibile 'tsp'
./tsp              # esegui il Genetic Algorithm
jupyter notebook  # apri Es09.1.ipynb
```

### Esercitazione 10 (MPI)

```bash
cd Esercitazione10
make               # compila l’eseguibile 'tsp_mpi'
mpirun -np 4 ./tsp_mpi  # esegui su 4 processi MPI oppure sostituire con il numero di processi desiderato 
jupyter notebook      # apri Es10.2.ipynb
```

### Esercitazione 11 e 12

Queste ultime due sono dedicate all’analisi dati con Jupyter:

```bash
cd Esercitazione11
jupyter notebook Es11.*.ipynb

cd ../Esercitazione12
jupyter notebook Es12.*.ipynb
```

---

## 🛠 Prerequisiti

* **make** e **g++** (comandi C++).
* **mpic++** e **mpirun** (Esercitazione10).
* **Jupyter Notebook** (per tutti i notebook).
* **Armadillo** (`libarmadillo-dev`) per Esercitazione4.

---



