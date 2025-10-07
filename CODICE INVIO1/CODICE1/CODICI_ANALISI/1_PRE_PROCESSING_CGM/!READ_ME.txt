DESCRIZIONE CODICE:

!! Ogni codice è commentato nel dettaglio !!

- Cleaning_and_gap_correction.m = esegue il pre-processing dei dati rimuovendo i campioni di calibrazione e sostituendo inizialmente i gap temporali con
			          il valore mediano a cavallo del gap temporale. Restituisce il segnale pulito e alcune informazioni utili relative ai 
			          gap temporli identificati. Rimuovendo i commenti nel codice si possono salvare agilmente i risultati.
- calibration_cleaning.m        = funzione che elimina i campioni relativi alla calibrazione del segnale CGM creando anche una nuova griglia temporale
- gap_correction.m 		= funzione che identifica e corregge il gap temporale con il valore mediano a cavallo del buco temporale. In uscita produce
				  una tablella contenente il segnale e alcune informazioni utili

FILE NECESSARI:

- Cleaning_and_gap_correction.m
  --> PD1.mat = file contenente per ciascun paziente una tabella ottenuta caricando manualmente e salvando singolarmente i dati per ogni paziente dal file
	        Excell relativo all'acquisizione in TIN del segnale glicemico. Utilizzare il toolbox 'Import Data' di Matlab. Operazione da eseguire 
		manualmente prima di lanciare i codici.