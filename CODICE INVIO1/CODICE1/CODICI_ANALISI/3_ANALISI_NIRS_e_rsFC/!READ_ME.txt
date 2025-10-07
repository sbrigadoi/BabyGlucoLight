DESCRIZIONE FILE

!!!! Ogni codice è commentato nel dettaglio !!!!

- CARTELLE: ogni cartella contiene i codici utilizzati per sincronizzare i dati NIRS e per effettuare la rsFC per ogni finestra da 5' identificata
	    nel primo ed ultimo intervallo di euglicemia.

	    !!! Dentro i file per la rsFC c'è una parte di codice commentata che contiene l'implementazione del metodo di correzione dei MA proposto
		da Yang et al.. Sebbene non sia stato utilizzato nella mia tesi (si è preferito il filtraggio Wavelet) potrebbe essere utile per futuri
		studi. Vedere tesi per tutte le informazioni e considerazioni su questo metodo !!!

		Yang M, Xia M, Zhang S, et al. Motion artifact correction for resting-state neonatal functional near-infrared spectroscopy through adaptive estimation of physiological oscillation denoising. Neurophotonics. 2022;9(4):045002. doi:10.1117/1.NPh.9.4.045002

- SINCRONIZZAZIONE NIRS-CGM

Data l'estrema complessità dei file NIRS, la sincronizzazione è stata effettuata manualmente per ogni paziente modificando e commentando opportunamente
il codice in base alla composizione del dataset.

File necessari:
 --> PDXX.mat 	    : per ogni paziente il file contenente il tracciato CGM (table: riferimento temporale formato datetime e valore numerico glicemia)
 --> idx_eu_int.mat : presente in questa cartella che contiene per ogni paziente l'indice di inizio e fine della prima ed ultima finestra di euglicemia
		      individuata dall'analisi CGM. Questo indice serve per estrarre le informazioni temporali della finestra di euglicemia.
 --> file acqusizione NIRS : per ogni paziente serve il database dell'acquisizione NIRS nel formato grezzo .nirs

Il codice restituisce:
 --> PDXX_synchro_result.mat: contenenti gli indici e informazioni utili sul primo ed ultimo intervallo di euglicemia sincronizzato
 --> PDXX_NIRS_first_eu_tot.mat / PDXX_NIRS_last_eu_tot.mat : contenente le porzioni di tracciato NIRS sincronizzate 

- ANALISI rsFC
Per ogni intervallo di euglicemia e per ogni paziente il codice esegue la rsFC (vedere capitolo materiali e metodi tesi)

File necessari: 
 --> PDXX_synchro_result.mat: in uscita da codice sincronizzazione
 --> PDXX_NIRS_first_eu_tot.mat / PDXX_NIRS_last_eu_tot.mat : in uscita da codice sincronizzazione
 --> Jacobiano e mesh ottenute dal codice "get_jac_infant_head.m" 
 --> Head Model estratto dal database Brigadoi et al. sulla base dell'età gestazionale del paziente analizzato
     Brigadoi S, Aljabar P, Kuklisova-Murgasova M, Arridge SR, Cooper RJ. A 4D neonatal head model for diffuse optical imaging of pre-term to term infants. Neuroimage. 2014 Oct 15;100:385-94. doi: 10.1016/j.neuroimage.2014.06.028. Epub 2014 Jun 18. PMID: 24954280.
 --> vol2gm.mat : per ogni età gestazionale, può essere creata all'interno del codice stesso togliendo i commenti nella sezione specifica

 --> Funzioni varie presenti nella cartella 'FUNCTION', il codice le carica in automatico all'occorrenza.
