DESCRIZIONE CODICE:

!! Ogni codice è commentato nel dettalgio !!

- pax_selection_MAIN.m                   = codice principale per effettuare l'analisi glicemica e la selezione dei pazienti (vedi capitolo tesi per criteri)
- pax_selection_glycemic_analysis.m      = funzione che effettua in automatico l'analisi glicemica, richiede le due funzioni:
	- pax_selection_euglycemia_analysis.m = funzione che effettua in automatico l'analisi degli intervalli di euglicemia
	- pax_selection_events_analysis.m     = funzione che effettua in automatico l'analisi degli eventi gicemici

- result_update_for_correlation_analysis = codice utilizzato per modificare manualmente alcuni risultati ottenuti nell'analisi glicemica
					   da utilizzare per effettuare l'analisi di correlazione finale con i dati rsFC NIRS (nei commenti 
					   iniziali del codice viene spiegato nel dettaglio cosa modifica questo codice)			

FILE NECESSARI

- pax_selection_MAIN.m 
  --> PD1.mat = structure con numero progressivo per ogni paziente del database contenente
		  PD.PD = table due colonne nella prima il riferimento temporale nel formato 'datetime', nella seconda i valori glicemici ripuliti
			e con al posto dei gap temporali il valore mediano. L'analisi glicemica effettua l'interpolazione di Fonda et. al (ev. si 
			possono salvare i risultati togliendo i commenti nel codice MAIN)
		  PD.infogap = table con tre colonne, nella prima l'indice di inizio del gap, nella seconda l'indice di fine e nella terza il numero
			     di campioni mancanti
		  PD.n_gap = numero di campioni mancanti

- result_update_for_correlation_analysis.m
  Entrambi i file necessari sono ottenuti in pax_selection_MAIN.m togliendo i commenti a fine codice
  --> pax_analysis_results.mat
  --> pax_event_report.mat


!!!!! OSSERVAZIONI IMPORTANTI SUL CODICE !!!!!

VENGONO CONTATI TUTTI GLI EVENTI, QUELLI CRITICI/INIZIANO/FINISCONO FUORI CGM SONO SEGNALATI

EVENTO INIZIALE MANTENUTO MA SEGNALATO 
EVENTO FINALE MANTENUTO MA SEGNALATO   
EVENTO CON GAP CRITICO SEGNALATO       
EUGLICEMIA CON GAP CRITICO DIVISO        (Diviso in due parti, una finesce prima dell'inizio del gap e una comincia dopo la fine del gap)

EVENTI IN OVERLAP = vengono uniti in un unico evento
EVENTI CHE INIZIANO NELLO STESSO INDICE = vengono uniti in un unico evento

Soglia GAP EU = 45' (è 25-esimo percentile della durata di tutti gli eventi)

OSS: C'è un intervallo di 5' tra la fine EU e inizio evento (stessa cosa con fine evento ed inizio EU) a causa di come è stato definito un evento.
     Va bene così, i risultati sono coerenti.
	 
