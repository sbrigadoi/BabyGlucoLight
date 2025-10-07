DESCRIZIONE FILE

- ANALISI_CORRELAZIONE.m
Effettua la correlazione tra le metriche glicemiche e la differenza di rsFC calcolata tra l'ultima e la prima finestra di euglicemia.
Necessita di:
 --> data.mat : contiene le metriche glicemiche sottoforma di tabella e le metriche di rsFC sottoforma di matrice dove il numero di righe sono
		il numero di pazienti mentre il numero di colonne sono il numero di coppie di ROI. Contiene anche label delle metriche per rappresentazione
		In questo file il PD14 ha NaN in corrispondenza di ROI escluse poichè non coperte da array good channels della finestra migliore (vedi
		capitolo materiali e metodi tesi)
 --> data_all_roi.mat : stessa composizione di data.mat ma con anche le ROI escluse del PD14

- Estrazione_metriche_CGM.m 
Estrae le metriche glicemiche in automatico
Necessita di:
 --> PDXX.mat : file contenente la tabella con i riferimenti temporali in formato datetime e valori glicemici per ogni paziente (segnale CGM)
 --> idx_eu_results.mat : file contenete una tabella con gli indici di inizio e fine del primo ed ultimo intervallo di euglicemia completo (file nella 
			  cartella "3_ANALISI_NIRS_e_rsFC")
 --> pax_analysis_results_for_cca.mat e pax_event_report_for_cca.mat : risultati della fase di analisi glicemica ottenuti dal codice
								       "result_update_for_correlation_analysis.m" nella cartella 
								       "2_ANALISI_CGM_e_SELEZIONE_PAX"

- Estrazione_metriche_rsFC.m
Estrae le metriche relative alla differenza di rsFC
Necessita di:
 --> 'PD_XX_FIRST_rsFC.mat' / 'PD_XX_LAST_rsFC.mat' : ottenute dai vari file rsFC.m nella cartella "2_ANALISI_CGM_e_SELEZIONE_PAX".
						      Ogni file contiene una tabella con il nome delle ROI, le relative coordinate 3D 
						      e un flag per identificare se è attiva o meno.
 --> Le matrici rsFC per HbO e HbR, matrici p_value HbO e HbR, valore p_value Bonferroni e un vettore con i nomi delle ROI per rappresentare i risultati

OSS. Per ottenere i risultati della correlazione basta avere a disposizione il file data.mat e il codice ANALISI_CORRELAZIONE.m