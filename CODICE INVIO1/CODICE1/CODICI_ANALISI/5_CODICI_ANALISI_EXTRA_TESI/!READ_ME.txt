In questa cartella ci sono i codici che ho utilizzato per creare le varie figure inserite nella tesi e per effettuare le analisi secondarie.

- codice_stampa_rsFC_Appendince.m
Stampa le matrici rsFC che ho inserito nell'appendice della tesi. 
Richiede:
 --> 'PD_XX_FIRST_rsFC.mat' e 'PD_XX_LAST_rsFC.mat' : risultati dell'analisi rsFC (ottenuti da codici cartella "2_ANALISI_CGM_e_SELEZIONE_PAX")

- controllo_prima_dopo_significativo.m
Rileva se le coppie di ROI hanno subito una variazione significativa nell'intervallo di tempo tra la prima ed ultima finestra di euglicemia.
Richiede:
 --> 'PD_XX_FIRST_rsFC.mat' e 'PD_XX_LAST_rsFC.mat' : risultati dell'analisi rsFC (ottenuti da codici cartella "2_ANALISI_CGM_e_SELEZIONE_PAX")

- controllo_rsFC_attivazioni.m
Conta il numero di ROI coppie di ROI intra- ed interlobo che hanno subito una variazione di valore superiore a 0.5 nel periodo di tempo tra la prima 
ed ultima finestra di euglicemia.
Richiede:
 --> 'PD_XX_FIRST_rsFC.mat' e 'PD_XX_LAST_rsFC.mat' : risultati dell'analisi rsFC (ottenuti da codici cartella "2_ANALISI_CGM_e_SELEZIONE_PAX")

- risultati_differenze_rsFC.m
Conta il numero di ROI che hanno subito una variazione di valore superiore a diverse soglie indicate, nel periodo di tempo tra la prima ed ultima finestra
di euglicemia.
Richiede:
 --> 'PD_XX_FIRST_rsFC.mat' e 'PD_XX_LAST_rsFC.mat' : risultati dell'analisi rsFC (ottenuti da codici cartella "2_ANALISI_CGM_e_SELEZIONE_PAX")

- Test_statistici_analisi_CGM.m
Effettua i test statistici tra varie metriche relative ad eventi ipo ed iperglicemici di diversa natura, vedi capitolo Risultati della tesi.