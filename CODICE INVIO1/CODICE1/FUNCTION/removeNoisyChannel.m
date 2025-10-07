function removeCh = removeNoisyChannel(nirsData,dRange,SNRrange)
% @nirsData = dati NIRS, righe time points colonne segnali per ogni canale
% per ogni lunghezza d'onda (x2)
% @dRange = [min max] valori per media segnale
% @SNRrange = soglia di SNR sotto la quale il segnale è considerato rumore
% @removeCh = vettore colonna 1 = segnale buono 0 = segnale da eliminare
% -------------------------------------------------------------------------
% La funzione prende in ingresso il segnale NIRS, il range di lavoro
% dell'intesità e il range di SNR accettabili. Restituisce un vettore di
% lunghezza numero di canali contenente 1 se segnale da tenere 0 se segnale
% da eliminare perchè fuori range intensità o sotto soglia SNR. Un segnale
% è da eliminare se almeno una lunghezza d'onda è bad channel.

nCh = size(nirsData,2)/2;

% Calcolo media, std e SNR dei dati
meanData = mean(nirsData);
stdData  = std(nirsData,[],1);
SNRData  = meanData./stdData;

% Inzialiazzo il vettore contenente gli indici dei segnali da eliminare con
% tutti zeri (0 = bad channel, 1 = good channel). Di conseguenza cerco 
% i segnali accettabili e gli metto ad 1 (cioè media dentro range 
% accettabile di intensità E SNR maggiore)
removeCh = zeros(nCh*2,1);
removeCh( meanData>dRange(1) & meanData<dRange(2) & SNRData>SNRrange) = 1;

% Se un canale ha 0 anche quello relativo all'altra lunghezza d'onda deve
% essere cancellato. 

% Creo la matrice tmp che ha sulla primma colonna il valore di removeCh
% corrispondente prima lunghezza d'onda e sulla seconda colonna il valore
% removeCh della seconda lunghezza d'onda. 
tmp1 = removeCh(1:nCh);
tmp2 = removeCh(nCh+1:end);
tmp = [tmp1 tmp2];
% Faccio somma per righe, se una riga ha 2 allora quella riga ha segnale da
% tenere, viceversa lascio zero. Poi per avere vettore in uscita con
% risultati uguali per ogni lunghezza d'onda faccio la copia della prima. 
removeCh = zeros(nCh*2,1);
removeCh(sum(tmp,2)==2) = 1;
removeCh(nCh+1:end) = removeCh(1:nCh);
end
