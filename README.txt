Eseguire il file main per avere tutti i parametri e gli script in ordine.
Sono presenti 3 stimatori: EKF è quello che riprende passo passo l'articolo di Bloesch, gli altri 2 (trackEKF e navEKF) 
utilizzano stimatori definiti da Matlab.
Le covarianze dei rumori di misura sono definite nel file measurement. Se il rumore dei bias viene portato a zero e si elimina
il passo di correzione dello stimatore il sistema funziona come previsto, mentre in fase di correzione lo stato diverge.