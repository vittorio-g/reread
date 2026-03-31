# Buone Idee per Migliorare il ReReReRe

Data: 2026-03-31. Idee emerse dalla discussione, ordinate per potenziale stimato.

---

## 1. Score di coerenza PER FATTORE, poi aggregati (priorita' alta)

### Il problema
Adesso il ReReReRe calcola UN singolo score di coerenza pesando tutte le coppie within-factor insieme. Questo diluisce il segnale quando il careless e' parziale (es. perde attenzione a meta' questionario) o quando un rispondente attento ha risposte estreme su un singolo fattore.

### L'idea
Calcolare uno score di coerenza separato per ogni fattore estratto dall'EFA, poi aggregare gli score fattoriali in un indice composito.

### Come funzionerebbe
1. L'EFA assegna ogni item a un fattore (gia' fatto in EFA-D)
2. Per ogni fattore con >=2 item: calcola `rowCor_weighted()` solo sulle coppie di quel fattore
3. Permutazione baseline per ogni fattore separatamente -> z-score per fattore
4. Aggregazione finale

### Perche' dovrebbe funzionare (giustificazione statistica)
- Con nF=10 e ipf=6, abbiamo **10 mini-test di coerenza indipendenti**
- Un rispondente careless avra' z basso su TUTTI i fattori -> la media dei z sara' bassa
- Un rispondente attento con opinioni estreme su 1 fattore avra' 1/10 basso ma 9/10 alti -> la media resta alta
- La proporzione di fattori "incoerenti" e' una statistica piu' robusta del singolo z globale
- E' come avere k test indipendenti — la reliability cresce con sqrt(k)

### Possibili varianti di aggregazione
- **Media z**: semplice, robusto, ma sensibile a fattori con pochi item
- **Mediana z**: piu' robusto agli outlier
- **Min z**: cattura il fattore peggiore — buono per carelessness parziale
- **% fattori con z < soglia**: binario, molto interpretabile ("quanti fattori sono incoerenti?")
- **Varianza z**: careless = bassa varianza + bassa media; attentivi = alta varianza + alta media

### Rischi
- Fattori con pochi item (ipf=3) producono z-score rumorosi per fattore
- Se un fattore ha solo 2 item, la "correlazione" per fattore e' un singolo prodotto — troppo rumoroso
- La parallel analysis potrebbe sotto/sovrastimare nF

---

## 2. Baseline random cross-factor (priorita' media-alta)

### Il problema
La permutation baseline attuale campiona coppie random tra TUTTI gli item, senza distinzione. Alcune coppie random cadono per caso within-factor, alzando il rand_mean e comprimendo lo z-score.

### L'idea
Ora che l'EFA identifica i fattori, possiamo campionare le coppie random SOLO tra item di fattori DIVERSI. Queste coppie sono garantite non informative.

### Perche' dovrebbe funzionare
- La baseline attuale include ~1/nF delle coppie come accidentalmente within-factor
- Con nF=5, il ~20% delle coppie random sono within-factor -> rand_mean inflazionato
- Eliminandole, rand_mean scende e rand_sd si restringe -> z-score piu' separati
- L'effetto e' proporzionalmente maggiore con pochi fattori (dove serve di piu')

### Rischi
- Con molti fattori (nF=30), l'effetto e' trascurabile (<3% coppie accidentalmente within-factor)
- Se l'EFA sbaglia l'assegnazione, le coppie "cross-factor" potrebbero includere paia effettivamente within-factor

---

## 4. EFA iterativa: pulisci e ri-stima (priorita' media)

### Il problema
L'EFA viene stimata sull'intero campione, inclusi i rispondenti careless. I careless aggiungono rumore alla matrice di correlazione, degradando la soluzione fattoriale e quindi la selezione delle coppie.

### L'idea
1. Primo passaggio: EFA sull'intero campione -> z-score -> flagga i sospetti careless (soglia liberale)
2. Rimuovi i flaggati dal campione
3. Secondo passaggio: ri-esegui EFA sul campione "pulito" -> nuova selezione coppie -> ri-calcola z-score su TUTTI i rispondenti usando la nuova struttura

### Perche' dovrebbe funzionare
- Analogo al principio dell'**M-estimator** in statistica robusta: stima i parametri escludendo i casi anomali, poi ri-valuta tutti con parametri puliti
- La struttura fattoriale senza careless e' piu' accurata -> coppie migliori -> segnale piu' forte
- Particolarmente utile con alta % careless (>25%)

### Rischi
- Se il primo passaggio ha troppi falsi positivi, il campione "pulito" e' biased
- Raddoppia il tempo di calcolo
- Rischio di circolarita' se si itera troppo
- Con pct_careless basso (<10%), il guadagno e' minimo

### Varianti
- **Single iteration**: un solo ciclo clean-refit (piu' sicuro)
- **Soft weighting**: pesare i rispondenti per 1-p(careless) nella stima dell'EFA invece di rimuoverli

---

## 5. Ensemble con altri metodi complementari (priorita' media, approccio diverso)

### Il problema
Il ReReReRe rileva solo **inconsistent carelessness** (random, mixed). Non rileva **consistent carelessness** (longstring, acquiescenza). Sono costrutti complementari, non alternativi. Lo abbiamo visto su PISA 2018: i careless acquiescenti avevano z-score ALTO nel ReReReRe.

### L'idea
Combinare ReReReRe con indicatori complementari:
- **IRV** (Intra-individual Response Variability): rileva low variability
- **LongString**: rileva sequenze identiche consecutive
- **Person-total correlation**: rileva pattern anomali rispetto al campione

### Strategie di combinazione
1. **OR logico**: flag se QUALUNQUE metodo flagga -> massima sensibilita'
2. **Score combinato**: media pesata degli z-score -> piu' smooth
3. **Two-step**: prima longstring/IRV (rapido), poi ReReReRe sui rimanenti

### Perche' dovrebbe funzionare
I metodi sono complementari per costruzione. L'unione dei detection set copre piu' careless respondents.

### Rischi
- Aumenta i falsi positivi (union = piu' flag)
- Perde la semplicita' di "una funzione, un risultato"
- Meglio come feature aggiuntiva / future work che come cambiamento al metodo core
