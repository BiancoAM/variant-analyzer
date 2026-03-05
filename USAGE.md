# Guida all'Uso - Variant Analyzer

## Installazione e Setup

### 1. Prerequisiti

- Python 3.9 o superiore
- pip (Python package manager)
- (Opzionale) PostgreSQL per database persistente

### 2. Installazione

```bash
cd variant-analyzer

# Installa le dipendenze
pip install -r requirements.txt
```

### 3. Configurazione

Crea un file `.env` nella root del progetto:

```bash
cp .env.example .env
```

Modifica `.env` con le tue credenziali:

```env
# Database (opzionale - per persistenza)
DATABASE_URL=postgresql://user:password@localhost:5432/variant_analyzer

# API Keys
NCBI_API_KEY=your_ncbi_api_key_here
NCBI_EMAIL=your_email@example.com

# Impostazioni App
DEBUG=True
LOG_LEVEL=INFO
```

**IMPORTANTE**: Per utilizzare PubMed, devi:
1. Configurare `NCBI_EMAIL` con una email valida
2. (Opzionale ma raccomandato) Ottenere una API key NCBI da https://www.ncbi.nlm.nih.gov/account/

### 4. Avvio del Server

```bash
cd backend
python -m api.main

# Oppure usa uvicorn direttamente:
uvicorn api.main:app --reload --host 0.0.0.0 --port 8000
```

Il server sarà disponibile su:
- **Frontend**: http://localhost:8000
- **API Docs**: http://localhost:8000/docs
- **API ReDoc**: http://localhost:8000/redoc

## Uso dell'Interfaccia Web

### Analisi di una Variante

1. Apri il browser su http://localhost:8000
2. Compila il form con i dati della variante:
   - **Gene**: Simbolo del gene (es. BRCA1)
   - **Cromosoma**: Numero o lettera del cromosoma (es. 17, X)
   - **Posizione**: Posizione genomica (es. 43044295)
   - **Allele Riferimento**: Allele di riferimento (es. G)
   - **Allele Alternativo**: Allele alternativo (es. A)
   - **HGVS Coding** (opzionale): Notazione HGVS coding (es. c.5266dup)
   - **HGVS Protein** (opzionale): Notazione HGVS protein (es. p.Gln1756fs)

3. Clicca su **Analizza Variante**

### Esempi Precaricati

Usa i pulsanti "Esempi rapidi" per caricare varianti di esempio:
- **Variante Missense**: BRCA1 p.Cys61Gly
- **Frameshift**: BRCA1 p.Gln1756fs
- **Splicing**: CFTR c.1585-1G>A

### Interpretazione dei Risultati

L'analisi restituisce:

#### 1. Classificazione ACMG/AMP
- **Pathogenic**: Variante patogenica
- **Likely Pathogenic**: Probabilmente patogenica
- **VUS**: Variante di significato incerto
- **Likely Benign**: Probabilmente benigna
- **Benign**: Benigna

#### 2. Predizioni In Silico

**Per Varianti Missense:**
- **SIFT**: Score ≤0.05 = dannosa
- **PolyPhen-2**: Score >0.85 = probabilmente dannosa
- **CADD**: PHRED >20 = top 1% dannose
- **REVEL**: Score >0.5 = probabilmente patogenica

**Per Varianti di Splicing:**
- **SpliceAI**: Delta >0.5 = alta confidenza alterazione splicing
- **MaxEntScan**: Cambiamento forza sito di splicing
- **SPIDEX**: Z-score >2 = impatto significativo

**Per Indel:**
- **Frameshift**: Tipicamente patogeniche (PVS1)
- **In-frame**: Impatto variabile, dipende da dimensione e posizione

#### 3. Letteratura Scientifica

Articoli da PubMed categorizzati per tipo:
- **Studi Funzionali**: Evidenze sperimentali di patogenicità
- **Case Report**: Casi clinici
- **Review**: Revisioni della letteratura

#### 4. Criteri ACMG

Evidenze secondo le linee guida ACMG/AMP 2015:
- **PVS1**: Variante null in gene LOF
- **PS3**: Studi funzionali supportano patogenicità
- **PM4**: Indel in-frame in regione non-repeat
- **PP3**: Predittori computazionali supportano patogenicità
- **BA1/BS1**: Alta frequenza in popolazione
- Etc.

## Uso dell'API REST

### Analisi Singola Variante

```bash
curl -X POST "http://localhost:8000/api/v1/analyze" \
  -H "Content-Type: application/json" \
  -d '{
    "chromosome": "17",
    "position": 43044295,
    "reference": "G",
    "alternate": "A",
    "gene": "BRCA1",
    "hgvs_c": "c.5266dup",
    "hgvs_p": "p.Gln1756fs"
  }'
```

### Analisi Batch

```bash
curl -X POST "http://localhost:8000/api/v1/analyze/batch" \
  -H "Content-Type: application/json" \
  -d '{
    "variants": [
      {
        "chromosome": "17",
        "position": 43044295,
        "reference": "G",
        "alternate": "A",
        "gene": "BRCA1"
      },
      {
        "chromosome": "7",
        "position": 117559590,
        "reference": "G",
        "alternate": "A",
        "gene": "CFTR"
      }
    ]
  }'
```

### Ricerca Letteratura

```bash
curl -X POST "http://localhost:8000/api/v1/literature/search" \
  -H "Content-Type: application/json" \
  -d '{
    "gene": "BRCA1",
    "variant_description": "p.Cys61Gly",
    "max_results": 20
  }'
```

### Lista Predittori Disponibili

```bash
curl "http://localhost:8000/api/v1/predictors"
```

### Lista Database Integrati

```bash
curl "http://localhost:8000/api/v1/databases"
```

## Tool di Predizione Integrati

### SNV/Missense
- **SIFT**: https://sift.bii.a-star.edu.sg/
- **PolyPhen-2**: http://genetics.bwh.harvard.edu/pph2/
- **CADD**: https://cadd.gs.washington.edu/
- **REVEL**: https://sites.google.com/site/revelgenomics/

### Splicing
- **SpliceAI**: https://github.com/Illumina/SpliceAI
- **MaxEntScan**: http://hollywood.mit.edu/burgelab/maxent/
- **SPIDEX**: Database pre-computato
- **dbscSNV**: Database pre-computato

### Database Popolazione
- **gnomAD**: https://gnomad.broadinstitute.org/
- **ClinVar**: https://www.ncbi.nlm.nih.gov/clinvar/

### Letteratura
- **PubMed**: Ricerca automatica via NCBI Entrez API

## Note Importanti

### Limitazioni Attuali

1. **API Esterne**: Molti predittori richiedono API keys separate o installazione locale
   - Alcuni tool (SpliceAI, SIFT, PolyPhen) sono meglio eseguiti localmente
   - CADD ha limitazioni di rate sulla API pubblica

2. **Database Popolazione**:
   - L'integrazione gnomAD/ClinVar richiede download database o API access
   - Attualmente mostra placeholder - da implementare con dati reali

3. **Studi Funzionali**:
   - La ricerca PubMed è automatica ma richiede interpretazione manuale
   - Non tutti gli articoli indicati sono necessariamente rilevanti

### Best Practices

1. **Sempre verificare** i risultati con database clinici (ClinVar, HGMD)
2. **Non usare** questo tool come unico strumento per diagnosi clinica
3. **Consultare** genetisti clinici per interpretazione finale
4. **Verificare** la qualità della variante chiamata nel file VCF originale
5. **Considerare** segregazione familiare e dati fenotipici

### Performance

- L'analisi di una singola variante richiede 5-15 secondi
- Il batch analysis è più efficiente per multiple varianti
- La ricerca PubMed può essere lenta per geni molto studiati
- Usa il caching per velocizzare analisi ripetute

## Troubleshooting

### Errore: "NCBI API Error"
- Verifica di aver configurato `NCBI_EMAIL` nel file `.env`
- Se hai molte richieste, ottieni una API key NCBI

### Errore: "Connection refused"
- Verifica che il server sia in esecuzione su porta 8000
- Controlla firewall e permessi di rete

### Risultati incompleti
- Alcuni predittori richiedono installazione locale
- Verifica connessione internet per API esterne
- Controlla i log del server per errori specifici

## Sviluppo Futuro

Funzionalità pianificate:
- [ ] Integrazione database gnomAD completa
- [ ] Download e parsing file VCF
- [ ] Annotazione automatica con VEP/SnpEff
- [ ] Integrazione SpliceAI locale
- [ ] Report PDF generabili
- [ ] Sistema di autenticazione utenti
- [ ] Database per storico analisi

## Riferimenti

- **ACMG/AMP Guidelines**: Richards et al. 2015, Genet Med
- **ClinGen**: https://clinicalgenome.org/
- **ACMG**: https://www.acmg.net/

## Supporto

Per domande, bug reports, o feature requests, apri una issue su GitHub.
