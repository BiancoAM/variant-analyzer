# Variant Analyzer - Progetto Completo

## Panoramica

Ho creato un **sistema completo di analisi per varianti genetiche** con le seguenti caratteristiche:

### ✅ Funzionalità Implementate

1. **Analisi Multi-Tipo di Varianti**
   - SNV/Missense
   - Varianti di splicing
   - Indel/Frameshift

2. **Predittori In Silico Integrati**
   - **Missense**: SIFT, PolyPhen-2, CADD, REVEL
   - **Splicing**: SpliceAI, MaxEntScan, SPIDEX, dbscSNV
   - **Indel**: Analisi frameshift/in-frame, CADD-Indel

3. **Integrazione Letteratura Scientifica**
   - Ricerca automatica su PubMed
   - Classificazione articoli (studi funzionali, case report, review)
   - Collegamento diretto agli articoli

4. **Classificazione ACMG/AMP**
   - Criteri ACMG automatici basati su predizioni
   - Interpretazione clinica
   - Raccomandazioni per analisi aggiuntive

5. **Collegamenti Database Pubblici**
   - gnomAD (frequenze popolazione)
   - ClinVar (significato clinico)
   - OMIM, GeneCards, GTEx
   - Ensembl, UCSC Genome Browser

6. **Interfaccia Web Completa**
   - Form intuitivo per input varianti
   - Visualizzazione risultati organizzata
   - Esempi precaricati
   - Design responsive

7. **API REST Completa**
   - Analisi singola variante
   - Analisi batch
   - Ricerca letteratura
   - Documentazione automatica (Swagger/OpenAPI)

## Struttura del Progetto

```
variant-analyzer/
├── backend/
│   ├── api/
│   │   ├── main.py              # API REST endpoints (FastAPI)
│   │   └── __init__.py
│   ├── models/
│   │   ├── variant.py           # Database models (SQLAlchemy)
│   │   └── __init__.py
│   ├── services/
│   │   ├── variant_service.py   # Orchestrazione analisi
│   │   ├── pubmed_service.py    # Integrazione PubMed
│   │   └── __init__.py
│   └── predictors/
│       ├── base_predictor.py         # Classe base predittori
│       ├── missense_predictors.py    # SIFT, PolyPhen, CADD, REVEL
│       ├── splicing_predictors.py    # SpliceAI, MaxEntScan, etc.
│       ├── indel_predictors.py       # Analisi indel
│       └── __init__.py
├── frontend/
│   ├── index.html               # Interfaccia web
│   ├── style.css                # Styling
│   └── app.js                   # JavaScript per API calls
├── database/                    # Database migrations (future)
├── tests/                       # Test suite (future)
├── requirements.txt             # Dipendenze Python
├── .env.example                 # Template configurazione
├── README.md                    # Documentazione principale
├── USAGE.md                     # Guida all'uso dettagliata
├── PROGETTO_COMPLETO.md         # Questo file
└── run.sh                       # Script avvio rapido

```

## Quick Start

### 1. Installazione

```bash
cd variant-analyzer

# Copia configurazione
cp .env.example .env

# Modifica .env con la tua email e API key NCBI (opzionale)
nano .env

# Installa dipendenze
pip install -r requirements.txt
```

### 2. Avvio

**Metodo 1: Script automatico**
```bash
./run.sh
```

**Metodo 2: Manuale**
```bash
cd backend
python -m api.main
```

### 3. Accesso

- **Frontend**: http://localhost:8000
- **API Docs**: http://localhost:8000/docs
- **Health Check**: http://localhost:8000/health

## Esempi di Uso

### Esempio 1: Variante Missense BRCA1

Input:
- Gene: `BRCA1`
- Cromosoma: `17`
- Posizione: `43045677`
- Riferimento: `G`
- Alternativo: `A`
- HGVS protein: `p.Cys61Gly`

Output:
- Predizioni da SIFT, PolyPhen-2, CADD, REVEL
- Consenso patogenicità
- Articoli PubMed rilevanti
- Classificazione ACMG
- Link a gnomAD, ClinVar, etc.

### Esempio 2: Frameshift BRCA1

Input:
- Gene: `BRCA1`
- Cromosoma: `17`
- Posizione: `43044295`
- Riferimento: `G`
- Alternativo: `GA`
- HGVS protein: `p.Gln1756fs`

Output:
- Classificazione frameshift (PVS1)
- Predizione loss-of-function
- Letteratura su frameshift BRCA1
- Classificazione "Likely Pathogenic"

### Esempio 3: Variante Splicing CFTR

Input:
- Gene: `CFTR`
- Cromosoma: `7`
- Posizione: `117559590`
- Riferimento: `G`
- Alternativo: `A`
- HGVS coding: `c.1585-1G>A`

Output:
- Predizioni SpliceAI, MaxEntScan
- Impatto su splicing
- Studi funzionali su CFTR
- Raccomandazioni per RNA studies

## Tool di Predizione - Collegamenti

### Per Varianti Missense

| Tool | Descrizione | Score | URL |
|------|-------------|-------|-----|
| **SIFT** | Predice effetto sostituzioni aminoacidiche | 0-1 (≤0.05 = dannosa) | https://sift.bii.a-star.edu.sg/ |
| **PolyPhen-2** | Impatto su struttura/funzione proteica | 0-1 (>0.85 = dannosa) | http://genetics.bwh.harvard.edu/pph2/ |
| **CADD** | Score integrato di deleteriousness | PHRED (>20 = top 1%) | https://cadd.gs.washington.edu/ |
| **REVEL** | Ensemble per varianti missense | 0-1 (>0.5 = patogenica) | https://sites.google.com/site/revelgenomics/ |

### Per Varianti di Splicing

| Tool | Descrizione | Score | URL |
|------|-------------|-------|-----|
| **SpliceAI** | Deep learning per alterazioni splicing | Δ 0-1 (>0.5 = alta confidenza) | https://github.com/Illumina/SpliceAI |
| **MaxEntScan** | Forza siti di splicing | Score relativo | http://hollywood.mit.edu/burgelab/maxent/ |
| **SPIDEX** | Predizione impatto splicing | Z-score (\|Z\|>2 = significativo) | Database |
| **dbscSNV** | Database SNV splice-altering | ada/rf score | Database |

### Database Popolazione e Clinici

| Database | Descrizione | URL |
|----------|-------------|-----|
| **gnomAD** | Frequenze alleliche (125k+ esomi) | https://gnomad.broadinstitute.org/ |
| **ClinVar** | Interpretazioni cliniche varianti | https://www.ncbi.nlm.nih.gov/clinvar/ |
| **COSMIC** | Mutazioni somatiche in cancro | https://cancer.sanger.ac.uk/cosmic |
| **HGMD** | Database mutazioni disease-causing | http://www.hgmd.cf.ac.uk/ |

### Risorse Geni

| Risorsa | Descrizione | URL |
|---------|-------------|-----|
| **OMIM** | Fenotipi e geni mendeliani | https://www.omim.org/ |
| **GeneCards** | Database integrato geni umani | https://www.genecards.org/ |
| **GTEx** | Espressione genica per tessuto | https://gtexportal.org/ |

## Caratteristiche Tecniche

### Backend (Python)
- **Framework**: FastAPI (async, alta performance)
- **ORM**: SQLAlchemy (per database futuro)
- **API PubMed**: BioPython Entrez
- **Validazione**: Pydantic
- **Caching**: In-memory con TTL

### Frontend (Vanilla JS)
- **No framework** - vanilla JavaScript ES6+
- **Responsive design** - CSS Grid/Flexbox
- **Async/await** - per chiamate API
- **Design moderno** - gradients, shadows, transitions

### API
- **REST completa** - endpoints per tutte le funzionalità
- **OpenAPI/Swagger** - documentazione automatica
- **CORS enabled** - per development
- **Error handling** - gestione errori robusta

## Criteri ACMG/AMP Implementati

### Evidenze Patogeniche
- **PVS1**: Variante null (frameshift, nonsense, splice)
- **PS3**: Studi funzionali supportano patogenicità
- **PM2**: Assente da database popolazione
- **PM4**: Indel in-frame in regione non-repeat
- **PP3**: Predittori computazionali supportano patogenicità

### Evidenze Benigne
- **BA1**: Alta frequenza allele (>5%)
- **BS1**: Frequenza allele elevata (>1%)
- **BP4**: Predittori computazionali suggeriscono benignità

### Combinazioni per Classificazione
- **Pathogenic**: PVS1 + (PS o 2×PM) o 2×PS + 1×PM
- **Likely Pathogenic**: PS + 1-2×PM o 1×PM + 2×PP
- **VUS**: Evidenza insufficiente
- **Likely Benign**: 1×BS + 1×BP
- **Benign**: BA1 stand-alone o 2×BS

## Sviluppi Futuri Consigliati

### Priorità Alta
1. **Integrazione gnomAD** - Download database locale o API access
2. **ClinVar API** - Recupero automatico significato clinico
3. **VEP/SnpEff** - Annotazione automatica varianti
4. **SpliceAI locale** - Installazione modello deep learning

### Priorità Media
5. **Database PostgreSQL** - Persistenza risultati
6. **Sistema autenticazione** - Utenti multipli
7. **Export PDF** - Report stampabili
8. **File VCF upload** - Analisi batch da file

### Priorità Bassa
9. **Visualizzazioni** - Grafici interattivi
10. **Integrazione IGV.js** - Viewer genomico
11. **API rate limiting** - Protezione sovraccarico
12. **Containerization** - Docker/Kubernetes

## Note Importanti

### ⚠️ Disclaimer
Questo software è un **tool di ricerca** e NON deve essere usato come unico strumento per diagnosi clinica. Sempre:
- Verificare risultati con genetisti clinici
- Consultare database clinici aggiornati
- Considerare segregazione familiare
- Valutare dati fenotipici del paziente

### 🔧 Configurazione NCBI
Per usare PubMed è **obbligatorio** configurare:
- `NCBI_EMAIL` nel file `.env`
- (Opzionale) `NCBI_API_KEY` per rate limits più alti

Ottieni API key: https://www.ncbi.nlm.nih.gov/account/

### 📊 Limiti Attuali
1. Alcuni predittori richiedono installazione locale
2. Database popolazione non ancora integrati (placeholder)
3. Rate limiting APIs esterne da gestire
4. Nessuna persistenza database (tutto in-memory)

## Supporto e Contributi

Per domande, bug, o nuove features:
1. Leggi `USAGE.md` per guide dettagliate
2. Consulta API docs su `/docs`
3. Apri issue su GitHub (se repository pubblico)

## Licenza e Citazioni

### Tool di Terze Parti
Quando usi questo software, cita i tool originali:
- **SIFT**: Kumar et al., Nat Protoc 2009
- **PolyPhen-2**: Adzhubei et al., Nat Methods 2010
- **CADD**: Kircher et al., Nat Genet 2014
- **SpliceAI**: Jaganathan et al., Cell 2019
- **ACMG Guidelines**: Richards et al., Genet Med 2015

### Questo Software
Variant Analyzer v1.0.0 - Sistema di Analisi Varianti Genetiche
Sviluppato con FastAPI, Python, e tecnologie open source.

---

**Buona analisi! 🧬**
