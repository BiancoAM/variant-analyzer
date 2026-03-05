# Variant Analyzer

Sistema completo per l'analisi di varianti genetiche con integrazione di tool di predizione in silico, PubMed e studi funzionali.

## Features

- **Analisi Multi-Tipo**: Supporto per SNV/Missense, Splicing, Indel/Frameshift
- **Predittori In Silico**: Integrazione con SIFT, PolyPhen-2, CADD, REVEL, SpliceAI, MaxEntScan
- **Letteratura Scientifica**: Ricerca automatica su PubMed per varianti specifiche
- **Studi Funzionali**: Collegamenti a evidenze funzionali di patogenicità
- **Database Pubblici**: Integrazione con ClinVar, gnomAD, COSMIC

## Tipi di Varianti Supportate

### SNV/Missense
- **SIFT**: Predice se una sostituzione aminoacidica influenza la funzione proteica
- **PolyPhen-2**: Predice l'impatto di varianti missense
- **REVEL**: Score ensemble per varianti missense
- **CADD**: Score di deleteriousness generale

### Splicing
- **SpliceAI**: Deep learning per predizione alterazioni splicing
- **MaxEntScan**: Predice forza siti di splicing
- **SPIDEX**: Percentile rank per impatto splicing

### Indel/Frameshift
- **CADD-Indel**: Score CADD per inserzioni/delezioni
- **Analisi frame**: Determina se in-frame o frameshift

## Struttura Progetto

```
variant-analyzer/
├── backend/
│   ├── api/              # REST API endpoints
│   ├── models/           # Database models
│   ├── services/         # Business logic
│   └── predictors/       # Integration con prediction tools
├── frontend/             # Web interface
├── database/             # Database migrations
└── tests/               # Test suite
```

## Setup

1. Installa le dipendenze:
```bash
pip install -r requirements.txt
```

2. Configura le variabili d'ambiente:
```bash
cp .env.example .env
# Modifica .env con le tue credenziali
```

3. Inizializza il database:
```bash
cd database
alembic upgrade head
```

4. Avvia il server:
```bash
cd backend
uvicorn api.main:app --reload
```

## API Keys Richieste

- **NCBI API Key**: Per accesso a PubMed (ottieni su https://www.ncbi.nlm.nih.gov/account/)
- Alcuni predittori richiedono registrazione/API keys separate

## Uso

L'applicazione web sarà disponibile su `http://localhost:8000`

API Documentation: `http://localhost:8000/docs`
