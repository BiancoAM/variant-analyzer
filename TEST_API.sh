#!/bin/bash

# Script per testare l'API Variant Analyzer

echo "🧬 Test API Variant Analyzer"
echo "=============================="
echo ""

API_URL="http://localhost:8000"

# Test 1: Health Check
echo "1️⃣  Test Health Check..."
curl -s "$API_URL/health" | python3 -m json.tool
echo ""
echo ""

# Test 2: Lista Predittori
echo "2️⃣  Test Lista Predittori..."
curl -s "$API_URL/api/v1/predictors" | python3 -m json.tool | head -50
echo "..."
echo ""
echo ""

# Test 3: Analisi Variante Missense
echo "3️⃣  Test Analisi Variante Missense (BRCA1 p.Cys61Gly)..."
curl -s -X POST "$API_URL/api/v1/analyze" \
  -H "Content-Type: application/json" \
  -d '{
    "chromosome": "17",
    "position": 43045677,
    "reference": "G",
    "alternate": "A",
    "gene": "BRCA1",
    "hgvs_p": "p.Cys61Gly"
  }' | python3 -m json.tool | head -100
echo "..."
echo ""
echo ""

# Test 4: Analisi Frameshift
echo "4️⃣  Test Analisi Frameshift (BRCA1 p.Gln1756fs)..."
curl -s -X POST "$API_URL/api/v1/analyze" \
  -H "Content-Type: application/json" \
  -d '{
    "chromosome": "17",
    "position": 43044295,
    "reference": "G",
    "alternate": "GA",
    "gene": "BRCA1",
    "hgvs_p": "p.Gln1756fs"
  }' | python3 -c "
import json
import sys
data = json.load(sys.stdin)
print('Tipo Variante:', data['variant_type'])
print('Classificazione:', data['acmg_classification']['classification'])
print('Evidenze Patogeniche:', len(data['acmg_classification']['pathogenic']))
print('Evidenze Benigne:', len(data['acmg_classification']['benign']))
if data['summary']['key_findings']:
    print('\nRisultati Principali:')
    for finding in data['summary']['key_findings']:
        print('  -', finding)
"
echo ""
echo ""

# Test 5: Ricerca Letteratura
echo "5️⃣  Test Ricerca Letteratura (BRCA1)..."
curl -s -X POST "$API_URL/api/v1/literature/search" \
  -H "Content-Type: application/json" \
  -d '{
    "gene": "BRCA1",
    "max_results": 5
  }' | python3 -c "
import json
import sys
data = json.load(sys.stdin)
print('Articoli Totali:', data['total_articles'])
print('Studi Funzionali:', data['categories']['functional_studies']['count'])
print('Case Reports:', data['categories']['case_reports']['count'])
print('Reviews:', data['categories']['reviews']['count'])
print('\nPrimi Articoli:')
for cat in ['functional_studies', 'case_reports']:
    articles = data['categories'][cat]['articles']
    if articles:
        print(f'\n{cat.upper()}:')
        for i, article in enumerate(articles[:2], 1):
            print(f'{i}. {article.get(\"title\", \"N/A\")} (PMID: {article.get(\"pmid\", \"N/A\")})')
"
echo ""
echo ""

# Test 6: Lista Database
echo "6️⃣  Test Lista Database..."
curl -s "$API_URL/api/v1/databases" | python3 -c "
import json
import sys
data = json.load(sys.stdin)
for category, dbs in data.items():
    print(f'\n{category.upper()}:')
    for name, info in dbs.items():
        print(f'  - {name}: {info[\"name\"]}')
"
echo ""
echo ""

echo "✅ Test completati!"
echo ""
echo "Per ulteriori test, apri il browser su:"
echo "  - Frontend: $API_URL"
echo "  - API Docs: $API_URL/docs"
