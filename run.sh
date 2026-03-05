#!/bin/bash

# Script per avviare Variant Analyzer

echo "🧬 Avvio Variant Analyzer..."
echo ""

# Check if .env exists
if [ ! -f .env ]; then
    echo "⚠️  File .env non trovato. Copiando da .env.example..."
    cp .env.example .env
    echo "✅ File .env creato. Modifica il file con le tue credenziali prima di continuare."
    echo ""
fi

# Check if virtual environment exists
if [ ! -d "venv" ]; then
    echo "📦 Creazione virtual environment..."
    python3 -m venv venv
    echo "✅ Virtual environment creato"
    echo ""
fi

# Activate virtual environment
echo "🔧 Attivazione virtual environment..."
source venv/bin/activate

# Install dependencies
echo "📥 Installazione dipendenze..."
pip install -r requirements.txt --quiet

echo ""
echo "✅ Setup completato!"
echo ""
echo "🚀 Avvio server su http://localhost:8000"
echo "   - Frontend: http://localhost:8000"
echo "   - API Docs: http://localhost:8000/docs"
echo ""
echo "Premi Ctrl+C per fermare il server"
echo ""

# Start server
cd backend
python -m api.main
