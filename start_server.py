#!/usr/bin/env python3
"""
Script per avviare il server Variant Analyzer
"""
import sys
import os

# Aggiungi la directory root al PYTHONPATH
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("backend.api.main:app", host="0.0.0.0", port=8000, reload=True)
