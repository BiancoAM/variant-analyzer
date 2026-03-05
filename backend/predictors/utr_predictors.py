from typing import Dict, Any, Optional
from .base_predictor import BasePredictor
import logging
import re

logger = logging.getLogger(__name__)


class RegulomeDBPredictor(BasePredictor):
    """
    RegulomeDB - Regulatory element database
    Scores variants in non-coding regions (UTRs, promoters, enhancers)
    based on functional genomic data.

    Score scale:
    - 1a: eQTL + TF binding + matched TF motif + matched DNase Footprint + DNase peak
    - 1b: eQTL + TF binding + any motif + DNase Footprint + DNase peak
    - 1c: eQTL + TF binding + DNase peak
    - 1d: eQTL + TF binding
    - 1e: eQTL + any motif + DNase peak
    - 1f: eQTL only
    - 2a: TF binding + matched TF motif + matched DNase Footprint + DNase peak
    - 2b: TF binding + any motif + DNase Footprint + DNase peak
    - 2c: TF binding + DNase peak
    - 3a: TF binding + any motif + DNase peak
    - 3b: TF binding + any motif
    - 4: TF binding
    - 5: DNase peak only
    - 6: Some data (protein binding microarray, chromatin state, etc.)
    - 7: No data

    Scores 1-2: Likely to affect binding (functional)
    Scores 3-4: Possible regulatory function
    Scores 5-7: Unlikely to affect binding

    API: https://regulomedb.org/regulome-search/
    """

    def __init__(self, api_url: Optional[str] = None):
        super().__init__(api_url or "https://regulomedb.org")

    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        cache_key = self._cache_key(**variant)
        cached = self._get_cached(cache_key)
        if cached:
            return cached

        chromosome = variant.get('chromosome', '').replace('chr', '')
        position = variant.get('position', '')

        result = {
            "tool": "RegulomeDB",
            "score": None,
            "score_category": None,
            "num_experiments": None,
            "interpretation": None,
            "url": self._generate_regulomedb_url(chromosome, position),
            "note": None
        }

        if chromosome and position:
            try:
                chr_formatted = f"chr{chromosome}"
                api_url = (
                    f"{self.api_url}/regulome-search/"
                    f"?regions={chr_formatted}:{position}-{position}"
                    f"&genome=GRCh38&format=json"
                )

                response = await self._make_request(api_url)

                if response and isinstance(response, dict):
                    # Parse RegulomeDB response
                    regions = response.get('regions', [])
                    if regions:
                        region_data = regions[0]
                        score_data = region_data.get('regulome_score', {})
                        raw_score = score_data.get('ranking', score_data.get('score'))
                        num_exp = score_data.get('num_experiments', region_data.get('features_count'))

                        result["score"] = str(raw_score) if raw_score else None
                        result["num_experiments"] = num_exp
                        result["score_category"] = self._interpret_score(str(raw_score) if raw_score else None)
                        result["interpretation"] = self._get_interpretation(str(raw_score) if raw_score else None)
                        result["note"] = f"Dati RegulomeDB v2 (GRCh38)"
                    else:
                        result["note"] = "Nessun dato regolatorio trovato per questa posizione"
                        result["interpretation"] = "Posizione non annotata in RegulomeDB"

            except Exception as e:
                logger.warning(f"RegulomeDB API call failed: {str(e)}")
                result["note"] = f"API non disponibile: {str(e)}"
                result["interpretation"] = "Verificare manualmente su RegulomeDB"

        logger.info(f"RegulomeDB predictor for variant: {self.format_variant_string(variant)}")

        self._set_cache(cache_key, result)
        return result

    def _generate_regulomedb_url(self, chromosome: str, position: Any) -> str:
        """Generate RegulomeDB direct link"""
        if chromosome and position:
            return (
                f"https://regulomedb.org/regulome-search/"
                f"?regions=chr{chromosome}:{position}-{position}&genome=GRCh38"
            )
        return "https://regulomedb.org/"

    def _interpret_score(self, score: Optional[str]) -> str:
        """Interpret RegulomeDB score category"""
        if not score:
            return "No data"
        if score.startswith('1'):
            return "Alta probabilità funzionale (eQTL + TF binding)"
        elif score.startswith('2'):
            return "Probabile funzione regolatoria (TF binding + DNase)"
        elif score.startswith('3') or score.startswith('4'):
            return "Possibile funzione regolatoria (TF binding)"
        elif score.startswith('5') or score.startswith('6'):
            return "Bassa probabilità funzionale"
        elif score == '7':
            return "Nessun dato funzionale"
        return f"Score: {score}"

    def _get_interpretation(self, score: Optional[str]) -> str:
        """Get clinical interpretation of RegulomeDB score"""
        if not score:
            return "Nessun dato disponibile in RegulomeDB"
        if score.startswith('1'):
            return (
                f"Score {score}: Alta evidenza di funzione regolatoria. "
                "Variante probabilmente funzionale - considerare studi funzionali."
            )
        elif score.startswith('2'):
            return (
                f"Score {score}: Buona evidenza di funzione regolatoria. "
                "Monitorare per potenziale impatto su espressione genica."
            )
        elif score.startswith('3') or score.startswith('4'):
            return (
                f"Score {score}: Evidenza moderata di funzione regolatoria. "
                "Possibile impatto su binding TF - richiede validazione."
            )
        elif score in ('5', '6', '7'):
            return f"Score {score}: Bassa evidenza funzionale - probabilmente benigna."
        return f"Score {score}"


class UTR5Predictor(BasePredictor):
    """
    5'UTR Variant Analyzer
    Analyzes variants in the 5' untranslated region for:
    - Upstream ORF (uORF) creation or disruption
    - Kozak sequence context changes (translational efficiency)
    - mRNA secondary structure alterations
    - Start codon (AUG) creation upstream of main ORF

    Clinical significance: 5'UTR variants can affect:
    - Translation initiation efficiency
    - mRNA stability
    - Protein expression levels

    Key tools:
    - UTRAnnotator (VEP plugin): https://github.com/ImperialCollegeLondon/UTRAnnotator
    - uORF DB: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4410766/
    """

    def __init__(self):
        super().__init__()

    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        cache_key = self._cache_key(**variant)
        cached = self._get_cached(cache_key)
        if cached:
            return cached

        hgvs_c = variant.get('hgvs_c', '')
        ref = variant.get('reference', '')
        alt = variant.get('alternate', '')

        # Detect potential uORF changes
        uorf_analysis = self._analyze_uorf_impact(hgvs_c, ref, alt)
        kozak_analysis = self._analyze_kozak_context(hgvs_c, ref, alt)

        result = {
            "tool": "UTR5Analyzer",
            "region": "5'UTR",
            "uorf_impact": uorf_analysis,
            "kozak_impact": kozak_analysis,
            "atg_creation": self._check_atg_creation(ref, alt),
            "interpretation": None,
            "url": "https://github.com/ImperialCollegeLondon/UTRAnnotator",
            "utrannotator_url": "https://www.ensembl.org/vep",
            "note": "Analisi in silico - validazione funzionale raccomandata"
        }

        result["interpretation"] = self._get_interpretation(uorf_analysis, kozak_analysis, ref, alt)

        logger.info(f"5'UTR predictor for variant: {self.format_variant_string(variant)}")

        self._set_cache(cache_key, result)
        return result

    def _analyze_uorf_impact(self, hgvs_c: str, ref: str, alt: str) -> Dict[str, Any]:
        """Analyze potential uORF creation/disruption"""
        impact = {
            "potential_uorf_creation": False,
            "potential_uorf_disruption": False,
            "note": ""
        }

        # ATG in alternate allele suggests potential uORF creation
        if 'ATG' in alt.upper() and 'ATG' not in ref.upper():
            impact["potential_uorf_creation"] = True
            impact["note"] = "Possibile creazione uORF (ATG nel allele alternativo)"
        elif 'ATG' in ref.upper() and 'ATG' not in alt.upper():
            impact["potential_uorf_disruption"] = True
            impact["note"] = "Possibile disruzione uORF esistente"

        return impact

    def _analyze_kozak_context(self, hgvs_c: str, ref: str, alt: str) -> Dict[str, Any]:
        """Analyze Kozak sequence context changes"""
        # Kozak consensus: (GCC)GCCRCCAUGG (R = purine)
        kozak_purines = {'A', 'G'}
        return {
            "potential_kozak_change": len(ref) == 1 and len(alt) == 1,
            "note": "Verificare contesto Kozak con sequenza +/-10bp dal sito di start"
        }

    def _check_atg_creation(self, ref: str, alt: str) -> bool:
        """Check if variant creates new ATG (potential uORF start)"""
        ref3 = ref.upper()
        alt3 = alt.upper()
        return 'ATG' in alt3 and 'ATG' not in ref3

    def _get_interpretation(
        self,
        uorf_analysis: Dict,
        kozak_analysis: Dict,
        ref: str,
        alt: str
    ) -> str:
        """Generate clinical interpretation for 5'UTR variant"""
        if uorf_analysis.get("potential_uorf_creation"):
            return (
                "Possibile creazione di nuovo uORF: può ridurre traduzione proteina principale. "
                "Studi funzionali (luciferasi reporter assay) raccomandati. "
                "Potenziale evidenza patogenica se confermata."
            )
        elif uorf_analysis.get("potential_uorf_disruption"):
            return (
                "Possibile disruzione uORF esistente: può alterare efficienza di traduzione. "
                "Impatto variabile - richiede analisi nel contesto del gene specifico."
            )
        else:
            return (
                "Variante 5'UTR: valutare impatto su struttura mRNA e traduzione. "
                "Utilizzare UTRAnnotator (plugin VEP) per annotazione completa."
            )


class UTR3Predictor(BasePredictor):
    """
    3'UTR Variant Analyzer
    Analyzes variants in the 3' untranslated region for:
    - miRNA binding site disruption/creation (seed region 2-8 nt)
    - Polyadenylation signal changes (AATAAA hexamer)
    - RNA-binding protein (RBP) motif alterations
    - mRNA stability elements (AU-rich elements, ARE)

    Clinical significance: 3'UTR variants can affect:
    - Post-transcriptional regulation
    - mRNA stability and half-life
    - Tissue-specific expression
    - Translation efficiency

    Key tools:
    - TargetScan: https://www.targetscan.org/
    - miRDB: https://mirdb.org/
    - RBPmap: http://rbpmap.technion.ac.il/
    - PolyASite: https://polyasite.unibas.ch/
    """

    def __init__(self):
        super().__init__()

    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        cache_key = self._cache_key(**variant)
        cached = self._get_cached(cache_key)
        if cached:
            return cached

        ref = variant.get('reference', '')
        alt = variant.get('alternate', '')
        gene = variant.get('gene', '')

        # Analyze potential miRNA site disruption
        mirna_analysis = self._analyze_mirna_impact(ref, alt)

        # Analyze polyadenylation signal
        polya_analysis = self._analyze_polya_signal(ref, alt)

        # Analyze ARE (AU-rich elements)
        are_analysis = self._analyze_are(ref, alt)

        result = {
            "tool": "UTR3Analyzer",
            "region": "3'UTR",
            "mirna_impact": mirna_analysis,
            "polyadenylation_impact": polya_analysis,
            "are_impact": are_analysis,
            "interpretation": None,
            "targetscan_url": f"https://www.targetscan.org/cgi-bin/targetscan/vert_80/targetscan.cgi?species=Human&gid={gene}" if gene else "https://www.targetscan.org/",
            "mirdb_url": f"https://mirdb.org/cgi-bin/search.cgi?searchType=miRNA&query={gene}" if gene else "https://mirdb.org/",
            "rbpmap_url": "http://rbpmap.technion.ac.il/",
            "polyasite_url": "https://polyasite.unibas.ch/",
            "note": "Analisi in silico - verificare con TargetScan, miRDB, RBPmap"
        }

        result["interpretation"] = self._get_interpretation(
            mirna_analysis, polya_analysis, are_analysis
        )

        logger.info(f"3'UTR predictor for variant: {self.format_variant_string(variant)}")

        self._set_cache(cache_key, result)
        return result

    def _analyze_mirna_impact(self, ref: str, alt: str) -> Dict[str, Any]:
        """Analyze potential miRNA binding site disruption"""
        return {
            "potential_site_disruption": len(ref) == 1 and len(alt) == 1,
            "potential_site_creation": len(ref) == 1 and len(alt) == 1,
            "note": (
                "Verificare con TargetScan e miRDB se la variante cade in un seed region "
                "di miRNA (posizioni 2-8 dal 5' del miRNA). "
                "Disruzione seed region: impatto significativo su regolazione post-trascrizionale."
            )
        }

    def _analyze_polya_signal(self, ref: str, alt: str) -> Dict[str, Any]:
        """Analyze polyadenylation signal (AATAAA) changes"""
        polya_signal = "AATAAA"
        polya_variants = ["ATTAAA", "AGTAAA", "TATAAA", "AATATA", "AATACA", "AATAGA"]

        ref_upper = ref.upper()
        alt_upper = alt.upper()

        # Check if variant is within known polyadenylation signal
        if ref_upper in polya_signal or alt_upper in polya_variants:
            return {
                "potential_polya_disruption": True,
                "note": "Possibile alterazione del segnale di poliadenilazione (AATAAA). Verificare con PolyASite."
            }

        return {
            "potential_polya_disruption": False,
            "note": "Verificare manualmente con PolyASite se la variante è vicina al segnale pA"
        }

    def _analyze_are(self, ref: str, alt: str) -> Dict[str, Any]:
        """Analyze AU-rich elements (ARE) - mRNA stability"""
        # ARE pentamer: AUUUA; nonamers: UUAUUUAUU
        return {
            "are_analysis_needed": True,
            "note": (
                "Valutare presenza di elementi AU-rich (ARE: ATTTA pentamer). "
                "Disruzione ARE può aumentare stabilità mRNA (effetto gain-of-function)."
            )
        }

    def _get_interpretation(
        self,
        mirna_analysis: Dict,
        polya_analysis: Dict,
        are_analysis: Dict
    ) -> str:
        """Generate clinical interpretation for 3'UTR variant"""
        findings = []

        if polya_analysis.get("potential_polya_disruption"):
            findings.append(
                "ATTENZIONE: Possibile alterazione segnale poliadenilazione - "
                "può alterare stabilità mRNA e livelli espressione."
            )

        findings.append(
            "Variante 3'UTR: verificare con TargetScan/miRDB se la posizione "
            "coincide con seed region di miRNA noti."
        )
        findings.append(
            "Valutare binding proteine RBP con RBPmap."
        )

        return " | ".join(findings) if findings else (
            "Variante 3'UTR: impatto su regolazione post-trascrizionale da valutare."
        )


class UTRPredictorAggregator:
    """
    Aggregates predictions for UTR variants (5'UTR and 3'UTR).
    Tools: RegulomeDB, UTR5Predictor/UTR3Predictor
    ACMG criteria relevant: BP4, PP3 (limitata evidenza per varianti UTR)
    """

    def __init__(self):
        self.regulomedb = RegulomeDBPredictor()
        self.utr5_predictor = UTR5Predictor()
        self.utr3_predictor = UTR3Predictor()

    async def predict_all(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        """Run all UTR predictors based on UTR type"""
        results = {}

        utr_type = variant.get('variant_type', '')

        # Always run RegulomeDB
        try:
            results["RegulomeDB"] = await self.regulomedb.predict(variant)
        except Exception as e:
            logger.error(f"Error in RegulomeDB predictor: {str(e)}")
            results["RegulomeDB"] = {"error": str(e), "tool": "RegulomeDB"}

        # Run type-specific predictor
        if utr_type == 'UTR_5':
            try:
                results["UTR5Analyzer"] = await self.utr5_predictor.predict(variant)
            except Exception as e:
                logger.error(f"Error in UTR5 predictor: {str(e)}")
                results["UTR5Analyzer"] = {"error": str(e), "tool": "UTR5Analyzer"}
        elif utr_type == 'UTR_3':
            try:
                results["UTR3Analyzer"] = await self.utr3_predictor.predict(variant)
            except Exception as e:
                logger.error(f"Error in UTR3 predictor: {str(e)}")
                results["UTR3Analyzer"] = {"error": str(e), "tool": "UTR3Analyzer"}

        results["aggregate"] = self._aggregate_predictions(results, utr_type)
        return results

    def _aggregate_predictions(self, results: Dict[str, Any], utr_type: str) -> Dict[str, Any]:
        """Provide overall interpretation for UTR variant"""
        interpretations = []
        acmg_criteria = []

        # RegulomeDB assessment
        regulome = results.get("RegulomeDB", {})
        regulome_score = regulome.get("score")

        if regulome_score:
            if regulome_score.startswith('1') or regulome_score.startswith('2'):
                interpretations.append(
                    f"RegulomeDB score {regulome_score}: alta evidenza di funzione regolatoria"
                )
                acmg_criteria.append("PP3 (supportante)")
            elif regulome_score in ('5', '6', '7'):
                interpretations.append(
                    f"RegulomeDB score {regulome_score}: bassa evidenza funzionale"
                )
                acmg_criteria.append("BP4 (supportante)")

        # UTR-specific findings
        if utr_type == 'UTR_5':
            utr5 = results.get("UTR5Analyzer", {})
            uorf = utr5.get("uorf_impact", {})
            if uorf.get("potential_uorf_creation"):
                interpretations.append("Possibile creazione uORF: può ridurre espressione proteica")
            elif uorf.get("potential_uorf_disruption"):
                interpretations.append("Possibile disruzione uORF: effetto variabile su traduzione")

        elif utr_type == 'UTR_3':
            utr3 = results.get("UTR3Analyzer", {})
            polya = utr3.get("polyadenylation_impact", {})
            if polya.get("potential_polya_disruption"):
                interpretations.append("Possibile disruzione segnale poliadenilazione - priorità alta")
                acmg_criteria.append("PP3 (supportante - se confermato)")

        if not interpretations:
            interpretations.append(
                "Variante UTR: analisi regolatoria manuale raccomandata"
            )

        region_label = "5'UTR" if utr_type == "UTR_5" else "3'UTR"

        return {
            "consensus": self._get_consensus(acmg_criteria, utr_type),
            "utr_region": region_label,
            "acmg_criteria": acmg_criteria,
            "details": interpretations,
            "recommendation": self._get_recommendation(utr_type)
        }

    def _get_consensus(self, acmg_criteria: list, utr_type: str) -> str:
        has_pp3 = any("PP3" in c for c in acmg_criteria)
        has_bp4 = any("BP4" in c for c in acmg_criteria)

        if has_pp3:
            return "Possibile impatto funzionale - studi di espressione raccomandati"
        elif has_bp4:
            return "Bassa evidenza di impatto funzionale"
        return "Variante UTR - interpretazione richiede dati funzionali"

    def _get_recommendation(self, utr_type: str) -> str:
        if utr_type == 'UTR_5':
            return (
                "5'UTR: "
                "1. Analizzare con UTRAnnotator (VEP plugin) per uORF e Kozak. "
                "2. Verificare con RegulomeDB per evidenza regolatoria. "
                "3. Considerare reporter assay se uORF creation/disruption identificato. "
                "4. Analisi RNA-seq per impatto su livelli espressione."
            )
        else:
            return (
                "3'UTR: "
                "1. Verificare seed region miRNA con TargetScan e miRDB. "
                "2. Analizzare binding RBP con RBPmap. "
                "3. Controllare segnali di poliadenilazione con PolyASite. "
                "4. Western blot / qPCR per valutare impatto su espressione."
            )
