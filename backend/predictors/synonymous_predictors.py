from typing import Dict, Any, Optional
from .base_predictor import BasePredictor
import logging

logger = logging.getLogger(__name__)


class MobiDetailsPredictor(BasePredictor):
    """
    MobiDetails - Comprehensive variant annotation and classification tool
    Provides ACMG/AMP classification, links to external databases, NMD prediction.
    Works best with HGVS nomenclature (NM accession + HGVS coding).
    API: https://mobidetails.iurc.montp.inserm.fr/MD/api/
    """

    def __init__(self, api_url: Optional[str] = None):
        super().__init__(api_url or "https://mobidetails.iurc.montp.inserm.fr/MD/api")

    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        cache_key = self._cache_key(**variant)
        cached = self._get_cached(cache_key)
        if cached:
            return cached

        transcript = variant.get('transcript', '')
        hgvs_c = variant.get('hgvs_c', '')
        gene = variant.get('gene', '')

        # Generate direct MobiDetails link
        mobidetails_url = self._generate_mobidetails_url(transcript, hgvs_c, gene)

        result = {
            "tool": "MobiDetails",
            "url": mobidetails_url,
            "api_url": None,
            "classification": None,
            "interpretation": None,
            "note": "Inserisci transcript NM + HGVS per annotazione completa"
        }

        # If transcript and hgvs_c are provided, try to call the API
        if transcript and hgvs_c and transcript.startswith('NM_'):
            try:
                import urllib.parse
                hgvs_encoded = urllib.parse.quote(hgvs_c, safe='')
                api_endpoint = f"{self.api_url}/variant/hgvs/{transcript}/{hgvs_encoded}/"

                response = await self._make_request(api_endpoint)

                if response and isinstance(response, dict):
                    # Extract key fields from MobiDetails response
                    result["api_url"] = api_endpoint
                    result["classification"] = response.get('acmg_classification', {}).get('acmg_class')
                    result["spliceai_score"] = response.get('spliceai', {}).get('max_score')
                    result["mpa_score"] = response.get('mpa_score')
                    result["interpretation"] = response.get('acmg_classification', {}).get('acmg_description')
                    result["note"] = "Dati recuperati da MobiDetails API"

            except Exception as e:
                logger.warning(f"MobiDetails API call failed: {str(e)} - using link only")

        logger.info(f"MobiDetails predictor for variant: {self.format_variant_string(variant)}")

        self._set_cache(cache_key, result)
        return result

    def _generate_mobidetails_url(self, transcript: str, hgvs_c: str, gene: str) -> str:
        """Generate MobiDetails direct link"""
        if transcript and hgvs_c:
            import urllib.parse
            query = urllib.parse.quote(f"{transcript}:{hgvs_c}")
            return f"https://mobidetails.iurc.montp.inserm.fr/MD/search/?query={query}"
        elif gene:
            import urllib.parse
            return f"https://mobidetails.iurc.montp.inserm.fr/MD/search/?query={urllib.parse.quote(gene)}"
        return "https://mobidetails.iurc.montp.inserm.fr/MD/"


class SILVAPredictor(BasePredictor):
    """
    SILVA - Synonymous variant Impact analysis
    Analyzes the impact of synonymous (silent) variants on:
    - Exonic Splicing Enhancers (ESE) and Silencers (ESS)
    - Codon usage bias
    - mRNA secondary structure changes
    - Potential cryptic splice site activation
    Reference: https://silva.bioinfo.cipf.es/
    """

    def __init__(self):
        super().__init__()

    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        cache_key = self._cache_key(**variant)
        cached = self._get_cached(cache_key)
        if cached:
            return cached

        hgvs_c = variant.get('hgvs_c', '')
        hgvs_p = variant.get('hgvs_p', '')
        gene = variant.get('gene', '')

        # Verify this is truly synonymous (p.(=) notation)
        is_synonymous = self._is_synonymous(hgvs_p)

        result = {
            "tool": "SILVA",
            "is_synonymous": is_synonymous,
            "ese_impact": None,
            "ess_impact": None,
            "codon_usage_change": None,
            "interpretation": None,
            "url": "https://silva.bioinfo.cipf.es/",
            "note": "Analisi manuale richiesta - inserire sequenza sul sito SILVA"
        }

        if is_synonymous:
            result["interpretation"] = (
                "Variante sinonima: valutare impatto su ESE/ESS e splicing. "
                "Utilizzare SILVA per analisi completa degli elementi regolatori dello splicing."
            )
        else:
            result["interpretation"] = (
                "Variante non sinonima - SILVA ottimizzato per varianti sinonime p.(=)"
            )

        logger.info(f"SILVA predictor for variant: {self.format_variant_string(variant)}")

        self._set_cache(cache_key, result)
        return result

    def _is_synonymous(self, hgvs_p: str) -> bool:
        """Check if variant is truly synonymous from HGVS protein notation"""
        if not hgvs_p:
            return False
        # p.(=) or p.= indicates synonymous
        return '=' in hgvs_p


class ESEfinderPredictor(BasePredictor):
    """
    ESEfinder - Exonic Splicing Enhancer prediction
    Identifies SR protein binding sites (SF2/ASF, SC35, SRp40, SRp55)
    that may be disrupted by synonymous variants.
    Critical for BP4/PP3 criteria for synonymous variants.
    Reference: http://krainer01.cshl.edu/cgi-bin/tools/ESE3/esefinder.cgi
    """

    def __init__(self):
        super().__init__()

    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        cache_key = self._cache_key(**variant)
        cached = self._get_cached(cache_key)
        if cached:
            return cached

        result = {
            "tool": "ESEfinder",
            "sr_proteins_affected": [],
            "ese_lost": None,
            "ess_created": None,
            "hexplorer_score": None,
            "interpretation": None,
            "url": "http://krainer01.cshl.edu/cgi-bin/tools/ESE3/esefinder.cgi",
            "hexplorer_url": "https://hexplorer.uni-due.de/",
            "note": "Analisi manuale richiesta - inserire sequenza flanking su ESEfinder/HEXplorer"
        }

        result["interpretation"] = (
            "Valutare impatto su siti ESE/ESS con ESEfinder e HEXplorer. "
            "Disruzione di ESE: evidenza BP4 (benigno) o PP3 (patogenico) per varianti sinonime."
        )

        logger.info(f"ESEfinder predictor for variant: {self.format_variant_string(variant)}")

        self._set_cache(cache_key, result)
        return result


class SynonymousPredictorAggregator:
    """
    Aggregates predictions for synonymous (silent) variants.
    Tools: MobiDetails, SILVA, ESEfinder
    ACMG criteria relevant: BP7, PP3, PS1
    """

    def __init__(self):
        self.predictors = {
            "MobiDetails": MobiDetailsPredictor(),
            "SILVA": SILVAPredictor(),
            "ESEfinder": ESEfinderPredictor()
        }

    async def predict_all(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        """Run all synonymous predictors and aggregate results"""
        results = {}

        for name, predictor in self.predictors.items():
            try:
                results[name] = await predictor.predict(variant)
            except Exception as e:
                logger.error(f"Error in {name} predictor: {str(e)}")
                results[name] = {"error": str(e), "tool": name}

        results["aggregate"] = self._aggregate_predictions(results)
        return results

    def _aggregate_predictions(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Provide overall interpretation for synonymous variant"""
        interpretations = []
        acmg_criteria = []

        silva = results.get("SILVA", {})
        mobidetails = results.get("MobiDetails", {})

        # Check if truly synonymous
        is_synonymous = silva.get("is_synonymous", False)

        if is_synonymous:
            interpretations.append(
                "Variante sinonima confermata (p.(=)): non modifica sequenza aminoacidica."
            )
            # BP7: synonymous variant with no predicted splice impact
            acmg_criteria.append("BP7 (potenziale - verificare splicing)")
        else:
            interpretations.append(
                "Verifica notazione HGVS proteica per confermare sinonimia."
            )

        # MobiDetails classification
        if mobidetails.get("classification"):
            interpretations.append(
                f"MobiDetails: {mobidetails['classification']}"
            )

        consensus = self._get_consensus(is_synonymous, acmg_criteria)

        return {
            "consensus": consensus,
            "is_synonymous": is_synonymous,
            "acmg_criteria": acmg_criteria,
            "details": interpretations,
            "recommendation": self._get_recommendation(is_synonymous)
        }

    def _get_consensus(self, is_synonymous: bool, acmg_criteria: list) -> str:
        if not is_synonymous:
            return "Verifica annotazione HGVS"
        if "BP7" in str(acmg_criteria):
            return "Probabile variante benigna (sinonima) - verificare effetto splicing"
        return "Variante sinonima - analisi splicing richiesta"

    def _get_recommendation(self, is_synonymous: bool) -> str:
        if is_synonymous:
            return (
                "1. Verificare impatto su splicing con SpliceAI e MaxEntScan. "
                "2. Analizzare ESE/ESS con ESEfinder/HEXplorer. "
                "3. Se vicino a sito di splicing (±10bp): RNA studies raccomandati. "
                "4. Considerare studi RNA se BP7 non applicabile."
            )
        return "Verificare notazione HGVS proteica prima di procedere."
