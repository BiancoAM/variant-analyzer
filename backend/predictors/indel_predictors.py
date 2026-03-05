from typing import Dict, Any
from .base_predictor import BasePredictor
import logging

logger = logging.getLogger(__name__)


class IndelAnalyzer(BasePredictor):
    """
    Analyzer for insertions and deletions
    Determines frame-shift vs in-frame and predicts impact
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

        # Calculate length change
        len_change = len(alt) - len(ref)
        is_insertion = len_change > 0
        is_deletion = len_change < 0

        # Determine if frameshift
        is_frameshift = abs(len_change) % 3 != 0

        result = {
            "tool": "IndelAnalyzer",
            "variant_type": self._classify_indel(is_insertion, is_deletion, is_frameshift),
            "length_change": len_change,
            "is_frameshift": is_frameshift,
            "is_inframe": not is_frameshift and len_change != 0,
            "predicted_effect": self._predict_effect(is_frameshift, len_change),
            "interpretation": self._get_interpretation(is_frameshift, len_change)
        }

        logger.info(f"Indel analysis for variant: {self.format_variant_string(variant)}")

        self._set_cache(cache_key, result)
        return result

    def _classify_indel(self, is_insertion: bool, is_deletion: bool, is_frameshift: bool) -> str:
        """Classify the type of indel"""
        if is_insertion:
            return "frameshift_insertion" if is_frameshift else "inframe_insertion"
        elif is_deletion:
            return "frameshift_deletion" if is_frameshift else "inframe_deletion"
        return "unknown"

    def _predict_effect(self, is_frameshift: bool, len_change: int) -> str:
        """Predict the functional effect of the indel"""
        if is_frameshift:
            return "Frameshift resulting in altered amino acid sequence and likely premature stop codon (PVS1)"
        elif abs(len_change) >= 9:  # Large inframe indel (≥3 amino acids)
            return "Large inframe indel likely to disrupt protein structure"
        else:
            return "Small inframe indel with variable impact depending on location"

    def _get_interpretation(self, is_frameshift: bool, len_change: int) -> str:
        """Provide clinical interpretation"""
        if is_frameshift:
            return "Likely pathogenic - frameshift variants typically result in loss of function (PVS1 criterion)"
        elif abs(len_change) >= 9:
            return "Potentially pathogenic - large inframe indels often disrupt protein function"
        else:
            return "Uncertain significance - small inframe indels require additional evidence"


class CADDIndelPredictor(BasePredictor):
    """
    CADD scores for indels
    Uses CADD API to score insertions and deletions
    """

    def __init__(self, api_url: str = "https://cadd.gs.washington.edu/api/v1.0"):
        super().__init__(api_url)

    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        cache_key = self._cache_key(**variant)
        cached = self._get_cached(cache_key)
        if cached:
            return cached

        result = {
            "tool": "CADD-Indel",
            "phred_score": None,
            "raw_score": None,
            "interpretation": None,
            "url": "https://cadd.gs.washington.edu/"
        }

        logger.info(f"CADD-Indel prediction for variant: {self.format_variant_string(variant)}")

        # Would call CADD API here
        # variant_vcf = self.format_variant_string(variant, "vcf")
        # response = await self._make_request(...)

        self._set_cache(cache_key, result)
        return result


class IndelPredictorAggregator:
    """Aggregates indel analysis with CADD scoring"""

    def __init__(self):
        self.analyzers = {
            "IndelAnalyzer": IndelAnalyzer(),
            "CADD": CADDIndelPredictor()
        }

    async def predict_all(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        """Run all indel predictors and aggregate results"""
        results = {}

        for name, predictor in self.analyzers.items():
            try:
                results[name] = await predictor.predict(variant)
            except Exception as e:
                logger.error(f"Error in {name}: {str(e)}")
                results[name] = {"error": str(e), "tool": name}

        # Add aggregate interpretation
        results["aggregate"] = self._aggregate_predictions(results)
        return results

    def _aggregate_predictions(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Provide overall interpretation"""
        indel_analysis = results.get("IndelAnalyzer", {})
        cadd_result = results.get("CADD", {})

        is_frameshift = indel_analysis.get("is_frameshift", False)
        cadd_phred = cadd_result.get("phred_score")

        interpretations = []

        # Frameshift variants
        if is_frameshift:
            interpretations.append("Frameshift variant - typically loss of function (PVS1)")
            consensus = "Likely pathogenic"
            acmg_criteria = ["PVS1"]
        else:
            # Inframe indel
            len_change = abs(indel_analysis.get("length_change", 0))
            if len_change >= 9:
                interpretations.append("Large inframe indel - likely disrupts protein structure")
                consensus = "Likely pathogenic"
                acmg_criteria = ["PM4"]
            else:
                interpretations.append("Small inframe indel - impact depends on location and conservation")
                consensus = "Uncertain significance"
                acmg_criteria = []

        # Add CADD interpretation
        if cadd_phred:
            if cadd_phred >= 30:
                interpretations.append(f"CADD score {cadd_phred:.1f} - highly deleterious (top 0.1%)")
                acmg_criteria.append("PP3")
            elif cadd_phred >= 20:
                interpretations.append(f"CADD score {cadd_phred:.1f} - deleterious (top 1%)")
                acmg_criteria.append("PP3")

        return {
            "consensus": consensus,
            "acmg_criteria": acmg_criteria,
            "details": interpretations,
            "recommendation": self._get_recommendation(is_frameshift, consensus)
        }

    def _get_recommendation(self, is_frameshift: bool, consensus: str) -> str:
        """Provide recommendation for further analysis"""
        if is_frameshift:
            return "Check for NMD escape and confirm with RNA studies if available. Consider segregation analysis."
        elif consensus == "Likely pathogenic":
            return "Consider functional studies and segregation analysis to support pathogenicity"
        else:
            return "Evaluate conservation, domain location, and segregation data. Functional studies recommended."
