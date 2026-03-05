from typing import Dict, Any, Optional
from .base_predictor import BasePredictor
import logging

logger = logging.getLogger(__name__)


class SIFTPredictor(BasePredictor):
    """
    SIFT (Sorting Intolerant From Tolerant)
    Predicts whether an amino acid substitution affects protein function
    Score range: 0-1 (≤0.05 = deleterious, >0.05 = tolerated)
    """

    def __init__(self, api_url: Optional[str] = None):
        super().__init__(api_url or "https://sift.bii.a-star.edu.sg/api")

    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        cache_key = self._cache_key(**variant)
        cached = self._get_cached(cache_key)
        if cached:
            return cached

        # Note: SIFT requires protein coordinates and amino acid change
        # In production, you'd need to annotate with VEP or similar first
        result = {
            "tool": "SIFT",
            "score": None,
            "prediction": None,
            "interpretation": None,
            "url": f"https://sift.bii.a-star.edu.sg/www/SIFT_seq_submit2.html"
        }

        # Placeholder for actual API call
        # In real implementation:
        # 1. Map genomic coordinates to protein coordinates
        # 2. Call SIFT API
        # 3. Parse response

        logger.info(f"SIFT prediction requested for variant: {self.format_variant_string(variant)}")

        # Example structure for when API is available:
        # result = await self._make_request(
        #     f"{self.api_url}/predict",
        #     method="POST",
        #     json={
        #         "sequence": protein_sequence,
        #         "position": aa_position,
        #         "substitution": aa_change
        #     }
        # )

        self._set_cache(cache_key, result)
        return result


class PolyPhen2Predictor(BasePredictor):
    """
    PolyPhen-2 (Polymorphism Phenotyping v2)
    Predicts the impact of amino acid substitutions on protein structure/function
    Score range: 0-1 (>0.85 = probably damaging, 0.15-0.85 = possibly damaging, <0.15 = benign)
    """

    def __init__(self, api_url: Optional[str] = None):
        super().__init__(api_url or "http://genetics.bwh.harvard.edu/ggi/pph2")

    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        cache_key = self._cache_key(**variant)
        cached = self._get_cached(cache_key)
        if cached:
            return cached

        result = {
            "tool": "PolyPhen-2",
            "score": None,
            "prediction": None,
            "interpretation": None,
            "url": "http://genetics.bwh.harvard.edu/pph2/"
        }

        logger.info(f"PolyPhen-2 prediction requested for variant: {self.format_variant_string(variant)}")

        self._set_cache(cache_key, result)
        return result


class CADDPredictor(BasePredictor):
    """
    CADD (Combined Annotation Dependent Depletion)
    Integrates multiple annotations into one metric (C-scores)
    Phred-scaled score: >20 = top 1% deleteriousness, >30 = top 0.1%
    """

    def __init__(self, api_url: Optional[str] = None):
        super().__init__(api_url or "https://cadd.gs.washington.edu/api/v1.0")

    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        cache_key = self._cache_key(**variant)
        cached = self._get_cached(cache_key)
        if cached:
            return cached

        result = {
            "tool": "CADD",
            "phred_score": None,
            "raw_score": None,
            "interpretation": None,
            "url": "https://cadd.gs.washington.edu/"
        }

        # CADD has a REST API that accepts VCF format
        variant_string = self.format_variant_string(variant, "vcf")

        logger.info(f"CADD prediction requested for variant: {self.format_variant_string(variant)}")

        # Example API call structure:
        # response = await self._make_request(
        #     f"{self.api_url}/score",
        #     method="POST",
        #     data={"variant": variant_string}
        # )

        self._set_cache(cache_key, result)
        return result


class REVELPredictor(BasePredictor):
    """
    REVEL (Rare Exome Variant Ensemble Learner)
    Ensemble method for predicting pathogenicity of missense variants
    Score range: 0-1 (>0.5 = likely pathogenic, >0.75 = pathogenic)
    """

    def __init__(self):
        super().__init__()
        # REVEL is typically accessed via downloaded database or VEP plugin

    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        cache_key = self._cache_key(**variant)
        cached = self._get_cached(cache_key)
        if cached:
            return cached

        result = {
            "tool": "REVEL",
            "score": None,
            "interpretation": None,
            "url": "https://sites.google.com/site/revelgenomics/"
        }

        logger.info(f"REVEL prediction requested for variant: {self.format_variant_string(variant)}")

        self._set_cache(cache_key, result)
        return result


class MissensePredictorAggregator:
    """Aggregates predictions from multiple missense predictors"""

    def __init__(self):
        self.predictors = {
            "SIFT": SIFTPredictor(),
            "PolyPhen2": PolyPhen2Predictor(),
            "CADD": CADDPredictor(),
            "REVEL": REVELPredictor()
        }

    async def predict_all(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        """Run all missense predictors and aggregate results"""
        results = {}

        for name, predictor in self.predictors.items():
            try:
                results[name] = await predictor.predict(variant)
            except Exception as e:
                logger.error(f"Error in {name} predictor: {str(e)}")
                results[name] = {"error": str(e), "tool": name}

        # Aggregate interpretation
        results["aggregate"] = self._aggregate_predictions(results)
        return results

    def _aggregate_predictions(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Provide overall interpretation based on multiple predictors"""
        pathogenic_count = 0
        total_count = 0

        # Simplified aggregation logic
        interpretations = []

        for tool, result in results.items():
            if "error" in result:
                continue

            total_count += 1

            # Check predictions based on tool-specific thresholds
            if tool == "SIFT" and result.get("score") is not None:
                if result["score"] <= 0.05:
                    pathogenic_count += 1
                    interpretations.append(f"SIFT: deleterious ({result['score']:.3f})")

            elif tool == "PolyPhen2" and result.get("score") is not None:
                if result["score"] >= 0.85:
                    pathogenic_count += 1
                    interpretations.append(f"PolyPhen-2: probably damaging ({result['score']:.3f})")

            elif tool == "CADD" and result.get("phred_score") is not None:
                if result["phred_score"] >= 20:
                    pathogenic_count += 1
                    interpretations.append(f"CADD: deleterious (PHRED {result['phred_score']:.1f})")

            elif tool == "REVEL" and result.get("score") is not None:
                if result["score"] >= 0.5:
                    pathogenic_count += 1
                    interpretations.append(f"REVEL: likely pathogenic ({result['score']:.3f})")

        if total_count == 0:
            consensus = "Insufficient data"
        elif pathogenic_count / total_count >= 0.75:
            consensus = "Strong evidence for pathogenicity"
        elif pathogenic_count / total_count >= 0.5:
            consensus = "Moderate evidence for pathogenicity"
        else:
            consensus = "Likely benign or uncertain"

        return {
            "consensus": consensus,
            "pathogenic_predictors": pathogenic_count,
            "total_predictors": total_count,
            "details": interpretations
        }
