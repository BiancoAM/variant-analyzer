from typing import Dict, Any, Optional
from .base_predictor import BasePredictor
import logging

logger = logging.getLogger(__name__)


class SpliceAIPredictor(BasePredictor):
    """
    SpliceAI - Deep learning-based splice prediction (Jaganathan et al., Cell 2019)
    Predicts splice altering variants with delta scores for:
    - Acceptor gain/loss
    - Donor gain/loss
    Delta score range: 0-1 (>0.5 = high confidence, 0.2-0.5 = moderate, <0.2 = low)

    Uses the SpliceAI Lookup API from the Broad Institute:
    https://spliceailookup-api.broadinstitute.org/spliceai/
    Supports both GRCh37 (hg19) and GRCh38 (hg38).
    """

    def __init__(self, api_url: Optional[str] = None, genome: str = "38"):
        super().__init__(api_url or "https://spliceailookup-api.broadinstitute.org")
        self.genome = genome  # "38" for GRCh38, "37" for GRCh37

    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        cache_key = self._cache_key(**variant)
        cached = self._get_cached(cache_key)
        if cached:
            return cached

        chromosome = variant.get('chromosome', '').replace('chr', '')
        position = variant.get('position', '')
        ref = variant.get('reference', '')
        alt = variant.get('alternate', '')

        result = {
            "tool": "SpliceAI",
            "delta_scores": {
                "acceptor_gain": None,
                "acceptor_loss": None,
                "donor_gain": None,
                "donor_loss": None
            },
            "max_delta": None,
            "interpretation": None,
            "affected_genes": [],
            "url": "https://spliceailookup.broadinstitute.org/",
            "lookup_url": None
        }

        if chromosome and position and ref and alt:
            # Construct variant string: chr-pos-ref-alt (no "chr" prefix for Broad API)
            variant_str = f"{chromosome}-{position}-{ref}-{alt}"
            lookup_url = (
                f"https://spliceailookup.broadinstitute.org/"
                f"#variant={variant_str}&hg={self.genome}&distance=500&mask=0"
            )
            result["lookup_url"] = lookup_url
            result["url"] = lookup_url

            try:
                api_endpoint = (
                    f"{self.api_url}/spliceai/"
                    f"?hg={self.genome}&variant={variant_str}"
                )

                response = await self._make_request(api_endpoint)

                if response and isinstance(response, dict):
                    scores_list = response.get('scores', [])

                    if scores_list:
                        # Use the first gene result (highest impact)
                        best_gene = None
                        max_delta = 0.0

                        for gene_result in scores_list:
                            delta = gene_result.get('delta_scores', [])
                            if isinstance(delta, list) and len(delta) >= 4:
                                gene_max = max(delta)
                            elif isinstance(delta, dict):
                                gene_max = max(delta.values()) if delta else 0
                            else:
                                gene_max = 0

                            if gene_max > max_delta:
                                max_delta = gene_max
                                best_gene = gene_result

                        if best_gene:
                            delta = best_gene.get('delta_scores', [])
                            if isinstance(delta, list) and len(delta) >= 4:
                                result["delta_scores"] = {
                                    "acceptor_gain": round(delta[0], 4),
                                    "acceptor_loss": round(delta[1], 4),
                                    "donor_gain": round(delta[2], 4),
                                    "donor_loss": round(delta[3], 4)
                                }
                            elif isinstance(delta, dict):
                                result["delta_scores"] = {
                                    "acceptor_gain": round(delta.get("acceptor_gain", 0), 4),
                                    "acceptor_loss": round(delta.get("acceptor_loss", 0), 4),
                                    "donor_gain": round(delta.get("donor_gain", 0), 4),
                                    "donor_loss": round(delta.get("donor_loss", 0), 4)
                                }

                            result["max_delta"] = round(max_delta, 4)
                            result["affected_genes"] = [
                                g.get("gene_name", g.get("gene_id", ""))
                                for g in scores_list
                                if g.get("gene_name") or g.get("gene_id")
                            ]
                            result["interpretation"] = self._interpret_score(max_delta)
                    else:
                        result["interpretation"] = "Nessun effetto di splicing predetto da SpliceAI"

            except Exception as e:
                logger.warning(f"SpliceAI Lookup API failed: {str(e)} - result has no scores")
                result["interpretation"] = f"API non disponibile: {str(e)}"

        logger.info(f"SpliceAI prediction for variant: {self.format_variant_string(variant)}")

        self._set_cache(cache_key, result)
        return result

    def _interpret_score(self, max_delta: float) -> str:
        """Interpret SpliceAI delta score"""
        if max_delta >= 0.5:
            return f"Alta confidenza effetto splicing (Δ={max_delta:.3f}) - PVS1/PS3 applicabile"
        elif max_delta >= 0.2:
            return f"Moderata confidenza effetto splicing (Δ={max_delta:.3f}) - RNA studies raccomandati"
        elif max_delta >= 0.1:
            return f"Bassa confidenza effetto splicing (Δ={max_delta:.3f}) - impatto incerto"
        else:
            return f"Nessun effetto di splicing significativo predetto (Δ={max_delta:.3f})"


class MaxEntScanPredictor(BasePredictor):
    """
    MaxEntScan - Maximum Entropy Model for splice site prediction
    Scores splice site strength for both reference and alternate alleles
    Higher score = stronger splice site
    """

    def __init__(self, api_url: Optional[str] = None):
        super().__init__(api_url)

    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        cache_key = self._cache_key(**variant)
        cached = self._get_cached(cache_key)
        if cached:
            return cached

        result = {
            "tool": "MaxEntScan",
            "ref_score": None,
            "alt_score": None,
            "score_change": None,
            "interpretation": None,
            "url": "http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html"
        }

        logger.info(f"MaxEntScan prediction requested for variant: {self.format_variant_string(variant)}")

        # MaxEntScan requires splice site sequence context
        # Would need to extract 9bp (donor) or 23bp (acceptor) sequence

        self._set_cache(cache_key, result)
        return result


class SPIDEXPredictor(BasePredictor):
    """
    SPIDEX - Deep learning prediction of splice disruption
    Provides percentile ranks for splicing impact
    |dpsi_zscore| > 2 suggests significant splicing impact
    """

    def __init__(self):
        super().__init__()
        # SPIDEX is accessed via pre-computed database

    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        cache_key = self._cache_key(**variant)
        cached = self._get_cached(cache_key)
        if cached:
            return cached

        result = {
            "tool": "SPIDEX",
            "dpsi_max_tissue": None,
            "dpsi_zscore": None,
            "percentile": None,
            "interpretation": None,
            "url": "http://www.openbioinformatics.org/annovar/spidex_download_form.php"
        }

        logger.info(f"SPIDEX prediction requested for variant: {self.format_variant_string(variant)}")

        self._set_cache(cache_key, result)
        return result


class DbscSNVPredictor(BasePredictor):
    """
    dbscSNV - Database of human splicing variants
    Pre-computed ensemble scores for splice-site SNVs
    ada_score and rf_score: higher = more likely to affect splicing
    """

    def __init__(self):
        super().__init__()

    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        cache_key = self._cache_key(**variant)
        cached = self._get_cached(cache_key)
        if cached:
            return cached

        result = {
            "tool": "dbscSNV",
            "ada_score": None,
            "rf_score": None,
            "interpretation": None,
            "url": "https://sites.google.com/site/jpopgen/dbNSFP"
        }

        logger.info(f"dbscSNV prediction requested for variant: {self.format_variant_string(variant)}")

        self._set_cache(cache_key, result)
        return result


class SplicingPredictorAggregator:
    """Aggregates predictions from multiple splicing predictors"""

    def __init__(self):
        self.predictors = {
            "SpliceAI": SpliceAIPredictor(),
            "MaxEntScan": MaxEntScanPredictor(),
            "SPIDEX": SPIDEXPredictor(),
            "dbscSNV": DbscSNVPredictor()
        }

    async def predict_all(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        """Run all splicing predictors and aggregate results"""
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
        splice_altering_count = 0
        total_count = 0
        interpretations = []

        for tool, result in results.items():
            if "error" in result:
                continue

            total_count += 1

            # Check predictions based on tool-specific thresholds
            if tool == "SpliceAI" and result.get("max_delta") is not None:
                if result["max_delta"] >= 0.5:
                    splice_altering_count += 1
                    interpretations.append(f"SpliceAI: high confidence splice alteration (Δ={result['max_delta']:.3f})")
                elif result["max_delta"] >= 0.2:
                    interpretations.append(f"SpliceAI: moderate confidence splice alteration (Δ={result['max_delta']:.3f})")

            elif tool == "MaxEntScan" and result.get("score_change") is not None:
                # Significant score change suggests splice site disruption
                if abs(result["score_change"]) >= 2:
                    splice_altering_count += 1
                    interpretations.append(f"MaxEntScan: significant score change ({result['score_change']:.2f})")

            elif tool == "SPIDEX" and result.get("dpsi_zscore") is not None:
                if abs(result["dpsi_zscore"]) >= 2:
                    splice_altering_count += 1
                    interpretations.append(f"SPIDEX: significant splicing impact (Z={result['dpsi_zscore']:.2f})")

            elif tool == "dbscSNV":
                ada = result.get("ada_score", 0)
                rf = result.get("rf_score", 0)
                if ada and rf and (ada >= 0.6 or rf >= 0.6):
                    splice_altering_count += 1
                    interpretations.append(f"dbscSNV: likely splice-altering (ada={ada:.3f}, rf={rf:.3f})")

        if total_count == 0:
            consensus = "Insufficient data"
        elif splice_altering_count >= 2:
            consensus = "High confidence splice-altering variant"
        elif splice_altering_count >= 1:
            consensus = "Moderate confidence splice-altering variant"
        else:
            consensus = "Unlikely to affect splicing"

        return {
            "consensus": consensus,
            "splice_altering_predictors": splice_altering_count,
            "total_predictors": total_count,
            "details": interpretations,
            "recommendation": self._get_recommendation(splice_altering_count, total_count)
        }

    def _get_recommendation(self, splice_count: int, total: int) -> str:
        """Provide clinical interpretation recommendation"""
        if total == 0:
            return "Run in silico predictions to assess splice impact"

        if splice_count >= 2:
            return "Consider RNA studies to validate splice-altering effect (PM4 evidence)"
        elif splice_count >= 1:
            return "Moderate evidence for splice impact - additional functional studies recommended"
        else:
            return "Low evidence for splice impact based on in silico predictions"
