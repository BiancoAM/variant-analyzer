from typing import Dict, Any, List
import logging
import re
from ..predictors import (
    MissensePredictorAggregator,
    SplicingPredictorAggregator,
    IndelPredictorAggregator,
    SynonymousPredictorAggregator,
    UTRPredictorAggregator,
    MobiDetailsPredictor
)
from .pubmed_service import PubMedService
from .myvariant_service import MyVariantService

logger = logging.getLogger(__name__)


class VariantAnalysisService:
    """
    Main orchestration service for variant analysis
    Coordinates prediction tools, database queries, and literature search
    """

    def __init__(self, pubmed_email: str, pubmed_api_key: str = None):
        """
        Initialize variant analysis service

        Args:
            pubmed_email: Email for PubMed API
            pubmed_api_key: Optional API key for PubMed
        """
        self.missense_predictor = MissensePredictorAggregator()
        self.splicing_predictor = SplicingPredictorAggregator()
        self.indel_predictor = IndelPredictorAggregator()
        self.synonymous_predictor = SynonymousPredictorAggregator()
        self.utr_predictor = UTRPredictorAggregator()
        self.mobidetails_predictor = MobiDetailsPredictor()
        self.pubmed_service = PubMedService(pubmed_email, pubmed_api_key)
        self.myvariant_service = MyVariantService()

    async def analyze_variant(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        """
        Perform comprehensive analysis of a variant

        Args:
            variant: Dictionary with variant information
                - chromosome: str
                - position: int
                - reference: str
                - alternate: str
                - gene: str (optional)
                - transcript: str (optional)
                - hgvs_c: str (optional)
                - hgvs_p: str (optional)

        Returns:
            Complete analysis results including:
            - Variant information
            - In silico predictions
            - Population frequencies
            - Literature evidence
            - ACMG classification support
        """
        logger.info(f"Starting analysis for variant: {self._format_variant(variant)}")

        # Determine variant type
        variant_type = self._determine_variant_type(variant)
        variant['variant_type'] = variant_type

        # Run appropriate predictors based on variant type
        prediction_results = await self._run_predictions(variant, variant_type)

        # Search literature
        literature = await self._search_literature(variant)

        # Get population data (would integrate with gnomAD, etc.)
        population_data = await self._get_population_data(variant)

        # Generate ACMG evidence
        acmg_evidence = self._generate_acmg_evidence(
            variant,
            prediction_results,
            literature,
            population_data
        )

        # Compile complete analysis
        analysis = {
            "variant": variant,
            "variant_type": variant_type,
            "predictions": prediction_results,
            "population_data": population_data,
            "literature": literature,
            "acmg_classification": acmg_evidence,
            "summary": self._generate_summary(
                variant_type,
                prediction_results,
                literature,
                acmg_evidence
            ),
            "tool_links": self._generate_tool_links(variant)
        }

        return analysis

    def _determine_variant_type(self, variant: Dict[str, Any]) -> str:
        """Determine the type of variant"""
        ref = variant.get('reference', '')
        alt = variant.get('alternate', '')
        hgvs_c = variant.get('hgvs_c') or ''
        hgvs_p = variant.get('hgvs_p') or ''

        # --- UTR variants (detected from HGVS coding notation) ---
        # 5'UTR: c.-NNN (negative offset before CDS start)
        if re.match(r'^c\.-\d+', hgvs_c):
            return 'UTR_5'
        # 3'UTR: c.*NNN (asterisk + offset after stop codon)
        if re.match(r'^c\.\*\d+', hgvs_c):
            return 'UTR_3'

        len_diff = len(alt) - len(ref)

        if len(ref) == 1 and len(alt) == 1:
            # Synonymous: p.(=) or p.= in HGVS protein notation
            if hgvs_p and ('=' in hgvs_p):
                return 'SNV_SYNONYMOUS'

            # Splicing: explicit splice notation or canonical ±1/±2 positions
            if 'splice' in hgvs_c.lower():
                return 'SPLICING'
            if re.search(r'[+-][12]', hgvs_c):
                return 'SPLICING'

            return 'SNV_MISSENSE'
        elif len_diff != 0:
            if len_diff % 3 == 0:
                return 'INDEL_INFRAME'
            else:
                return 'INDEL_FRAMESHIFT'
        else:
            return 'COMPLEX'

    async def _run_predictions(
        self,
        variant: Dict[str, Any],
        variant_type: str
    ) -> Dict[str, Any]:
        """Run appropriate prediction tools based on variant type"""
        results = {}

        # Always run CADD for general deleteriousness
        try:
            from ..predictors.missense_predictors import CADDPredictor
            cadd = CADDPredictor()
            results['CADD'] = await cadd.predict(variant)
        except Exception as e:
            logger.error(f"Error running CADD: {str(e)}")
            results['CADD'] = {"error": str(e)}

        # Type-specific predictors
        if variant_type in ['SNV_MISSENSE', 'COMPLEX']:
            try:
                results['missense'] = await self.missense_predictor.predict_all(variant)
            except Exception as e:
                logger.error(f"Error running missense predictors: {str(e)}")
                results['missense'] = {"error": str(e)}

        # Synonymous variants
        if variant_type == 'SNV_SYNONYMOUS':
            try:
                results['synonymous'] = await self.synonymous_predictor.predict_all(variant)
            except Exception as e:
                logger.error(f"Error running synonymous predictors: {str(e)}")
                results['synonymous'] = {"error": str(e)}

        # UTR variants (5' or 3')
        if variant_type in ['UTR_5', 'UTR_3']:
            try:
                results['utr'] = await self.utr_predictor.predict_all(variant)
            except Exception as e:
                logger.error(f"Error running UTR predictors: {str(e)}")
                results['utr'] = {"error": str(e)}

        # Always check splicing impact (except pure UTR3/UTR5 without splice proximity)
        # For synonymous and missense SNVs, splicing check is mandatory
        if variant_type not in ['UTR_3']:
            try:
                results['splicing'] = await self.splicing_predictor.predict_all(variant)
            except Exception as e:
                logger.error(f"Error running splicing predictors: {str(e)}")
                results['splicing'] = {"error": str(e)}

        # Indel-specific analysis
        if 'INDEL' in variant_type:
            try:
                results['indel'] = await self.indel_predictor.predict_all(variant)
            except Exception as e:
                logger.error(f"Error running indel predictors: {str(e)}")
                results['indel'] = {"error": str(e)}

        return results

    async def _search_literature(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        """Search for relevant literature"""
        try:
            articles = await self.pubmed_service.search_variant(variant, max_results=20)

            # Categorize articles by evidence type
            functional_studies = [
                a for a in articles
                if a.get('evidence_type') == 'functional_study'
            ]
            case_reports = [
                a for a in articles
                if a.get('evidence_type') == 'case_report'
            ]
            reviews = [
                a for a in articles
                if a.get('evidence_type') == 'review'
            ]
            other = [
                a for a in articles
                if a.get('evidence_type') not in ['functional_study', 'case_report', 'review']
            ]

            return {
                "total_articles": len(articles),
                "functional_studies": functional_studies[:5],
                "case_reports": case_reports[:5],
                "reviews": reviews[:3],
                "other": other[:5],
                "search_url": self._generate_pubmed_search_url(variant)
            }

        except Exception as e:
            logger.error(f"Error searching literature: {str(e)}")
            return {
                "error": str(e),
                "search_url": self._generate_pubmed_search_url(variant)
            }

    async def _get_population_data(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        """
        Get population frequency data from MyVariant.info (gnomAD, dbSNP, ClinVar)
        """
        try:
            # Query MyVariant for real data
            myvariant_data = await self.myvariant_service.get_variant_annotations(variant)

            gnomad = myvariant_data.get('gnomad', {})
            clinvar = myvariant_data.get('clinvar', {})
            dbsnp = myvariant_data.get('dbsnp', {})
            cadd = myvariant_data.get('cadd', {})
            predictions = myvariant_data.get('predictions', {})

            return {
                # dbSNP
                "dbsnp_rsid": dbsnp.get('rsid'),
                "dbsnp_url": dbsnp.get('url'),

                # gnomAD
                "gnomad_genome_af": gnomad.get('genome', {}).get('af'),
                "gnomad_genome_af_popmax": gnomad.get('genome', {}).get('af_popmax'),
                "gnomad_exome_af": gnomad.get('exome', {}).get('af'),
                "gnomad_exome_af_popmax": gnomad.get('exome', {}).get('af_popmax'),
                "gnomad_max_af": gnomad.get('max_af'),
                "gnomad_url": self._generate_gnomad_url(variant),

                # ClinVar
                "clinvar_significance": clinvar.get('significance'),
                "clinvar_review_status": clinvar.get('review_status'),
                "clinvar_url": clinvar.get('url') or self._generate_clinvar_url(variant),

                # CADD from MyVariant
                "cadd_phred": cadd.get('phred'),
                "cadd_raw": cadd.get('raw'),
                "cadd_interpretation": cadd.get('interpretation'),

                # Predictions from dbNSFP
                "sift_score": predictions.get('sift', {}).get('score'),
                "sift_prediction": predictions.get('sift', {}).get('prediction'),
                "polyphen2_score": predictions.get('polyphen2', {}).get('score'),
                "polyphen2_prediction": predictions.get('polyphen2', {}).get('prediction'),
                "revel_score": predictions.get('revel', {}).get('score'),
                "vest4_score": predictions.get('vest4', {}).get('score'),

                # Metadata
                "data_source": "MyVariant.info (gnomAD v3, dbSNP, ClinVar, dbNSFP)",
                "found": myvariant_data.get('found', False)
            }
        except Exception as e:
            logger.error(f"Error getting population data from MyVariant: {str(e)}")
            return {
                "gnomad_max_af": None,
                "gnomad_url": self._generate_gnomad_url(variant),
                "clinvar_significance": None,
                "clinvar_url": self._generate_clinvar_url(variant),
                "dbsnp_rsid": None,
                "error": str(e),
                "found": False
            }

    def _generate_acmg_evidence(
        self,
        variant: Dict[str, Any],
        predictions: Dict[str, Any],
        literature: Dict[str, Any],
        population: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Generate ACMG/AMP classification evidence

        Following ACMG/AMP 2015 guidelines for variant interpretation
        """
        evidence = {
            "pathogenic": [],
            "benign": [],
            "classification": "Uncertain Significance (VUS)",
            "confidence": "Low"
        }

        # PVS1 - Null variant (frameshift, nonsense, splice site)
        if variant.get('variant_type') == 'INDEL_FRAMESHIFT':
            evidence["pathogenic"].append({
                "code": "PVS1",
                "description": "Null variant (frameshift) in a gene where LOF is a known mechanism",
                "strength": "Very Strong"
            })

        # PP3 / BP4 - Multiple in silico predictors
        missense_results = predictions.get('missense', {})
        if 'aggregate' in missense_results:
            aggregate = missense_results['aggregate']
            pathogenic_count = aggregate.get('pathogenic_predictors', 0)
            total_count = aggregate.get('total_predictors', 0)

            if total_count > 0:
                if pathogenic_count / total_count >= 0.75:
                    evidence["pathogenic"].append({
                        "code": "PP3",
                        "description": f"Multiple computational predictors support pathogenicity ({pathogenic_count}/{total_count})",
                        "strength": "Supporting"
                    })
                elif pathogenic_count / total_count <= 0.25:
                    evidence["benign"].append({
                        "code": "BP4",
                        "description": f"Multiple computational predictors suggest benign impact ({total_count-pathogenic_count}/{total_count})",
                        "strength": "Supporting"
                    })

        # PM4 - Protein length change (inframe indel)
        if variant.get('variant_type') == 'INDEL_INFRAME':
            indel_results = predictions.get('indel', {})
            indel_analysis = indel_results.get('IndelAnalyzer', {})
            if abs(indel_analysis.get('length_change', 0)) >= 9:
                evidence["pathogenic"].append({
                    "code": "PM4",
                    "description": "Large inframe indel in non-repeat region",
                    "strength": "Moderate"
                })

        # PS3 / BS3 - Functional studies
        functional_studies = literature.get('functional_studies', [])
        if functional_studies:
            evidence["pathogenic"].append({
                "code": "PS3",
                "description": f"Functional studies available ({len(functional_studies)} found) - review required",
                "strength": "Strong (if supportive)",
                "studies": [
                    {
                        "pmid": study.get('pmid'),
                        "title": study.get('title'),
                        "url": study.get('pubmed_url')
                    }
                    for study in functional_studies[:3]
                ]
            })

        # BP7 - Synonymous variant with no predicted splicing impact
        if variant.get('variant_type') == 'SNV_SYNONYMOUS':
            synonymous_results = predictions.get('synonymous', {})
            aggregate = synonymous_results.get('aggregate', {})
            splicing_results = predictions.get('splicing', {})
            splicing_agg = splicing_results.get('aggregate', {})

            splice_count = splicing_agg.get('splice_altering_predictors', 0) if splicing_agg else 0

            if splice_count == 0:
                evidence["benign"].append({
                    "code": "BP7",
                    "description": (
                        "Variante sinonima senza impatto predetto su splicing "
                        "(SpliceAI e predittori in silico non suggeriscono alterazione)"
                    ),
                    "strength": "Supporting"
                })
            else:
                evidence["pathogenic"].append({
                    "code": "PP3",
                    "description": (
                        f"Variante sinonima con possibile impatto su splicing "
                        f"({splice_count} predittori suggeriscono alterazione) - RNA studies raccomandati"
                    ),
                    "strength": "Supporting"
                })

        # UTR variants - limited ACMG criteria, rely on functional evidence
        if variant.get('variant_type') in ['UTR_5', 'UTR_3']:
            utr_results = predictions.get('utr', {})
            utr_agg = utr_results.get('aggregate', {})
            utr_criteria = utr_agg.get('acmg_criteria', []) if utr_agg else []

            for criterion in utr_criteria:
                code = criterion.split(' ')[0]
                if code.startswith('PP'):
                    evidence["pathogenic"].append({
                        "code": code,
                        "description": (
                            f"Variante UTR con evidenza regolatoria (RegulomeDB): {criterion}"
                        ),
                        "strength": "Supporting"
                    })
                elif code.startswith('BP'):
                    evidence["benign"].append({
                        "code": code,
                        "description": (
                            f"Variante UTR con bassa evidenza regolatoria (RegulomeDB): {criterion}"
                        ),
                        "strength": "Supporting"
                    })

        # PP2 - Missense in gene where missense is common mechanism
        # Would need gene-specific knowledge

        # BA1 / BS1 / PM2 - Population frequency
        gnomad_af = population.get('gnomad_max_af')
        if gnomad_af is not None:
            if gnomad_af > 0.05:
                evidence["benign"].append({
                    "code": "BA1",
                    "description": f"High allele frequency in population (AF={gnomad_af:.6f}, gnomAD)",
                    "strength": "Stand-alone"
                })
            elif gnomad_af > 0.01:
                evidence["benign"].append({
                    "code": "BS1",
                    "description": f"Elevated allele frequency (AF={gnomad_af:.6f}, gnomAD)",
                    "strength": "Strong"
                })
            elif gnomad_af == 0:
                evidence["pathogenic"].append({
                    "code": "PM2",
                    "description": "Absent from gnomAD population databases (supports rarity)",
                    "strength": "Moderate"
                })
        # If data found but AF is None, variant likely not in gnomAD
        elif population.get('found'):
            evidence["pathogenic"].append({
                "code": "PM2",
                "description": "Not found in gnomAD population databases (supports rarity)",
                "strength": "Moderate"
            })

        # Calculate overall classification
        evidence["classification"] = self._calculate_classification(evidence)

        return evidence

    def _calculate_classification(self, evidence: Dict[str, Any]) -> str:
        """
        Calculate ACMG classification based on evidence

        Simplified logic - full implementation would follow ACMG combination rules
        """
        pathogenic = evidence.get("pathogenic", [])
        benign = evidence.get("benign", [])

        # Check for stand-alone evidence
        for ev in benign:
            if ev["strength"] == "Stand-alone":
                return "Benign"

        for ev in pathogenic:
            if ev["strength"] == "Very Strong":
                if any(e["strength"] in ["Strong", "Moderate"] for e in pathogenic):
                    return "Pathogenic"
                return "Likely Pathogenic"

        # Count evidence by strength
        path_strong = sum(1 for e in pathogenic if e["strength"] == "Strong")
        path_moderate = sum(1 for e in pathogenic if e["strength"] == "Moderate")
        path_supporting = sum(1 for e in pathogenic if e["strength"] == "Supporting")

        ben_strong = sum(1 for e in benign if e["strength"] == "Strong")
        ben_supporting = sum(1 for e in benign if e["strength"] == "Supporting")

        # Simplified classification rules
        if path_strong >= 2 or (path_strong >= 1 and path_moderate >= 2):
            return "Pathogenic"
        elif path_strong >= 1 and path_moderate >= 1:
            return "Likely Pathogenic"
        elif ben_strong >= 2:
            return "Benign"
        elif ben_strong >= 1 and ben_supporting >= 1:
            return "Likely Benign"

        return "Uncertain Significance (VUS)"

    def _generate_summary(
        self,
        variant_type: str,
        predictions: Dict[str, Any],
        literature: Dict[str, Any],
        acmg: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Generate human-readable summary of analysis"""
        summary = {
            "classification": acmg.get("classification", "Unknown"),
            "key_findings": [],
            "recommendations": []
        }

        # Add key findings based on predictions
        if 'missense' in predictions:
            missense_consensus = predictions['missense'].get('aggregate', {}).get('consensus', '')
            if missense_consensus:
                summary["key_findings"].append(f"Missense prediction: {missense_consensus}")

        if 'synonymous' in predictions:
            syn_consensus = predictions['synonymous'].get('aggregate', {}).get('consensus', '')
            if syn_consensus:
                summary["key_findings"].append(f"Synonymous prediction: {syn_consensus}")

        if 'utr' in predictions:
            utr_consensus = predictions['utr'].get('aggregate', {}).get('consensus', '')
            utr_region = predictions['utr'].get('aggregate', {}).get('utr_region', 'UTR')
            if utr_consensus:
                summary["key_findings"].append(f"{utr_region} prediction: {utr_consensus}")

        if 'splicing' in predictions:
            splicing_consensus = predictions['splicing'].get('aggregate', {}).get('consensus', '')
            if splicing_consensus:
                summary["key_findings"].append(f"Splicing prediction: {splicing_consensus}")

        # Literature findings
        functional_count = len(literature.get('functional_studies', []))
        if functional_count > 0:
            summary["key_findings"].append(f"{functional_count} functional studies found in literature")

        # Recommendations
        if functional_count == 0 and acmg.get("classification") == "Uncertain Significance (VUS)":
            summary["recommendations"].append("Consider functional studies to establish pathogenicity")

        if variant_type in ['INDEL_FRAMESHIFT', 'INDEL_INFRAME']:
            summary["recommendations"].append("Consider RNA studies to confirm variant effect")

        if variant_type == 'SNV_SYNONYMOUS':
            summary["recommendations"].append(
                "Analizzare ESE/ESS con ESEfinder e HEXplorer (impatto su splicing)"
            )
            summary["recommendations"].append(
                "Se vicino a sito di splicing (±10bp): RNA studies raccomandati"
            )

        if variant_type == 'UTR_5':
            summary["recommendations"].append(
                "5'UTR: verificare creazione/disruzione uORF con UTRAnnotator (VEP)"
            )
            summary["recommendations"].append(
                "Reporter assay (luciferasi) per valutare impatto su traduzione"
            )

        if variant_type == 'UTR_3':
            summary["recommendations"].append(
                "3'UTR: verificare seed region miRNA con TargetScan e miRDB"
            )
            summary["recommendations"].append(
                "Analizzare binding RBP con RBPmap e segnali poliadenilazione con PolyASite"
            )

        summary["recommendations"].append("Review segregation data if family members available")
        summary["recommendations"].append("Check for additional cases in clinical databases (ClinVar, HGMD)")

        return summary

    def _generate_tool_links(self, variant: Dict[str, Any]) -> Dict[str, str]:
        """Generate URLs to external prediction tools and databases"""
        chr = variant.get('chromosome', '').replace('chr', '')
        pos = variant.get('position', '')
        ref = variant.get('reference', '')
        alt = variant.get('alternate', '')
        gene = variant.get('gene', '')

        links = {
            "gnomAD": self._generate_gnomad_url(variant),
            "ClinVar": self._generate_clinvar_url(variant),
            "UCSC_Genome_Browser": f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr{chr}:{pos}-{pos}",
            "Ensembl": f"https://www.ensembl.org/Homo_sapiens/Variation/Explore?v={chr}-{pos}-{ref}-{alt}",
            "PubMed_Search": self._generate_pubmed_search_url(variant),
            "CADD": "https://cadd.gs.washington.edu/",
            "SpliceAI_Lookup": f"https://spliceailookup.broadinstitute.org/#variant={chr}-{pos}-{ref}-{alt}&hg=38",
            "MobiDetails": self._generate_mobidetails_url(variant),
        }

        if gene:
            links["GeneCards"] = f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}"
            links["OMIM"] = f"https://www.omim.org/search?search={gene}"
            links["GTEx"] = f"https://gtexportal.org/home/gene/{gene}"

        # UTR-specific links
        variant_type = variant.get('variant_type', '')
        if variant_type == 'UTR_3':
            links["TargetScan"] = f"https://www.targetscan.org/cgi-bin/targetscan/vert_80/targetscan.cgi?species=Human&gid={gene}" if gene else "https://www.targetscan.org/"
            links["miRDB"] = "https://mirdb.org/"
            links["RBPmap"] = "http://rbpmap.technion.ac.il/"
            links["PolyASite"] = "https://polyasite.unibas.ch/"
        elif variant_type == 'UTR_5':
            links["UTRAnnotator_VEP"] = "https://www.ensembl.org/vep"
            links["RegulomeDB"] = f"https://regulomedb.org/regulome-search/?regions=chr{chr}:{pos}-{pos}&genome=GRCh38"
        elif variant_type == 'SNV_SYNONYMOUS':
            links["ESEfinder"] = "http://krainer01.cshl.edu/cgi-bin/tools/ESE3/esefinder.cgi"
            links["HEXplorer"] = "https://hexplorer.uni-due.de/"
            links["SILVA"] = "https://silva.bioinfo.cipf.es/"

        return links

    def _generate_mobidetails_url(self, variant: Dict[str, Any]) -> str:
        """Generate MobiDetails direct link"""
        import urllib.parse
        transcript = variant.get('transcript', '')
        hgvs_c = variant.get('hgvs_c', '')
        gene = variant.get('gene', '')

        if transcript and hgvs_c and transcript.startswith('NM_'):
            query = urllib.parse.quote(f"{transcript}:{hgvs_c}")
            return f"https://mobidetails.iurc.montp.inserm.fr/MD/search/?query={query}"
        elif gene:
            return f"https://mobidetails.iurc.montp.inserm.fr/MD/search/?query={urllib.parse.quote(gene)}"
        return "https://mobidetails.iurc.montp.inserm.fr/MD/"

    def _generate_gnomad_url(self, variant: Dict[str, Any]) -> str:
        """Generate gnomAD URL for variant"""
        chr = variant.get('chromosome', '').replace('chr', '')
        pos = variant.get('position', '')
        ref = variant.get('reference', '')
        alt = variant.get('alternate', '')
        return f"https://gnomad.broadinstitute.org/variant/{chr}-{pos}-{ref}-{alt}"

    def _generate_clinvar_url(self, variant: Dict[str, Any]) -> str:
        """Generate ClinVar search URL"""
        chr = variant.get('chromosome', '').replace('chr', '')
        pos = variant.get('position', '')
        return f"https://www.ncbi.nlm.nih.gov/clinvar/?term={chr}[chr]+AND+{pos}[chrpos37]"

    def _generate_pubmed_search_url(self, variant: Dict[str, Any]) -> str:
        """Generate PubMed search URL for variant"""
        gene = variant.get('gene', '')
        if gene:
            return f"https://pubmed.ncbi.nlm.nih.gov/?term={gene}+variant"
        return "https://pubmed.ncbi.nlm.nih.gov/"

    def _format_variant(self, variant: Dict[str, Any]) -> str:
        """Format variant as readable string"""
        chr = variant.get('chromosome', '')
        pos = variant.get('position', '')
        ref = variant.get('reference', '')
        alt = variant.get('alternate', '')
        gene = variant.get('gene', '')

        if gene:
            return f"{gene} {chr}:{pos}{ref}>{alt}"
        return f"{chr}:{pos}{ref}>{alt}"
