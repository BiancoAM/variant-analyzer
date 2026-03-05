"""
MyVariant.info Service
Provides access to gnomAD, dbSNP, CADD, SIFT, PolyPhen and more through a single API
"""
from typing import Dict, Any, Optional, Tuple
import logging
import myvariant
import asyncio
from pyliftover import LiftOver

logger = logging.getLogger(__name__)


class MyVariantService:
    """
    Service for querying MyVariant.info API
    Aggregates data from: gnomAD, dbSNP, ClinVar, CADD, dbNSFP (SIFT, PolyPhen, etc.)
    """

    def __init__(self):
        self.mv = myvariant.MyVariantInfo()
        # Initialize liftover for hg38 -> hg19 conversion
        try:
            self.lo = LiftOver('hg38', 'hg19')
            logger.info("LiftOver hg38->hg19 initialized successfully")
        except Exception as e:
            logger.warning(f"Could not initialize LiftOver: {e}. Will attempt queries with original coordinates only.")
            self.lo = None

    def format_hgvs(self, variant: Dict[str, Any]) -> str:
        """
        Format variant as HGVS genomic notation for MyVariant query
        Returns: chr17:g.43045677G>A
        """
        chr = variant.get('chromosome', '').replace('chr', '')
        pos = variant.get('position')
        ref = variant.get('reference', '')
        alt = variant.get('alternate', '')

        # Handle different variant types
        if len(ref) == 1 and len(alt) == 1:
            # SNV
            return f"chr{chr}:g.{pos}{ref}>{alt}"
        elif len(ref) > len(alt):
            # Deletion
            if len(alt) == 1:
                # Simple deletion
                start = pos + 1
                end = pos + len(ref) - 1
                return f"chr{chr}:g.{start}_{end}del"
            else:
                # Delins
                start = pos + 1
                end = pos + len(ref) - 1
                return f"chr{chr}:g.{start}_{end}delins{alt[1:]}"
        elif len(alt) > len(ref):
            # Insertion
            if len(ref) == 1:
                # Simple insertion
                return f"chr{chr}:g.{pos}_{pos+1}ins{alt[1:]}"
            else:
                # Delins
                start = pos + 1
                end = pos + len(ref) - 1
                return f"chr{chr}:g.{start}_{end}delins{alt[1:]}"
        else:
            # Complex
            return f"chr{chr}:g.{pos}{ref}>{alt}"

    def _convert_to_hg19(self, variant: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Convert variant coordinates from hg38 to hg19 using pyliftover
        Returns: New variant dict with hg19 coordinates, or None if conversion fails
        """
        if not self.lo:
            return None

        try:
            chr = variant.get('chromosome', '').replace('chr', '')
            pos = variant.get('position')

            if not chr or not pos:
                logger.warning("Missing chromosome or position for conversion")
                return None

            # LiftOver expects 'chrN' format and returns 0-based coordinates
            # We need to convert back to 1-based
            chr_with_prefix = f'chr{chr}'

            # Convert coordinate (position is 1-based, convert to 0-based for liftover)
            converted = self.lo.convert_coordinate(chr_with_prefix, pos - 1)

            if not converted or len(converted) == 0:
                logger.warning(f"Could not convert {chr_with_prefix}:{pos} from hg38 to hg19")
                return None

            # LiftOver returns list of tuples: [(target_chr, target_pos, strand, conversion_chain_score)]
            target_chr = converted[0][0].replace('chr', '')
            target_pos = int(converted[0][1]) + 1  # Convert back to 1-based

            logger.info(f"Converted hg38 chr{chr}:{pos} -> hg19 chr{target_chr}:{target_pos}")

            # Create new variant dict with hg19 coordinates
            variant_hg19 = variant.copy()
            variant_hg19['chromosome'] = target_chr
            variant_hg19['position'] = target_pos

            return variant_hg19

        except Exception as e:
            logger.error(f"Error converting coordinates from hg38 to hg19: {e}")
            return None

    async def get_variant_annotations(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        """
        Query MyVariant.info for comprehensive variant annotations
        Automatically tries hg38 coordinates first, then converts to hg19 if needed

        Returns dict with:
        - dbsnp: rsID and other dbSNP data
        - gnomad: population frequencies
        - cadd: CADD scores
        - dbnsfp: SIFT, PolyPhen, REVEL, etc.
        - clinvar: clinical significance
        """
        try:
            # Query MyVariant with specific fields
            fields = [
                'dbsnp.rsid',
                'dbsnp.alleles',
                'gnomad_genome.af',
                'gnomad_genome.af_popmax',
                'gnomad_exome.af',
                'gnomad_exome.af_popmax',
                'cadd.phred',
                'cadd.raw',
                'dbnsfp.sift.score',
                'dbnsfp.sift.pred',
                'dbnsfp.polyphen2.hdiv.score',
                'dbnsfp.polyphen2.hdiv.pred',
                'dbnsfp.revel.score',
                'dbnsfp.vest4.score',
                'dbnsfp.mutationtaster.score',
                'dbnsfp.mutationtaster.pred',
                'clinvar.rcv',
                'clinvar.clinical_significance',
                'clinvar.review_status',
                'vcf.position',
                'vcf.ref',
                'vcf.alt'
            ]

            # Try with original coordinates first (assume hg38)
            hgvs_id = self.format_hgvs(variant)
            logger.info(f"Querying MyVariant for: {hgvs_id} (trying as hg38 first)")

            # Run synchronous MyVariant call in executor
            loop = asyncio.get_event_loop()
            result = await loop.run_in_executor(
                None,
                lambda: self.mv.getvariant(hgvs_id, fields=','.join(fields))
            )

            # If not found and we have liftover, try converting to hg19
            if not result and self.lo:
                logger.info(f"Not found with hg38 coordinates, converting to hg19...")
                variant_hg19 = self._convert_to_hg19(variant)

                if variant_hg19:
                    hgvs_id_hg19 = self.format_hgvs(variant_hg19)
                    logger.info(f"Querying MyVariant with hg19 coordinates: {hgvs_id_hg19}")

                    result = await loop.run_in_executor(
                        None,
                        lambda: self.mv.getvariant(hgvs_id_hg19, fields=','.join(fields))
                    )

                    if result:
                        logger.info(f"✅ Found data with hg19 coordinates!")
                else:
                    logger.warning("Could not convert coordinates to hg19")

            if not result:
                logger.warning(f"No data found for {hgvs_id} (tried both hg38 and hg19)")
                return self._empty_result()

            # Parse and structure the result
            annotations = self._parse_myvariant_result(result)
            return annotations

        except Exception as e:
            logger.error(f"Error querying MyVariant: {str(e)}")
            return self._empty_result()

    def _parse_myvariant_result(self, result: Dict[str, Any]) -> Dict[str, Any]:
        """Parse MyVariant.info result into structured format"""

        annotations = {
            'found': True,
            'dbsnp': self._parse_dbsnp(result.get('dbsnp', {})),
            'gnomad': self._parse_gnomad(result),
            'cadd': self._parse_cadd(result.get('cadd', {})),
            'predictions': self._parse_predictions(result.get('dbnsfp', {})),
            'clinvar': self._parse_clinvar(result.get('clinvar', {})),
            'raw': result  # Keep raw data for debugging
        }

        return annotations

    def _parse_dbsnp(self, dbsnp_data: Dict[str, Any]) -> Dict[str, Any]:
        """Parse dbSNP data"""
        if not dbsnp_data:
            return {'rsid': None, 'alleles': []}

        # Handle both single rsid and list
        rsid = dbsnp_data.get('rsid')
        if isinstance(rsid, list):
            rsid = rsid[0] if rsid else None

        return {
            'rsid': rsid,
            'alleles': dbsnp_data.get('alleles', []),
            'url': f"https://www.ncbi.nlm.nih.gov/snp/{rsid}" if rsid else None
        }

    def _parse_gnomad(self, result: Dict[str, Any]) -> Dict[str, Any]:
        """Parse gnomAD population frequency data"""
        gnomad = {
            'genome': {
                'af': None,
                'af_popmax': None,
                'ac': None,
                'an': None
            },
            'exome': {
                'af': None,
                'af_popmax': None,
                'ac': None,
                'an': None
            }
        }

        # Genome data
        if 'gnomad_genome' in result:
            genome = result['gnomad_genome']
            if isinstance(genome, dict):
                gnomad['genome']['af'] = genome.get('af', {}).get('af') if isinstance(genome.get('af'), dict) else genome.get('af')
                gnomad['genome']['af_popmax'] = genome.get('af_popmax')
                gnomad['genome']['ac'] = genome.get('ac')
                gnomad['genome']['an'] = genome.get('an')

        # Exome data
        if 'gnomad_exome' in result:
            exome = result['gnomad_exome']
            if isinstance(exome, dict):
                gnomad['exome']['af'] = exome.get('af', {}).get('af') if isinstance(exome.get('af'), dict) else exome.get('af')
                gnomad['exome']['af_popmax'] = exome.get('af_popmax')
                gnomad['exome']['ac'] = exome.get('ac')
                gnomad['exome']['an'] = exome.get('an')

        # Use the higher AF from genome or exome
        af_genome = gnomad['genome']['af'] or 0
        af_exome = gnomad['exome']['af'] or 0
        gnomad['max_af'] = max(af_genome, af_exome) if (af_genome or af_exome) else None

        return gnomad

    def _parse_cadd(self, cadd_data: Dict[str, Any]) -> Dict[str, Any]:
        """Parse CADD scores"""
        if not cadd_data:
            return {'phred': None, 'raw': None}

        return {
            'phred': cadd_data.get('phred'),
            'raw': cadd_data.get('raw'),
            'interpretation': self._interpret_cadd(cadd_data.get('phred'))
        }

    def _interpret_cadd(self, phred_score: Optional[float]) -> Optional[str]:
        """Interpret CADD PHRED score"""
        if phred_score is None:
            return None
        if phred_score >= 30:
            return "Highly deleterious (top 0.1%)"
        elif phred_score >= 20:
            return "Deleterious (top 1%)"
        elif phred_score >= 10:
            return "Moderately deleterious (top 10%)"
        else:
            return "Likely benign"

    def _parse_predictions(self, dbnsfp_data: Dict[str, Any]) -> Dict[str, Any]:
        """Parse prediction scores from dbNSFP"""
        predictions = {
            'sift': {'score': None, 'prediction': None},
            'polyphen2': {'score': None, 'prediction': None},
            'revel': {'score': None},
            'vest4': {'score': None},
            'mutationtaster': {'score': None, 'prediction': None}
        }

        if not dbnsfp_data:
            return predictions

        # SIFT
        if 'sift' in dbnsfp_data:
            sift = dbnsfp_data['sift']
            if isinstance(sift, dict):
                predictions['sift']['score'] = self._get_first_value(sift.get('score'))
                predictions['sift']['prediction'] = self._get_first_value(sift.get('pred'))

        # PolyPhen2
        if 'polyphen2' in dbnsfp_data and 'hdiv' in dbnsfp_data['polyphen2']:
            polyphen = dbnsfp_data['polyphen2']['hdiv']
            if isinstance(polyphen, dict):
                predictions['polyphen2']['score'] = self._get_first_value(polyphen.get('score'))
                predictions['polyphen2']['prediction'] = self._get_first_value(polyphen.get('pred'))

        # REVEL
        if 'revel' in dbnsfp_data:
            revel = dbnsfp_data['revel']
            if isinstance(revel, dict):
                predictions['revel']['score'] = self._get_first_value(revel.get('score'))
            elif isinstance(revel, list):
                predictions['revel']['score'] = self._get_first_value(revel)

        # VEST4
        if 'vest4' in dbnsfp_data:
            vest4 = dbnsfp_data['vest4']
            if isinstance(vest4, dict):
                predictions['vest4']['score'] = self._get_first_value(vest4.get('score'))

        # MutationTaster
        if 'mutationtaster' in dbnsfp_data:
            mt = dbnsfp_data['mutationtaster']
            if isinstance(mt, dict):
                predictions['mutationtaster']['score'] = self._get_first_value(mt.get('score'))
                predictions['mutationtaster']['prediction'] = self._get_first_value(mt.get('pred'))

        return predictions

    def _parse_clinvar(self, clinvar_data: Dict[str, Any]) -> Dict[str, Any]:
        """Parse ClinVar data"""
        if not clinvar_data:
            return {
                'significance': None,
                'review_status': None,
                'rcv': None
            }

        # Handle list of submissions
        if isinstance(clinvar_data, list):
            clinvar_data = clinvar_data[0] if clinvar_data else {}

        significance = clinvar_data.get('clinical_significance')
        if isinstance(significance, list):
            significance = ', '.join(significance)

        return {
            'significance': significance,
            'review_status': clinvar_data.get('review_status'),
            'rcv': clinvar_data.get('rcv'),
            'url': f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_data.get('variation_id')}/" if clinvar_data.get('variation_id') else None
        }

    def _get_first_value(self, value):
        """Get first value from list or return value"""
        if isinstance(value, list):
            return value[0] if value else None
        return value

    def _empty_result(self) -> Dict[str, Any]:
        """Return empty result structure"""
        return {
            'found': False,
            'dbsnp': {'rsid': None, 'alleles': []},
            'gnomad': {
                'genome': {'af': None, 'af_popmax': None},
                'exome': {'af': None, 'af_popmax': None},
                'max_af': None
            },
            'cadd': {'phred': None, 'raw': None},
            'predictions': {
                'sift': {'score': None, 'prediction': None},
                'polyphen2': {'score': None, 'prediction': None},
                'revel': {'score': None}
            },
            'clinvar': {'significance': None, 'review_status': None}
        }
