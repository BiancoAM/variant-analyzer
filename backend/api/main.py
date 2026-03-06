from fastapi import FastAPI, HTTPException, Depends
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
import logging
import os
from pathlib import Path
from dotenv import load_dotenv

from services.variant_service import VariantAnalysisService

# Load environment variables
load_dotenv()

# Configure logging
logging.basicConfig(
    level=os.getenv('LOG_LEVEL', 'INFO'),
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Initialize FastAPI app
app = FastAPI(
    title="Variant Analyzer API",
    description="Comprehensive genetic variant analysis with in silico predictions and literature search",
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc"
)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, specify actual origins
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files for frontend
frontend_path = Path(__file__).parent.parent.parent / "frontend"
if frontend_path.exists():
    app.mount("/static", StaticFiles(directory=str(frontend_path)), name="static")


# Pydantic models for request/response
class VariantInput(BaseModel):
    """Input model for variant analysis"""
    chromosome: str = Field(..., description="Chromosome (e.g., '1', 'chr1', 'X')")
    position: int = Field(..., description="Genomic position", gt=0)
    reference: str = Field(..., description="Reference allele")
    alternate: str = Field(..., description="Alternate allele")
    gene: Optional[str] = Field(None, description="Gene symbol (optional)")
    transcript: Optional[str] = Field(None, description="Transcript ID (optional)")
    hgvs_c: Optional[str] = Field(None, description="HGVS coding notation (optional)")
    hgvs_p: Optional[str] = Field(None, description="HGVS protein notation (optional)")

    class Config:
        schema_extra = {
            "example": {
                "chromosome": "17",
                "position": 43044295,
                "reference": "G",
                "alternate": "A",
                "gene": "BRCA1",
                "hgvs_c": "c.5266dup",
                "hgvs_p": "p.Gln1756fs"
            }
        }


class VariantBatchInput(BaseModel):
    """Input model for batch variant analysis"""
    variants: List[VariantInput] = Field(..., description="List of variants to analyze")
    max_variants: int = Field(100, description="Maximum variants per batch", le=100)


class LiteratureSearchInput(BaseModel):
    """Input for literature search"""
    gene: str = Field(..., description="Gene symbol")
    variant_description: Optional[str] = Field(None, description="Variant description")
    max_results: int = Field(20, description="Maximum results", le=50)


# Dependency for service injection
def get_analysis_service() -> VariantAnalysisService:
    """Dependency to get analysis service instance"""
    pubmed_email = os.getenv('NCBI_EMAIL', 'user@example.com')
    pubmed_api_key = os.getenv('NCBI_API_KEY')

    return VariantAnalysisService(pubmed_email, pubmed_api_key)


@app.get("/")
async def root():
    """Serve the frontend application"""
    frontend_file = Path(__file__).parent.parent.parent / "frontend" / "index.html"
    if frontend_file.exists():
        return FileResponse(str(frontend_file))
    return {
        "message": "Variant Analyzer API",
        "version": "1.0.0",
        "docs": "/docs",
        "endpoints": {
            "analyze_variant": "/api/v1/analyze",
            "batch_analyze": "/api/v1/analyze/batch",
            "search_literature": "/api/v1/literature/search",
            "health": "/health"
        }
    }


@app.get("/style.css")
async def get_css():
    """Serve CSS file"""
    css_file = Path(__file__).parent.parent.parent / "frontend" / "style.css"
    if css_file.exists():
        return FileResponse(str(css_file), media_type="text/css")
    raise HTTPException(status_code=404, detail="CSS file not found")


@app.get("/app.js")
async def get_js():
    """Serve JavaScript file"""
    js_file = Path(__file__).parent.parent.parent / "frontend" / "app.js"
    if js_file.exists():
        return FileResponse(str(js_file), media_type="application/javascript")
    raise HTTPException(status_code=404, detail="JS file not found")


@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "service": "variant-analyzer",
        "version": "1.0.0"
    }


@app.post("/api/v1/analyze", response_model=Dict[str, Any])
async def analyze_variant(
    variant: VariantInput,
    service: VariantAnalysisService = Depends(get_analysis_service)
):
    """
    Analyze a single genetic variant

    Performs comprehensive analysis including:
    - In silico predictions (SIFT, PolyPhen, CADD, SpliceAI, etc.)
    - Population frequency data
    - Literature search (PubMed)
    - ACMG/AMP classification evidence

    Returns:
        Complete analysis results with predictions, literature, and classification
    """
    try:
        logger.info(f"Analyzing variant: {variant.gene} {variant.chromosome}:{variant.position}")

        # Convert to dict
        variant_dict = variant.dict()

        # Perform analysis
        result = await service.analyze_variant(variant_dict)

        logger.info(f"Analysis complete: {result.get('acmg_classification', {}).get('classification', 'Unknown')}")

        return result

    except Exception as e:
        logger.error(f"Error analyzing variant: {str(e)}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Error analyzing variant: {str(e)}")


@app.post("/api/v1/analyze/batch")
async def analyze_variants_batch(
    batch: VariantBatchInput,
    service: VariantAnalysisService = Depends(get_analysis_service)
):
    """
    Analyze multiple variants in batch

    Args:
        batch: List of variants to analyze

    Returns:
        List of analysis results for each variant
    """
    try:
        if len(batch.variants) > batch.max_variants:
            raise HTTPException(
                status_code=400,
                detail=f"Too many variants. Maximum {batch.max_variants} per batch."
            )

        logger.info(f"Batch analysis: {len(batch.variants)} variants")

        results = []
        errors = []

        for idx, variant in enumerate(batch.variants):
            try:
                variant_dict = variant.dict()
                result = await service.analyze_variant(variant_dict)
                results.append({
                    "index": idx,
                    "variant": variant_dict,
                    "analysis": result,
                    "status": "success"
                })
            except Exception as e:
                logger.error(f"Error analyzing variant {idx}: {str(e)}")
                errors.append({
                    "index": idx,
                    "variant": variant.dict(),
                    "error": str(e),
                    "status": "error"
                })

        return {
            "total_variants": len(batch.variants),
            "successful": len(results),
            "failed": len(errors),
            "results": results,
            "errors": errors
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error in batch analysis: {str(e)}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Error in batch analysis: {str(e)}")


@app.post("/api/v1/literature/search")
async def search_literature(
    search_params: LiteratureSearchInput,
    service: VariantAnalysisService = Depends(get_analysis_service)
):
    """
    Search PubMed for literature related to a gene/variant

    Args:
        search_params: Gene and variant information

    Returns:
        Categorized literature results from PubMed
    """
    try:
        logger.info(f"Literature search for: {search_params.gene}")

        # Create variant dict for search
        variant = {
            "gene": search_params.gene,
            "hgvs_p": search_params.variant_description or ""
        }

        # Search literature
        articles = await service.pubmed_service.search_variant(
            variant,
            max_results=search_params.max_results
        )

        # Categorize articles
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
            "gene": search_params.gene,
            "variant_description": search_params.variant_description,
            "total_articles": len(articles),
            "categories": {
                "functional_studies": {
                    "count": len(functional_studies),
                    "articles": functional_studies
                },
                "case_reports": {
                    "count": len(case_reports),
                    "articles": case_reports
                },
                "reviews": {
                    "count": len(reviews),
                    "articles": reviews
                },
                "other": {
                    "count": len(other),
                    "articles": other
                }
            },
            "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/?term={search_params.gene}+variant"
        }

    except Exception as e:
        logger.error(f"Error searching literature: {str(e)}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Error searching literature: {str(e)}")


@app.get("/api/v1/predictors")
async def list_predictors():
    """
    List all available prediction tools and their descriptions

    Returns:
        Dictionary of available predictors organized by variant type
    """
    return {
        "missense_predictors": {
            "SIFT": {
                "name": "SIFT",
                "description": "Sorting Intolerant From Tolerant - predicts amino acid substitution effects",
                "score_range": "0-1",
                "interpretation": "≤0.05 = deleterious, >0.05 = tolerated",
                "url": "https://sift.bii.a-star.edu.sg/"
            },
            "PolyPhen-2": {
                "name": "PolyPhen-2",
                "description": "Polymorphism Phenotyping v2 - predicts impact on protein structure/function",
                "score_range": "0-1",
                "interpretation": ">0.85 = probably damaging, 0.15-0.85 = possibly damaging, <0.15 = benign",
                "url": "http://genetics.bwh.harvard.edu/pph2/"
            },
            "CADD": {
                "name": "CADD",
                "description": "Combined Annotation Dependent Depletion - integrates multiple annotations",
                "score_range": "Phred-scaled",
                "interpretation": ">20 = top 1% deleteriousness, >30 = top 0.1%",
                "url": "https://cadd.gs.washington.edu/"
            },
            "REVEL": {
                "name": "REVEL",
                "description": "Rare Exome Variant Ensemble Learner - ensemble method for missense variants",
                "score_range": "0-1",
                "interpretation": ">0.5 = likely pathogenic, >0.75 = pathogenic",
                "url": "https://sites.google.com/site/revelgenomics/"
            }
        },
        "splicing_predictors": {
            "SpliceAI": {
                "name": "SpliceAI",
                "description": "Deep learning-based splice prediction",
                "score_range": "Delta scores 0-1",
                "interpretation": ">0.5 = high confidence, 0.2-0.5 = moderate",
                "url": "https://github.com/Illumina/SpliceAI"
            },
            "MaxEntScan": {
                "name": "MaxEntScan",
                "description": "Maximum Entropy Model for splice site strength",
                "score_range": "Varies",
                "interpretation": "Higher score = stronger splice site",
                "url": "http://hollywood.mit.edu/burgelab/maxent/"
            },
            "SPIDEX": {
                "name": "SPIDEX",
                "description": "Deep learning prediction of splice disruption",
                "score_range": "Z-scores",
                "interpretation": "|dpsi_zscore| > 2 suggests significant impact",
                "url": "http://www.openbioinformatics.org/annovar/spidex_download_form.php"
            }
        },
        "indel_analyzers": {
            "IndelAnalyzer": {
                "name": "Indel Analyzer",
                "description": "Classifies insertions/deletions as frameshift or in-frame",
                "interpretation": "Frameshift = likely pathogenic (PVS1), Large inframe = potentially pathogenic (PM4)"
            },
            "CADD-Indel": {
                "name": "CADD for Indels",
                "description": "CADD scoring for insertions and deletions",
                "score_range": "Phred-scaled",
                "interpretation": ">20 = top 1%, >30 = top 0.1%",
                "url": "https://cadd.gs.washington.edu/"
            }
        },
        "synonymous_predictors": {
            "MobiDetails": {
                "name": "MobiDetails",
                "description": "Comprehensive variant annotation with ACMG/AMP classification. Requires NM accession + HGVS coding.",
                "url": "https://mobidetails.iurc.montp.inserm.fr/MD/",
                "api": "https://mobidetails.iurc.montp.inserm.fr/MD/api/",
                "applies_to": "All coding variants"
            },
            "SILVA": {
                "name": "SILVA",
                "description": "Synonymous variant impact analysis: ESE/ESS motifs, codon usage, mRNA structure",
                "url": "https://silva.bioinfo.cipf.es/",
                "applies_to": "SNV_SYNONYMOUS (p.(=))"
            },
            "ESEfinder": {
                "name": "ESEfinder / HEXplorer",
                "description": "Predicts SR protein binding sites (ESE). Disruption may affect splicing of synonymous variants.",
                "url": "http://krainer01.cshl.edu/cgi-bin/tools/ESE3/esefinder.cgi",
                "hexplorer_url": "https://hexplorer.uni-due.de/",
                "applies_to": "SNV_SYNONYMOUS"
            }
        },
        "utr_predictors": {
            "RegulomeDB": {
                "name": "RegulomeDB",
                "description": "Scores regulatory variants based on functional genomic data (TF binding, DNase, eQTL)",
                "score_range": "1a-7 (1a = most functional, 7 = no data)",
                "interpretation": "Score 1-2 = likely functional; 5-7 = unlikely functional",
                "url": "https://regulomedb.org/",
                "applies_to": "UTR_5, UTR_3, non-coding"
            },
            "UTR5Analyzer": {
                "name": "5'UTR Analyzer (UTRAnnotator)",
                "description": "Detects upstream ORF creation/disruption, Kozak sequence changes, translation impact",
                "url": "https://github.com/ImperialCollegeLondon/UTRAnnotator",
                "vep_plugin": "https://www.ensembl.org/vep",
                "applies_to": "UTR_5 (c.-NNN variants)"
            },
            "UTR3Analyzer": {
                "name": "3'UTR Analyzer",
                "description": "Analyzes miRNA seed region disruption, polyadenylation signals, RBP motifs, ARE elements",
                "url_targetscan": "https://www.targetscan.org/",
                "url_mirdb": "https://mirdb.org/",
                "url_rbpmap": "http://rbpmap.technion.ac.il/",
                "url_polyasite": "https://polyasite.unibas.ch/",
                "applies_to": "UTR_3 (c.*NNN variants)"
            }
        },
        "variant_types_supported": {
            "SNV_MISSENSE": "Single nucleotide variant affecting protein sequence",
            "SNV_SYNONYMOUS": "Synonymous (silent) variant - p.(=) notation",
            "SPLICING": "Variant affecting canonical splice sites (±1, ±2)",
            "INDEL_FRAMESHIFT": "Frameshift insertion/deletion",
            "INDEL_INFRAME": "In-frame insertion/deletion",
            "UTR_5": "5'UTR variant (HGVS: c.-NNN)",
            "UTR_3": "3'UTR variant (HGVS: c.*NNN)",
            "COMPLEX": "Complex substitution"
        }
    }


@app.get("/api/v1/databases")
async def list_databases():
    """
    List external databases and resources integrated or referenced

    Returns:
        Dictionary of databases with descriptions and URLs
    """
    return {
        "population_databases": {
            "gnomAD": {
                "name": "Genome Aggregation Database",
                "description": "Population allele frequencies from 125,748 exomes and 15,708 genomes",
                "url": "https://gnomad.broadinstitute.org/"
            },
            "ExAC": {
                "name": "Exome Aggregation Consortium",
                "description": "Predecessor to gnomAD with 60,706 exomes",
                "url": "http://exac.broadinstitute.org/"
            }
        },
        "clinical_databases": {
            "ClinVar": {
                "name": "ClinVar",
                "description": "Public archive of variant interpretations and clinical significance",
                "url": "https://www.ncbi.nlm.nih.gov/clinvar/"
            },
            "HGMD": {
                "name": "Human Gene Mutation Database",
                "description": "Comprehensive collection of disease-causing mutations",
                "url": "http://www.hgmd.cf.ac.uk/"
            },
            "COSMIC": {
                "name": "Catalogue Of Somatic Mutations In Cancer",
                "description": "Somatic mutation database for cancer genetics",
                "url": "https://cancer.sanger.ac.uk/cosmic"
            }
        },
        "literature": {
            "PubMed": {
                "name": "PubMed",
                "description": "NCBI database of biomedical literature",
                "url": "https://pubmed.ncbi.nlm.nih.gov/"
            }
        },
        "gene_resources": {
            "OMIM": {
                "name": "Online Mendelian Inheritance in Man",
                "description": "Comprehensive knowledge base of human genes and genetic phenotypes",
                "url": "https://www.omim.org/"
            },
            "GeneCards": {
                "name": "GeneCards",
                "description": "Integrated database of human genes",
                "url": "https://www.genecards.org/"
            },
            "GTEx": {
                "name": "Genotype-Tissue Expression",
                "description": "Gene expression data across tissues",
                "url": "https://gtexportal.org/"
            }
        }
    }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000, log_level="info")
