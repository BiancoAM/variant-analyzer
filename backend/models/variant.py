from sqlalchemy import Column, Integer, String, Float, DateTime, JSON, Text
from sqlalchemy.ext.declarative import declarative_base
from datetime import datetime

Base = declarative_base()


class Variant(Base):
    __tablename__ = 'variants'

    id = Column(Integer, primary_key=True, index=True)
    chromosome = Column(String(10), nullable=False, index=True)
    position = Column(Integer, nullable=False, index=True)
    reference = Column(String(500), nullable=False)
    alternate = Column(String(500), nullable=False)
    gene = Column(String(50), index=True)
    transcript = Column(String(50))
    hgvs_c = Column(String(200))  # HGVS coding
    hgvs_p = Column(String(200))  # HGVS protein
    variant_type = Column(String(50))  # SNV, INDEL, SPLICING, etc.
    created_at = Column(DateTime, default=datetime.utcnow)


class PredictionResult(Base):
    __tablename__ = 'prediction_results'

    id = Column(Integer, primary_key=True, index=True)
    variant_id = Column(Integer, index=True)

    # SNV/Missense Predictors
    sift_score = Column(Float)
    sift_prediction = Column(String(50))
    polyphen2_score = Column(Float)
    polyphen2_prediction = Column(String(50))
    cadd_phred = Column(Float)
    revel_score = Column(Float)

    # Splicing Predictors
    spliceai_score = Column(Float)
    spliceai_delta_score = Column(JSON)  # Store all delta scores
    maxentscan_ref = Column(Float)
    maxentscan_alt = Column(Float)

    # Population Frequencies
    gnomad_af = Column(Float)
    gnomad_af_popmax = Column(Float)

    # Clinical Annotations
    clinvar_significance = Column(String(100))
    clinvar_review_status = Column(String(100))

    # Raw results storage
    raw_results = Column(JSON)

    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)


class LiteratureEvidence(Base):
    __tablename__ = 'literature_evidence'

    id = Column(Integer, primary_key=True, index=True)
    variant_id = Column(Integer, index=True)

    # PubMed Info
    pmid = Column(String(20), index=True)
    title = Column(Text)
    abstract = Column(Text)
    authors = Column(Text)
    publication_date = Column(String(50))
    journal = Column(String(200))

    # Evidence Classification
    evidence_type = Column(String(50))  # functional_study, case_report, etc.
    evidence_level = Column(String(20))  # strong, moderate, supporting

    # Study details
    study_type = Column(String(100))  # in vitro, in vivo, clinical, etc.
    pathogenicity_support = Column(String(50))  # pathogenic, benign, uncertain

    added_at = Column(DateTime, default=datetime.utcnow)
