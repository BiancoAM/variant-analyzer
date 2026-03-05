from typing import Dict, Any, List, Optional
from Bio import Entrez
import asyncio
import logging
from datetime import datetime
import re

logger = logging.getLogger(__name__)


class PubMedService:
    """
    Service for searching and retrieving literature from PubMed
    Uses NCBI Entrez API via BioPython
    """

    def __init__(self, email: str, api_key: Optional[str] = None):
        """
        Initialize PubMed service

        Args:
            email: Email for NCBI Entrez (required by NCBI)
            api_key: NCBI API key (optional, increases rate limits)
        """
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.cache = {}

    async def search_variant(
        self,
        variant: Dict[str, Any],
        max_results: int = 20
    ) -> List[Dict[str, Any]]:
        """
        Search PubMed for articles related to a specific variant

        Args:
            variant: Variant information dictionary
            max_results: Maximum number of results to return

        Returns:
            List of article dictionaries with PubMed data
        """
        # Build search query
        queries = self._build_search_queries(variant)

        all_articles = []
        seen_pmids = set()

        for query in queries:
            try:
                articles = await self._search_pubmed(query, max_results)

                # Deduplicate
                for article in articles:
                    pmid = article.get('pmid')
                    if pmid and pmid not in seen_pmids:
                        seen_pmids.add(pmid)
                        all_articles.append(article)

            except Exception as e:
                logger.error(f"Error searching PubMed with query '{query}': {str(e)}")

        # Sort by relevance and date
        all_articles.sort(
            key=lambda x: (
                self._calculate_relevance(x, variant),
                x.get('pub_date', '')
            ),
            reverse=True
        )

        return all_articles[:max_results]

    def _build_search_queries(self, variant: Dict[str, Any]) -> List[str]:
        """
        Build multiple search queries for the variant

        Returns queries ordered by specificity (most specific first)
        """
        gene = variant.get('gene') or ''
        hgvs_p = variant.get('hgvs_p') or ''
        hgvs_c = variant.get('hgvs_c') or ''

        queries = []

        # Most specific: gene + protein change
        if gene and hgvs_p:
            protein_change = self._extract_protein_change(hgvs_p)
            if protein_change:
                queries.append(f'{gene}[Gene] AND {protein_change}')
                queries.append(f'{gene} AND {protein_change} AND (pathogenic OR mutation OR variant)')

        # Gene + coding change
        if gene and hgvs_c:
            coding_change = self._extract_coding_change(hgvs_c)
            if coding_change:
                queries.append(f'{gene}[Gene] AND {coding_change}')

        # Gene + general variant terms
        if gene:
            queries.append(
                f'{gene}[Gene] AND (pathogenic OR pathogenicity OR "functional study" OR '
                f'"in vitro" OR "in vivo" OR "patient-derived")'
            )

        # Genomic coordinates
        chr_pos = f"chr{variant.get('chromosome', '')}:{variant.get('position', '')}"
        if chr_pos != "chr:":
            queries.append(f'"{chr_pos}" AND variant')

        return queries

    def _extract_protein_change(self, hgvs_p: str) -> Optional[str]:
        """Extract simple protein change notation (e.g., p.Arg123His -> R123H)"""
        # Match pattern like p.Arg123His
        match = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs_p)
        if match:
            aa_from = self._three_to_one(match.group(1))
            position = match.group(2)
            aa_to = self._three_to_one(match.group(3))
            return f'{aa_from}{position}{aa_to}'
        return None

    def _extract_coding_change(self, hgvs_c: str) -> Optional[str]:
        """Extract coding change notation"""
        # Simple extraction - in production would need more robust parsing
        match = re.search(r'c\.(.+)', hgvs_c)
        if match:
            return match.group(1)
        return None

    def _three_to_one(self, aa_three: str) -> str:
        """Convert three-letter amino acid code to one-letter"""
        aa_map = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
            'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
            'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
            'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
            'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
            'Ter': '*', 'Sec': 'U'
        }
        return aa_map.get(aa_three, aa_three)

    async def _search_pubmed(self, query: str, max_results: int) -> List[Dict[str, Any]]:
        """
        Execute PubMed search and fetch article details

        Args:
            query: PubMed search query
            max_results: Maximum results to return

        Returns:
            List of article dictionaries
        """
        # Run Entrez operations in thread pool (they're blocking)
        loop = asyncio.get_event_loop()

        # Search for PMIDs
        search_handle = await loop.run_in_executor(
            None,
            lambda: Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=max_results,
                sort="relevance"
            )
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()

        pmids = search_results.get('IdList', [])

        if not pmids:
            return []

        # Fetch article details
        fetch_handle = await loop.run_in_executor(
            None,
            lambda: Entrez.efetch(
                db="pubmed",
                id=','.join(pmids),
                rettype="medline",
                retmode="xml"
            )
        )
        records = Entrez.read(fetch_handle)
        fetch_handle.close()

        # Parse articles
        articles = []
        for record in records.get('PubmedArticle', []):
            article = self._parse_article(record)
            if article:
                articles.append(article)

        return articles

    def _parse_article(self, record: Dict) -> Optional[Dict[str, Any]]:
        """Parse PubMed XML record into structured dictionary"""
        try:
            medline_citation = record.get('MedlineCitation', {})
            article = medline_citation.get('Article', {})

            pmid = str(medline_citation.get('PMID', ''))

            # Title
            title = article.get('ArticleTitle', '')

            # Abstract
            abstract_texts = article.get('Abstract', {}).get('AbstractText', [])
            abstract = ' '.join(str(text) for text in abstract_texts) if abstract_texts else ''

            # Authors
            author_list = article.get('AuthorList', [])
            authors = []
            for author in author_list[:10]:  # Limit to first 10 authors
                last_name = author.get('LastName', '')
                initials = author.get('Initials', '')
                if last_name:
                    authors.append(f"{last_name} {initials}".strip())

            authors_str = ', '.join(authors)
            if len(author_list) > 10:
                authors_str += ' et al.'

            # Journal
            journal = article.get('Journal', {})
            journal_title = journal.get('Title', '')

            # Publication date
            pub_date = journal.get('JournalIssue', {}).get('PubDate', {})
            year = pub_date.get('Year', '')
            month = pub_date.get('Month', '')
            pub_date_str = f"{year}-{month}" if year and month else year

            # Keywords and MeSH terms
            keywords = []
            keyword_list = medline_citation.get('KeywordList', [])
            for kw_list in keyword_list:
                keywords.extend([str(kw) for kw in kw_list])

            mesh_terms = []
            mesh_list = medline_citation.get('MeshHeadingList', [])
            for mesh in mesh_list:
                descriptor = mesh.get('DescriptorName', '')
                if descriptor:
                    mesh_terms.append(str(descriptor))

            # Classify evidence type
            evidence_type = self._classify_evidence(title, abstract, keywords, mesh_terms)

            return {
                'pmid': pmid,
                'title': title,
                'abstract': abstract,
                'authors': authors_str,
                'journal': journal_title,
                'pub_date': pub_date_str,
                'keywords': keywords,
                'mesh_terms': mesh_terms,
                'evidence_type': evidence_type,
                'pubmed_url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                'relevance_score': 0  # Will be calculated separately
            }

        except Exception as e:
            logger.error(f"Error parsing article: {str(e)}")
            return None

    def _classify_evidence(
        self,
        title: str,
        abstract: str,
        keywords: List[str],
        mesh_terms: List[str]
    ) -> str:
        """
        Classify the type of evidence provided by the article

        Returns:
            Evidence type: functional_study, case_report, population_study, review, etc.
        """
        text = f"{title} {abstract}".lower()
        all_terms = ' '.join(keywords + mesh_terms).lower()

        # Functional studies
        functional_keywords = [
            'functional study', 'in vitro', 'in vivo', 'cell line',
            'luciferase', 'western blot', 'expression analysis',
            'protein function', 'enzyme activity', 'functional assay',
            'mutagenesis', 'transfection', 'minigene'
        ]
        if any(kw in text or kw in all_terms for kw in functional_keywords):
            return 'functional_study'

        # Case reports / clinical studies
        case_keywords = [
            'case report', 'patient', 'clinical', 'phenotype',
            'case series', 'proband', 'affected individual'
        ]
        if any(kw in text or kw in all_terms for kw in case_keywords):
            return 'case_report'

        # Population / cohort studies
        population_keywords = [
            'cohort', 'population', 'gwas', 'exome sequencing',
            'prevalence', 'frequency', 'screening'
        ]
        if any(kw in text or kw in all_terms for kw in population_keywords):
            return 'population_study'

        # Reviews
        if 'review' in text or 'Review' in mesh_terms:
            return 'review'

        # Structural studies
        structural_keywords = [
            'crystal structure', 'protein structure', 'molecular dynamics',
            'structural analysis', 'nmr'
        ]
        if any(kw in text or kw in all_terms for kw in structural_keywords):
            return 'structural_study'

        return 'other'

    def _calculate_relevance(self, article: Dict[str, Any], variant: Dict[str, Any]) -> float:
        """
        Calculate relevance score for article based on variant

        Returns:
            Relevance score (0-1, higher = more relevant)
        """
        score = 0.0

        text = f"{article.get('title', '')} {article.get('abstract', '')}".lower()
        gene = (variant.get('gene') or '').lower()
        hgvs_p = (variant.get('hgvs_p') or '').lower()

        # Gene mentioned
        if gene and gene in text:
            score += 0.3

        # Specific variant mentioned
        if hgvs_p:
            protein_change = self._extract_protein_change(hgvs_p)
            if protein_change and protein_change.lower() in text:
                score += 0.5

        # Evidence type bonuses
        evidence_type = article.get('evidence_type', '')
        if evidence_type == 'functional_study':
            score += 0.2
        elif evidence_type == 'case_report':
            score += 0.15

        # Pathogenicity keywords
        pathogenicity_keywords = [
            'pathogenic', 'pathogenicity', 'deleterious', 'disease-causing',
            'damaging', 'loss of function', 'gain of function'
        ]
        if any(kw in text for kw in pathogenicity_keywords):
            score += 0.1

        return min(score, 1.0)

    async def get_article_details(self, pmid: str) -> Optional[Dict[str, Any]]:
        """
        Fetch detailed information for a specific PMID

        Args:
            pmid: PubMed ID

        Returns:
            Article dictionary or None if not found
        """
        loop = asyncio.get_event_loop()

        try:
            fetch_handle = await loop.run_in_executor(
                None,
                lambda: Entrez.efetch(
                    db="pubmed",
                    id=pmid,
                    rettype="medline",
                    retmode="xml"
                )
            )
            records = Entrez.read(fetch_handle)
            fetch_handle.close()

            if records.get('PubmedArticle'):
                return self._parse_article(records['PubmedArticle'][0])

        except Exception as e:
            logger.error(f"Error fetching article {pmid}: {str(e)}")

        return None

    async def search_functional_studies(
        self,
        gene: str,
        max_results: int = 10
    ) -> List[Dict[str, Any]]:
        """
        Search specifically for functional studies of a gene

        Args:
            gene: Gene symbol
            max_results: Maximum results to return

        Returns:
            List of functional study articles
        """
        query = (
            f'{gene}[Gene] AND ('
            f'"functional study"[Title/Abstract] OR '
            f'"in vitro"[Title/Abstract] OR '
            f'"in vivo"[Title/Abstract] OR '
            f'"functional analysis"[Title/Abstract] OR '
            f'"functional characterization"[Title/Abstract]'
            f')'
        )

        articles = await self._search_pubmed(query, max_results)

        # Filter to only functional studies
        functional_articles = [
            article for article in articles
            if article.get('evidence_type') == 'functional_study'
        ]

        return functional_articles
