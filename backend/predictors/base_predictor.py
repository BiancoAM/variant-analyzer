from abc import ABC, abstractmethod
from typing import Dict, Any, Optional
import aiohttp
import asyncio
from datetime import datetime, timedelta


class CacheEntry:
    def __init__(self, data: Any, ttl: int = 3600):
        self.data = data
        self.timestamp = datetime.utcnow()
        self.ttl = ttl

    def is_valid(self) -> bool:
        return datetime.utcnow() - self.timestamp < timedelta(seconds=self.ttl)


class BasePredictor(ABC):
    """Base class for all variant prediction tools"""

    def __init__(self, api_url: Optional[str] = None, cache_ttl: int = 3600):
        self.api_url = api_url
        self.cache: Dict[str, CacheEntry] = {}
        self.cache_ttl = cache_ttl

    def _cache_key(self, **kwargs) -> str:
        """Generate cache key from parameters"""
        return f"{self.__class__.__name__}:" + ":".join(f"{k}={v}" for k, v in sorted(kwargs.items()))

    def _get_cached(self, key: str) -> Optional[Any]:
        """Get cached result if valid"""
        if key in self.cache:
            entry = self.cache[key]
            if entry.is_valid():
                return entry.data
            else:
                del self.cache[key]
        return None

    def _set_cache(self, key: str, data: Any):
        """Store result in cache"""
        self.cache[key] = CacheEntry(data, self.cache_ttl)

    @abstractmethod
    async def predict(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        """
        Predict effect of variant

        Args:
            variant: Dictionary containing variant information
                - chromosome: str
                - position: int
                - reference: str
                - alternate: str
                - gene: str (optional)
                - transcript: str (optional)

        Returns:
            Dictionary with prediction results
        """
        pass

    async def _make_request(self, url: str, method: str = "GET", **kwargs) -> Dict[str, Any]:
        """Make async HTTP request with error handling"""
        try:
            async with aiohttp.ClientSession() as session:
                if method == "GET":
                    async with session.get(url, **kwargs) as response:
                        response.raise_for_status()
                        return await response.json()
                elif method == "POST":
                    async with session.post(url, **kwargs) as response:
                        response.raise_for_status()
                        return await response.json()
        except aiohttp.ClientError as e:
            return {"error": str(e), "success": False}
        except Exception as e:
            return {"error": f"Unexpected error: {str(e)}", "success": False}

    def format_variant_string(self, variant: Dict[str, Any], format_type: str = "standard") -> str:
        """Format variant in various string representations"""
        chr = variant.get('chromosome', '').replace('chr', '')
        pos = variant.get('position')
        ref = variant.get('reference')
        alt = variant.get('alternate')

        if format_type == "standard":
            return f"{chr}-{pos}-{ref}-{alt}"
        elif format_type == "vcf":
            return f"{chr}\t{pos}\t.\t{ref}\t{alt}"
        elif format_type == "hgvs":
            # Simplified HGVS - in real scenario would need proper HGVS construction
            if len(ref) == 1 and len(alt) == 1:
                return f"chr{chr}:g.{pos}{ref}>{alt}"
            else:
                return f"chr{chr}:g.{pos}_{pos+len(ref)-1}delins{alt}"
        return f"{chr}:{pos}{ref}>{alt}"
