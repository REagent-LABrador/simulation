"""Force ligand_filter's default source offline: no paperclip, no retries."""
import ligand_filter as lf
class _Offline(lf.ChemCompSource):
    def _fetch_batch(self, batch, *, attempts=3):
        for c in batch:
            self._cache.setdefault(c, None)
