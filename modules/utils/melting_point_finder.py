import re
from dataclasses import dataclass
from typing import Optional, Dict, List

import requests


@dataclass
class MeltingPointResult:
    melting_temp: Optional[float]  # Kelvin
    source: str
    confidence: float
    raw: Optional[str] = None
    found_but_no_mp: bool = False
    cas: Optional[str] = None


class MeltingPointFinder:
    PUBCHEM_CID_URL = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{}/cids/JSON"
    )
    PUBCHEM_RECORD_URL = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON"
    )
    PUBCHEM_SYN_URL = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/synonyms/JSON"
    )

    NIST_BASE_URL = "https://webbook.nist.gov/cgi/cbook.cgi"

    # ---------------- PubChem: InChIKey → CID ----------------

    def _pubchem_cid_from_inchikey(self, inchikey: str) -> Optional[int]:
        try:
            r = requests.get(self.PUBCHEM_CID_URL.format(inchikey), timeout=10)
            data = r.json()
            cids = data.get("IdentifierList", {}).get("CID", [])
            return cids[0] if cids else None
        except Exception:
            return None

    # ---------------- PubChem: CID → full record ----------------

    def _pubchem_record(self, cid: int) -> Optional[dict]:
        try:
            r = requests.get(self.PUBCHEM_RECORD_URL.format(cid), timeout=10)
            return r.json()
        except Exception:
            return None

    # ---------------- PubChem: CAS extraction ----------------

    def _extract_cas_from_synonyms(self, cid: int) -> Optional[str]:
        try:
            r = requests.get(self.PUBCHEM_SYN_URL.format(cid), timeout=10)
            info = r.json().get("InformationList", {}).get("Information", [])
            if not info:
                return None

            syns = info[0].get("Synonym", [])
            cas_pattern = re.compile(r"\b\d{2,7}-\d{2}-\d\b")

            for s in syns:
                m = cas_pattern.search(s)
                if m:
                    return m.group(0)

        except Exception:
            pass

        return None

    def _extract_cas_from_record(self, record: dict) -> Optional[str]:
        cas_pattern = re.compile(r"\b\d{2,7}-\d{2}-\d\b")

        def walk(node):
            if isinstance(node, dict):
                if "RegistryNumber" in node:
                    rn = node["RegistryNumber"]
                    if isinstance(rn, str):
                        m = cas_pattern.search(rn)
                        if m:
                            return m.group(0)

                if "StringWithMarkup" in node:
                    for entry in node["StringWithMarkup"]:
                        text = entry.get("String", "")
                        m = cas_pattern.search(text)
                        if m:
                            return m.group(0)

                for v in node.values():
                    out = walk(v)
                    if out:
                        return out

            elif isinstance(node, list):
                for item in node:
                    out = walk(item)
                    if out:
                        return out

            return None

        return walk(record)

    # ---------------- PubChem: melting point parsing ----------------

    def _parse_temperature(self, text: str) -> Optional[float]:
        m = re.search(r"-?\d+(\.\d+)?", text)
        if not m:
            return None

        value = float(m.group(0))

        text_lower = text.lower()

        if "°c" in text_lower or " c" in text_lower:
            return value + 273.15

        if "°f" in text_lower or " f" in text_lower:
            return (value - 32) * 5.0 / 9.0 + 273.15

        if " k" in text_lower or text_lower.endswith("k"):
            return value

        # Default: assume Celsius
        return value + 273.15

    def _collect_all_pubchem_mps(self, record) -> List[float]:
        temps = []

        def walk(node):
            if isinstance(node, dict):
                if node.get("TOCHeading") == "Melting Point":
                    for info in node.get("Information", []):
                        value = info.get("Value", {})
                        strings = value.get("StringWithMarkup", [])
                        for s in strings:
                            t = self._parse_temperature(s.get("String", ""))
                            if t is not None:
                                temps.append(t)

                for v in node.values():
                    walk(v)

            elif isinstance(node, list):
                for item in node:
                    walk(item)

        walk(record)
        return temps

    def _choose_best_mp(self, temps: List[float]) -> Optional[float]:
        if not temps:
            return None

        # Filter to realistic melting points
        filtered = [t for t in temps if 50 < t < 400]
        if filtered:
            return min(filtered)

        return min(temps)

    # ---------------- PubChem wrapper ----------------

    def _pubchem_get_mp_and_cas(self, inchikey: str):
        cid = self._pubchem_cid_from_inchikey(inchikey)
        if cid is None:
            return None, None, "no_cid"

        record = self._pubchem_record(cid)
        if not record:
            return None, None, "no_record"

        # CAS from record or synonyms
        cas = self._extract_cas_from_record(record)
        if cas is None:
            cas = self._extract_cas_from_synonyms(cid)

        # Melting points
        temps = self._collect_all_pubchem_mps(record)
        mp_k = self._choose_best_mp(temps)

        return mp_k, cas, "pubchem"

    # ---------------- NIST ----------------

    def _nist_get_mp(self, cas: str) -> MeltingPointResult:
        if not cas:
            return MeltingPointResult(None, "nist", 0.0, "no_cas")

        params = {"ID": cas, "Units": "SI", "cTP": "on"}
        try:
            r = requests.get(self.NIST_BASE_URL, params=params, timeout=10)
            if not r.ok:
                return MeltingPointResult(None, "nist", 0.0, "http_error")

            html = r.text

            m = re.search(r"Tfus.*?<td[^>]*>(.*?)</td>", html, re.DOTALL | re.IGNORECASE)
            if not m:
                return MeltingPointResult(None, "nist", 0.0, "no_tfus", found_but_no_mp=True)

            cell = m.group(1)
            num = re.search(r"-?\d+(\.\d+)?", cell)
            if not num:
                return MeltingPointResult(None, "nist", 0.0, "no_number", found_but_no_mp=True)

            mp_k = float(num.group(0))
            return MeltingPointResult(mp_k, "nist", 1.0, "nist")

        except Exception as e:
            return MeltingPointResult(None, "nist", 0.0, str(e))

    # ---------------- Public API ----------------

    def get_all_sources(self, inchikey: str) -> Dict[str, MeltingPointResult]:
        inchikey = inchikey.strip().upper()

        # PubChem
        mp_pubchem, cas, raw = self._pubchem_get_mp_and_cas(inchikey)
        if mp_pubchem is None:
            pubchem_res = MeltingPointResult(None, "pubchem", 0.0, raw, found_but_no_mp=True)
        else:
            pubchem_res = MeltingPointResult(mp_pubchem, "pubchem", 0.6, raw)

        pubchem_res.cas = cas

        # NIST
        nist_res = self._nist_get_mp(cas)
        nist_res.cas = cas

        return {
            "pubchem": pubchem_res,
            "nist": nist_res,
        }

    def get_best(self, inchikey: str) -> MeltingPointResult:
        all_res = self.get_all_sources(inchikey)
        return max(all_res.values(), key=lambda r: r.confidence)


if __name__ == "__main__":
    finder = MeltingPointFinder()

    test_inchikey = "BAPJBEWLBFYGME-UHFFFAOYSA-N"  # methyl acrylate

    results = finder.get_all_sources(test_inchikey)
    print(f"All sources for {test_inchikey}:")
    for name, res in results.items():
        print(
            f"{name:8s} -> mp={res.melting_temp}, cas={res.cas}, "
            f"conf={res.confidence}, raw={res.raw}, found_but_no_mp={res.found_but_no_mp}"
        )

    best = finder.get_best(test_inchikey)
    print("\nBest:")
    print(best)
