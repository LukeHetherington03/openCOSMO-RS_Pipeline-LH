"""
tests/utils/test_melting_point_finder.py

Unit tests for MeltingPointFinder.

Run with:
    cd /home/lunet/cglh4/openCOSMO-RS_Pipeline-LH
    pytest tests/utils/test_melting_point_finder.py -v

Network tests are skipped by default. To include them:
    pytest tests/utils/test_melting_point_finder.py -v -m network
"""

from unittest.mock import MagicMock, patch

import pytest

from modules.utils.melting_point_finder import MeltingPointFinder, MeltingPointResult


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture
def finder():
    return MeltingPointFinder()


# ─────────────────────────────────────────────────────────────────────────────
# _parse_temperature — unit conversion
# ─────────────────────────────────────────────────────────────────────────────

class TestParseTemperature:
    def test_celsius_symbol(self, finder):
        assert abs(finder._parse_temperature("25 °C") - 298.15) < 0.01

    def test_celsius_letter(self, finder):
        assert abs(finder._parse_temperature("100 c") - 373.15) < 0.01

    def test_fahrenheit_symbol(self, finder):
        # 77 °F → 25 °C → 298.15 K
        assert abs(finder._parse_temperature("77 °F") - 298.15) < 0.01

    def test_fahrenheit_letter(self, finder):
        # 32 °F → 0 °C → 273.15 K
        assert abs(finder._parse_temperature("32 f") - 273.15) < 0.01

    def test_kelvin_suffix(self, finder):
        assert abs(finder._parse_temperature("298 K") - 298.0) < 0.01

    def test_kelvin_endswith(self, finder):
        assert abs(finder._parse_temperature("300k") - 300.0) < 0.01

    def test_negative_celsius(self, finder):
        assert abs(finder._parse_temperature("-10 °C") - 263.15) < 0.01

    def test_ambiguous_bare_number_defaults_celsius(self, finder):
        # No unit → assumed Celsius
        result = finder._parse_temperature("20")
        assert abs(result - 293.15) < 0.01

    def test_no_number_returns_none(self, finder):
        assert finder._parse_temperature("no number here") is None

    def test_decimal(self, finder):
        assert abs(finder._parse_temperature("36.6 °C") - 309.75) < 0.01


# ─────────────────────────────────────────────────────────────────────────────
# _choose_best_mp — source selection logic
# ─────────────────────────────────────────────────────────────────────────────

class TestChooseBestMp:
    def test_empty_returns_none(self, finder):
        assert finder._choose_best_mp([]) is None

    def test_picks_value_in_range(self, finder):
        # 50 < T < 400 K filter
        temps = [10.0, 200.0, 500.0]
        assert finder._choose_best_mp(temps) == 200.0

    def test_picks_minimum_in_range(self, finder):
        temps = [150.0, 250.0, 350.0]
        assert finder._choose_best_mp(temps) == 150.0

    def test_all_outside_range_picks_minimum(self, finder):
        temps = [10.0, 20.0, 600.0]
        assert finder._choose_best_mp(temps) == 10.0

    def test_single_value_in_range(self, finder):
        assert finder._choose_best_mp([300.0]) == 300.0

    def test_single_value_outside_range(self, finder):
        assert finder._choose_best_mp([5.0]) == 5.0


# ─────────────────────────────────────────────────────────────────────────────
# PubChem path (mocked HTTP)
# ─────────────────────────────────────────────────────────────────────────────

# Minimal PUG-View JSON stub that contains one melting point entry
_PUBCHEM_RECORD_STUB = {
    "Record": {
        "Section": [
            {
                "TOCHeading": "Chemical and Physical Properties",
                "Section": [
                    {
                        "TOCHeading": "Melting Point",
                        "Information": [
                            {
                                "Value": {
                                    "StringWithMarkup": [
                                        {"String": "135 °C"}
                                    ]
                                }
                            }
                        ],
                    }
                ],
            }
        ]
    }
}

# Aspirin CID
_CID_STUB = {"IdentifierList": {"CID": [2244]}}

# Synonyms stub with a CAS number
_SYN_STUB = {
    "InformationList": {
        "Information": [{"Synonym": ["50-78-2", "aspirin", "acetylsalicylic acid"]}]
    }
}


def _mock_requests_get(url, **kwargs):
    resp = MagicMock()
    resp.ok = True
    if "cids" in url:
        resp.json.return_value = _CID_STUB
    elif "pug_view" in url:
        resp.json.return_value = _PUBCHEM_RECORD_STUB
    elif "synonyms" in url:
        resp.json.return_value = _SYN_STUB
    else:
        resp.json.return_value = {}
    return resp


class TestPubChemMocked:
    @patch("modules.utils.melting_point_finder.requests.get", side_effect=_mock_requests_get)
    def test_extracts_mp_and_converts_to_kelvin(self, mock_get, finder):
        # 135 °C → 408.15 K
        result = finder.get_best("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        assert result.melting_temp is not None
        assert abs(result.melting_temp - 408.15) < 0.5

    @patch("modules.utils.melting_point_finder.requests.get", side_effect=_mock_requests_get)
    def test_source_is_pubchem_or_nist(self, mock_get, finder):
        result = finder.get_best("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        assert result.source in ("pubchem", "nist")

    @patch("modules.utils.melting_point_finder.requests.get", side_effect=_mock_requests_get)
    def test_result_is_melting_point_result_instance(self, mock_get, finder):
        result = finder.get_best("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        assert isinstance(result, MeltingPointResult)


# ─────────────────────────────────────────────────────────────────────────────
# Fallback order
# ─────────────────────────────────────────────────────────────────────────────

class TestFallbackOrder:
    @patch("modules.utils.melting_point_finder.requests.get")
    def test_pubchem_no_mp_returns_none_melting_temp(self, mock_get, finder):
        """PubChem has the compound but no melting point listed."""
        empty_record = {"Record": {"Section": []}}

        def side_effect(url, **kwargs):
            resp = MagicMock()
            resp.ok = True
            if "cids" in url:
                resp.json.return_value = _CID_STUB
            elif "pug_view" in url:
                resp.json.return_value = empty_record
            elif "synonyms" in url:
                resp.json.return_value = _SYN_STUB
            elif "nist" in url or "webbook" in url:
                resp.text = "<html>no tfus here</html>"
            else:
                resp.json.return_value = {}
            return resp

        mock_get.side_effect = side_effect
        result = finder.get_best("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        assert result.melting_temp is None

    @patch("modules.utils.melting_point_finder.requests.get")
    def test_both_sources_fail_returns_none(self, mock_get, finder):
        """CID lookup fails → both sources return None."""
        def side_effect(url, **kwargs):
            resp = MagicMock()
            resp.ok = True
            resp.json.return_value = {}  # no CID found
            return resp

        mock_get.side_effect = side_effect
        result = finder.get_best("AAAAAAAAAA-UHFFFAOYSA-N")
        assert result.melting_temp is None


# ─────────────────────────────────────────────────────────────────────────────
# NIST path (mocked HTTP)
# ─────────────────────────────────────────────────────────────────────────────

class TestNistMocked:
    def test_nist_no_cas_returns_none(self, finder):
        result = finder._nist_get_mp("")
        assert result.melting_temp is None

    @patch("modules.utils.melting_point_finder.requests.get")
    def test_nist_parses_tfus_from_html(self, mock_get, finder):
        html = """
        <html><body>
        <table>
        <tr><td>Tfus</td><td>293.15</td></tr>
        </table>
        </body></html>
        """
        resp = MagicMock()
        resp.ok = True
        resp.text = html
        mock_get.return_value = resp

        result = finder._nist_get_mp("50-78-2")
        assert result.melting_temp is not None
        assert abs(result.melting_temp - 293.15) < 0.1


# ─────────────────────────────────────────────────────────────────────────────
# Live / network tests (skipped by default)
# Run with: pytest tests/utils/test_melting_point_finder.py -v -m network
# ─────────────────────────────────────────────────────────────────────────────

# Molecules confirmed to have MP data in PubChem (verified live).
# Reference values in Kelvin from acr_t22.csv (Chemical Book).
# PubChem reports these in °F — _parse_temperature converts correctly.
_LIVE_MOLECULES_WITH_PUBCHEM_MP = [
    ("ethyl methacrylate",   "SUPCQIBBMFXVTL-UHFFFAOYSA-N", 198.15),  # PubChem: -103 °F
    ("methyl acrylate",      "BAPJBEWLBFYGME-UHFFFAOYSA-N", 198.15),  # PubChem: -105.7 °F
    ("ethyl acrylate",       "JIGUQPWFLRLWPJ-UHFFFAOYSA-N", 202.15),  # PubChem: -96.2 °F
    ("methyl methacrylate",  "VVQNEPGJFQJSBK-UHFFFAOYSA-N", 225.15),  # PubChem: -54 °F
    ("TPO",                  "VFHVQBAGLAREND-UHFFFAOYSA-N", 366.15),  # PubChem: 93 °C
]

# Molecules in PubChem but with no MP data in their records (13 of 18).
# Included so the test suite documents the gap rather than silently skipping.
_LIVE_MOLECULES_NO_PUBCHEM_MP = [
    ("2-(dimethylamino)ethyl acrylate",              "DPBJAVGHACCNRL-UHFFFAOYSA-N"),
    ("2-(2-oxoimidazolidin-1-yl)ethyl methacrylate", "PFPUZMSQZJFLBK-UHFFFAOYSA-N"),
    ("tert-butyl methacrylate",                      "SJMYWORNLPSJQO-UHFFFAOYSA-N"),
    ("pentabromobenzyl acrylate",                    "GRKDVZMVHOLESV-UHFFFAOYSA-N"),
    ("8-methylnonyl methacrylate",                   "COCLLEMEIJQBAG-UHFFFAOYSA-N"),
    ("2-hydroxypropyl acrylate",                     "GWZMWHWAWHPNHN-UHFFFAOYSA-N"),
    ("2,2,2-trifluoroethyl methacrylate",            "QTKPMCIBUROOGY-UHFFFAOYSA-N"),
    ("phenyl methacrylate",                          "QIWKUEJZZCOPFV-UHFFFAOYSA-N"),
    ("2-methylbutyl acrylate",                       "NCTBYWFEJFTVEL-UHFFFAOYSA-N"),
    ("methyl 3-phenylacrylate",                      "CCRCUPLGCSFEDV-BQYQJAHWSA-N"),
    ("ethyl 3-phenylacrylate",                       "KBEBGUQPQBELIU-CMDGGOBGSA-N"),
    ("triazine triacrylate",                         "YIJYFLXQHDOQGW-UHFFFAOYSA-N"),
    ("Irgacure651",                                  "OQCIXXVMXBXDSC-UHFFFAOYSA-N"),
]

_TOLERANCE_K = 5.0  # PubChem NTP values vs Chemical Book may differ by a few K


@pytest.mark.network
@pytest.mark.parametrize(
    "name,inchikey,ref_mp_k",
    _LIVE_MOLECULES_WITH_PUBCHEM_MP,
    ids=[m[0] for m in _LIVE_MOLECULES_WITH_PUBCHEM_MP],
)
def test_live_melting_point(finder, name, inchikey, ref_mp_k):
    """
    Query MeltingPointFinder against real PubChem and check the returned
    melting point is within ±25 K of the reference value from acr_t22.csv.
    Only molecules confirmed to have MP data in PubChem are tested here.
    """
    result = finder.get_best(inchikey)
    assert result.melting_temp is not None, (
        f"{name}: no melting point found (source={result.source}, raw={result.raw})"
    )
    assert abs(result.melting_temp - ref_mp_k) <= _TOLERANCE_K, (
        f"{name}: got {result.melting_temp:.2f} K, expected {ref_mp_k:.2f} K "
        f"(±{_TOLERANCE_K} K), source={result.source}"
    )


@pytest.mark.network
@pytest.mark.parametrize(
    "name,inchikey",
    _LIVE_MOLECULES_NO_PUBCHEM_MP,
    ids=[m[0] for m in _LIVE_MOLECULES_NO_PUBCHEM_MP],
)
def test_live_no_pubchem_mp(finder, name, inchikey):
    """
    These molecules are in PubChem but have no MP in their records.
    Assert melting_temp is None so we detect if PubChem adds data later.
    """
    result = finder.get_best(inchikey)
    assert result.melting_temp is None, (
        f"{name}: expected no MP from PubChem but got {result.melting_temp} K "
        f"(source={result.source}) — PubChem may have added data, update test list"
    )
