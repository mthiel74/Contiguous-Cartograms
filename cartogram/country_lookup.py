"""Resolve arbitrary country keys to a canonical identifier.

Users hand `world_cartogram` dicts like ``{"USA": 1.0, "UK": 2.0,
"United States": 3.0}`` and expect them to work even when the
underlying GeoDataFrame keys on a specific column. This module
normalises an input key (ISO-3166-1 alpha-3 / alpha-2, ISO numeric,
Natural Earth ADM0_A3, Wolfram-style canonical name, common English
display name) to a canonical sovereign 3-letter code.

The table here is intentionally small and covers the cases users hit
in practice -- full coverage is tracked against the ISO-3166-1 list.
Additions welcome: add a line to ``_ALIASES`` below.
"""

from __future__ import annotations

from typing import Dict, Iterable, Mapping, Optional

# --- canonical target identifier ------------------------------------------
# Every entry below is a many->one map from a user-supplied key to the
# ISO-3166-1 alpha-3 code (matching Natural Earth ISO_A3 / ADM0_A3 /
# SOV_A3 for most countries). The canonical form is ALWAYS uppercase
# alpha-3 with no whitespace.

_ALIASES: Dict[str, str] = {}


def _add(canonical: str, *keys: str) -> None:
    for k in keys:
        _ALIASES[_normalise(k)] = canonical


def _normalise(key: str) -> str:
    return "".join(key.strip().upper().split())


# United States -----------------------------------------------------------
_add("USA",
     "USA", "US", "U.S.", "U.S.A.", "U.S.A",
     "840",                       # ISO numeric
     "UnitedStates",              # Wolfram-style
     "United States",
     "United States of America",
     "USofA", "America")

# United Kingdom ----------------------------------------------------------
_add("GBR",
     "GBR", "GB", "UK", "U.K.",
     "826",
     "UnitedKingdom",
     "United Kingdom",
     "Great Britain", "Britain",
     "UnitedKingdomofGreatBritain",
     "England")   # conventional even though GB covers all four nations

# Common European / G20 members ------------------------------------------
_add("DEU", "DEU", "DE", "276", "Germany")
_add("FRA", "FRA", "FR", "250", "France", "FrenchRepublic")
_add("ITA", "ITA", "IT", "380", "Italy", "ItalianRepublic")
_add("ESP", "ESP", "ES", "724", "Spain")
_add("PRT", "PRT", "PT", "620", "Portugal")
_add("NLD", "NLD", "NL", "528", "Netherlands", "Holland")
_add("BEL", "BEL", "BE", "056", "Belgium")
_add("CHE", "CHE", "CH", "756", "Switzerland")
_add("AUT", "AUT", "AT", "040", "Austria")
_add("IRL", "IRL", "IE", "372", "Ireland")
_add("NOR", "NOR", "NO", "578", "Norway")
_add("SWE", "SWE", "SE", "752", "Sweden")
_add("DNK", "DNK", "DK", "208", "Denmark")
_add("FIN", "FIN", "FI", "246", "Finland")
_add("POL", "POL", "PL", "616", "Poland")
_add("CZE", "CZE", "CZ", "203", "Czechia", "CzechRepublic")
_add("HUN", "HUN", "HU", "348", "Hungary")
_add("GRC", "GRC", "GR", "300", "Greece")
_add("ROU", "ROU", "RO", "642", "Romania")
_add("RUS", "RUS", "RU", "643", "Russia", "RussianFederation")
_add("UKR", "UKR", "UA", "804", "Ukraine")
_add("TUR", "TUR", "TR", "792", "Turkey", "Türkiye", "Turkiye")

# Major Asian economies ---------------------------------------------------
_add("CHN", "CHN", "CN", "156", "China", "PeoplesRepublicofChina")
_add("JPN", "JPN", "JP", "392", "Japan")
_add("KOR", "KOR", "KR", "410",
     "SouthKorea", "Korea, South", "RepublicofKorea", "Korea")
_add("PRK", "PRK", "KP", "408",
     "NorthKorea", "Korea, North",
     "DemocraticPeoplesRepublicofKorea")
_add("IND", "IND", "IN", "356", "India")
_add("IDN", "IDN", "ID", "360", "Indonesia")
_add("VNM", "VNM", "VN", "704", "Vietnam", "VietNam")
_add("THA", "THA", "TH", "764", "Thailand")
_add("PHL", "PHL", "PH", "608", "Philippines")
_add("PAK", "PAK", "PK", "586", "Pakistan")
_add("BGD", "BGD", "BD", "050", "Bangladesh")
_add("MYS", "MYS", "MY", "458", "Malaysia")
_add("SGP", "SGP", "SG", "702", "Singapore")
_add("SAU", "SAU", "SA", "682", "SaudiArabia", "Saudi Arabia")
_add("ARE", "ARE", "AE", "784", "UnitedArabEmirates", "UAE")
_add("ISR", "ISR", "IL", "376", "Israel")
_add("IRN", "IRN", "IR", "364", "Iran", "IslamicRepublicofIran")
_add("IRQ", "IRQ", "IQ", "368", "Iraq")

# Americas ----------------------------------------------------------------
_add("CAN", "CAN", "CA", "124", "Canada")
_add("MEX", "MEX", "MX", "484", "Mexico")
_add("BRA", "BRA", "BR", "076", "Brazil")
_add("ARG", "ARG", "AR", "032", "Argentina")
_add("CHL", "CHL", "CL", "152", "Chile")
_add("COL", "COL", "CO", "170", "Colombia")
_add("VEN", "VEN", "VE", "862", "Venezuela",
     "BolivarianRepublicofVenezuela")
_add("PER", "PER", "PE", "604", "Peru")
_add("CUB", "CUB", "CU", "192", "Cuba")

# Oceania -----------------------------------------------------------------
_add("AUS", "AUS", "AU", "036", "Australia")
_add("NZL", "NZL", "NZ", "554", "NewZealand", "New Zealand")

# Africa (biggest by population / GDP) -----------------------------------
_add("NGA", "NGA", "NG", "566", "Nigeria")
_add("EGY", "EGY", "EG", "818", "Egypt", "ArabRepublicofEgypt")
_add("ETH", "ETH", "ET", "231", "Ethiopia")
_add("COD", "COD", "CD", "180",
     "DRC", "DemocraticRepublicoftheCongo",
     "CongoDRC", "DRCongo")
_add("ZAF", "ZAF", "ZA", "710", "SouthAfrica", "South Africa")
_add("KEN", "KEN", "KE", "404", "Kenya")
_add("TZA", "TZA", "TZ", "834",
     "Tanzania", "UnitedRepublicofTanzania")
_add("MAR", "MAR", "MA", "504", "Morocco")
_add("DZA", "DZA", "DZ", "012", "Algeria")


def canonical(key: str) -> Optional[str]:
    """Return the canonical ISO-3 for ``key``, or ``None`` if unknown.

    Lookup is case-insensitive and whitespace-insensitive. Examples::

        canonical("us")                -> "USA"
        canonical("United States")     -> "USA"
        canonical("UnitedStates")      -> "USA"   # Wolfram-style
        canonical("840")               -> "USA"
        canonical("UK")                -> "GBR"
        canonical("  Great Britain ")  -> "GBR"
    """
    if not isinstance(key, str):
        return None
    return _ALIASES.get(_normalise(key))


def normalise_mapping(
    values: Mapping[str, float],
    *,
    default_to_key: bool = True,
) -> Dict[str, float]:
    """Return ``{ISO3: value}`` from a user-supplied mapping.

    Keys that cannot be resolved fall back to their uppercase, stripped
    form if ``default_to_key`` is True (so unknown but already-canonical
    codes still work) or are dropped otherwise.
    """
    out: Dict[str, float] = {}
    for k, v in values.items():
        c = canonical(k)
        if c is None:
            if default_to_key:
                c = _normalise(k)
            else:
                continue
        out[c] = float(v)
    return out


def known_aliases() -> Iterable[str]:
    """All aliases the module currently recognises (handy for docs)."""
    return _ALIASES.keys()
