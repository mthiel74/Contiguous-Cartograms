"""Public-API fetchers for per-country indicator data.

This module turns one-line calls such as

    fetch("worldbank", "SP.POP.TOTL", year=2022)

into a ``{iso3: value}`` dictionary that can be fed directly to
:meth:`cartogram.WorldMap.with_values`.

Two sources are built in:

* **World Bank Open Data** — ``https://api.worldbank.org/v2/...``. No
  key required. Coverage: population, GDP, CO₂, health, education,
  environment, trade, finance, … thousands of indicators.
* **WHO Global Health Observatory** — ``https://ghoapi.azureedge.net``.
  No key required. Coverage: life expectancy, mortality, risk
  factors, health-workforce metrics, vaccine coverage, …

Responses are cached on disk under ``~/.cache/cartogram/data/`` so the
GUI can refresh cartograms quickly and offline after first fetch.

A curated catalogue of ~20 popular indicators lives in
:data:`INDICATOR_CATALOGUE` and drives the GUI dropdown.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

CACHE_DIR = Path.home() / ".cache" / "cartogram" / "data"
HTTP_TIMEOUT = 30  # seconds


# ----------------------------------------------------------------------
@dataclass(frozen=True)
class Indicator:
    source: str          # 'worldbank' or 'who'
    code: str            # provider's indicator code
    label: str           # human-readable name
    unit: str = ""       # displayed alongside the value
    log_recommended: bool = False  # true for heavy-tailed quantities
    default_year: Optional[int] = None


#: Curated catalogue used by the GUI. Extend freely.
INDICATOR_CATALOGUE: list[Indicator] = [
    # --- World Bank ---
    Indicator("worldbank", "SP.POP.TOTL",       "Population (total)",               "people", False, 2022),
    Indicator("worldbank", "NY.GDP.MKTP.CD",    "GDP (current US$)",                "USD",    True,  2022),
    Indicator("worldbank", "NY.GDP.PCAP.CD",    "GDP per capita (current US$)",     "USD",    True,  2022),
    Indicator("worldbank", "SP.DYN.LE00.IN",    "Life expectancy at birth",         "years",  False, 2022),
    Indicator("worldbank", "SP.URB.TOTL",       "Urban population",                 "people", True,  2022),
    Indicator("worldbank", "EN.ATM.CO2E.KT",    "CO₂ emissions",                    "kt",     True,  2020),
    Indicator("worldbank", "EN.ATM.CO2E.PC",    "CO₂ emissions per capita",         "t",      True,  2020),
    Indicator("worldbank", "MS.MIL.XPND.CD",    "Military expenditure",             "USD",    True,  2022),
    Indicator("worldbank", "IT.NET.USER.ZS",    "Internet users",                   "% pop.", False, 2022),
    Indicator("worldbank", "SP.DYN.IMRT.IN",    "Infant mortality rate",            "per 1k", True,  2022),
    Indicator("worldbank", "AG.LND.FRST.K2",    "Forest area",                      "km²",    True,  2021),
    Indicator("worldbank", "EG.USE.ELEC.KH.PC", "Electric power consumption / cap", "kWh",    True,  2014),
    Indicator("worldbank", "ST.INT.ARVL",       "Intl. tourist arrivals",           "people", True,  2020),
    Indicator("worldbank", "SH.DYN.MORT",       "Under-5 mortality rate",           "per 1k", True,  2022),
    Indicator("worldbank", "SE.XPD.TOTL.GD.ZS", "Government education expenditure", "% GDP",  False, 2022),
    # --- WHO GHO ---
    Indicator("who",       "WHOSIS_000001",     "Life expectancy at birth (WHO)",   "years",  False, 2021),
    Indicator("who",       "WHOSIS_000015",     "Healthy life expectancy",          "years",  False, 2021),
    Indicator("who",       "NCD_BMI_30A",       "Adult obesity rate",               "% 18+",  False, 2016),
    Indicator("who",       "SA_0000001462",     "Alcohol consumption / capita",     "L pure", True,  2019),
    Indicator("who",       "MDG_0000000007",    "Malaria incidence",                "per 1k", True,  2021),
]


def catalogue_labels() -> list[str]:
    return [f"{i.label} — {i.source}:{i.code}" for i in INDICATOR_CATALOGUE]


def indicator_by_label(label: str) -> Indicator:
    for i in INDICATOR_CATALOGUE:
        if f"{i.label} — {i.source}:{i.code}" == label:
            return i
    raise KeyError(label)


# ----------------------------------------------------------------------
def _cache_path(key: str) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    h = hashlib.sha256(key.encode("utf-8")).hexdigest()[:20]
    return CACHE_DIR / f"{h}.json"


def _cached_get(url: str, ttl: float = 7 * 24 * 3600) -> str:
    """GET ``url``, caching the body on disk for ``ttl`` seconds."""
    cache = _cache_path(url)
    if cache.exists() and (time.time() - cache.stat().st_mtime) < ttl:
        return cache.read_text(encoding="utf-8")
    req = urllib.request.Request(
        url, headers={"User-Agent": "cartogram/0.1 (+https://github.com)"}
    )
    with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT) as resp:
        body = resp.read().decode("utf-8")
    cache.write_text(body, encoding="utf-8")
    return body


# ----------------------------------------------------------------------
def fetch_worldbank(code: str, year: Optional[int] = None) -> Dict[str, float]:
    """Return ``{ISO3: value}`` for a World Bank indicator.

    If ``year`` is None, the latest year with any non-null values is used.
    """
    date_clause = f"&date={year}" if year else ""
    url = (
        f"https://api.worldbank.org/v2/country/all/indicator/{code}"
        f"?format=json&per_page=20000{date_clause}"
    )
    body = _cached_get(url)
    payload = json.loads(body)
    if not isinstance(payload, list) or len(payload) < 2 or payload[1] is None:
        raise RuntimeError(
            f"World Bank indicator {code!r} returned no rows (year={year})"
        )

    rows = payload[1]
    if year is None:
        # Find the most recent year with the densest coverage.
        by_year: dict[str, int] = {}
        for row in rows:
            if row.get("value") is None:
                continue
            by_year[row["date"]] = by_year.get(row["date"], 0) + 1
        if not by_year:
            raise RuntimeError(f"no values found for {code!r}")
        best_year = max(by_year.items(), key=lambda kv: (kv[1], kv[0]))[0]
        rows = [r for r in rows if r["date"] == best_year]

    out: Dict[str, float] = {}
    for row in rows:
        v = row.get("value")
        if v is None:
            continue
        country = row.get("countryiso3code") or row.get("country", {}).get("id")
        if not country or len(country) != 3:
            continue
        try:
            out[country] = float(v)
        except (TypeError, ValueError):
            continue
    if not out:
        raise RuntimeError(
            f"World Bank indicator {code!r} fetched but yielded no ISO3 values"
        )
    return out


# ----------------------------------------------------------------------
def fetch_who(code: str, year: Optional[int] = None) -> Dict[str, float]:
    """Return ``{ISO3: value}`` for a WHO GHO indicator.

    WHO does not always have ISO3 in the ``SpatialDim`` field — it uses
    ``SpatialDimType == "COUNTRY"`` and the country code in
    ``SpatialDim``. We filter accordingly. If multiple rows per country
    exist (age/sex breakdowns), we take the aggregate row (usually
    ``Dim1 == "BTSX"`` and no age filter).
    """
    url = f"https://ghoapi.azureedge.net/api/{code}"
    body = _cached_get(url)
    payload = json.loads(body)
    rows = payload.get("value") or []
    if not rows:
        raise RuntimeError(f"WHO indicator {code!r} returned no rows")

    # Keep only country-level rows.
    rows = [r for r in rows if r.get("SpatialDimType") == "COUNTRY"]
    # Prefer the total (both-sex, all-ages) disaggregation.
    preferred = [
        r for r in rows
        if (r.get("Dim1") in (None, "BTSX", "TOTAL", "")) and
           (r.get("Dim2") in (None, "TOTAL", ""))
    ]
    if preferred:
        rows = preferred

    if year is None:
        years = [r.get("TimeDim") for r in rows if r.get("TimeDim") is not None]
        if not years:
            raise RuntimeError(f"no yearly rows for WHO indicator {code!r}")
        year = max(years)

    rows = [r for r in rows if r.get("TimeDim") == year]
    out: Dict[str, float] = {}
    for r in rows:
        iso = r.get("SpatialDim")
        if not iso or len(iso) != 3:
            continue
        val = r.get("NumericValue")
        if val is None:
            continue
        try:
            out[iso] = float(val)
        except (TypeError, ValueError):
            continue

    if not out:
        raise RuntimeError(
            f"WHO indicator {code!r} yielded no ISO3 values for year {year}"
        )
    return out


# ----------------------------------------------------------------------
def fetch(
    source: str, code: str, year: Optional[int] = None
) -> Dict[str, float]:
    """Dispatch helper used by the CLI / GUI."""
    source = source.lower()
    if source == "worldbank":
        return fetch_worldbank(code, year)
    if source == "who":
        return fetch_who(code, year)
    raise ValueError(f"unknown data source {source!r}")
