"""Tests for the country-key normaliser."""

from cartogram.country_lookup import canonical, normalise_mapping


def test_canonical_basic():
    assert canonical("USA") == "USA"
    assert canonical("us") == "USA"
    assert canonical("U.S.") == "USA"
    assert canonical("United States") == "USA"
    assert canonical("united states of america") == "USA"
    assert canonical("UnitedStates") == "USA"          # Wolfram-style
    assert canonical("840") == "USA"                    # ISO numeric


def test_canonical_uk_variants():
    assert canonical("UK") == "GBR"
    assert canonical("Great Britain") == "GBR"
    assert canonical("UnitedKingdom") == "GBR"
    assert canonical("826") == "GBR"


def test_canonical_unknown():
    assert canonical("Atlantis") is None
    assert canonical("") is None
    assert canonical("ZZZ") is None


def test_canonical_non_string():
    assert canonical(None) is None      # type: ignore[arg-type]
    assert canonical(840) is None       # type: ignore[arg-type]


def test_normalise_mapping_mixed_keys():
    out = normalise_mapping({
        "us": 1.0,
        "UnitedKingdom": 2.0,
        "France": 3.0,
        "JPN": 4.0,
        "deutschland": 5.0,             # unknown -> falls back to "DEUTSCHLAND"
    })
    assert out["USA"] == 1.0
    assert out["GBR"] == 2.0
    assert out["FRA"] == 3.0
    assert out["JPN"] == 4.0
    # Fallback: uppercased stripped form so a user-supplied ISO3 that
    # the table doesn't know still round-trips.
    assert out["DEUTSCHLAND"] == 5.0


def test_normalise_mapping_drop_unknown():
    out = normalise_mapping(
        {"us": 1.0, "atlantis": 2.0},
        default_to_key=False,
    )
    assert out == {"USA": 1.0}
