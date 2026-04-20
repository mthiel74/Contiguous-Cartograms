"""Tests for CartogramConfig."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest


def test_defaults():
    from cartogram import CartogramConfig

    cfg = CartogramConfig()
    assert cfg.metric == "population"
    assert cfg.performance_goal == "speed"
    assert cfg.grid_size == (256, 512)
    assert cfg.missing_countries == "hide"


def test_from_dict_unknown_key_raises():
    from cartogram import CartogramConfig

    with pytest.raises(ValueError, match="unknown"):
        CartogramConfig.from_dict({"metric": "gdp", "bogus_key": 1})


def test_from_dict_normalises_tuples():
    from cartogram import CartogramConfig

    cfg = CartogramConfig.from_dict({
        "grid_size": [128, 256],
        "bbox": [-180.0, -60.0, 180.0, 85.0],
    })
    assert cfg.grid_size == (128, 256)
    assert cfg.bbox == (-180.0, -60.0, 180.0, 85.0)


def test_replace_and_to_dict():
    from cartogram import CartogramConfig

    a = CartogramConfig(metric="gdp", performance_goal="quality")
    b = a.replace(grid_size=(512, 1024))
    assert a.grid_size == (256, 512)          # original untouched
    assert b.grid_size == (512, 1024)
    assert b.performance_goal == "quality"

    d = b.to_dict()
    assert d["metric"] == "gdp"
    assert d["grid_size"] == (512, 1024)


def test_from_json(tmp_path: Path):
    from cartogram import CartogramConfig

    p = tmp_path / "cfg.json"
    p.write_text(json.dumps({
        "metric": "gdp_per_capita",
        "grid_size": [256, 512],
        "performance_goal": "quality",
        "min_floor_fraction": 0.1,
    }))
    cfg = CartogramConfig.from_file(p)
    assert cfg.metric == "gdp_per_capita"
    assert cfg.grid_size == (256, 512)
    assert cfg.performance_goal == "quality"
    assert cfg.min_floor_fraction == 0.1


def test_from_toml(tmp_path: Path):
    # tomllib is stdlib on 3.11+
    if sys.version_info < (3, 11):
        pytest.importorskip("tomli")
    from cartogram import CartogramConfig

    p = tmp_path / "cfg.toml"
    p.write_text(
        'metric = "gdp"\n'
        'grid_size = [128, 256]\n'
        'performance_goal = "speed"\n'
    )
    cfg = CartogramConfig.from_file(p)
    assert cfg.metric == "gdp"
    assert cfg.grid_size == (128, 256)


def test_from_yaml(tmp_path: Path):
    pytest.importorskip("yaml")
    from cartogram import CartogramConfig

    p = tmp_path / "cfg.yml"
    p.write_text(
        "metric: gdp\n"
        "grid_size: [128, 256]\n"
        "performance_goal: quality\n"
        "missing_countries: grey\n"
    )
    cfg = CartogramConfig.from_file(p)
    assert cfg.metric == "gdp"
    assert cfg.performance_goal == "quality"
    assert cfg.missing_countries == "grey"


def test_from_file_unknown_suffix(tmp_path: Path):
    from cartogram import CartogramConfig

    p = tmp_path / "cfg.txt"
    p.write_text("metric = gdp")
    with pytest.raises(ValueError, match="unknown config format"):
        CartogramConfig.from_file(p)
