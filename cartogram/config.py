"""Reproducible config for :func:`world_cartogram` runs.

A single :class:`CartogramConfig` dataclass captures every knob that
affects the output. Use it to:

* document defaults in one place rather than spread across
  ``Cartogram.__init__``, ``world_cartogram``, the CLI, and the GUI;
* run a cartogram from a YAML/TOML/JSON file so experiments are
  reproducible from a single committed artefact;
* build variants with ``replace(cfg, grid_size=(512, 1024))``.

All fields are optional; the defaults match the ones built into
:func:`world_cartogram`.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields, replace
from pathlib import Path
from typing import Any, Dict, Literal, Mapping, Optional, Tuple, Union


Metric = Union[str, Mapping[str, float]]


@dataclass
class CartogramConfig:
    """Every option :func:`world_cartogram` cares about, in one place.

    Callable metrics (``Function[row] -> float``) are not supported
    here, since they don't serialise; pass them directly to
    ``world_cartogram`` instead. Dicts, column names, and paths work
    fine.
    """

    metric: Metric = "population"
    label: Optional[str] = None

    # Loader
    resolution: str = "110m"
    projection: str = "EPSG:8857"
    drop_antarctica: bool = True
    merge_sovereign: bool = False

    # Rasterisation
    grid_size: Tuple[int, int] = (256, 512)
    bbox: Optional[Tuple[float, float, float, float]] = None

    # Diffusion / preprocessing
    mean_floor: float = 0.02
    blur_sigma: float = 1.5
    sea_density: Union[str, float, None] = "auto"
    min_floor_fraction: Optional[float] = None
    rho_ceil_multiplier: Optional[float] = None

    # Solver
    performance_goal: Literal["speed", "quality"] = "speed"
    snapshots: int = 60
    tol: float = 1e-3

    # Output
    render: bool = True
    cmap: str = "plasma"
    missing_countries: Literal["hide", "grey"] = "hide"

    # ------------------------------------------------------------------
    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "CartogramConfig":
        """Construct from a plain mapping; unknown keys raise."""
        known = {f.name for f in fields(cls)}
        extra = set(data) - known
        if extra:
            raise ValueError(
                f"unknown CartogramConfig fields: {sorted(extra)}"
            )
        kwargs: Dict[str, Any] = dict(data)
        # Normalise tuple-shaped fields from the list form that YAML /
        # JSON use natively.
        if "grid_size" in kwargs and kwargs["grid_size"] is not None:
            kwargs["grid_size"] = tuple(kwargs["grid_size"])  # type: ignore[arg-type]
        if "bbox" in kwargs and kwargs["bbox"] is not None:
            kwargs["bbox"] = tuple(kwargs["bbox"])  # type: ignore[arg-type]
        return cls(**kwargs)

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> "CartogramConfig":
        """Load config from a YAML file. Requires PyYAML."""
        import yaml
        p = Path(path)
        with p.open("r", encoding="utf-8") as fh:
            data = yaml.safe_load(fh) or {}
        if not isinstance(data, Mapping):
            raise ValueError(
                f"{path}: top-level YAML must be a mapping, got {type(data).__name__}"
            )
        return cls.from_dict(data)

    @classmethod
    def from_toml(cls, path: Union[str, Path]) -> "CartogramConfig":
        """Load config from a TOML file. Uses stdlib ``tomllib`` (py>=3.11)."""
        try:
            import tomllib  # type: ignore[attr-defined]
        except ImportError:  # pragma: no cover
            import tomli as tomllib  # type: ignore[no-redef]
        p = Path(path)
        with p.open("rb") as fh:
            data = tomllib.load(fh)
        return cls.from_dict(data)

    @classmethod
    def from_json(cls, path: Union[str, Path]) -> "CartogramConfig":
        """Load config from a JSON file."""
        import json
        p = Path(path)
        with p.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
        if not isinstance(data, Mapping):
            raise ValueError(
                f"{path}: top-level JSON must be an object, got {type(data).__name__}"
            )
        return cls.from_dict(data)

    @classmethod
    def from_file(cls, path: Union[str, Path]) -> "CartogramConfig":
        """Dispatch on file extension: .yaml/.yml, .toml, or .json."""
        p = Path(path)
        suffix = p.suffix.lower()
        if suffix in (".yaml", ".yml"):
            return cls.from_yaml(p)
        if suffix == ".toml":
            return cls.from_toml(p)
        if suffix == ".json":
            return cls.from_json(p)
        raise ValueError(
            f"{path}: unknown config format {suffix!r}; expected .yaml/.yml, .toml, or .json"
        )

    # ------------------------------------------------------------------
    def replace(self, **changes: Any) -> "CartogramConfig":
        """Return a copy with fields overridden."""
        return replace(self, **changes)

    def to_dict(self) -> Dict[str, Any]:
        """Plain-dict form suitable for logging / JSON serialisation."""
        return asdict(self)


def run_config(cfg: CartogramConfig):
    """Shortcut: ``run_config(cfg)`` is ``world_cartogram(**cfg_as_kwargs)``.

    Lives here so tests don't have to import the full GIS stack just
    to construct a config.
    """
    from .one_liner import world_cartogram

    kwargs = cfg.to_dict()
    metric = kwargs.pop("metric")
    return world_cartogram(metric, **kwargs)
