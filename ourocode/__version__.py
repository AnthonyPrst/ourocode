from pathlib import Path

try:
    import tomllib
except ImportError:
    try:
        import tomli as tomllib
    except ImportError:
        tomllib = None

if tomllib is not None:
    _pyproject_path = Path(__file__).parent.parent / "pyproject.toml"
    with open(_pyproject_path, "rb") as f:
        _data = tomllib.load(f)
    __version__ = _data["project"]["version"]
else:
    # Fallback si tomllib/tomli non disponible (Python < 3.11 sans tomli)
    __version__ = "1.9.0"
