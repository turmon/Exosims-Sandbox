"""__init__: Allow loading of modules containing - by substituting _

Some of the python utilities here were written with - in their names.
They are command-line utilities, so for actual use, this is OK. However,
tooling (like mkdocs) imports them to examine docstrings, and the
import of "drm-tabulate" can fail. This shim makes a mapping to translate
hyphen-scripts <-> underscore_scripts.

Synopsis:
  `from util import drm_tabulate`
is in effect translated to:
  `from util import drm-tabulate`

by the mapping defined here. 

Module-level `__getattr__` (PEP 562, Python 3.7+) is only invoked
when normal attribute lookup fails, so non-hyphenated imports are
completely unaffected by this change. Hyphenated scripts are loaded
on first access and cached in `sys.modules` for subsequent imports.
The one caveat:
  `from pkg import my_script` 
triggers `__getattr__`, but 
  `import pkg; pkg.my_script` 
on a second call hits `sys.modules` directly so we
get lazy + cached behavior for free.

"""

import importlib.util, pathlib, sys

_HYPHEN_MAP = {
    p.stem.replace("-", "_"): p
    for p in pathlib.Path(__file__).parent.glob("*-*.py")
}

def __getattr__(name: str):
    if name not in _HYPHEN_MAP:
        raise AttributeError(f"module 'pkg' has no attribute {name!r}")
    path = _HYPHEN_MAP[name]
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[f"pkg.{name}"] = mod
    spec.loader.exec_module(mod)
    return mod
