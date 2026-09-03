#!/usr/bin/env python
r'''
PlanetNames.py -- Display names for the customizable "earthlike" planet class.

The `earthlike` group of `config-reduce.json` re-defines, numerically, which
planets are counted as "earthlike" (see `PlanetBins.py`). Once that class is
re-defined, calling it an "Earth" in plot titles, axis labels, and table
headers is wrong -- the SAG23 scenarios, for instance, re-define it to be a
Sub-Neptune population.

This module holds the *display names* for that class. They live in the same
`earthlike` group, and all of them are optional:

```
"earthlike": {
    "Earth_Rp_lo": 1.00, "Earth_Rp_hi": 3.50,     <- the numbers (PlanetBins.py)
    "name":        "Sub-Neptune",                 <- the names (this module)
    "name_plural": "Sub-Neptunes",
    "name_adj":    "Sub-Neptune-like",
    "name_short":  "SubNep",
    "symbol":      "\\mathrm{SN}"
    }
```

Omitted entries are derived from `name`; if `name` itself is absent, the
Earth-based defaults are used, so un-customized scenarios are unaffected.

Can run as an executable, to verify customization, like so:

```
$ PlanetNames.py [DIR]
```

where `DIR` is an optional directory to look for a `config-reduce.json`
customization file. A useful idiom, paralleling the one in `PlanetBins.py`:

```
$ diff <(util/reduce_drm_tools/PlanetNames.py sims/SCENARIO) <(util/reduce_drm_tools/PlanetNames.py)
```

which will output any differences from the customization.
'''

# turmon sep 2026

import os
import re
import sys
import argparse
from pathlib import Path

# idiom for imports from "."
try:
    from . import utils
except ImportError:
    import utils


# group within config-reduce.json holding these names
CONFIG_GROUP = 'earthlike'

# the recognized name-keys within that group
#   (PlanetBins.RpLBins skips these when customizing the numeric binner)
CONFIG_KEYS = ('name', 'name_plural', 'name_adj', 'name_short', 'symbol')

# prefix used when relaying names through the reduce_info dict
INFO_PREFIX = 'planet_'

# the historical names -- used when nothing is customized
DEFAULTS = dict(
    name        = 'Earth',
    name_plural = 'Earths',
    name_adj    = 'Earthlike',
    name_short  = 'Earth',
    symbol      = r'\oplus',
    )


class PlanetNames:
    r'''Display names for the (customizable) earthlike planet class.

    Attributes:
      + `name`   -- singular noun, e.g. "Sub-Neptune"
      + `plural` -- plural/count noun, e.g. "Sub-Neptunes"
      + `adj`    -- class-membership adjective, e.g. "Sub-Neptune-like"
      + `short`  -- compact form for tick labels and glued compounds, e.g. "SubNep"
      + `symbol` -- LaTeX subscript for the occurrence rate, e.g. r"\mathrm{SN}"
    '''

    def __init__(self, name=None, name_plural=None, name_adj=None,
                 name_short=None, symbol=None):
        if name:
            # a class name was given: derive any un-given forms from it
            self.name   = name
            self.plural = name_plural or (name + 's')
            self.adj    = name_adj    or (name if ' ' in name else name + '-like')
            self.short  = name_short  or name
            self.symbol = symbol      or self._derive_symbol(self.short)
        else:
            # no class name: fall back to Earth, form by form
            self.name   = DEFAULTS['name']
            self.plural = name_plural or DEFAULTS['name_plural']
            self.adj    = name_adj    or DEFAULTS['name_adj']
            self.short  = name_short  or DEFAULTS['name_short']
            self.symbol = symbol      or DEFAULTS['symbol']

    @staticmethod
    def _derive_symbol(short):
        r'''Make a LaTeX subscript from the short name (upright, no whitespace).'''
        squeezed = re.sub(r'[\s\-_]+', '', short)
        return r'\mathrm{%s}' % squeezed if squeezed else DEFAULTS['symbol']

    @classmethod
    def _from_keymap(cls, keymap):
        r'''Build from a dict keyed by the CONFIG_KEYS. Non-strings are ignored.'''
        kwargs = {}
        for key in CONFIG_KEYS:
            value = keymap.get(key, None)
            # tolerate NaN/None/numbers arriving from CSV or hand-edited JSON
            if isinstance(value, str) and value.strip():
                kwargs[key] = value.strip()
        return cls(**kwargs)

    @classmethod
    def from_config(cls, config):
        r'''Build from a loaded config-reduce dict (which may be None or empty).'''
        if not config:
            return cls()
        group = config.get(CONFIG_GROUP, None)
        if not isinstance(group, dict):
            return cls()
        return cls._from_keymap(group)

    @classmethod
    def from_dir(cls, dirname, log_origin=None):
        r'''Build from the config-reduce.json reachable from dirname.

        Uses utils.load_reduce_config(), so the .fam/.exp parent-lookup rule
        stays in one place. These names are display-only, so an unreadable or
        malformed config yields the defaults with a warning, rather than an
        exception that would abort a plotting or indexing run.
        '''
        try:
            config = utils.load_reduce_config(Path(dirname), log_origin=log_origin)
        except Exception as e:
            print(f'{log_origin or "PlanetNames.py"}: Warning: '
                  f'Could not load reduction config in {dirname} ({e}). '
                  f'Using default planet names.', file=sys.stderr)
            return cls()
        return cls.from_config(config)

    @classmethod
    def from_reduce_info(cls, reduce_info):
        r'''Build from the reduce_info dict handed to each plot function.'''
        if not reduce_info:
            return cls()
        keymap = {key: reduce_info.get(INFO_PREFIX + key, None) for key in CONFIG_KEYS}
        return cls._from_keymap(keymap)

    def to_reduce_info(self):
        r'''Export as plain strings, for merging into the reduce_info dict.

        Plain strings so the names survive pickling into multiprocessing
        workers -- unlike a class-level customization such as RpLBins.'''
        return {
            INFO_PREFIX + 'name':        self.name,
            INFO_PREFIX + 'name_plural': self.plural,
            INFO_PREFIX + 'name_adj':    self.adj,
            INFO_PREFIX + 'name_short':  self.short,
            INFO_PREFIX + 'symbol':      self.symbol,
            }

    def mapping(self):
        r'''Export as a str.format() mapping, for templated label strings.'''
        return dict(
            planet        = self.name,
            planets       = self.plural,
            planet_adj    = self.adj,
            planet_short  = self.short,
            planet_symbol = self.symbol,
            eta           = self.eta,
            eta_html      = self.eta_html,
            )

    @property
    def eta(self):
        r'''The occurrence rate as LaTeX, e.g. r"\eta_{\oplus}".'''
        return r'\eta_{%s}' % self.symbol

    @property
    def eta_html(self):
        r'''The occurrence rate as HTML, e.g. "eta<sub>Earth</sub>".'''
        return f'eta<sub>{self.short}</sub>'

    def is_default(self):
        r'''True if these are the historical Earth-based names.'''
        return (self.name   == DEFAULTS['name']   and
                self.plural == DEFAULTS['name_plural'] and
                self.adj    == DEFAULTS['name_adj'] and
                self.short  == DEFAULTS['name_short'] and
                self.symbol == DEFAULTS['symbol'])

    def show(self):
        r'''Print the names we resolved to.'''
        print(f'Group: {CONFIG_GROUP}')
        for label, value in (('name', self.name), ('name_plural', self.plural),
                             ('name_adj', self.adj), ('name_short', self.short),
                             ('symbol', self.symbol)):
            print(f'  attribute: {label}')
            print(f'    {value}')
        print(f'  derived: eta')
        print(f'    {self.eta}')

    def __repr__(self):
        return (f'{self.__class__.__name__}(name={self.name!r}, '
                f'name_plural={self.plural!r}, name_adj={self.adj!r}, '
                f'name_short={self.short!r}, symbol={self.symbol!r})')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Show planet-class display names.",
                                     epilog='')
    parser.add_argument('dir', metavar='DIR', type=str, nargs='?', default='',
                            help='directory to initialize from')
    args = parser.parse_args()
    args.progname = os.path.basename(sys.argv[0])

    if args.dir:
        p = Path(args.dir)
        if not p.is_dir():
            print(f'{args.progname}: Directory {args.dir} is not a readable directory',
                  file=sys.stderr)
            sys.exit(1)
        config = utils.load_reduce_config(p, log_origin=args.progname)
        if config is None:
            print(f'{args.progname}: No customization found! Empty customization used.')
            config = {}
        else:
            print(f'{args.progname}: Loaded customization from file.')
    else:
        config = {}
        print(f'{args.progname}: No customization given, using defaults.')

    names = PlanetNames.from_config(config)
    print(f'{args.progname}: Planet-name summary:')
    names.show()
    print(f'{args.progname}: Done.')
    sys.exit(0)
