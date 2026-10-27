"""GFE format linter - shared by GFE_format_linter.ipynb and the web page (index.html).

Checks whether a data file conforms to the Global Fireball Exchange (GFE)
standard, version 1.2 (https://github.com/Global-Fireball-Exchange-standard).

Entry point (no file-system access, no prompts):
    check_file(text, filename='') -> str          the linter report
"""
import os
import re
import warnings
import numbers
from datetime import datetime
import numpy as np
from astropy.table import Table

# For every item: mand = is it mandatory; kind = expected value type;
#   dec = minimum number of decimal places the standard asks for (0 = none);
#   lo / hi = inclusive valid range for numeric values (None = unbounded).
# kind is one of:
#   'float'   floating-point number
#   'num'     integer or floating-point number
#   'int'     integer
#   'str'     free text
#   'isodate' ISO-8601 string  YYYY-MM-DDThh:mm:ss.sss  (>= 3 decimals)
#   'bool'    True / False

META_SPEC = {
    # --- the three mandatory metadata items -------------------------------- #
    'obs_latitude':  dict(mand=True,  kind='float', dec=6, lo=-90.0,  hi=90.0),
    'obs_longitude': dict(mand=True,  kind='float', dec=6, lo=-180.0, hi=180.0),
    'obs_elevation': dict(mand=True,  kind='num',   dec=0, lo=-500.0, hi=9000.0),
    # --- optional metadata ------------------------------------------------- #
    'origin':        dict(mand=False, kind='str'),
    'location':      dict(mand=False, kind='str'),
    'telescope':     dict(mand=False, kind='str'),
    'camera_id':     dict(mand=False, kind='str'),
    'observer':      dict(mand=False, kind='str'),
    'comment':       dict(mand=False, kind='str', nobraces=True),
    'instrument':    dict(mand=False, kind='str'),
    'lens':          dict(mand=False, kind='str'),
    'cx':            dict(mand=False, kind='int', lo=1),
    'cy':            dict(mand=False, kind='int', lo=1),
    'photometric_band': dict(mand=False, kind='str'),
    'image_file':    dict(mand=False, kind='str'),
    'isodate_start_obs': dict(mand=False, kind='isodate'),
    'isodate_calib': dict(mand=False, kind='isodate'),
    'exposure_time': dict(mand=False, kind='float', dec=2, lo=0.0),
    'astrometry_number_stars': dict(mand=False, kind='int', lo=0),
    'mag_label':     dict(mand=False, kind='str'),
    'no_frags':      dict(mand=False, kind='int', lo=1),
    'obs_az':        dict(mand=False, kind='float', dec=3, lo=0.0,    hi=360.0),
    'obs_ev':        dict(mand=False, kind='float', dec=3, lo=-90.0,  hi=90.0),
    'obs_rot':       dict(mand=False, kind='float', dec=3, lo=-360.0, hi=360.0),
    'fov_horiz':     dict(mand=False, kind='float', dec=3, lo=0.0,    hi=360.0),
    'fov_vert':      dict(mand=False, kind='float', dec=3, lo=0.0,    hi=360.0),
}

# Point-observation columns.  These are the BASE names; a fragment number and/or
# a trailing "V" (velocity-pick) may be appended - see Appendix 4.
COL_SPEC = {
    'datetime': dict(mand=True,  kind='isodate'),
    'ra':       dict(mand=True,  kind='float', dec=6, lo=0.0,   hi=360.0),
    'dec':      dict(mand=True,  kind='float', dec=6, lo=-90.0, hi=90.0),
    'azimuth':  dict(mand=True,  kind='float', dec=6, lo=0.0,   hi=360.0),
    'altitude': dict(mand=True,  kind='float', dec=6, lo=-90.0, hi=90.0),
    'x_image':  dict(mand=False, kind='float'),
    'y_image':  dict(mand=False, kind='float'),
    'err_minus_altitude': dict(mand=False, kind='float'),
    'err_plus_altitude':  dict(mand=False, kind='float'),
    'err_minus_azimuth':  dict(mand=False, kind='float'),
    'err_plus_azimuth':   dict(mand=False, kind='float'),
    'saturated_pixels':   dict(mand=False, kind='bool'),
    'integrated_pixel_value': dict(mand=False, kind='int'),
    'err_minus_mag':      dict(mand=False, kind='float'),
    'err_plus_mag':       dict(mand=False, kind='float'),
}

# The light-curve column is named by the mag_label metadata item.  These are the
# names the standard documents; the actual name from mag_label is added at run
# time.  Its value is floating point.
MAG_COL_BASES = {'mag', 'FLUX_AUTO', 'no_mag_data'}

# The mandatory per-fragment position columns (datetime is shared, not repeated).
PER_FRAGMENT_MANDATORY = ['ra', 'dec', 'azimuth', 'altitude']

ISO_RE = re.compile(r'^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3,}$')


# --------------------------------------------------------------------------- #
# Report collector
class Report:
    def __init__(self, name):
        self.name = name
        self.errors = []
        self.warnings = []
        self.not_in_standard = []
        self.info = []

    def error(self, msg):   self.errors.append(msg)
    def warn(self, msg):    self.warnings.append(msg)
    def extra(self, msg):   self.not_in_standard.append(msg)
    def note(self, msg):    self.info.append(msg)

    def render(self):
        line = '=' * 72
        out = [line, 'GFE format linter  -  ' + os.path.basename(str(self.name)), line]
        out += ['  . ' + m for m in self.info]
        out.append('')

        out.append('EXCEPTIONS - ERRORS (break conformance to the GFE standard): %d'
                   % len(self.errors))
        out += ['  [ERROR]   ' + m for m in self.errors]
        out.append('')

        out.append('EXCEPTIONS - WARNINGS (conforms, but departs from the standard): %d'
                   % len(self.warnings))
        out += ['  [WARNING] ' + m for m in self.warnings]
        out.append('')

        out.append('DATA ITEMS PRESENT BUT NOT DEFINED IN THE GFE STANDARD: %d'
                   % len(self.not_in_standard))
        out += ['  [EXTRA]   ' + m for m in self.not_in_standard]
        out.append('')

        out.append(line)
        if self.errors:
            out.append('VERDICT:  DOES NOT CONFORM to the GFE standard  '
                       '(%d error(s), %d warning(s)).'
                       % (len(self.errors), len(self.warnings)))
        else:
            out.append('VERDICT:  CONFORMS to the GFE standard  '
                       '(0 errors, %d warning(s)).' % len(self.warnings))
        out.append(line)
        return '\n'.join(out)


# --------------------------------------------------------------------------- #
# Small helpers
def is_bool(v):
    return isinstance(v, bool)

def is_int(v):
    return isinstance(v, numbers.Integral) and not isinstance(v, bool)

def is_real(v):
    return isinstance(v, numbers.Real) and not isinstance(v, bool)

def count_decimals(raw):
    """Number of decimal places written in a raw numeric string (0 if none)."""
    if raw is None:
        return None
    s = str(raw).strip().strip("'\"")
    if '.' not in s:
        return 0
    frac = s.split('.', 1)[1]
    frac = re.split('[eE]', frac)[0]
    return len(frac)

def valid_iso(s):
    if not isinstance(s, str) or not ISO_RE.match(s):
        return False
    try:
        datetime.strptime(s, '%Y-%m-%dT%H:%M:%S.%f')
        return True
    except ValueError:
        return False

def parse_raw(text):
    """Return (raw_meta, raw_firstrow): dicts of as-written string values.

    Best effort - used only for the 'number of decimal places' checks, which the
    parsed numeric values cannot answer. Never fatal.
    """
    raw_meta, raw_firstrow = {}, {}
    known_meta = set(META_SPEC)
    header_names = None
    for ln in text.splitlines():
        if ln.startswith('#'):
            m = re.search(r'-\s*\{\s*([A-Za-z_][\w]*)\s*:\s*(.*?)\s*\}\s*$', ln)
            if m and m.group(1) in known_meta:
                raw_meta[m.group(1)] = m.group(2)
        else:
            cells = ln.split(',')
            if header_names is None:
                header_names = cells            # first non-# line = column names
            else:
                for name, val in zip(header_names, cells):   # first data row
                    raw_firstrow[name] = val
                break
    return raw_meta, raw_firstrow


# --------------------------------------------------------------------------- #
# Value checks and column-name decomposition
def check_scalar(name, value, spec, raw, report, where):
    """Check one metadata value against its spec entry."""
    kind = spec['kind']
    ok_type = True

    if kind == 'str':
        if not isinstance(value, str):
            report.error('%s "%s" should be free text but is %s'
                         % (where, name, type(value).__name__))
            ok_type = False
        elif spec.get('nobraces') and ('{' in value or '}' in value):
            report.error('%s "%s" must not contain the characters { or }' % (where, name))
    elif kind == 'int':
        if not is_int(value):
            report.error('%s "%s" should be an integer but is %r (%s)'
                         % (where, name, value, type(value).__name__))
            ok_type = False
    elif kind in ('float', 'num'):
        if not (is_real(value) or is_int(value)):
            report.error('%s "%s" should be a number but is %r (%s)'
                         % (where, name, value, type(value).__name__))
            ok_type = False
    elif kind == 'isodate':
        if not valid_iso(value):
            report.error('%s "%s" should be a single-quoted ISO-8601 date '
                         'YYYY-MM-DDThh:mm:ss.sss (>=3 decimals); got %r'
                         % (where, name, value))
            ok_type = False

    if not ok_type:
        return

    if kind in ('float', 'num', 'int'):
        lo, hi = spec.get('lo'), spec.get('hi')
        if lo is not None and value < lo:
            report.error('%s "%s" = %s is below the valid minimum of %s'
                         % (where, name, value, lo))
        if hi is not None and value > hi:
            report.error('%s "%s" = %s is above the valid maximum of %s'
                         % (where, name, value, hi))

    need = spec.get('dec', 0)
    if need and kind in ('float', 'num'):
        got = count_decimals(raw)
        if got is not None and got < need:
            msg = ('%s "%s" is written with %d decimal place(s); the standard '
                   'asks for at least %d' % (where, name, got, need))
            report.warn(msg)


def dtype_kind(col):
    try:
        return col.dtype.kind
    except Exception:
        return '?'


def check_column(colname, base, col, spec, raw, report):
    """Check one data column against the spec for its base name."""
    kind = spec['kind']
    k = dtype_kind(col)

    if kind == 'isodate':
        if k not in ('U', 'S', 'O'):
            report.error('column "%s" should hold ISO-8601 date strings but has '
                         'numeric/other type (%s)' % (colname, col.dtype))
            return
        bad = [str(v) for v in col if not valid_iso(str(v))]
        if bad:
            report.error('column "%s" has %d value(s) that are not valid ISO-8601 '
                         'dates (e.g. %r)' % (colname, len(bad), bad[0]))
        return

    if kind == 'bool':
        if k == 'b':
            return
        vals = set(str(v).strip() for v in col)
        if not vals <= {'True', 'False', 'true', 'false'}:
            report.error('column "%s" should be True/False but has values %s'
                         % (colname, sorted(vals)[:4]))
        return

    if kind == 'int':
        if k not in ('i', 'u'):
            report.error('column "%s" should be integer but has type %s'
                         % (colname, col.dtype))
        return

    # numeric (float)
    if k not in ('f', 'i', 'u'):
        report.error('column "%s" should be a floating-point number but has type %s'
                     % (colname, col.dtype))
        return
    if k in ('i', 'u'):
        msg = ('column "%s" is stored as integer; the standard specifies '
               'floating point' % colname)
        report.warn(msg)

    lo, hi = spec.get('lo'), spec.get('hi')
    if lo is not None or hi is not None:
        try:
            arr = np.asarray(col, dtype=float)
            finite = arr[np.isfinite(arr)]
            if finite.size:
                if lo is not None and finite.min() < lo:
                    report.error('column "%s" has value(s) below the valid minimum '
                                 'of %s (min = %s)' % (colname, lo, finite.min()))
                if hi is not None and finite.max() > hi:
                    report.error('column "%s" has value(s) above the valid maximum '
                                 'of %s (max = %s)' % (colname, hi, finite.max()))
        except Exception:
            pass

    need = spec.get('dec', 0)
    if need:
        got = count_decimals(raw)
        if got is not None and got < need:
            msg = ('column "%s" first value is written with %d decimal place(s); '
                   'the standard asks for at least %d' % (colname, got, need))
            report.warn(msg)


def decompose(colname, known_bases):
    """Return (base, frag, is_velocity) for a column name.

    'ra'   -> ('ra', 0, False)      'ra2'  -> ('ra', 2, False)
    'ra0'  -> ('ra', 0, False)      'ra2V' -> ('ra', 2, True)
    Only treats a trailing number as a fragment index if the remaining stem is a
    recognised base name, so names like 'no_mag_data' are left intact.
    """
    is_v = colname.endswith('V')
    core = colname[:-1] if is_v else colname
    m = re.match(r'^(.*?)(\d+)$', core)
    if m and m.group(1) in known_bases:
        return m.group(1), int(m.group(2)), is_v
    return core, 0, is_v


# --------------------------------------------------------------------------- #
# Metadata and column checks
def check_metadata(meta, raw_meta, report):
    # mandatory present?
    for name, spec in META_SPEC.items():
        if spec['mand'] and name not in meta:
            report.error('mandatory metadata item "%s" is missing' % name)

    # validate each present item; anything not in the standard is reported the
    # same way, whatever its name.
    for name in list(meta.keys()):
        if name in META_SPEC:
            check_scalar(name, meta[name], META_SPEC[name], raw_meta.get(name),
                         report, 'metadata')
        else:
            report.extra('metadata key "%s"  (not defined in the GFE standard)' % name)


def check_columns(table, raw_firstrow, report):
    colnames = list(table.colnames)
    colset = set(colnames)

    # magnitude column bases (documented names plus the one named by mag_label)
    mag_bases = set(MAG_COL_BASES)
    mag_label = table.meta.get('mag_label')
    if isinstance(mag_label, str) and mag_label:
        mag_bases.add(mag_label)
        if mag_label not in colset:
            report.warn('metadata "mag_label" = "%s" but there is no column with '
                        'that name' % mag_label)

    known_bases = set(COL_SPEC) | mag_bases

    if 'datetime' not in colset:
        report.error('mandatory column "datetime" is missing')

    nf = table.meta.get('no_frags', 1)
    if not is_int(nf) or nf < 1:
        nf = 1
    for n in range(nf):
        for base in PER_FRAGMENT_MANDATORY:
            cands = [base + str(n)]
            if n == 0:
                cands.append(base)      # fragment 0 may omit the "0"
            if not any(c in colset for c in cands):
                report.error('mandatory column "%s" (fragment %d) is missing'
                             % (base if n == 0 else base + str(n), n))

    for colname in colnames:
        base, frag, is_v = decompose(colname, known_bases)
        if base in COL_SPEC:
            check_column(colname, base, table[colname], COL_SPEC[base],
                         raw_firstrow.get(colname), report)
        elif base in mag_bases:
            check_column(colname, base, table[colname],
                         dict(kind='float'), raw_firstrow.get(colname), report)
        else:
            report.extra('column "%s"  (not defined in the GFE standard)' % colname)


# --------------------------------------------------------------------------- #
# Reading and the top-level routines
def read_table(text, report):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            table = Table.read(text, format='ascii.ecsv')
        except Exception:
            return None

        for w in caught:
            report.note('astropy read warning: %s' % w.message)
    return table


def lint(text, filename=''):
    """Lint the text of a GFE .ecsv file and return the Report."""
    report = Report(filename or '(unnamed file)')
    if filename:
        report.note('file: %s' % filename)
        if not filename.lower().endswith('.ecsv'):
            report.error('the file extension is not ".ecsv"')

    table = read_table(text, report)
    if table is None:
        report.error('File cannot be read because it is not an astropy table in ECSV format')
        return report

    report.note('read OK: %d column(s), %d data row(s)'
                % (len(table.colnames), len(table)))
    raw_meta, raw_firstrow = parse_raw(text)

    check_metadata(table.meta, raw_meta, report)
    check_columns(table, raw_firstrow, report)
    return report


def check_file(text, filename=''):
    """Lint the text of a GFE .ecsv file and return the report as a string."""
    return lint(text, filename).render()

