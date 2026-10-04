"""Development-only Python 2 semantic adapter for the pinned TAMA scripts.
Prefer a real CPython 2.7 oracle when available. This fallback executes the
unaltered source AST with seed-0 64-bit dict ordering, integer division,
Python 2 float formatting/rounding and text subprocess streams. Validated
against the existing CPython 2 gmap golden outputs; it is not a general
Python 2 interpreter. Requires Python 3.9-3.12 and Biopython. The unused
pysam import is stubbed and fails if any attribute is used.
Usage: python3 oracle/python2_compat.py ../tama/tama_collapse.py [args...]
"""
import ast, sys, runpy, decimal, functools, builtins
from pathlib import Path
from collections.abc import MutableMapping
from lib2to3.refactor import RefactoringTool, get_fixers_from_package
MASK = (1 << 64) - 1
EMPTY = object()
DUMMY = object()

def phash(key):
    if isinstance(key, int):
        return key if key != -1 else -2
    if isinstance(key, str):
        b = key.encode()
        h = b[0] << 7 if b else 0
        for v in b:
            h = (h * 1000003 ^ v) & MASK
        h ^= len(b)
        h = h if h < 1 << 63 else h - (1 << 64)
        return h if h != -1 else -2
    raise TypeError(type(key))

class D(MutableMapping):

    def __init__(self, items=()):
        self.t = [EMPTY] * 8
        self.used = 0
        self.fill = 0
        for k, v in items:
            self[k] = v

    def slot(self, k):
        h = phash(k)
        p = h & MASK
        i = p & len(self.t) - 1
        free = None
        while True:
            v = self.t[i]
            if v is EMPTY:
                return i if free is None else free
            if v is DUMMY:
                if free is None:
                    free = i
            elif v[0] == k:
                return i
            i = i * 5 + p + 1 & len(self.t) - 1
            p >>= 5

    def __getitem__(self, k):
        v = self.t[self.slot(k)]
        if v is EMPTY or v is DUMMY:
            raise KeyError(k)
        return v[1]

    def __setitem__(self, k, v):
        i = self.slot(k)
        old = self.t[i]
        if old is EMPTY:
            self.fill += 1
        if old is EMPTY or old is DUMMY:
            self.used += 1
        self.t[i] = (k, v)
        if self.fill * 3 >= len(self.t) * 2:
            size = 8
            minimum = self.used * (2 if self.used > 50000 else 4)
            while size <= minimum:
                size *= 2
            old = self.t
            self.t = [EMPTY] * size
            self.used = 0
            self.fill = 0
            for v in old:
                if v is not EMPTY and v is not DUMMY:
                    i = self.slot(v[0])
                    self.t[i] = v
                    self.used += 1
                    self.fill += 1

    def __delitem__(self, k):
        i = self.slot(k)
        self[k]
        self.t[i] = DUMMY
        self.used -= 1

    def __iter__(self):
        for v in self.t:
            if v is not EMPTY and v is not DUMMY:
                yield v[0]

    def __len__(self):
        return self.used

    def copy(self):
        return D(self.items())

def div(a, b):
    return a // b if isinstance(a, int) and isinstance(b, int) else a / b

def string(x=''):
    if isinstance(x, float):
        s = format(x, '.12g')
        if '.' not in s and 'e' not in s and (s not in ('nan', 'inf', '-inf')):
            s += '.0'
        return s
    return str(x)

def rounding(x, n=0):
    return float(decimal.Decimal.from_float(float(x)).quantize(decimal.Decimal(1).scaleb(-n), rounding=decimal.ROUND_HALF_UP))

class Transform(ast.NodeTransformer):

    def visit_Dict(self, n):
        self.generic_visit(n)
        return ast.copy_location(ast.Call(ast.Name('D', ast.Load()), [ast.List([ast.Tuple([k, v], ast.Load()) for k, v in zip(n.keys, n.values)], ast.Load())], []), n)

    def visit_BinOp(self, n):
        self.generic_visit(n)
        if isinstance(n.op, ast.Div):
            return ast.copy_location(ast.Call(ast.Name('div', ast.Load()), [n.left, n.right], []), n)
        return n

    def visit_Call(self, n):
        self.generic_visit(n)
        if isinstance(n.func, ast.Name) and n.func.id == 'str':
            n.func.id = 'string'
        if isinstance(n.func, ast.Name) and n.func.id == 'Popen':
            n.keywords.append(ast.keyword(arg='universal_newlines', value=ast.Constant(True)))
        return n
import types
try:
    import pysam
except ImportError:
    stub = types.ModuleType('pysam')

    def unavailable(name):
        raise RuntimeError('unused pysam dependency accessed: ' + name)
    stub.__getattr__ = unavailable
    sys.modules['pysam'] = stub
source = Path(sys.argv[1])
sys.argv = sys.argv[1:]
sys.path.insert(0, str(source.parent))
import hashlib, marshal
import tempfile
cache = Path(tempfile.gettempdir()) / ('braid-oracle-' + hashlib.sha256(source.read_bytes() + Path(__file__).read_bytes()).hexdigest() + '.marshal')
if cache.exists():
    code = marshal.loads(cache.read_bytes())
else:
    s = str(RefactoringTool([f for f in get_fixers_from_package('lib2to3.fixes') if not f.endswith('fix_import')]).refactor_string(source.read_text(), str(source)))
    t = ast.fix_missing_locations(Transform().visit(ast.parse(s)))
    code = compile(t, str(source), 'exec')
    cache.write_bytes(marshal.dumps(code))
ns = {'__name__': '__main__', '__file__': str(source), 'D': D, 'div': div, 'string': string, 'round': rounding}
exec(code, ns)
