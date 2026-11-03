# Transcribes the pipe stem tables out of McMillan (2010) into data/pipe_stems_mcmillan2010.csv.
# Input is the text layer: pdftotext -layout mcmillan2010.pdf mcm_layout.txt
# Columns in the PDF tables are positional, so counts are assigned to the bore
# size whose header centre is nearest, then checked against the published Binford date.

import re, csv, sys

lines = open('data/mcmillan2010_pdf_text.txt').read().split('\n')
bores = ['9/64"','8/64"','7/64"','6/64"','5/64"','4/64"']

def col_centres(h):
    pos = {}
    for b in bores:
        i = h.find(b)
        if i >= 0:
            pos[b] = i + len(b)/2.0   # centre of the header label
    pos['Sum'] = h.find('Sum') + 1.5
    return pos

rows = []
i = 0
while i < len(lines):
    h = lines[i]
    if '9/64"' in h and 'Sum' in h and 'Context' in h:
        pos = col_centres(h)
        j = i + 1
        block = []
        while j < len(lines) and not re.search(r'Table 4\.\d+:', lines[j]):
            if lines[j].strip():
                block.append(lines[j])
            j += 1
        caption = lines[j] if j < len(lines) else ''
        for ln in block:
            # context label = text before the date column
            m = re.search(r'\s{2,}((?:MCD: )?\d{4}[^\s]*(?:\s*-\s*\d{4}[^\s]*)?|\d{4}\'s[^\s]*)\s', ln)
            if not m:
                continue
            ctx = ln[:m.start()].strip()
            date = m.group(1).strip()
            counts = {b: '' for b in bores}
            total = None
            for nm in re.finditer(r'\d+', ln[m.end():]):
                val = int(nm.group())
                centre = m.end() + (nm.start() + nm.end())/2.0
                # nearest header column centre
                b = min(pos, key=lambda k: abs(pos[k] - centre))
                if b == 'Sum':
                    total = val
                else:
                    counts[b] = val
            rows.append(dict(caption=caption.strip(), context=ctx, date=date,
                             **{b: counts[b] for b in bores}, sum=total))
        i = j
    i += 1

bad = []
for r in rows:
    s = sum(v for v in (r[b] for b in bores) if v != '')
    if r['sum'] != s:
        bad.append((r['caption'][:40], r['context'], s, r['sum']))

with open('data/pipe_counts_raw.csv','w',newline='') as f:
    w = csv.DictWriter(f, fieldnames=['caption','context','date']+bores+['sum'])
    w.writeheader(); w.writerows(rows)
print('rows:', len(rows))
print('sum mismatches:', len(bad))
for b in bad: print('  ', b)

import re, csv
lines = open('data/mcmillan2010_pdf_text.txt').read().split('\n')[4726:4815]
state = None
rows = []
pending = None
for ln in lines:
    m = re.match(r'\s+(Maryland|Virginia|North Carolina|South Carolina) Formula Results', ln)
    if m:
        state = m.group(1); continue
    if 'Site & Context' in ln or not ln.strip() or ln.strip().isdigit() or ln.startswith('APPENDIX'):
        continue
    # row: name, n, date, mean, binford, diff, hanson, diff, hd, diff
    m = re.match(r'\s+(\S.*?)\s{2,}(\d+)\s+(MCD: \d{4}|\S+)\s+(\d{4})\s+(\d{4})\s+([+-]?\d+)\s+(\d{4})\s+([+-]?\d+)\s+(\d{4})\s+([+-]?\d+)\s*$', ln)
    if m:
        g = m.groups()
        rows.append(dict(state=state, context=g[0].strip(), n=int(g[1]), date=g[2],
                         mean_date=int(g[3]), binford=int(g[4]), hanson=int(g[6]), hd=int(g[8])))
        continue
    # split row (Abell's Trash Pit): name alone, then numbers alone
    if re.match(r'^\s+\D+$', ln) and 'cont' not in ln:
        pending = ln.strip(); continue
    m = re.match(r'\s+(\d+)\s+(\S+)\s+(\d{4})\s+(\d{4})\s+([+-]?\d+)\s+(\d{4})\s+([+-]?\d+)\s+(\d{4})\s+([+-]?\d+)\s*$', ln)
    if m and pending:
        g = m.groups()
        rows.append(dict(state=state, context=pending, n=int(g[0]), date=g[1],
                         mean_date=int(g[2]), binford=int(g[3]), hanson=int(g[5]), hd=int(g[7])))
        pending = None
        continue
    print('UNPARSED:', repr(ln))

with open('data/appendixA.csv','w',newline='') as f:
    w = csv.DictWriter(f, fieldnames=['state','context','n','date','mean_date','binford','hanson','hd'])
    w.writeheader(); w.writerows(rows)
print('appendix rows:', len(rows))

import csv, re
bores = ['9/64"','8/64"','7/64"','6/64"','5/64"','4/64"']
sizes = [9,8,7,6,5,4]
cnt = list(csv.DictReader(open('data/pipe_counts_raw.csv')))
app = list(csv.DictReader(open('data/appendixA.csv')))

for c in cnt:
    c['n'] = sum(int(c[b]) for b in bores if c[b])
    c['mean_bore'] = sum(int(c[b])*s for b, s in zip(bores, sizes) if c[b]) / c['n']
    c['binford_calc'] = 1931.85 - 38.26 * c['mean_bore']
    m = re.match(r"Table 4\.\d+: (.+?) Pipe Stem Counts?", c['caption'])
    name = m.group(1)
    m2 = re.match(r"(.+?) \((.+?)\)$", name)
    c['site'], c['site_number'] = (m2.group(1), m2.group(2)) if m2 else (name, '')

# join on sample size + the Binford date the counts imply
used, out = set(), []
for c in cnt:
    cand = [(abs(c['binford_calc'] - int(a['binford'])), i)
            for i, a in enumerate(app) if i not in used and int(a['n']) == c['n']]
    d, i = min(cand) if cand else (None, None)
    if i is not None:
        used.add(i)
    out.append((c, app[i] if i is not None else None, d))

def bounds(s):
    s = s.strip()
    if s.startswith('MCD'):
        y = int(re.search(r'\d{4}', s).group()); return y, y, 'mean ceramic date'
    m = re.match(r"^(\d{4})'s-(\d{4})'s$", s)
    if m: return int(m.group(1)), int(m.group(2)) + 10, 'decades'
    m = re.match(r"^(\d{4})'s$", s)
    if m: return int(m.group(1)), int(m.group(1)) + 10, 'decade'
    m = re.match(r"^(\d{4})-(\d{4})'s$", s)
    if m: return int(m.group(1)), int(m.group(2)) + 10, 'range to decade'
    m = re.match(r"^(\d{4})-(\d{4})$", s)
    if m: return int(m.group(1)), int(m.group(2)), 'range'
    m = re.match(r"^(\d{4})$", s)
    if m: return int(m.group(1)), int(m.group(1)), 'point'
    return None, None, 'UNPARSED'

cols = (['state', 'site', 'site_number', 'context', 'date_raw', 'date_start', 'date_end',
         'date_type', 'n'] + ['b' + b[0] for b in bores] +
        ['mean_bore', 'binford_calc', 'binford_pub', 'hanson_pub', 'hd_pub', 'mean_date_pub',
         'binford_check'])
rows, probs = [], []
for c, a, d in out:
    s, e, t = bounds(c['date'])
    r = dict(state=a['state'] if a else '', site=c['site'], site_number=c['site_number'],
             context=c['context'], date_raw=c['date'], date_start=s, date_end=e, date_type=t,
             n=c['n'], mean_bore=round(c['mean_bore'], 4),
             binford_calc=round(c['binford_calc'], 1),
             binford_pub=a['binford'] if a else '', hanson_pub=a['hanson'] if a else '',
             hd_pub=a['hd'] if a else '', mean_date_pub=a['mean_date'] if a else '',
             binford_check='ok' if d is not None and d <= 1.5 else 'DIFFERS')
    for b in bores:
        r['b' + b[0]] = c[b] if c[b] else 0
    rows.append(r)
    if a and (s is not None) and abs((s + e) / 2 - int(a['mean_date'])) > 0.5:
        probs.append((c['site'], c['context'], c['date'], (s + e) / 2, a['mean_date']))

with open('data/pipe_stems_mcmillan2010.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader(); w.writerows(rows)
print('rows:', len(rows), '| unparsed dates:', sum(r['date_type'] == 'UNPARSED' for r in rows))
print('binford_check DIFFERS:', sum(r['binford_check'] == 'DIFFERS' for r in rows))
print('midpoint vs published mean date mismatches:', len(probs))
for p in probs: print('   ', p)
from collections import Counter
print(Counter(r['date_type'] for r in rows))
