"""Flat, self-contained LaTeX source for journal submission.

Run from the project root after `make paper`:
    python3 paper/py/make_submission_source.py OUTDIR

The submission system parses the uploaded .tex for the title, authors and abstract,
so the source must not depend on generated macros or \\input files. This script
  1. inlines every \\input (sections, tables, numbers.tex) and the bibliography (.bbl);
  2. replaces each \\Num... macro with its literal value and drops the definitions;
  3. renames figures to their number in the manuscript (Figure1.pdf, ...), flat;
  4. copies the MDPI class folder (Definitions/), which the template expects as is.
The result compiles with a single `pdflatex` run pair and no BibTeX.
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper"
out = Path(sys.argv[1]).resolve()
out.mkdir(parents=True, exist_ok=True)

# 1. inline \input files and the bibliography; drop comments
tex = subprocess.run(["latexpand", "--empty-comments", "--expand-bbl", "main.bbl", "main.tex"],
                     cwd=PAPER, check=True, capture_output=True, text=True).stdout

# 2. literal values for the number macros (definitions come from numbers.tex)
defs = {}
for line in (PAPER / "numbers.tex").read_text().splitlines():
    m = re.match(r"\\newcommand\{\\(Num[A-Za-z]+)\}\{(.*)\\xspace\}$", line)
    if m:
        defs[m.group(1)] = m.group(2)
tex = "\n".join(l for l in tex.splitlines()
                if not re.match(r"\\newcommand\{\\Num[A-Za-z]+\}\{.*\\xspace\}$", l))


def literal(m):
    name = m.group(1)
    if name not in defs:                     # \NumMissing and similar helpers stay as they are
        return m.group(0)
    return defs[name]


tex = re.sub(r"\\(Num[A-Za-z]+)(?![A-Za-z])(\{\})?", literal, tex)
left = sorted(set(re.findall(r"\\(Num[A-Za-z]+)", tex)) - {"NumMissing"})
if left:
    sys.exit(f"unresolved number macros: {left}")

# 3. figures, numbered in order of appearance
figs = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{figures/([^}]+)\}", tex)
for k, f in enumerate(dict.fromkeys(figs), start=1):
    src = PAPER / "figures" / (f if f.endswith(".pdf") else f + ".pdf")
    shutil.copy(src, out / f"Figure{k}.pdf")
    tex = tex.replace("{figures/" + f + "}", "{Figure%d.pdf}" % k)

# 4. class files
shutil.copytree(PAPER / "Definitions", out / "Definitions", dirs_exist_ok=True)
(out / "manuscript.tex").write_text(tex)
print(f"{out / 'manuscript.tex'}: {len(defs)} macros replaced, {len(set(figs))} figures")
