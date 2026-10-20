#!/usr/bin/env python3
"""Rename a paper PDF to Author_Year, store docling full text next to it, and update the literature index.

Usage:
    register_paper.py --pdf PATH --author SURNAME --year YYYY --docling-md PATH [--literature-dir DIR]
                      [--meta-json JSON] [--dry-run]

Prints the assigned key (e.g. Schreiner_2023 or Schreiner_2023b) and the output paths as JSON.
Only the standard library is used so it runs with any python3.
"""
import argparse
import html
import json
import re
import shutil
import sys
import unicodedata
from pathlib import Path


def slug_surname(name):
    """'Bullón Tarrasó' -> 'BullonTarraso', "O'Neil-Smith" -> 'ONeilSmith', 'van der Berg' -> 'vanderBerg'."""
    ascii_name = unicodedata.normalize('NFKD', name).encode('ascii', 'ignore').decode()
    parts = re.split(r"[\s\-']+", ascii_name.strip())
    out = ''.join(p[:1].upper() + p[1:] if i == 0 or p[:1].isupper() else p for i, p in enumerate(parts) if p)
    return re.sub(r'[^A-Za-z0-9]', '', out)


def load_index(path):
    if path.exists():
        return json.loads(path.read_text())
    return {'papers': []}


def assign_key(base, index, pdf_path, lit_dir):
    """Return base, or base+'b'/'c'... if another paper already owns it. Re-registering the same file keeps its key."""
    taken = {p['key']: p for p in index['papers']}
    for suffix in [''] + [chr(c) for c in range(ord('b'), ord('z') + 1)]:
        key = base + suffix
        owner = taken.get(key)
        if owner is None and not (lit_dir / f'{key}.pdf').exists():
            return key
        if owner is not None and (lit_dir / owner['pdf']).resolve() == pdf_path.resolve():
            return key
        if owner is None and (lit_dir / f'{key}.pdf').resolve() == pdf_path.resolve():
            return key
    sys.exit(f'Too many papers with key {base}')


def clean_markdown(text):
    """Undo docling's escaping so the full text reads naturally."""
    text = html.unescape(text)
    text = re.sub(r'[\x00-\x08\x0b-\x1f\x7f]', '', text)  # stray control chars (e.g. from PDF ligatures) make grep see binary
    text = text.replace('\\_', '_')
    text = re.sub(r'^(1234567890\(\):,;\s*)+', '', text)  # Nature-style margin noise at the top
    text = re.sub(r'(?m)^\s*\d{1,4}\s*$', '', text)  # manuscript line numbers / page numbers on their own line
    text = re.sub(r'\n{3,}', '\n\n', text)
    return text.strip() + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pdf', required=True)
    ap.add_argument('--author', required=True, help='first-author surname as printed in the paper')
    ap.add_argument('--year', required=True)
    ap.add_argument('--docling-md', required=True, help='markdown saved by docling (save_docling_document)')
    ap.add_argument('--literature-dir', default=None, help='defaults to the folder containing the PDF')
    ap.add_argument('--meta-json', default='{}', help='JSON object merged into the index entry (title, authors, doi, ...)')
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()

    pdf = Path(a.pdf).resolve()
    lit = Path(a.literature_dir).resolve() if a.literature_dir else pdf.parent
    notes = lit / 'notes'
    index_path = lit / 'index.json'
    index = load_index(index_path)

    key = assign_key(f'{slug_surname(a.author)}_{a.year}', index, pdf, lit)
    new_pdf = lit / f'{key}.pdf'
    fulltext = notes / f'{key}.fulltext.md'
    summary = notes / f'{key}.md'

    result = {'key': key, 'pdf': str(new_pdf), 'fulltext_md': str(fulltext), 'summary_md': str(summary),
              'original_filename': pdf.name}
    if a.dry_run:
        print(json.dumps(result, indent=2))
        return

    notes.mkdir(parents=True, exist_ok=True)
    if pdf != new_pdf:
        shutil.move(str(pdf), str(new_pdf))
    fulltext.write_text(clean_markdown(Path(a.docling_md).read_text()))

    entry = next((p for p in index['papers'] if p['key'] == key), None)
    if entry is None:
        entry = {'key': key, 'original_filename': pdf.name}
        index['papers'].append(entry)
    entry.update(json.loads(a.meta_json))
    entry.update({'key': key, 'pdf': new_pdf.name, 'fulltext_md': f'notes/{fulltext.name}',
                  'summary_md': f'notes/{summary.name}'})
    index['papers'].sort(key=lambda p: p['key'].lower())
    index_path.write_text(json.dumps(index, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
