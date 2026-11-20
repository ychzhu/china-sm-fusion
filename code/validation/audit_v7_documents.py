"""Structural and rendered-output audit for the final V7 manuscript package."""
from pathlib import Path
import json
import re

import pymupdf
from docx import Document

ROOT = Path(__file__).resolve().parents[3]
V7 = ROOT / "paper" / "delivery_v7_final"
SOURCE = V7 / "source"


def document_text(doc):
    chunks = [p.text for p in doc.paragraphs]
    for table in doc.tables:
        for row in table.rows:
            chunks.extend(cell.text for cell in row.cells)
    return "\n".join(chunks)


def main():
    checks = {}
    bib_keys = set(re.findall(r"@\w+\{([^,]+),", (SOURCE / "references.bib").read_text(encoding="utf-8")))
    main_md = (SOURCE / "manuscript_word.md").read_text(encoding="utf-8")
    used_keys = {key for key in bib_keys if re.search(rf"@{re.escape(key)}\b", main_md)}
    checks["bibliography_entries"] = len(bib_keys)
    checks["bibliography_has_30_entries"] = len(bib_keys) == 30
    checks["cited_bibliography_keys"] = len(used_keys)
    checks["all_30_bibliography_entries_cited"] = used_keys == bib_keys and len(bib_keys) == 30
    checks["main_source_figures"] = len(re.findall(r"^!\[\]\(\.\./figures/fig\d", main_md, flags=re.M))
    si_md = (SOURCE / "supplement_word.md").read_text(encoding="utf-8")
    checks["si_source_figures"] = len(re.findall(r"^!\[\]\(\.\./figures/figS", si_md, flags=re.M))
    checks["source_placeholders_absent"] = "[[FIGURE" not in main_md + si_md
    checks["abstract_doi_not_duplicated"] = "018rl; https://doi.org/10.57760/sciencedb.018rl" not in main_md
    checks["si_has_no_empty_references_heading"] = "# References" not in si_md

    for name, expected_tables, expected_figures in (("manuscript", 6, 9), ("supplement", 11, 5)):
        doc = Document(V7 / f"{name}.docx")
        text = document_text(doc)
        checks[f"{name}_docx_tables"] = len(doc.tables)
        checks[f"{name}_docx_inline_figures"] = len(doc.inline_shapes)
        checks[f"{name}_docx_table_count_pass"] = len(doc.tables) == expected_tables
        checks[f"{name}_docx_figure_count_pass"] = len(doc.inline_shapes) == expected_figures
        checks[f"{name}_unresolved_citations_absent"] = "[@" not in text
        checks[f"{name}_placeholders_absent"] = "[[FIGURE" not in text
        with pymupdf.open(V7 / f"{name}.pdf") as pdf:
            checks[f"{name}_pdf_pages"] = len(pdf)
            checks[f"{name}_pdf_a4_pages"] = all(abs(page.rect.width - 595.276) < 1 and abs(page.rect.height - 841.89) < 1 for page in pdf)
            pdf_text = "\n".join(page.get_text() for page in pdf)
            if name == "manuscript":
                main_pdf_text = pdf_text
            checks[f"{name}_pdf_unresolved_citations_absent"] = "[@" not in pdf_text
            checks[f"{name}_pdf_placeholders_absent"] = "[[FIGURE" not in pdf_text
    checks["pdf_abstract_doi_not_duplicated"] = "018rl; https://doi.org/10.57760/sciencedb.018rl" not in main_pdf_text

    bool_checks = [value for value in checks.values() if isinstance(value, bool)]
    report = {"status": "PASS" if all(bool_checks) else "FAIL", "checks": checks}
    out = V7 / "WORD_PDF_QA.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
