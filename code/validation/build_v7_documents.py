"""Build V7 DOCX with md-to-word and render both deliverables to PDF with Word."""
from pathlib import Path
import json
import subprocess
import sys

import win32com.client
from docx import Document

PYTHON = Path(sys.executable)
ENGINE = PYTHON.parent / "Scripts" / "md-to-word.exe"
ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "paper" / "delivery_v7_final" / "source"
DELIVERY = SOURCE.parent
LOGS = ROOT / "analysis" / "v7_work" / "document_build"
LOGS.mkdir(parents=True, exist_ok=True)

if not ENGINE.is_file():
    raise FileNotFoundError(ENGINE)

builds = []
for stem in ("manuscript", "supplement"):
    source = SOURCE / f"{stem}_word.md"
    output = DELIVERY / f"{stem}.docx"
    command = [str(ENGINE), source.name, "-r", "format.txt", "-o", str(output), "--json"]
    result = subprocess.run(command, cwd=SOURCE, capture_output=True, encoding="utf-8", errors="strict", timeout=300)
    (LOGS / f"{stem}_stdout.json").write_text(result.stdout, encoding="utf-8")
    (LOGS / f"{stem}_stderr.txt").write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(f"{stem} md-to-word failed ({result.returncode}): {result.stdout}\n{result.stderr}")
    payload = json.loads(result.stdout)
    if not payload.get("ok") or not output.is_file():
        raise RuntimeError(f"{stem} export did not pass its gates: {payload}")
    builds.append({"document": stem, "engine": payload, "docx": str(output)})

# Keep each image with its own caption, but do not chain a caption to the next image.
for item in builds:
    docx = Path(item["docx"])
    document = Document(docx)
    adjusted = 0
    for paragraph in document.paragraphs:
        if paragraph.style.name == "Image Caption":
            paragraph.paragraph_format.keep_with_next = False
            adjusted += 1
    document.save(docx)
    item["image_captions_unchained"] = adjusted

word = win32com.client.DispatchEx("Word.Application")
word.Visible = False
word.DisplayAlerts = 0
try:
    for item in builds:
        docx = Path(item["docx"])
        pdf = docx.with_suffix(".pdf")
        document = word.Documents.Open(str(docx), ReadOnly=True)
        try:
            document.ExportAsFixedFormat(str(pdf), 17)
        finally:
            document.Close(False)
        if not pdf.is_file() or pdf.stat().st_size == 0:
            raise RuntimeError(f"PDF export failed: {pdf}")
        item["pdf"] = str(pdf)
finally:
    word.Quit()

summary = {"status": "PASS", "documents": builds}
(LOGS / "build_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps({"status": "PASS", "outputs": [{"docx": b["docx"], "pdf": b["pdf"]} for b in builds]}, ensure_ascii=False, indent=2))
