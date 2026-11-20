"""Prepare deterministic md-to-word sources from the authoritative V6 Markdown."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[3]
V6 = ROOT / "paper" / "delivery_v6_final"
SOURCE = V6 / "source"

REPLACEMENTS = {
    "(Reichle et al., 2004; Liu et al., 2011; Gruber et al., 2019)": "[@reichle2004; @liu2011; @gruber2019]",
    "(Blyverket et al., 2019; Gruber et al., 2019; Cheng et al., 2023)": "[@blyverket2019; @gruber2019; @cheng2023]",
    "(Blyverket et al., 2019; Gruber et al., 2019)": "[@blyverket2019; @gruber2019]",
    "(Cheng et al., 2023; Song et al., 2022)": "[@cheng2023; @song2022]",
    "(Li et al., 2022; Meng et al., 2021)": "[@li2022; @meng2021]",
    "(Muñoz-Sabater et al., 2021; Rodell et al., 2004)": "[@munozsabater2021; @rodell2004]",
    "(Hersbach et al., 2020; Seo et al., 2021)": "[@hersbach2020; @seo2021]",
    "(Stoffelen, 1998; Gruber et al., 2019)": "[@stoffelen1998; @gruber2019]",
    "(TCA; Stoffelen, 1998; McColl et al., 2014)": "(TCA) [@stoffelen1998; @mccoll2014]",
    "(Li et al., 2022; Song et al., 2022; Meng et al., 2021; Cheng et al., 2023)": "[@li2022; @song2022; @meng2021; @cheng2023]",
    "(Gruber et al., 2019; Muñoz-Sabater et al., 2021)": "[@gruber2019; @munozsabater2021]",
    "(Blyverket et al., 2019)": "[@blyverket2019]",
    "(Brodzik et al., 2012)": "[@brodzik2012]",
    "(Cheng et al., 2023)": "[@cheng2023]",
    "(Colliander et al., 2017)": "[@colliander2017]",
    "(Didan, 2021)": "[@didan2021]",
    "(Entekhabi et al., 2010)": "[@entekhabi2010]",
    "(Gruber et al., 2019)": "[@gruber2019]",
    "(Hersbach et al., 2020)": "[@hersbach2020]",
    "(Li et al., 2022)": "[@li2022]",
    "(Liu et al., 2011)": "[@liu2011]",
    "(McColl et al., 2014)": "[@mccoll2014]",
    "(Meng et al., 2021)": "[@meng2021]",
    "(Muñoz-Sabater et al., 2021)": "[@munozsabater2021]",
    "(O’Neill et al., 2023)": "[@oneill2023]",
    "(Reichle et al., 2004)": "[@reichle2004]",
    "(Rodell et al., 2004)": "[@rodell2004]",
    "(Seo et al., 2021)": "[@seo2021]",
    "(Song et al., 2022)": "[@song2022]",
    "(Stoffelen, 1998)": "[@stoffelen1998]",
    "(Entekhabi et al., 2010; O'Neill et al., 2023)": "[@entekhabi2010; @oneill2023]",
    "(Entekhabi et al., 2010; O’Neill et al., 2023)": "[@entekhabi2010; @oneill2023]",
    "(Hersbach et al., 2020; Rodell et al., 2004)": "[@hersbach2020; @rodell2004]",
    "(Blyverket et al., 2019; Seo et al., 2021)": "[@blyverket2019; @seo2021]",
    "(Hersbach et al., 2020; Muñoz-Sabater et al., 2021)": "[@hersbach2020; @munozsabater2021]",
    "(O'Neill et al., 2023)": "[@oneill2023]",
    "Song et al. (2022)": "@song2022",
    "Meng et al. (2021)": "@meng2021",
    "Cheng et al. (2023)": "@cheng2023",
    "Reichle and Koster (2004)": "@reichle2004",
    "Didan (2021)": "@didan2021",
    "Entekhabi et al. (2010)": "@entekhabi2010",
    "Gruber et al. (2019)": "@gruber2019",
    "Hersbach et al. (2020)": "@hersbach2020",
    "Li et al. (2022)": "@li2022",
    "Liu et al. (2011)": "@liu2011",
    "Muñoz-Sabater et al. (2021)": "@munozsabater2021",
    "O’Neill et al. (2023)": "@oneill2023",
    "O'Neill et al. (2023)": "@oneill2023",
    "Reichle et al. (2004)": "@reichle2004",
}

FIGURES = {
    "**Figure 1.": "fig01_station_map.png",
    "**Figure 2.": "fig02_workflow.png",
    "**Figure 3.": "fig03_fusion_maps.png",
    "**Figure 4.": "fig04_scatter_density.png",
    "**Figure 5.": "fig05_monthly.png",
    "**Figure 6.": "fig06_seasonal.png",
    "**Figure 7.": "fig07_density.png",
    "**Figure 8.": "fig08_ndvi.png",
    "**Figure 9.": "fig09_yearly.png",
}

SI_FIGURES = {
    "[[FIGURE 1]]": "figS1_coverage.png",
    "[[FIGURE 2]]": "figS2_within_cell.png",
    "[[FIGURE 3]]": "figS3_hl_correlation.png",
    "[[FIGURE 4]]": "figS4_aerr_reliability.png",
    "[[FIGURE 5]]": "figS5_reference.png",
}


def convert_citations(text: str) -> str:
    for old, new in sorted(REPLACEMENTS.items(), key=lambda item: len(item[0]), reverse=True):
        text = text.replace(old, new)
    return text


def style_common(text: str) -> str:
    lines = text.splitlines()
    if lines and lines[0].startswith("# "):
        lines[0] = "%% 文档标题 " + lines[0][2:]
    styled = []
    for line in lines:
        numbered = re.match(r"^#{2,6} (\d+(?:\.\d+)*)\s", line)
        if numbered:
            depth = numbered.group(1).count(".") + 1
            line = "#" * depth + " " + line.split(" ", 1)[1]
        elif line.startswith("## "):
            line = line[1:]
        if re.match(r"^(?:\*\*)?Table (?:S?\d+)\.", line):
            if re.match(r"^(?:\*\*)?Table 1\.", line):
                styled.append("{{< pagebreak >}}")
                styled.append("")
            line = "%% 表题 " + line
        elif re.match(r"^(?:\*\*)?Figure S\d+\.", line):
            line = "%% 图片题注 " + line
        styled.append(line)
    return "\n".join(styled)


def main() -> None:
    manuscript = (SOURCE / "manuscript.md").read_text(encoding="utf-8")
    body, references_and_figures = manuscript.split("## References", 1)
    figure_text = references_and_figures.split("## Figure captions", 1)[1]
    body = style_common(convert_citations(body).rstrip())
    figure_text = convert_citations(figure_text).strip()
    for caption, filename in FIGURES.items():
        image = f"![](../figures/{filename}){{width=15.5cm}}"
        if figure_text.count(caption) != 1:
            raise RuntimeError(f"Expected one {caption} caption")
        if not (V6 / "figures" / filename).is_file():
            raise FileNotFoundError(filename)
        figure_text = figure_text.replace(caption, image + "\n\n%% 图片题注 " + caption)
    manuscript_export = (
        "---\nbibliography: references.bib\ncsl: copernicus-publications.csl\n---\n\n"
        + body
        + "\n\n# References\n\n# Figure captions\n\n"
        + figure_text
        + "\n"
    )
    (SOURCE / "manuscript_word.md").write_text(manuscript_export, encoding="utf-8")

    supplement = convert_citations((SOURCE / "supplement.md").read_text(encoding="utf-8"))
    for placeholder, filename in SI_FIGURES.items():
        if supplement.count(placeholder) != 1:
            raise RuntimeError(f"Expected one {placeholder}")
        if not (V6 / "figures" / filename).is_file():
            raise FileNotFoundError(filename)
        supplement = supplement.replace(placeholder, f"![](../figures/{filename}){{width=15.5cm}}")
    supplement = style_common(supplement)
    (SOURCE / "supplement_word.md").write_text(supplement, encoding="utf-8")

    for name, text in (("manuscript", manuscript_export), ("supplement", supplement)):
        author_names = "Blyverket|Brodzik|Cheng|Colliander|Didan|Entekhabi|Gruber|Hersbach|Li|Liu|McColl|Meng|Muñoz-Sabater|O’Neill|O'Neill|Reichle|Rodell|Seo|Song|Stoffelen"
        leftovers = re.findall(rf"(?:{author_names})(?: et al\.| and Koster)?,? \(?(?:19|20)\d{{2}}\)?", text)
        leftovers = [item for item in leftovers if "Figure" not in item and "Table" not in item]
        if leftovers:
            raise RuntimeError(f"Unconverted author-year citations in {name}: {leftovers}")
        if "[[" in text or "]]" in text:
            raise RuntimeError(f"Unresolved placeholder in {name}")
    print("Prepared manuscript_word.md and supplement_word.md")
    print(f"Main figures inserted: {len(FIGURES)}; SI figures inserted: {len(SI_FIGURES)}")
    print("References: 19-entry DOI-backed Copernicus bibliography in the main manuscript")


if __name__ == "__main__":
    main()
