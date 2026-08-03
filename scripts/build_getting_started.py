#!/usr/bin/env python
"""Regenerate README.md in full from the getting started guide,
doc/html/GettingStarted.html.

doc/html/GettingStarted.html is exported from doc/GettingStarted.mlx (in
Matlab: Live Editor > Save > Export to HTML), and is regenerated locally
whenever the live script changes. This script does not run Matlab or execute
anything -- it parses that already-rendered HTML, extracts headings, text,
lists, code and their text/figure outputs in document order, saves each
figure into doc/getting-started/*.png, and rewrites README.md entirely: a
small fixed header (the logo image and File Exchange badge -- the one thing
that can't live inside a Live Script, which has no way to embed a clickable
image) followed by the live script's own content verbatim, one level down
from its title (title -> "#", used by the header above; everything else ->
"##"/"###"). README.md has no other hand-maintained content of its own:
every section lives in the live script, so there is exactly one place to
edit each piece of documentation, never both.

Run this after re-exporting the live script to HTML:

    python scripts/build_getting_started.py

CI (.github/workflows/regen-docs.yml) also runs this on every push that
touches doc/html/GettingStarted.html, committing the result back to main.
"""
from __future__ import annotations

import base64
import re
import sys
from pathlib import Path

from bs4 import BeautifulSoup, Tag

REPO = Path(__file__).resolve().parent.parent
HTML_PATH = REPO / "doc" / "html" / "GettingStarted.html"
IMAGES_DIR = REPO / "doc" / "getting-started"
README_PATH = REPO / "README.md"

# Images are embedded via absolute raw.githubusercontent.com URLs (matching
# the logo at the top of the README and the convention used in gx2-py),
# since a repo-relative path wouldn't resolve on renderers with no base URL
# to resolve against (e.g. a File Exchange-synced description).
RAW_BASE = "https://raw.githubusercontent.com/abhranildas/gx2-matlab/main"

# The only hand-maintained content in README.md. A Live Script has no way to
# embed a *clickable* image (only plain text can be a hyperlink), so the logo
# and File Exchange badge -- both linked images -- can't be authored inside
# doc/GettingStarted.mlx itself and live here instead. Everything else in
# README.md comes from the live script.
README_HEADER = (
    "<!-- This file is generated in full from doc/GettingStarted.mlx (via its "
    "HTML export, doc/html/GettingStarted.html) -- do not edit by hand, other "
    "than the logo/badge line below; regenerate with "
    "`python scripts/build_getting_started.py` -->\n\n"
    '<p align="center">\n'
    f'  <img src="{RAW_BASE}/gx2_icon.png" alt="gx2" width="260">\n'
    "</p>\n\n"
    "# Generalized chi-square distribution "
    "[![View Generalized chi-square distribution on File Exchange]"
    "(https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)]"
    "(https://www.mathworks.com/matlabcentral/fileexchange/85028-generalized-chi-square-distribution)\n\n"
)

# Live-script headings render as h1 (title -- skipped, since README's own
# header above already covers it) / h2 / h3, kept at the same level: h2
# sections sit directly under the README's title, exactly as they do under
# the live script's.
HEADING_MARKDOWN = {"h2": "##", "h3": "###"}

VALUE_CLASSES = {
    "embeddedOutputsVariableMatrixElement",
    "embeddedOutputsVariableStringElement",
    "embeddedOutputsVariableElement",
    "embeddedOutputsTextElement",
}
DIAGNOSTIC_CLASSES = {"embeddedOutputsWarningElement", "embeddedOutputsErrorElement"}


def slugify(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


def latexify(tag: Tag) -> str:
    """Replace Matlab's rendered-equation <img> spans with their source LaTeX,
    inline-styled (bold/italic/monospace) spans with markdown syntax, and
    hyperlinks with markdown link syntax."""
    for eq in tag.select("span[texencoding]"):
        eq.replace_with(f"${eq['texencoding']}$")
    for span in tag.select("span[style]"):
        text = span.get_text()
        style = span["style"]
        if not text or not isinstance(style, str):
            continue
        if "monospace" in style:
            text = f"`{text}`"
        if "font-weight: bold" in style:
            text = f"**{text}**"
        if "font-style: italic" in style:
            text = f"_{text}_"
        span.replace_with(text)
    for a in tag.select("a[href]"):
        a.replace_with(f"[{a.get_text()}]({a['href']})")
    return tag.get_text()


def flatten_scaling_factors(el: Tag) -> None:
    """Collapse Matlab's '10<sup>4</sup> x' scale-factor markup into one '10^4 x' run."""
    for sf in el.select(".veScalingFactor"):
        base = next((c for c in sf.contents if isinstance(c, str)), "10").strip()
        sup = sf.find("sup")
        exp = sup.get_text().strip() if sup else ""
        mult = sf.select_one(".multiply")
        symbol = mult.get_text().strip() if mult else "×"
        sf.replace_with(f"{base}^{exp} {symbol}")


def clean_block_text(el: Tag) -> str:
    el = el.__copy__()
    for layer in el.select(".doNotExport"):
        layer.decompose()
    for summ in el.select(".veVariableValueSummary"):
        summ.decompose()
    flatten_scaling_factors(el)
    text = el.get_text("\n")
    text = re.sub(r"[ \t]+\n", "\n", text)
    text = re.sub(r"\n{2,}", "\n\n", text).strip("\n")
    return text


def figure_bytes(out_div: Tag) -> bytes | None:
    img = out_div.select_one("img.figureImage")
    if img is None:
        return None
    return base64.b64decode(img["src"].split(",", 1)[1])


def render_output(out_div: Tag) -> str | None:
    """Render a non-figure output div (variable/text/warning/error) as text.

    Collapses consecutive duplicate warning/error blocks, which Matlab's HTML
    export otherwise repeats once per internal occurrence of the same message.
    """
    blocks = []
    last_key = None
    for wrapper in out_div.select(".eoOutputWrapper"):
        classes = set(wrapper.get("class", []))
        if classes & DIAGNOSTIC_CLASSES:
            key_kind = "diag"
        elif classes & VALUE_CLASSES:
            key_kind = "val"
        else:
            continue
        text = clean_block_text(wrapper)
        key = (key_kind, text)
        if key == last_key:
            continue
        last_key = key
        if text:
            blocks.append(text)
    return "\n\n".join(blocks) if blocks else None


def process_codeblock(
    block: Tag, heading_slug: str, fig_counter: list[int], used_files: set[str]
) -> str:
    lines = []
    pending: list[str] = []
    first_comment = [None]

    def flush_code():
        if pending:
            code = "".join(pending).rstrip("\n")
            lines.append("```matlab")
            lines.append(code)
            lines.append("```")
            pending.clear()

    for wrapper in block.find_all("div", class_="inlineWrapper", recursive=False):
        children = wrapper.find_all("div", recursive=False)
        if not children:
            continue
        code_container = children[0]
        out_container = children[1] if len(children) > 1 else None

        line_text = latexify(code_container)
        for raw in line_text.split("\n"):
            stripped = raw.strip()
            if stripped.startswith("%") and first_comment[0] is None:
                first_comment[0] = stripped.lstrip("%").strip()
        pending.append(line_text if line_text.endswith("\n") else line_text + "\n")

        if out_container is not None:
            flush_code()
            png = figure_bytes(out_container)
            if png is not None:
                fig_counter[0] += 1
                stem = f"{heading_slug}-{fig_counter[0]}" if heading_slug else f"fig-{fig_counter[0]}"
                filename = f"{stem}.png"
                used_files.add(filename)
                (IMAGES_DIR / filename).write_bytes(png)
                alt = first_comment[0] or f"Plot output {fig_counter[0]}"
                lines.append(f"![{alt}]({RAW_BASE}/doc/getting-started/{filename})")
            else:
                text = render_output(out_container)
                if text:
                    lines.append("```")
                    lines.append(text)
                    lines.append("```")
            first_comment[0] = None
    flush_code()
    return "\n".join(lines)


def build_readme_markdown(soup: BeautifulSoup) -> tuple[str, set[str]]:
    content = soup.select_one("div.rtcContent")
    out_lines = []
    fig_counter = [0]
    used_files: set[str] = set()
    heading_slug = ""
    # the live script's own title (h1) is already covered by the README's
    # own header above -- skip emitting it, but still start capturing every
    # paragraph/list/code block that follows it.
    started = False
    for child in content.children:
        if not isinstance(child, Tag):
            continue
        name = child.name
        if name == "h1":
            started = True
            continue
        if name in ("h2", "h3"):
            started = True
            heading_slug = slugify(latexify(child).strip())
            out_lines.append(f"{HEADING_MARKDOWN[name]} {latexify(child).strip()}")
            out_lines.append("")
            continue
        if not started:
            continue
        classes = child.get("class") or []
        if "CodeBlock" in classes:
            out_lines.append(process_codeblock(child, heading_slug, fig_counter, used_files))
            out_lines.append("")
        elif name == "div":
            # Matlab's HTML export assigns each distinct paragraph style a
            # per-document "Sn" class number, which shifts whenever content
            # is added/removed elsewhere -- so match any non-CodeBlock div
            # as a plain text paragraph, rather than specific Sn numbers.
            text = latexify(child).strip()
            if text:
                out_lines.append(text)
                out_lines.append("")
        elif name in ("ol", "ul"):
            for i, li in enumerate(child.find_all("li"), start=1):
                marker = f"{i}." if name == "ol" else "-"
                out_lines.append(f"{marker} {latexify(li).strip()}")
            out_lines.append("")
        elif name == "br":
            continue
        else:
            print(f"warning: unhandled top-level element <{name}> class={classes}", file=sys.stderr)

    return "\n".join(out_lines).rstrip() + "\n", used_files


def write_readme(body_markdown: str) -> None:
    README_PATH.write_text(README_HEADER + body_markdown, encoding="utf-8")


def main() -> None:
    IMAGES_DIR.mkdir(exist_ok=True)
    # start from a clean folder each run, so a renamed/removed plot leaves no
    # stale, no-longer-referenced image behind
    for png in IMAGES_DIR.glob("*.png"):
        png.unlink()

    soup = BeautifulSoup(HTML_PATH.read_text(encoding="utf-8"), "html.parser")
    body_markdown, used_files = build_readme_markdown(soup)
    write_readme(body_markdown)
    print(
        f"Updated {README_PATH.relative_to(REPO)} and {len(used_files)} "
        f"image(s) in {IMAGES_DIR.relative_to(REPO)}/"
    )


if __name__ == "__main__":
    main()
