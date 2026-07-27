#!/usr/bin/env python
"""Sync the getting-started example section of README.md from the getting
started guide, doc/html/GettingStarted.html.

doc/html/GettingStarted.html is exported from doc/GettingStarted.mlx (in
Matlab: Live Editor > Save > Export to HTML), and is regenerated locally
whenever the live script changes. This script does not run Matlab or execute
anything -- it parses that already-rendered HTML, extracts headings, code,
and their text/figure outputs in document order, saves each figure into
getting-started/*.png, and rewrites the auto-generated block of README.md's
"Examples" section to match.

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
IMAGES_DIR = REPO / "getting-started"
README_PATH = REPO / "README.md"

# Images are embedded via absolute raw.githubusercontent.com URLs (matching
# the logo at the top of the README and the convention used in gx2-py),
# since a repo-relative path wouldn't resolve on renderers with no base URL
# to resolve against (e.g. a File Exchange-synced description).
RAW_BASE = "https://raw.githubusercontent.com/abhranildas/gx2-matlab/main"

BEGIN_MARKER = (
    "<!-- BEGIN GENERATED: getting-started "
    "(do not edit by hand; regenerate with `python scripts/build_getting_started.py`) -->"
)
END_MARKER = "<!-- END GENERATED: getting-started -->"

# Live-script headings render as h1 (title, skipped -- README's own header
# covers it) / h2 / h3; demoted by one level to sit under README's own
# "## Examples" heading.
HEADING_MARKDOWN = {"h2": "###", "h3": "####"}

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
    """Replace Matlab's rendered-equation <img> spans with their source LaTeX."""
    for eq in tag.select("span[texencoding]"):
        eq.replace_with(f"${eq['texencoding']}$")
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
                lines.append(f"![{alt}]({RAW_BASE}/getting-started/{filename})")
            else:
                text = render_output(out_container)
                if text:
                    lines.append("```")
                    lines.append(text)
                    lines.append("```")
            first_comment[0] = None
    flush_code()
    return "\n".join(lines)


def build_examples_markdown(soup: BeautifulSoup) -> tuple[str, set[str]]:
    content = soup.select_one("div.rtcContent")
    out_lines = []
    fig_counter = [0]
    used_files: set[str] = set()
    heading_slug = ""
    # the guide's h1 intro (title, `doc ...` list, author/citation) is
    # already covered by the README's own header -- start after it, at the
    # first h2.
    started = False
    for child in content.children:
        if not isinstance(child, Tag):
            continue
        name = child.name
        if name in ("h1", "h2", "h3"):
            if name == "h1":
                started = False
                continue
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
        elif name == "div" and ("S1" in classes or "S5" in classes):
            text = latexify(child).strip()
            if text:
                out_lines.append(text)
                out_lines.append("")
        elif name == "ol":
            for li in child.find_all("li"):
                out_lines.append(f"- {latexify(li).strip()}")
            out_lines.append("")
        elif name == "br":
            continue
        else:
            print(f"warning: unhandled top-level element <{name}> class={classes}", file=sys.stderr)

    return "\n".join(out_lines).rstrip() + "\n", used_files


def update_readme(examples_markdown: str) -> None:
    readme = README_PATH.read_text(encoding="utf-8")
    pattern = re.compile(re.escape(BEGIN_MARKER) + r".*?" + re.escape(END_MARKER), re.DOTALL)
    if pattern.search(readme):
        replacement = f"{BEGIN_MARKER}\n\n{examples_markdown}\n{END_MARKER}"
        readme = pattern.sub(lambda _: replacement, readme)
    else:
        anchor = "The following are the worked examples from it.\n"
        if anchor not in readme:
            raise SystemExit(
                "README.md is missing the anchor line to insert the generated block after:\n"
                f"{anchor!r}"
            )
        block = f"\n{BEGIN_MARKER}\n\n{examples_markdown}\n{END_MARKER}\n"
        readme = readme.replace(anchor, anchor + block, 1)
    README_PATH.write_text(readme, encoding="utf-8")


def main() -> None:
    IMAGES_DIR.mkdir(exist_ok=True)
    # start from a clean folder each run, so a renamed/removed plot leaves no
    # stale, no-longer-referenced image behind
    for png in IMAGES_DIR.glob("*.png"):
        png.unlink()

    soup = BeautifulSoup(HTML_PATH.read_text(encoding="utf-8"), "html.parser")
    examples_markdown, used_files = build_examples_markdown(soup)
    update_readme(examples_markdown)
    print(
        f"Updated {README_PATH.relative_to(REPO)} and {len(used_files)} "
        f"image(s) in {IMAGES_DIR.relative_to(REPO)}/"
    )


if __name__ == "__main__":
    main()
