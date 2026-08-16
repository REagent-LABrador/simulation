"""Strip an mmCIF down to the categories `StructureContext` actually reads.

`ligand_filter._MMCIF_CONTEXT_CATEGORIES` is six header categories. `_atom_site`
is not among them and is ~95% of a typical entry file, so the context fixtures
do not need coordinates at all.

That matters here for a specific reason: `scripts/deploy.ts` zips each
`.claude/skills/<dir>/` **whole, with no exclusions**, so anything parked beside
a skill ships to the Skills API. Full entry mmCIFs for the seven context cases
are ~3.6 MB; trimmed they are ~0.3 MB, and the harness produces identical
verdicts either way (verify with `test_v2.py` — the `context` block must stay
9/9).

    python3 trim_cif.py <in.cif> <out.cif>
    python3 trim_cif.py --check <trimmed.cif>     # list categories retained

Keeps the `data_` header, every wanted category (loop_ and key-value forms), and
drops everything else. Provenance stays readable: the retained `_struct_conn`
rows are the covalent-linkage evidence the polymer-conjugate rules turn on, so a
reader can still audit a verdict from the file.
"""
import pathlib
import sys

WANTED = (
    "_entity.",
    "_entity_poly.",
    "_entity_poly_seq.",
    "_struct_asym.",
    "_struct_conn.",
    "_struct_ref.",
)


def _is_tag(line):
    return line.startswith("_")


def _category_of(tag_line):
    tag = tag_line.split()[0]
    return tag[: tag.index(".") + 1] if "." in tag else tag


def trim(text):
    lines = text.splitlines()
    out, i, n = [], 0, len(lines)
    while i < n:
        line = lines[i]
        stripped = line.strip()

        if stripped.startswith("data_") or stripped == "#":
            out.append(line)
            i += 1
            continue

        if stripped == "loop_":
            # Collect the tag block, decide once, then take/skip the rows.
            j = i + 1
            tags = []
            while j < n and lines[j].strip().startswith("_"):
                tags.append(lines[j])
                j += 1
            keep = bool(tags) and _category_of(tags[0].strip()) in WANTED
            body = []
            while j < n:
                s = lines[j].strip()
                if s == "loop_" or (s.startswith("_") and not s.startswith("_ ")):
                    break
                if s == "#":
                    break
                body.append(lines[j])
                j += 1
            if keep:
                out.append(line)
                out.extend(tags)
                out.extend(body)
                out.append("#")
            i = j
            continue

        if _is_tag(stripped):
            cat = _category_of(stripped)
            block = [line]
            j = i + 1
            # A value may run onto following lines (semicolon text field).
            if j < n and lines[j].startswith(";"):
                block.append(lines[j])
                j += 1
                while j < n and not lines[j].startswith(";"):
                    block.append(lines[j])
                    j += 1
                if j < n:
                    block.append(lines[j])
                    j += 1
            if cat in WANTED:
                out.extend(block)
            i = j
            continue

        i += 1

    return "\n".join(out) + "\n"


def main(argv):
    if len(argv) == 3 and argv[1] == "--check":
        text = pathlib.Path(argv[2]).read_text(errors="replace")
        cats = sorted({
            _category_of(ln.strip())
            for ln in text.splitlines() if ln.startswith("_")
        })
        print(f"{argv[2]}: {len(text):,} bytes")
        for c in cats:
            print(f"  {c}{'' if c in WANTED else '   <- NOT a context category'}")
        return 0
    if len(argv) != 3:
        print(__doc__)
        return 2
    src, dst = pathlib.Path(argv[1]), pathlib.Path(argv[2])
    before = src.read_text(errors="replace")
    after = trim(before)
    dst.write_text(after)
    print(f"{src.name}: {len(before):,} -> {len(after):,} bytes "
          f"({100 * len(after) / len(before):.1f}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
