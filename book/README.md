# Documentation site builder

The site is rendered with [bookdown](https://bookdown.org/) from `docs/*.md` and `CHANGELOG.md`; those files are the only source. Do not edit the rendered HTML.

| File | Purpose |
|---|---|
| `prepare_book.R` | Chapter order; copies the docs into `_build/`, adds heading IDs, rewrites cross-document links, fails on broken links |
| `render_book.R` | Runs `prepare_book.R`, renders the gitbook into `_book/`, checks every internal link in the HTML |
| `index.Rmd` | Home page |
| `_bookdown.yml`, `_output.yml`, `style.css` | bookdown settings and theme |
| `publish_pages.sh` | Fallback: render locally and force-push to the `bookdown` branch |

Build locally with `Rscript book/render_book.R` (needs R, bookdown and pandoc) and open `book/_book/index.html`.

## Publishing

`.github/workflows/docs.yml` renders the book on every push to `main` that changes `docs/` or `book/`, and on release tags, then force-pushes `book/_book/` to the orphan branch `bookdown`. Enable GitHub Pages once: **Settings → Pages → Build and deployment → Deploy from a branch → `bookdown` / `(root)`**. The site is then served at <https://sbresnahan.github.io/placenta-xqtl-federated/>.

To add a chapter, write `docs/<NAME>.md` with one level-1 heading and add it to the `chapters` table in `prepare_book.R`. Link to other pages as you would on GitHub (`[Memory](HPC.md#memory)`); the builder converts the links.
