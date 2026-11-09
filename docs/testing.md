# Testing the YTAB app

Thank you for trying this. You need **Docker** and about **3 GB of disk**; everything else — R, Shiny, Python and the
analysis tools — is inside the image, so nothing is installed on your machine and nothing is left behind except the
clone and the image.

## Run it

```bash
git clone --depth 1 https://github.com/jayelazuno/ytab.git
cd ytab
docker run --rm -p 3838:3838 -v "$PWD:/opt/ytab" ghcr.io/jayelazuno/ytab-app:v0.1
```

Then open <http://localhost:3838>. Stop it with Ctrl-C.

The first run downloads the image (a few minutes). After that it starts in seconds. On an Apple Silicon Mac Docker
runs it under emulation, so the first page may take a little longer to appear; that is expected and does not affect
what you see.

`--depth 1` clones only the current state rather than the full history: about 800 MB instead of 3 GB. The clone
carries demo projects, so there is data to look at straight away.

## What to look at

The app opens on a landing page listing the projects in your clone. Everything is read-only browsing unless you start
an analysis, and nothing you do affects anyone else's copy.

Worth exercising:

- **Quality Control** — the library table and its flags; each subtab's plots (Reads & Mapping, Insertions, Genes,
  Saturation); hover a flagged cell for the reason.
- **Site Normalization** — pick a run, compare methods, try the Top sites subtab.
- **Treated vs Control** — the three methods, their options and the plots; try the gene labels and the downloads.
- **Condition Comparisons**, **Essentiality**, **Comparative View**, **Gene Explorer**.
- Plot customization: the style controls, and whether downloads match what is on screen.

## Telling us what you found

Please open an issue at <https://github.com/jayelazuno/ytab/issues> with:

1. the tab and the exact option you chose,
2. what you expected and what happened,
3. a screenshot if it is a layout or label problem,
4. the text in the terminal where Docker is running, if something failed.

Small things are worth reporting: a confusing label, a plot that is hard to read, an option whose effect is unclear.
The interface is what we are testing.

## If something goes wrong

**"the YTAB repository is not mounted at /opt/ytab"** — run the command from inside your clone; `$PWD` has to be the
`ytab` directory.

**Port 3838 already in use** — pick another: `-p 8080:3838`, then open <http://localhost:8080>.

**The page does not load** — give it a minute on the first run, then check the terminal for errors.

**Docker says permission denied or cannot connect** — Docker Desktop is probably not running; start it and try again.

## Without Docker

Possible, but slower to set up and less predictable, since package versions are resolved on your machine rather than
fixed in the image:

```bash
micromamba create -y -f environment.local.yml   # or: conda env create -f environment.local.yml
micromamba activate ytab-local
Rscript app/shiny/run_app.R
```

If you take this route, run both R and Python from the `ytab-local` environment: the app calls Python for some
actions, and a mismatch there produces errors part-way through a click rather than at startup.
