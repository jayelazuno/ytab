# Testing the YTAB app

Thank you for trying this. You need **Docker** and about **8 GB of free disk**: the image is 5.6 GB and the clone
is 2.4 GB. Everything else — R, Shiny, Python and the analysis tools — is inside the image, so nothing is installed
on your machine and nothing is left behind except the clone and the image.

## Run it

You need a container runtime. Any of these works, and the launcher uses whichever it finds:

| | |
|---|---|
| [Docker Desktop](https://docs.docker.com/get-started/get-docker/) | the usual choice; free for personal, educational and small-business use, but larger organizations need a paid subscription — worth checking with your institution |
| [Podman Desktop](https://podman-desktop.io) | open source, no subscription |
| [Colima](https://github.com/abiosoft/colima) | open source, no subscription; `brew install colima docker` then `colima start` |

Install one and start it. Then download YTAB from
the [releases page](https://github.com/jayelazuno/ytab/releases) — "Source code (zip)", about 0.8 GB — and unzip
it. No git needed. Then double-click the launcher for your computer:

- **macOS**: `YTAB.command`
- **Windows**: `YTAB.bat`

It checks Docker, downloads YTAB the first time, starts the app and opens your browser. Leave the window it opens
alone while you work; closing it stops YTAB. Your results are saved in the YTAB folder under `output/`.

On macOS a downloaded file may refuse to open. If double-clicking does nothing, open Terminal in the YTAB folder
and run `chmod +x YTAB.command`, then try again.

## Analysing your own reads

YTAB can analyse any library of its own: the clone carries genomes and bowtie2 indexes for eight species
(*C. glabrata*, *C. albicans*, *C. auris*, *S. cerevisiae*, *S. pombe*, *C. neoformans*, *Y. lipolytica*,
*K. phaffii*) and ten library designs across Hermes, Ac/Ds and piggyBac.

**Start from FASTQ › Browse local folders** opens on **My files**, your home folder, so you can navigate to
wherever your reads are kept. The launcher shares it with YTAB read-only: YTAB only ever reads your FASTQs, and
everything it writes goes to the YTAB folder.

To share a narrower folder instead of your whole home folder, set `YTAB_FILES` first:

```bash
YTAB_FILES=~/sequencing open YTAB.command      # macOS
```

Running the command yourself, hand it whatever folder you like and it appears as **Mounted data**:

```bash
docker run --rm -p 3838:3838 -v "$PWD:/opt/ytab" -v /path/to/your/fastqs:/data:ro \
  ghcr.io/jayelazuno/ytab-app:v0.1
```

(This is a Docker matter: a container only sees folders handed to it, which is why the app has to be told where
your files are. Running YTAB without Docker, it simply sees them.)

If you would rather type the commands yourself:

```bash
git clone --depth 1 https://github.com/jayelazuno/ytab.git
cd ytab
docker run --rm -p 3838:3838 -v "$PWD:/opt/ytab" ghcr.io/jayelazuno/ytab-app:v0.1
```

Then open <http://localhost:3838>. Stop it with Ctrl-C.

The first run downloads the image: about four minutes on a fast connection. After that it starts in seconds, and the
first page takes 20-30 seconds to appear while R loads the app. On an Apple Silicon Mac Docker runs the image under
emulation, so that first page is slower still; that is expected and does not affect what you see.

`--depth 1` clones only the current state rather than the full history: about 0.9 GB to download, 2.4 GB on disk.
The clone carries one worked dataset — Gale et al. 2020, the screen the tutorial follows — already processed, so
there is data to look at straight away. It also carries a sample of that study's reads (200,000 pairs per library),
so **Run › Mapping** can be tried without supplying your own. The other published datasets are downloaded per
dataset from the Datasets page rather than shipped with the app.

## What to look at

The app opens on a landing page listing the projects in your clone. Everything is read-only browsing unless you start
an analysis, and nothing you do affects anyone else's copy.

Worth exercising:

- **Quality Control** — the library table and its flags; each subtab's plots (Reads & Mapping, Insertions & Genes,
  Saturation & Bias, Library Agreement); hover a flagged cell for the reason. **Build QC table** recomputes it, which
  is how newer metrics such as the reads from the transposon donor appear in a stored run.
- **Normalization** — pick a run, compare methods, try the Top sites subtab.
- **Condition Comparisons** — treated vs control with its three methods, their options and the plots; try the gene
  labels and the downloads; then the dose-series and genetic-interaction views.
- **Essentiality**, **Gene Explorer**, **Screen Planning**.
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

**Port 3838 already in use** — the launcher moves to the next free port by itself and prints the address. Running
the command by hand, pick another: `-p 8080:3838`, then open <http://localhost:8080>.

**"Start from FASTQ" cannot see your sequencing folder** — the container sees only what is mounted into it, which is
the clone. To try the FASTQ workflow on your own reads, mount their folder as well and browse to `/data`:

```bash
docker run --rm -p 3838:3838 -v "$PWD:/opt/ytab" -v /path/to/your/fastqs:/data:ro \
  ghcr.io/jayelazuno/ytab-app:v0.1
```

**The page does not load** — give it a minute on the first run, then check the terminal for errors.

**Docker says permission denied or cannot connect** — the Docker engine is not running. Start Docker Desktop, or
`colima start` if you use Colima, and try again.

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
