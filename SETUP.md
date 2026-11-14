# Setup: packages and install commands

What you need to run `cp_visualizer.py` and
`Airfoil_polyline_discretization/airfoil_discretize.py`. Commands are for macOS
(Terminal, zsh). Windows differences are noted where they exist.

## What is needed

Python 3.8 or newer, plus these packages (also listed in `requirements.txt`):

| Package | `cp_visualizer.py` | `airfoil_discretize.py` | Used for |
|---|---|---|---|
| numpy | yes | yes | arrays, interpolation, integration |
| matplotlib | yes | yes | the pressure plot and the outline check plot |
| scipy | no | yes | the spline through the airfoil points |
| ezdxf | no | yes | writing the DXF file |

Everything else the scripts import (`argparse`, `csv`, `math`, `os`, `re`, `sys`)
comes with Python.

Tested with Python 3.13, numpy 2.5, matplotlib 3.11, scipy 1.18 and ezdxf 1.4. The
minimum versions in `requirements.txt` are estimates and have not been tested.

Not Python packages, and not installed by the steps below:

- **XFOIL** produces the `.cp` files that `cp_visualizer.py` reads. The example
  `e423_a12.cp` is already in the repo, so XFOIL is only needed for new cases.
- **MATLAB or GNU Octave** is only for `airfoil_discretize.m`, which needs no
  toolboxes or extra packages.

## Step 0: check Python

```bash
python3 --version
```

If that prints 3.8 or higher, go to Option A or B. If the command is not found,
install Python with Homebrew:

```bash
brew install python
```

## Option A (recommended): a virtual environment

A virtual environment is a private folder of packages for this project only, so
nothing else on your computer is affected. Run these from the repo folder.

Create it (once):

```bash
python3 -m venv .venv
```

Turn it on (every new Terminal window):

```bash
source .venv/bin/activate
```

Install the packages (once):

```bash
python3 -m pip install -r requirements.txt
```

Your prompt shows `(.venv)` while it is on. To turn it off:

```bash
deactivate
```

On Windows, the turn-on command is `.venv\Scripts\activate` and `python3` is
usually just `python`.

## Option B: install globally

No environment to switch on, but the packages are shared with every other Python
project on the computer.

```bash
python3 -m pip install --user -r requirements.txt
```

If this stops with `error: externally-managed-environment`, your Python came from
Homebrew, which blocks global installs. Use Option A instead.

To install only what `cp_visualizer.py` needs:

```bash
python3 -m pip install --user numpy matplotlib
```

## Check that it worked

```bash
python3 -c "import numpy, matplotlib, scipy, ezdxf; print('all packages found')"
```

Then run each tool from the repo folder:

```bash
python3 cp_visualizer.py e423.dat e423_a12.cp --V 34 --spacing 3 --thickness 0.125 --buckets 6 --bucket-spacing adaptive --alpha 12 --save out/e423 --no-show
```

```bash
cd Airfoil_polyline_discretization
```

```bash
python3 airfoil_discretize.py e423.dat --chord 25.5 --wall 0.3 --n 40 --spacing cosine --units in
```

The first ends with `Done: 0 warning(s), 0 error(s).` and the second prints the
leading edge at (0, 0) and the trailing edge at x = 25.5.

## If something goes wrong

| Message | Cause | Fix |
|---|---|---|
| `ModuleNotFoundError: No module named 'scipy'` (or another package) | packages not installed, or the virtual environment is off | `source .venv/bin/activate`, then install again |
| `error: externally-managed-environment` | global install blocked by Homebrew Python | use Option A |
| `command not found: python3` | Python is not installed | `brew install python` |
| `command not found: pip` | pip is not on the path | use `python3 -m pip` as written above |
