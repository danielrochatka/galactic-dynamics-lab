# Windows setup

## Recommended approach

For Windows users, the recommended way to build and run Galactic Dynamics Lab is **WSL (Windows Subsystem for Linux)**.

This repository has one authoritative simulator runtime in `engine/` and currently uses a Linux/macOS-style build flow (`cd engine && make`). WSL matches that workflow directly and avoids Windows-native toolchain drift.

> Status note: native Windows/MSVC support is not the documented path in this repository at this time. Use WSL for the supported workflow.

## 1) Install WSL

From **PowerShell (Administrator)**:

```powershell
wsl --install
```

Then:

1. Reboot if prompted.
2. Launch Ubuntu (or your selected distro) from the Start menu.
3. Create your Linux username/password when prompted.

If WSL is already installed, open your WSL shell and continue.

## 2) Install dependencies in WSL

Inside the WSL terminal:

```bash
sudo apt update
sudo apt install -y build-essential python3 python3-pip git
```

These cover the current build/run/plot workflow:
- `build-essential` (includes `make`, `g++`)
- `python3` + `python3-pip` (tooling/plot scripts)
- `git` (clone/update)

## 3) Clone the repo

Inside WSL:

```bash
git clone https://github.com/danielrochatka/galactic-dynamics-lab.git
cd galactic-dynamics-lab
```

## 4) Build the engine

From repo root in WSL:

```bash
cd engine
make
```

This builds the C++ simulator runtime in `engine/`.

## 5) Run a config

From `engine/`, run using a root config file:

```bash
./galaxy_sim --config=../configs/example.cfg
```

Outputs are written under root `outputs/<run_id>/`.

## 6) Plot results

From repo root:

```bash
cd ..
python3 plot_cpp_run.py outputs/<run_id>
```

Replace `<run_id>` with the folder created by your run (for example, check `outputs/`).

## 7) Run tests

Engine tests:

```bash
cd engine && make test
```

Python tooling tests (from repo root):

```bash
pytest -q python_tests
```

If you use the project helper script, you can also run:

```bash
./run_tests.sh
```

## 8) Common Windows pitfalls

- **Running commands in PowerShell instead of WSL**: run build/run/test commands inside your WSL Linux shell.
- **Path confusion (`C:\...` vs `/home/...`)**: keep the command context clear; WSL paths are Linux-style.
- **Missing `make` / `g++`**: install `build-essential` in WSL.
- **Python environment mismatch**: use WSL `python3`/`pip3`, not a separate Windows Python install, for repo tooling commands.
- **Slow I/O on mounted Windows paths**: for better performance, keep the repo in your WSL Linux home directory (for example, `/home/<user>/...`) rather than under `/mnt/c/...`.

## Operational surfaces (restructured repo)

- Authoritative runtime: `engine/`
- Authoritative configs: root `configs/`
- Authoritative outputs: root `outputs/`
- Python: tooling only (plotting/analysis/tests), not a simulation runtime.
