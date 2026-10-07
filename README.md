# simulSPC

C++ simulation tool for single-molecule confocal microscopy, Fluorescence Correlation Spectroscopy (FCS), and Time-Correlated Single Photon Counting (TCSPC).

It generates time-tagged photon records from simulated diffusing emitters in solution or in a microfluidic channel.

---

## Features

- **Diffusion:** 3D Brownian motion (translation) and rotational diffusion using quaternions.
- **Flow:** Optional laminar flow (rectangular Poiseuille profile) along the microchannel.
- **Optics:** Confocal observation volume modeled as a 3D Gaussian profile on a 2D lookup grid.
- **Photophysics:** Absorption based on dipole orientation, quantum yield, and photon collection efficiency.
- **TCSPC:** Photon arrival times with exponential lifetime decay and Gaussian instrument response (IRF).
- **Polarization:** Two detection channels (parallel and perpendicular) to simulate time-resolved anisotropy.
- **Detector:** Non-paralyzable dead-time applied per channel.
- **Polydispersity:** Support for multiple particle species with size and lifetime distributions.
- **Parallel processing:** Multi-threaded particle updates using OpenMP.
- **Output:** Raw 64-bit binary file (`.ptn`).

---

## Repository Structure

```text
.
├── CMakeLists.txt
├── config/             # Configuration files (.cfg)
├── include/            # C++ headers
├── src/                # C++ source files
└── tests/              # Unit tests

```

---

## Build Instructions

### Requirements

* CMake (>= 3.15)
* C++17 compiler with OpenMP support (GCC, Clang, or MSVC)

### Linux

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

```

Run tests:

```bash
ctest --test-dir build --output-on-failure

```

Executables are generated in `build/` (`build/simulSPC` and test binaries).

### Windows

From the **x64 Native Tools Command Prompt for VS**:

```cmd
cmake -B build
cmake --build build --config Release

```

The executable is generated in `build\Release\simulSPC.exe`.

---

## Usage

Run the program with a configuration file and an output file path:

```bash
# Linux
./build/simulSPC config/parameters.cfg output.ptn

# Windows
.\build\Release\simulSPC.exe config\parameters.cfg output.ptn

```

If no arguments are given, default paths are `config/parameters.cfg` and `photon_list.ptn`.

---

## Binary Output Format

The output file is a raw binary array of 64-bit unsigned integers (`uint64`):

* **Bits 0–39:** Macrotime (clock tick index)
* **Bits 40–55:** Microtime (TCSPC channel, 0 to 4095)
* **Bits 56–63:** Detector channel (0: parallel, 1: perpendicular)

Output data files can be loaded and analyzed using [pySPC](https://github.com/MLoum/pySPC).

```

```
