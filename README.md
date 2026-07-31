# Heterozygous EGFR Spheroid 2D Cross Section Agent-Based Model

## Model Description

This model simulates the 2D cross-sectional growth of a cancer spheroid with heterozygous EGFR mutation using an agent-based modelling (ABM) approach implemented in Chaste. Each cell is represented as an individual agent with defined rules governing proliferation, adhesion, and spatial organisation. The model captures emergent spheroid morphology arising from protein (EGFR, TGFα, and HIF) dynamics in presence of heterozygous EGFR for studying tumour growth dynamics in a controlled in-silico environment. To run this simulation, you have follow the instructions below. You may also skip the Docker introduction if you are already familiar with the interface:

---

## Table of Contents

1. [Why Docker?](#why-docker)
2. [Installing Docker](#installing-docker)
3. [Installing Chaste via Docker](#installing-chaste-via-docker)
4. [Container Directory Structure](#container-directory-structure)
5. [Working with Projects](#working-with-projects)
6. [Using VS Code](#using-vs-code-alternative-method)
7. [Test the Project](#test-the-project)
7. [Further Resources](#further-resources)

---

## Why Docker?

Scientific software like Chaste has complex dependency chains - specific versions of compilers, numerical libraries, and system packages must all be present and correctly configured for the code to build and run. Getting this right manually is time-consuming, fragile, and notoriously difficult to reproduce across different machines or operating systems.

**Docker solves this by packaging the entire software environment into a portable, isolated container.** A container behaves like a lightweight virtual machine: it includes the operating system layer, all dependencies, and the application itself, but shares the host machine's kernel for efficiency. This means:

- **Reproducibility** - The same container runs identically on any machine (Linux, macOS, Windows), eliminating "works on my machine" problems.
- **No installation conflicts** - Chaste and its dependencies are fully isolated from your host system. You won't break anything already installed.
- **Fast onboarding** - Instead of spending hours configuring a build environment, you can be running simulations within minutes.
- **Version control for environments** - Docker images are tagged and versioned, so you can pin your project to a specific Chaste release and guarantee that collaborators or reviewers use the exact same software stack.
- **Cleaner host system** - When you're done, simply remove the container. Nothing lingers on your host.

For computational biology projects like this one - where simulation correctness depends heavily on the precise software environment - Docker is not just convenient; it is essential for scientific rigour.

---

## Installing Docker

### Linux [Recommended for Chaste]

Install [Docker Desktop for Linux](https://docs.docker.com/desktop/install/linux-install/). On Linux, Docker containers share all available RAM and processing cores with the host by default - no additional configuration is needed.

### macOS

1. Install [Docker Desktop for Mac](https://docs.docker.com/desktop/install/mac-install/) - select the **Intel** or **Apple Silicon** installer to match your chip.
2. **Increase resource limits.** Docker Desktop restricts RAM and CPU by default. Go to **Settings → Resources** and allocate at least **4 GB RAM** and half your available CPU cores, then click **Apply & Restart**.
3. Open Terminal and proceed with the commands in the [Installing Chaste via Docker](#installing-chaste-via-docker) section.

> **Apple Silicon (M1/M2/M3):** The Chaste image is built for `amd64`. It runs under Rosetta 2 emulation - functionally correct but slower. Add `--platform linux/amd64` to your `docker run` command if you see a platform mismatch warning.
### Windows

Docker on Windows uses WSL2 (Windows Subsystem for Linux 2) as its backend. **Requires Windows 10 64-bit version 22H2 (build 19045) or Windows 11.**

1. **Install WSL2 and Ubuntu.** Open PowerShell as Administrator and run:
   ```powershell
   wsl --install -d ubuntu
   ```
   Restart your machine when prompted. This installs WSL2 and the Ubuntu distribution.

2. **Install [Docker Desktop for Windows](https://docs.docker.com/desktop/install/windows-install/).** During setup, ensure **"Use WSL 2 instead of Hyper-V"** is selected when prompted.

3. **Configure Docker Desktop:**
   - **Open Docker Desktop** → **Settings → General** and confirm **"Use the WSL 2 based engine"** is checked.
   - Docker automatically enables integration with your default WSL distribution (Ubuntu), so `docker` commands should work in the Ubuntu terminal straight away. If not, go to [Docker WSL Integration Guide](https://docs.docker.com/desktop/features/wsl/) and follow their steps.
   > **Note:** The Resources → Advanced panel (for RAM/CPU limits) is not available in WSL2 mode - Windows manages those resources. Use `.wslconfig` instead (see step 4).

4. **Limit WSL2 memory usage** (recommended). By default WSL2 can consume a large portion of your system RAM. To cap it, create a file at `C:\Users\<YourUsername>\.wslconfig` with the following content:
   ```
   [wsl2]
   memory=8GB
   processors=4
   ```
   Adjust values to suit your machine, then restart WSL2 with `wsl --shutdown`.

5. **Launch the Ubuntu app** from the Start menu. All Docker commands below should be run from this Ubuntu terminal. Note: Windows drives are automatically accessible inside WSL2 at `/mnt/c/`, `/mnt/d/` etc. - no additional file sharing configuration is needed.

---

## Installing Chaste via Docker

Once Docker is installed and running, open a terminal (Linux/WSL2 Ubuntu) and run the following command [(see detailed instruction here)](https://github.com/Chaste/chaste-docker):

```bash
docker run --init -it --rm -v chaste_data:/home/chaste chaste/release
```

You should see a prompt like this:

```
chaste@301291afbedf:~$
```

This is a bash shell inside an isolated Docker container. Chaste is fully compiled and ready to use now - no further build steps are required.

**What each flag does:**

| Flag | Purpose |
|------|---------|
| `--init` | Ensures signals (e.g. Ctrl+C) are handled correctly inside the container |
| `-it` | Runs the container interactively with a terminal |
| `--rm` | Automatically removes the container when you exit (your data in the [volume](https://docs.docker.com/engine/storage/volumes/) is preserved) |
| `-v chaste_data:/home/chaste` | Mounts a persistent Docker volume so your work survives container restarts |

### Using a Specific Release

To use a particular Chaste version (recommended for reproducibility), specify a tag:

```bash
docker run --init -it --rm -v chaste_data:/home/chaste chaste/release:2024.2
```

Available tags are listed at: https://hub.docker.com/r/chaste/release/tags

### Using the Development Branch (not recommended)

To use the latest development code from the `develop` branch:

```bash
docker run --init -it --rm -v chaste_data:/home/chaste chaste/develop
```

> ⚠️ The development branch may contain unstable or untested changes. Use a tagged release for any serious simulation work.

---

## Container Directory Structure

Once inside the container, you will be in `/home/chaste` with the following layout:

```
/home/chaste/
├── build/      # Precompiled Chaste binaries and libraries
├── projects/   # Symlink to /home/chaste/src/projects - your custom    projects go here
├── src/        # Chaste source code
└── output/     # Output folder used by the testing framework
```

---

## Working with Projects

### Create a New Project

```bash
new_project.sh
```

This creates a project template in `~/projects` with the correct directory structure and CMake configuration.

### Build a Project

```bash
build_project.sh <TestMyProject> c
```

Replace `<TestMyProject>` with your project name. The `c` argument tells the build system to reconfigure CMake - only needed when new source files have been added. Build output and simulation results will appear in `~/output`.

---

## Using VS Code (Alternative Method)

If you prefer working in VS Code rather than a terminal:

1. Clone the [Chaste repository](https://github.com/Chaste/Chaste).
2. Open the cloned folder in VS Code.
3. Install the [Remote Development extension pack](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack) when prompted.
4. Click **"Reopen in Container"** when prompted. VS Code will automatically pull the Docker image and configure the environment.

This gives you a full IDE experience while running code inside the container.

---

## Test the Project

Clone or copy this project into `~/projects` inside the container, then build it:

```bash
cd ~/build
cmake .
make TestPC9spheroid2D
ctest -V -R TestPC9spheroid2D
```

### Output Files

After a successful run, two folders will appear in `~/output`:

- `cellData/` - contains `.dat` files with per-cell simulation data
- `results_from_time_0/` - contains the mesh output files:
  - `mesh_results_0.vtu` - cell mesh geometry (visualise this)
  - `voronoi_results_0.vtu` - Voronoi tessellation of the cell population

These `.vtu` files can be visualised in [ParaView](https://www.paraview.org/).

### Modifying the Simulation

For your first run, no changes are needed. Once you have a baseline result, the most impactful parameter to adjust is the **cell cycle duration**, which controls how fast cells proliferate and directly affects spheroid growth rate. This can be found and modified in `TestPC9spheroid2D.hpp`.

---


## Further Resources

- Chaste Docker repository: https://github.com/Chaste/chaste-docker
- Chaste user tutorials: https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials
- Chaste documentation: https://chaste.github.io/docs/
- Docker documentation: https://docs.docker.com/
