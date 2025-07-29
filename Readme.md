# GRAPE-Tauola Simulation Setup

## Setup Instructions

### 1. Clone the Repository

```bash
git clone <repository-url>
cd <repository-directory>
```

**Note:** All commands below should be run from within the main directory of the repository.

---

### 2. Switch to C Shell

```bash
csh
```

---

## Step 1: Compile the Programs

### Build HEPMC3 and TAUOLA

```bash
chmod +x build_hep.sh
./build_hep.sh
```
- This step may take some time but only needs to be done once.

### Build GRAPE
First, go into grape-dilepton_adapted.1k/src/set_grape_spring and change the file paths so they point to the same files in *your* directory (I will automate this down the line)! Then, return to the highest directory and run:
```bash
chmod +x build.sh
./build.sh
```

- GRAPE executables will be located in:
  ```
  ./grape-dilepton/build
  ```


---

## Step 2: Run GRAPE

1. Edit the `grape.cards` file to configure run parameters.

2. Run GRAPE:

```bash
chmod +x run_grape.sh
./run_grape.sh
```

---

## Step 3: Run TAUOLA

```bash
chmod +x run_tauola.sh
./run_tauola.sh
```

---

## Output

- Output is printed to the terminal.
- Partially processed files will be found in:
  ```
  ./grape/build
  ```
