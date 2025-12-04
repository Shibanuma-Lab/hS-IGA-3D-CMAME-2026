# 🚀 Installation Guide

## For New Computers

This guide provides a complete installation workflow for setting up the S-IGA Circular Crack solver on a new computer.

---

## 📋 Prerequisites

Before starting, ensure you have:
- Ubuntu/Debian Linux or WSL (Windows Subsystem for Linux)
- Git installed
- Sudo privileges
- Internet connection

---

## 🔄 Complete Installation Workflow

### Step 1: Clone the Repository

```bash
git clone <repository-url>
cd S-IGA-circular-crack-in-3D-solid-linear
```

Verify directory structure:
```bash
ls -d sfem_linear circular_crack
```

You should see both directories exist.

---

### Step 2: Run Automated Setup

**Make the script executable and run:**

```bash
chmod +x setup.sh
./setup.sh
```

The script will automatically handle:

#### 2.1 Python 3.10 Installation

The script checks if Python 3.10 is installed. If not, it will:
- Install build dependencies
- Download Python 3.10.6 source
- Compile and install Python 3.10
- Verify installation

**Manual verification:**
```bash
python3.10 --version
# Should output: Python 3.10.6 (or similar)
```

#### 2.2 Pipenv Installation

The script checks if Pipenv is available. If not, it will:
- Install Pipenv via pip
- Add `$HOME/.local/bin` to PATH in `~/.bashrc`
- Update current session's PATH

**Manual verification:**
```bash
pipenv --version
# Should output: pipenv, version XX.XX.XX
```

**Note:** If Pipenv is not found after installation, run:
```bash
source ~/.bashrc
```

#### 2.3 Virtual Environment Setup

The script will:
- Create a Pipenv virtual environment using Python 3.10
- Install dependencies from `Pipfile`:
  - numpy >= 1.20.0
  - scipy >= 1.7.0
  - logzero >= 1.7.0

**Manual verification:**
```bash
pipenv run python3 -c "import numpy, scipy, logzero; print('✅ All packages installed')"
```

#### 2.4 Directory Structure Check

The script ensures required directories exist:
- Creates `circular_crack/logs/` if missing
- Verifies `sfem_linear/bin/sfem_linear` solver binary

---

## ✅ Verify Installation

After `setup.sh` completes, verify everything is working:

### Test 1: Check Python Environment

```bash
pipenv shell
python3 --version  # Should show Python 3.10.x
python3 -c "import numpy; print(f'NumPy: {numpy.__version__}')"
python3 -c "import scipy; print(f'SciPy: {scipy.__version__}')"
python3 -c "import logzero; print(f'Logzero: {logzero.__version__}')"
exit  # Exit the shell
```

### Test 2: Run the Solver Help

```bash
cd circular_crack
pipenv run python3 main.py --help
```

You should see the command-line interface help message.

### Test 3: Check Logs Directory

```bash
ls -la circular_crack/logs/
```

The directory should exist (even if empty initially).

---

## 🎮 Usage After Installation

### Method 1: Using Pipenv Shell (Interactive)

```bash
cd S-IGA-circular-crack-in-3D-solid-linear
pipenv shell  # Enter virtual environment

cd circular_crack
python3 main.py --help
python3 main.py [options]  # Run your simulation

exit  # Exit virtual environment when done
```

### Method 2: Direct Run (Non-interactive)

```bash
cd S-IGA-circular-crack-in-3D-solid-linear/circular_crack
pipenv run python3 main.py [options]
```

This method doesn't require entering/exiting the shell.

---

## 🔧 Manual Installation (Alternative)

If automated setup fails, follow manual steps:

### 1. Install Python 3.10

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y build-essential zlib1g-dev libncurses5-dev \
                    libgdbm-dev libnss3-dev libssl-dev libsqlite3-dev \
                    libreadline-dev libffi-dev libbz2-dev wget

cd /tmp
wget https://www.python.org/ftp/python/3.10.6/Python-3.10.6.tgz
tar -xf Python-3.10.6.tgz
cd Python-3.10.6
./configure --enable-optimizations
make -j$(nproc)
sudo make altinstall

python3.10 --version
```

### 2. Install Pipenv

```bash
python3.10 -m pip install --user pipenv

# Add to PATH
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

pipenv --version
```

### 3. Setup Environment

```bash
cd S-IGA-circular-crack-in-3D-solid-linear

# Create virtual environment
pipenv --python /usr/local/bin/python3.10

# Install dependencies
pipenv install

# Create logs directory
mkdir -p circular_crack/logs
```

---

## ⚠️ Troubleshooting

### Issue: Python 3.10 installation takes very long

**Solution:** The compilation step `make -j$(nproc)` uses all CPU cores. This is normal and can take 5-15 minutes depending on your system.

### Issue: Pipenv not found after installation

**Solution:** Update your PATH:
```bash
export PATH="$HOME/.local/bin:$PATH"
source ~/.bashrc
```

### Issue: Permission denied errors during Python installation

**Solution:** Ensure you're using `sudo` for system installation:
```bash
sudo make altinstall
```

### Issue: ModuleNotFoundError when running main.py

**Solution:** Ensure you're running inside the Pipenv environment:
```bash
pipenv shell
# or
pipenv run python3 main.py
```

### Issue: Pipfile requires Python 3.10 but you have a different version

**Solution:** Specify Python 3.10 explicitly:
```bash
pipenv --python /usr/local/bin/python3.10
```

### Issue: Logs directory error

**Solution:** Create it manually:
```bash
mkdir -p circular_crack/logs
```

---

## 📦 What Gets Installed

### Python Packages (via Pipenv)

| Package | Version | Purpose |
|---------|---------|---------|
| numpy | ≥ 1.20.0 | Numerical computing |
| scipy | ≥ 1.7.0 | Scientific computing (MATLAB file I/O) |
| logzero | ≥ 1.7.0 | Logging framework |

### System Packages (during Python build)

- build-essential (gcc, g++, make)
- Various development libraries (zlib, ssl, sqlite, etc.)

---

## 🔄 Updating the Installation

To update Python dependencies:

```bash
cd S-IGA-circular-crack-in-3D-solid-linear

# Update Pipfile if needed (edit manually)
# Then sync environment
pipenv update
```

To rebuild the Fortran solver (if needed):

```bash
cd sfem_linear
make clean
make
```

---

## 🗑️ Uninstallation

To remove the virtual environment:

```bash
cd S-IGA-circular-crack-in-3D-solid-linear
pipenv --rm
```

To remove Python 3.10 (if you want):

```bash
sudo rm /usr/local/bin/python3.10
sudo rm -rf /usr/local/lib/python3.10
```

---

## 📞 Support

If you encounter issues not covered here:

1. Check that directory structure is correct (`sfem_linear/` and `circular_crack/` exist)
2. Verify Python 3.10 is installed: `python3.10 --version`
3. Verify Pipenv is installed: `pipenv --version`
4. Check virtual environment: `pipenv run python3 -c "import sys; print(sys.version)"`
5. Review the setup script output for error messages

---

## 🎉 You're Ready!

Once installation is complete, you can:

- Run simulations with custom parameters
- Process results and calculate J-integrals
- Visualize crack propagation
- Extend the framework for your research

Happy computing! 🚀
