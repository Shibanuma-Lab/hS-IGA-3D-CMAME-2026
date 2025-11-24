# 🚀 Quick Reference Card

## New Computer Setup (3 Commands)

```bash
git clone <repository-url>
cd S-IGA-circular-crack-in-3D-solid-linear
./setup.sh
```

---

## Daily Usage

### Start Working
```bash
cd S-IGA-circular-crack-in-3D-solid-linear/circular_crack
pipenv shell
python3 main.py [options]
```

### Quick Run (Without Shell)
```bash
cd S-IGA-circular-crack-in-3D-solid-linear/circular_crack
pipenv run python3 main.py [options]
```

---

## What setup.sh Does

| Step | Action | Check |
|------|--------|-------|
| 1️⃣ | Verify directories (sfem_linear, circular_crack) | ✅ |
| 2️⃣ | Install Python 3.10 (if needed) | `python3.10 --version` |
| 3️⃣ | Install Pipenv (if needed) | `pipenv --version` |
| 4️⃣ | Create virtual environment | `pipenv --python 3.10` |
| 5️⃣ | Install numpy, scipy, logzero | `pipenv install` |
| 6️⃣ | Create logs directory | `mkdir -p circular_crack/logs` |

---

## Key Files

- **setup.sh** - One-command automated installer
- **Pipfile** - Python 3.10 dependencies (numpy, scipy, logzero)
- **INSTALL.md** - Detailed installation guide
- **README.md** - Project documentation

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| Pipenv not found | `source ~/.bashrc` |
| ModuleNotFoundError | Use `pipenv shell` or `pipenv run` |
| Python 3.10 not found | Run setup.sh again |
| Permission denied | `chmod +x setup.sh` |

---

## Important Notes

⚠️ **Python Version**: Must use Python 3.10 (specified in Pipfile)

⚠️ **Logs Directory**: Required for program execution (auto-created by setup.sh)

⚠️ **Virtual Environment**: Always use `pipenv shell` or `pipenv run` to execute

✅ **One-Time Setup**: Run setup.sh only once per computer

---

## Commands Cheat Sheet

```bash
# Check installation
python3.10 --version
pipenv --version
pipenv run python3 -c "import numpy, scipy, logzero"

# Enter/exit environment
pipenv shell     # Enter
exit            # Exit

# Run program
pipenv run python3 main.py --help

# Update dependencies (if Pipfile changes)
pipenv update

# Remove environment (if you want to start fresh)
pipenv --rm
```
