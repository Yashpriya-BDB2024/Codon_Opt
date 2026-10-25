#!/bin/bash
# ============================================================
# DNA Chisel — Environment Setup
# License: MIT (commercial use allowed)
# Algorithm: Sliding-window, deterministic, multi-parameter
# Best for: CHO, mammalian — CAI + GC + restriction sites
#           + hairpins simultaneously
#
# For CentOS 7 with Python 3.12 at /usr/local/bin/python3.12
# No sudo required
# Usage: bash setup_dnachisel.sh
# ============================================================

set -e

PYTHON=/usr/local/bin/python3.12
TOOL_DIR="$HOME/codon_optimization/dnachisel"
ENV_DIR="$TOOL_DIR/venv"

echo "================================================"
echo "  DNA Chisel — Setup"
echo "  License: MIT | Commercial use: YES"
echo "================================================"
echo ""

# ── Check Python ─────────────────────────────────────────────
echo "[1/4] Checking Python 3.12..."
if [ ! -f "$PYTHON" ]; then
    echo "ERROR: $PYTHON not found. Run check_system.sh first."
    exit 1
fi
$PYTHON --version
echo "✓ Python 3.12 confirmed"
echo ""

# ── Create directories and virtual environment ───────────────
echo "[2/4] Creating environment at $TOOL_DIR ..."
mkdir -p "$TOOL_DIR"
if [ ! -d "$ENV_DIR" ]; then
    $PYTHON -m venv "$ENV_DIR"
fi
"$ENV_DIR/bin/pip" install --upgrade pip --quiet
echo "✓ Virtual environment ready"
echo ""

# ── Install DNA Chisel and dependencies ──────────────────────
echo "[3/4] Installing DNA Chisel..."
"$ENV_DIR/bin/pip" install "numpy<2.0" --quiet
"$ENV_DIR/bin/pip" install dnachisel --quiet
"$ENV_DIR/bin/pip" install biopython --quiet
"$ENV_DIR/bin/pip" install python-codon-tables --quiet
"$ENV_DIR/bin/pip" install pandas matplotlib seaborn --quiet
"$ENV_DIR/bin/pip" install jupyter notebook --quiet
echo "✓ All packages installed"
echo ""

# ── Verify ───────────────────────────────────────────────────
echo "[4/4] Verifying installation..."
"$ENV_DIR/bin/python3" -c "
import dnachisel, Bio, python_codon_tables, pandas, matplotlib
print(f'  ✓ dnachisel          {dnachisel.__version__}')
print(f'  ✓ biopython          {Bio.__version__}')
print(f'  ✓ python-codon-tables OK')
print(f'  ✓ pandas             {pandas.__version__}')
print(f'  ✓ matplotlib         {matplotlib.__version__}')
print()
print('  All packages verified successfully!')
"
echo ""

echo "================================================"
echo "  Setup Complete!"
echo "================================================"
echo ""
echo "  Activate:  source $ENV_DIR/bin/activate"
echo "  Run test:  python3 $TOOL_DIR/test_dnachisel.py"
echo "  Jupyter:   jupyter notebook"
echo ""
