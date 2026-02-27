#!/bin/bash
set -e

echo "=== DLNM Mission Init ==="

# Install R packages if missing
echo "Checking R packages..."
Rscript -e '
pkgs <- c("rextendr", "mixmeta", "testthat", "pkgload", "bench", "data.table", "arrow")
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  cat("Installing:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, repos = "https://cloud.r-project.org", quiet = TRUE)
} else {
  cat("All R packages present\n")
}
'

# Install Python packages in venv (PEP 668 requires venv on externally-managed Python)
echo "Checking Python venv..."
if [ ! -d ".venv" ]; then
  echo "Creating Python venv..."
  python3 -m venv .venv
fi
.venv/bin/python3 -c "import matplotlib, numpy, pandas" 2>/dev/null || {
  echo "Installing Python packages in venv..."
  .venv/bin/pip install matplotlib numpy pandas --quiet
}

# Verify Rust toolchain
echo "Checking Rust toolchain..."
rustc --version
cargo --version

# Verify benchmark data exists
echo "Checking benchmark data..."
for f in benchmarks/data/scale_10mb.rds benchmarks/data/scale_100mb.rds benchmarks/data/scale_1gb.rds; do
  if [ -f "$f" ]; then
    echo "  OK: $f"
  else
    echo "  MISSING: $f"
  fi
done

if [ -d "benchmarks/data/scale_10gb" ]; then
  echo "  OK: benchmarks/data/scale_10gb/ ($(ls benchmarks/data/scale_10gb/*.parquet 2>/dev/null | wc -l) files)"
else
  echo "  MISSING: benchmarks/data/scale_10gb/"
fi

# Load dlnm package to verify it works
echo "Verifying dlnm loads from source..."
Rscript -e 'pkgload::load_all(".", quiet=TRUE); cat("dlnm loaded OK\n")'

echo "=== Init complete ==="
