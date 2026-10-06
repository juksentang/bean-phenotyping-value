# ============================================================================
# Install R dependencies for Layer 1 + Layer 2
# ============================================================================

pkgs <- c("AlphaSimR", "yaml", "data.table")
installed <- installed.packages()[, "Package"]
to_install <- pkgs[!pkgs %in% installed]

if (length(to_install) > 0) {
  cat("Installing:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install, repos = "https://mirrors.163.com/cran/")
} else {
  cat("All packages already installed.\n")
}

# Verify
for (pkg in pkgs) {
  cat(pkg, ":", requireNamespace(pkg, quietly = TRUE), "\n")
}
