packages <- c(
  "cran::wrMisc",
  "drieslab/Giotto@v3.3.2",
  "xuranw/MuSiC@f21fe67f5670d5e9fca0ad7550abaae3423eb59c",
  "YingMa0107/CARD@2d64b91abb5cdd0c7f576b1c5d4727c84e7c93a0",
  "saezlab/DOT",
  "federicomarini/quantiseqr@1bc02b04350d190ed29de3f9c179f92ad70337af",
  "saezlab/OmnipathR@v4.0.4",
  "omnideconv/spacedeconv@main"
)

for (package in packages) {
  pak::pkg_install(package, dependencies = FALSE, upgrade = FALSE)
}
