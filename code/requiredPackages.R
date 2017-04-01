requiredPackages <- c("jsonlite", "pracma", "nloptr")
for (name in requiredPackages) {
  if (!require(name, character.only = TRUE)) {
    install.packages(name, repos="https://cran.rstudio.com")
  }
}

## remove temp vars
rm(name, requiredPackages)
