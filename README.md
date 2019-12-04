NanoStringClustR will be an R package for evaluating NanoString nCounter normalization and
is under development! 

To obtain the latest version from github, install devtools in R and use the following:

```{r}
library(devtools)
devtools::install_github("MarthaCooper/NanoStringClustR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
```

To run NanoStringClustR, load the library and follow the examples in the vignette.

```{r}
library(NanoStringClustR)
browseVignettes("NanoStringClustR")
```

For more details, contact Martha Cooper: martha.cooper92@gmail.com 

Version 0.1.1

