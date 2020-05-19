
library("devtools")
library(roxygen2)
setwd("/Users/madjl/code")
#create("bedlm")

setwd("bedlm")
document()
check()

## Make distributable packages
setwd("~/code")
system("rm -f bedlm.tar.gz")
system("tar --exclude=.git -czvf bedlm.tar.gz bedlm")
install.packages("bedlm.tar.gz",repos = NULL, type="source")
