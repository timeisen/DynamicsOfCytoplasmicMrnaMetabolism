## Compiles the global model param files.
## TJE 2019 07 28.
args <- commandArgs(trailingOnly = TRUE)
library(tidyverse)
m <- read_tsv(args[1],col_names = FALSE)

mCat <- group_by(m, X1) %>%
	 mutate(n = length(X1),
			bootstrap_no = median(X1),
			si = median(X3),
			loc	= median(X4),
			sc	= median(X5),
			resid = median(X6)) %>%
	ungroup() %>%
	select(n, bootstrap_no, si, loc, sc, resid) %>%
	unique() %>%
	print() %>% 
	select(-n)

write_tsv(mCat, args[2])
