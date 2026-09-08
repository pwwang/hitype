SHELL := /bin/bash

readme:
	Rscript <(echo "devtools::build_readme()")

docs:
	Rscript <(echo "devtools::document(); pkgdown::build_site()")

test:
	Rscript <(echo "devtools::test()")

.PHONY: readme docs test
