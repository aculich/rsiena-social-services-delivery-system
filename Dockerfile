FROM aculich/rockyter

RUN conda install --yes r-igraph

USER root

RUN apt-get update && apt-get install --yes tcl

USER $NB_USER
RUN Rscript -e "install.packages('RSiena', dependencies = TRUE, repos='http://cran.us.r-project.org')"
