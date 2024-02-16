# get shiny serves plus tidyverse packages image
FROM rocker/shiny:3.6.3

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libxml2-dev \
    libjpeg-dev
    
# install R packages required 
# (change it dependeing on the packages you need)
RUN R -e "install.packages('shiny', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(pkgs = 'DESeq2', update = T, ask = F)"
RUN R -e "install.packages('shinyjs', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('DT', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggplot2', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggrepel', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('viridis', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('plotly', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('stringr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('htmlwidgets', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinyWidgets', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinycssloaders', repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('shinydashboard', repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('lubridate', repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('magrittr', repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('glue', repos='http://cran.rstudio.com/')"

# copy the app to the image
COPY app.R /srv/shiny-server/
COPY R /srv/shiny-server/R
COPY data /srv/shiny-server/data
COPY js /srv/shiny-server/js
COPY www /srv/shiny-server/www

# select port
EXPOSE 3838

# allow permission
RUN sudo chown -R shiny:shiny /srv/shiny-server

# run app
CMD ["/usr/bin/shiny-server.sh"]

