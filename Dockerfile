FROM rocker/tidyverse:latest

# Установка дополнительных пакетов, включая plotly, htmlwidgets, taxize, writexl и philentropy
RUN R -e "install.packages(c('plotly', 'htmlwidgets', 'taxize', 'writexl', 'philentropy'), repos='https://cloud.r-project.org')"

COPY compare_abundance.R convert_bracken.R convert_gtdbtk.R convert_kraken.R convert_metaphlan.R convert_truth.R /usr/local/bin/

RUN chmod +x /usr/local/bin/*.R

# Создаем удобные алиасы без .R
RUN for f in /usr/local/bin/*.R; do ln -s \"$f\" \"${f%.R}\"; done

ENTRYPOINT ["/bin/bash"]
