FROM rocker/tidyverse:latest

# Установка дополнительных пакетов (добавил vegan)
RUN R -e "install.packages(c('plotly', 'htmlwidgets', 'taxize', 'writexl', 'philentropy', 'vegan'), repos='https://cloud.r-project.org')"

# Копируем R-скрипты из папки scripts в /usr/local/bin
COPY scripts/ /usr/local/bin/

# Делаем все скрипты исполняемыми
RUN chmod +x /usr/local/bin/*.R

# Создаем алиасы без расширения .R
RUN for f in /usr/local/bin/*.R; do ln -s "$f" "${f%.R}"; done

# По умолчанию запускаем bash
ENTRYPOINT ["/bin/bash"]
