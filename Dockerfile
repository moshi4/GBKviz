FROM python:3.9-slim

ENV GBKVIZ_VERSION=1.1.0

# Install MUMmer
RUN apt-get update && \
    apt-get install -y mummer

# Install GBKviz & Clear dependencies cache
RUN pip install -U pip && \
    pip install gbkviz==${GBKVIZ_VERSION} --no-cache-dir

# Launch GBKviz on port 8501
CMD ["gbkviz_webapp"]
EXPOSE 8501
