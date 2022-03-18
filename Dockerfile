FROM python:3.9-slim

# Install MUMmer
RUN apt-get update && \
    apt-get install -y mummer

# Install GBKviz & Clear dependencies cache
RUN pip install -U pip && \
    pip install git+git://github.com/moshi4/GBKviz.git --no-cache-dir

# Launch GBKviz on port 8501
CMD ["gbkviz_webapp"]
EXPOSE 8501
