# Custom container for nf-xenium-seg-search Python scripts.
# Provides: numpy, pandas, scipy, scanpy, tifffile, geopandas, matplotlib,
#           plotly, pyyaml, pyarrow for all bin/ scripts.
FROM python:3.11-slim

LABEL maintainer="Altos Labs"
LABEL org.opencontainers.image.source="https://github.com/altos-labs/nf-xenium-seg-search"

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    numpy==1.26.4 \
    pandas==2.2.2 \
    scipy==1.13.1 \
    scanpy==1.10.1 \
    tifffile==2024.2.12 \
    "geopandas==0.14.4" \
    shapely==2.0.4 \
    matplotlib==3.8.4 \
    plotly==5.22.0 \
    kaleido==0.2.1 \
    pyyaml==6.0.1 \
    pyarrow==15.0.2 \
    scikit-image==0.23.2 \
    h5py==3.11.0 \
    anndata==0.10.7

# Copy bin scripts so they're available on PATH inside the container
COPY bin/ /usr/local/bin/
RUN chmod +x /usr/local/bin/*.py
