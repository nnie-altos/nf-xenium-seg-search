# Custom container for nf-xenium-seg-search Python scripts.
# Provides: numpy, pandas, scipy, scanpy, tifffile, geopandas, matplotlib,
#           plotly, pyyaml, pyarrow, spatialdata for all bin/ scripts.
FROM python:3.11-slim

LABEL maintainer="Altos Labs"
LABEL org.opencontainers.image.source="https://github.com/altos-labs/nf-xenium-seg-search"

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir uv

RUN uv pip install --system --no-cache \
    numpy==2.0.2 \
    pandas==2.2.2 \
    scipy==1.14.1 \
    scanpy==1.10.1 \
    tifffile==2025.3.30 \
    "geopandas==0.14.4" \
    shapely==2.0.4 \
    matplotlib==3.8.4 \
    plotly==5.22.0 \
    kaleido==0.2.1 \
    pyyaml==6.0.1 \
    pyarrow==17.0.0 \
    scikit-image==0.24.0 \
    h5py==3.11.0 \
    anndata==0.10.7 \
    statsmodels==0.14.4 \
    patsy==0.5.6 \
    spatialdata==0.7.2 \
    dask==2025.5.1 \
    xarray==2025.3.1 \
    zarr==3.0.10

# Pre-warm numba cache in a world-writable directory so non-root container
# users (Nextflow passes -u $(id -u):$(id -g)) can read/write the cache.
ENV NUMBA_CACHE_DIR=/opt/numba_cache
RUN mkdir -p /opt/numba_cache && \
    python -c "import scanpy" && \
    chmod -R 777 /opt/numba_cache

# Copy bin scripts so they're available on PATH inside the container
COPY bin/ /usr/local/bin/
RUN chmod +x /usr/local/bin/*.py
