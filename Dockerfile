FROM dolfinx/dolfinx:stable@sha256:2ae4bfbc0d9077268880faf04c72750528bee986c94ab223a2c159969bd56fa8

WORKDIR /workspace

RUN python -m pip install --no-cache-dir \
    gmsh==4.15.2 \
    matplotlib==3.10.9

COPY . /workspace

ENV PYTHONPATH=/workspace/src:${PYTHONPATH}
CMD ["bash"]
