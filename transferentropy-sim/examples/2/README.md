# TE-Causality - Example #2

This example expands the first one by using shared volumes, and is intended for the Docker installation path.

The output is redirected to a `/output` directory in the container, which is mirrored to the `output` folder on the host.

## Setup

First let's create the output directory on the host machine:
```
mkdir /tmp/te-causality-output
```

Then we can run the reconstruction:
```
docker run -it \
  --mount "type=bind,source=/tmp/te-causality-output,target=/output" \
  te-causality bash -c "./te-extended examples/2/control.txt"
```
