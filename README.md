Foam-Extend CPL socket
======================

This repository assembles a Docker image which has what's needed to run coupled
[Foam-Extend 5.0](https://git.code.sf.net/p/foam-extend/foam-extend-4.1) (`nextRelease` branch) - 
[LAMMPS](https://github.com/lammps/lammps.git)
(Back from [this commit](https://github.com/lammps/lammps/commit/6354777d098deafc18a600877d00dbfcd8ce15c3) in the stable branch)
through the [CPL-library](https://github.com/Crompulence/cpl-library).

Everything is compiled on Ubuntu 20.04 LTS with Gcc-9 and MPICH2.
The image is basically used to develop the Foam-Extend socket for CPL, so you'll have to compile it separately.

Important notes about the Docker image
======================================

- The image is MPI-ready and you can use it on a Docker Swarm (or `docker-machine`) to spawn a cluster of containers.
- The default user is named `openfoam` with UID 1001 (Take a look at the `Dockerfile` for more info)
- `/cpl-library` is where the source code of the CPL library
- `~/foam` holds the patched Foam-Extend installation
- `~/lammps` holds the patched LAMMPS installation


How to use the Docker image
===========================

1. Install Docker if you don't have it yet: `bash <(curl -s https://get.docker.com/)`
2. Pull the image `docker pull foamscience/fe4-mpich-cpl-lammps`
3. Create a temporary container to play around with the library `docker run --rm -it foamscience/fe4-mpich-cpl-lammps bash`
4. You can then clone this repo there and compile with `make` after sourcing `SOURCEME.sh` file

`examples/LAMMPS-OPENFOAM` from this repo is an attempt to adapt `/cpl-library/examples/LAMMPS-OPENFOAM` (the original example)
to work with the example solver (`icoFoam`).

Installation of the Foam-Extend socket
======================================

Please skim through the Docker files for installation instructions on Ubuntu 20.04 LTS
(you can ignore the SSH-related commands in there).

I don't want to run things in Docker containers
===============================================

Unfortunately, we need to override some code in the `Pstreams` part of the `foam` library. If you do this
on your regular installation, you'll probably break it.

You can always dedicate a Foam-Extend installation to this task, but then you'll have to apply the patches present in `DockerFiles`
directory.

> If your Foam-Extend is not compiled with MPICH, you can get away with pre-loading the patched `libfoam.so` shared
> library (you can get it from the Docker image) to prevent the original one from loading with 
> `LD_PRELOAD=/path/to/patched/libfoam.so <your-solver>`
>
> However that is highly discouraged, please use the docker container instead!

License
=======

CPL_APP_OpenFOAM is released under the GNU GPL v3 license. Details are found in
the file LICENSE that is included with the release.
