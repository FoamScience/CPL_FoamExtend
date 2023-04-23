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
- The default user is named `openfoam` with UID 1001 (Take a look at the `DockerFiles/01-fe4-mpich-cpl.Dockerfile` for more info)
- `/cpl-library` is where the source code of the CPL library lives
- `~/foam/foam-extend-5.0/` holds the patched Foam-Extend installation
- `~/lammps` holds the patched LAMMPS installation


How to use the Docker image
===========================

1. Install Docker if you don't have it yet: `bash <(curl -s https://get.docker.com/)`
2. Create a temporary container to play around with the library.
   - Basically: `docker run --rm -v /tmp/data:/home/openfoam/data -u $(id -u) -it foamscience/fe4-mpich-cpl-lammps:latest bash`
4. You can then clone this repo there and compile with `make` after sourcing the `SOURCEME.sh` file

`examples/LAMMPS-OPENFOAM` from this repo is an attempt to adapt `/cpl-library/examples/LAMMPS-OPENFOAM` (the original example)
to work with the example solver (`icoFoam`).

The docker image is designed to develop the FoamExtend-LAMMPS socket, typically you'd want to create a permanent `cpl`
container, which has access to your hosts's display:
```
docker run --name cpl --net host -e "DISPLAY" \
    -v /tmp/.X11-unix:/tmp/.X11-unix -v /tmp/data:/home/openfoam/data \
    -it foamscience/fe4-mpich-cpl-lammps:latest
```

On your host, run `xauth list`, you'll something like this:
```bash
MACHINE_NAME/unix:  ...  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
```

Get into the container, and run this command to give permission to the container so we can attach to the display from there:
```bash
docker exec -it cpl bash
# :1 here is the display; can be different for you, run this on the host: echo $DISPLAY 
# Replace MACHINE_NAME and the hash with values from the host
xauth add :1  MACHINE_NAME/unix  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
```

Now if you install ParaView and run it, you should see it pop up on your display.

> All changes to any other components (CPL library itself, FE4 source code or LAMMPS's code) should be applied with Git patches;
> It eases Docker image builds

Installation of the Foam-Extend socket
======================================

Please skim through the Docker files for installation instructions on Ubuntu 20.04 LTS
(you can ignore the SSH-related commands in there).

I don't want to run things in Docker containers
===============================================

Unfortunately, we need to override some code in the `Pstreams` part of the `foam` library. If you do this
on your regular installation, you'll probably break it.

You can always dedicate a Foam-Extend installation to CPL-related tasks. You'll just have to apply the patches present in `DockerFiles`
directory.

License
=======

CPL_APP_OpenFOAM is released under the GNU GPL v3 license. Details are found in
the file LICENSE that is included with the release.
