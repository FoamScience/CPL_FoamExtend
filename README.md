
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


Installation
============

Please skim through the Docker files for installation instructions on Ubuntu 20.04 LTS (you can ignore the SSH-related commands in there).

License
=======

CPL_APP_OpenFOAM is released under the GNU GPL v3 license. Details are found in
the file LICENSE that is included with the release.
