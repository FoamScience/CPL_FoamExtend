# Build this image:  docker build -t fe4-mpi .
#
# Builds Foam-extend 4.1 nextRelease branch with MPICH
# from Ubuntu 20.04 LTS (focal) in Opt mode

# Base image
FROM cpl-lib-local

# Main user
ENV USER openfoam

# Disable apt prompts and set user home
ENV DEBIAN_FRONTEND=noninteractive \
    HOME=/home/${USER} 

# Install requirements
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends sudo apt-utils vim && \
    apt-get install -y --no-install-recommends openssh-server python-dev \
        gfortran binutils && \
    apt-get install -y git-core build-essential binutils-dev cmake flex libfl-dev \
        zlib1g-dev libncurses5-dev curl bison \
        libxt-dev rpm mercurial graphviz gcc-7 g++-7 && \
    apt-get remove libopenmpi-dev openmpi-bin openmpi-common openmpi-doc && \
    apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# ------------------------------------------------------------
# SSH mess
# ------------------------------------------------------------

RUN mkdir /var/run/sshd
RUN echo 'root:${USER}' | chpasswd
RUN sed -i 's/PermitRootLogin without-password/PermitRootLogin yes/' /etc/ssh/sshd_config
# SSH login fix. Otherwise user is kicked off right after login
RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd
ENV NOTVISIBLE "in users profile"
RUN echo "export VISIBLE=now" >> /etc/profile

# ------------------------------------------------------------
# Add an 'openfoam' user with root access
# ------------------------------------------------------------

RUN adduser --disabled-password --gecos "" ${USER} -u 1001 && \
    echo "ALL ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers

# ------------------------------------------------------------
# Set-Up SSH with a dummy key
# ------------------------------------------------------------

ENV SSHDIR ${HOME}/.ssh/

RUN mkdir -p ${SSHDIR}

ADD ssh/config ${SSHDIR}/config
ADD ssh/id_rsa.fe4 ${SSHDIR}/id_rsa
ADD ssh/id_rsa.fe4.pub ${SSHDIR}/id_rsa.pub
ADD ssh/id_rsa.fe4.pub ${SSHDIR}/authorized_keys

RUN chmod -R 600 ${SSHDIR}* && \
    chown -R ${USER}:${USER} ${SSHDIR}

# ------------------------------------------------------------
# Get, Patch and Compile foam-extend-5.0
# ------------------------------------------------------------

USER openfoam

ENV FOAM_REPO_URL https://git.code.sf.net/p/foam-extend/foam-extend-4.1
ENV FOAM_BRANCH   nextRelease
ENV FOAM_VNAME    foam-extend-5.0

RUN mkdir -p ${HOME}/foam
WORKDIR ${HOME}/foam
RUN git clone --depth 1 --single-branch --branch ${FOAM_BRANCH} ${FOAM_REPO_URL} ${FOAM_VNAME}
RUN git config --global user.email "you@example.com" && \
    git config --global user.name "Your Name"

WORKDIR ${HOME}/foam/${FOAM_VNAME}
COPY 0001-compile-on-Ubuntu-20.04-with-MPICH.patch .
RUN git am 0001-compile-on-Ubuntu-20.04-with-MPICH.patch
SHELL ["/bin/bash", "-c"]
ENV MPI_ARCH_PATH /home/${USER}/foam/${FOAM_VNAME}/ThirdParty/mpich-1.2.4/platforms/linux64GccDPInt32Opt
RUN mkdir -p ${MPI_ARCH_PATH}
RUN ln -s /usr/lib/mpich ${MPI_ARCH_PATH}/include
RUN bash -c "source etc/bashrc; ./Allwmake.firstInstall"
RUN bash -c "source etc/bashrc; ./Allwmake"
RUN echo 'source ~/foam/${FOAM_VNAME}/etc/bashrc' >> ${HOME}/.bashrc

# ------------------------------------------------------------
# Get, Patch and Compile lammps
# ------------------------------------------------------------

ENV LAMMPS_COMMIT "6354777d098deafc18a600877d00dbfcd8ce15c3"
RUN git clone -b stable https://github.com/lammps/lammps.git ${HOME}/lammps && \
    git clone https://github.com/Crompulence/CPL_APP_LAMMPS-DEV.git ${HOME}/CPL_APP_LAMMPS-DEV

WORKDIR ${HOME}/lammps
RUN git checkout ${LAMMPS_COMMIT}
WORKDIR ${HOME}/CPL_APP_LAMMPS-DEV
RUN echo "${HOME}/lammps" > ${HOME}/CPL_APP_LAMMPS-DEV/CODE_INST_DIR && \
    echo granular >> config/lammps_packages.in && \
    cd config && \
    sh ./enable-packages.sh make && \
    cd ../ && \
    make patch-lammps

RUN make -j $(nproc)
ENV PATH="${HOME}/CPL_APP_LAMMPS-DEV/bin:${PATH}"

# ------------------------------------------------------------
# Patch Pstreams in Foam Extend
# ------------------------------------------------------------

WORKDIR ${HOME}/foam/${FOAM_VNAME}/src
COPY 0002-Fix-MPI_init-called-twice.patch .
RUN git am 0002-Fix-MPI_init-called-twice.patch
RUN bash -c "source ${HOME}/foam/${FOAM_VNAME}/etc/bashrc; ./Allwmake"

# ------------------------------------------------------------
# Configure MPI and set ownership on ~/data
# ------------------------------------------------------------

USER root

RUN rm -fr ${HOME}/.openmpi && mkdir -p ${HOME}/.openmpi
ADD default-mca-params.conf ${HOME}/.openmpi/mca-params.conf
RUN chown -R ${USER}:${USER} ${HOME}/.openmpi
RUN mkdir ${HOME}/data
RUN chown -R ${USER}:${USER} ${HOME}/data

# ------------------------------------------------------------
# Final preparations
# ------------------------------------------------------------

WORKDIR ${HOME}/data

USER root
EXPOSE 22
RUN /usr/bin/ssh-keygen -A 
USER ${USER}

LABEL maintainer="Mohammed Elwardi Fadeli <fadeli.elwardi@tu-darmstadt.de>"
# Make sure SSH has the keys
CMD bash -c 'sudo /usr/sbin/sshd -d'
