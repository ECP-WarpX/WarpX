ARG BUILDPLATFORM

ARG WARPXCOMPUTE=OMP

FROM 167833485543.dkr.ecr.us-west-2.amazonaws.com/simteam/ubuntu-mirror:focal-${BUILDPLATFORM}

ENV WARPXCOMPUTE=${WARPXCOMPUTE}

ARG MPI_IMPL=ompi

# The '--login' is the way to get mpich to load happily, which in
# turn is necessary for a successful warp pinstall (and any mpi command
# actually) - NOT SURE HOW OPENMPI WORKS WITH THIS.

COPY WarpX/medocker/docker_shared/ubuntu_prep.sh /root
RUN bash --login /root/ubuntu_prep.sh

COPY WarpX/medocker/docker_shared/ubuntu_packages.sh /root
RUN bash --login /root/ubuntu_packages.sh

COPY WarpX/medocker/docker_shared/ubuntu_${MPI_IMPL}.sh /root
RUN bash --login /root/ubuntu_${MPI_IMPL}.sh

# Switch to the new user
USER me_user

COPY WarpX/medocker/docker_shared/python_packages.sh /root
# create directory so that it's added to $PATH before login
RUN mkdir -p /home/me_user/.local/bin
#RUN mkdir -p /root/.local/bin
RUN bash --login /root/python_packages.sh

# Copy in files
COPY --chown=me_user:me_user WarpX /merunset/WarpX
COPY --chown=me_user:me_user warpx-data /merunset/warpx-data
COPY --chown=me_user:me_user WarpX/medocker/docker_shared/run_warpx.sh /merunset/run/

WORKDIR /merunset

COPY WarpX/medocker/docker_shared/build_me_packages.sh /root
RUN /bin/bash --login /root/build_me_packages.sh

# This is needed to get around a bug in openmpi that would print copious error messages
# Unfortunately, this turns off CMA and uses shared memory for communication.
# An alternative is to do "docker run --cap-add SYS_PTRACE ...", which keeps CMA.
# ENV OMPI_MCA_btl_vader_single_copy_mechanism none

# DIRECTORY is used to select S3 directory
ENV DIRECTORY=null
# BUCKET is used to select S3 bucket
ENV BUCKET=null

WORKDIR /merunset/run
ENTRYPOINT ["/bin/bash", "--login", "-c", "source /etc/profile && ./run_warpx.sh $DIRECTORY $BUCKET"]
