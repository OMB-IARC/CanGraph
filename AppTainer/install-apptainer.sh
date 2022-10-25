#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2022 Pablo Marcos <software@loreak.org>
#
# SPDX-License-Identifier: MIT

# *********************** #
# Install APPTainer
# *********************** #

# This is the a simplified script that installs AppTainer v1.1.2
# on any debian-based OS (basically, one using BaSH and APT), using
# and installing Go v1.19 and the GoLangCI v1.43.0

# This script was taken from:
# https://github.com/apptainer/apptainer/blob/main/INSTALL.md

# Another option is to use a custom repo; see:
# https://docs.sylabs.io/guides/3.0/user-guide/installation.html#install-the-debian-ubuntu-package-using-apt

# *********************** #

# First, we update the repos and install some dependencies
sudo apt-get update && \
sudo apt-get install -y build-essential libseccomp-dev pkg-config \
    uidmap squashfs-tools squashfuse fuse2fs fuse-overlayfs fakeroot \
    cryptsetup curl wget git

# Then, we void the Go Folder and re-install the desired Go version
sudo rm -r /usr/local/go
export VERSION=1.19 OS=linux ARCH=amd64
wget -O /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz https://dl.google.com/go/go${VERSION}.${OS}-${ARCH}.tar.gz && \
sudo tar -C /usr/local -xzf /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz

# We appropriate the env vars accordingly using source and bashrc
echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
source ~/.bashrc

# We will also need to install the GoLangCI
curl -sfL https://install.goreleaser.com/github.com/golangci/golangci-lint.sh |
sh -s -- -b $(go env GOPATH)/bin v1.43.0

# Then, we clone the git repo for APPTainer and checkout to the v3.6.3 branch
mkdir -p ${GOPATH}/src/github.com/apptainer && \
cd ${GOPATH}/src/github.com/apptainer && \
rm -rf apptainer/
git clone https://github.com/apptainer/apptainer.git && \
cd apptainer
git checkout v1.1.2

# Finally, we install apptainer using the included files
./mconfig && \
cd ./builddir && \
make && \
sudo make install

# And, we can check it
apptainer --version
