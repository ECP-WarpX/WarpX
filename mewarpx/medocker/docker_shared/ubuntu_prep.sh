#!/bin/bash

set -e
set -x

# Create a new user
useradd --create-home me_user
adduser me_user sudo

mkdir /merunset

chown me_user:me_user /merunset
chmod go+rX /root
