FROM ubuntu:20.04

ARG USER_ID
ARG GROUP_ID
ARG DEBIAN_FRONTEND=noninteractive

RUN apk update && apk add \
  make execline doxygen graphviz bash g++ git

RUN addgroup --gid $GROUP_ID user
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user

USER root
RUN mkdir /math
RUN chown $USER_ID:$GROUP_ID /math
USER user
