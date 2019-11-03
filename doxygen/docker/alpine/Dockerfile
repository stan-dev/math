FROM nnadeau/docker-doxygen

RUN apk update && apk add \
  make execline

RUN mkdir math
COPY . ./math
WORKDIR ./math
