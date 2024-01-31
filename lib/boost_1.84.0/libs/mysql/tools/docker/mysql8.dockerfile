#
# Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

FROM mysql:8.1.0

ENV MYSQL_ALLOW_EMPTY_PASSWORD=1
ENV MYSQL_ROOT_PASSWORD=

COPY tools/docker/mysql_entrypoint.sh /
COPY tools/docker/ssl.cnf /etc/mysql/conf.d/
COPY tools/ssl/*.pem /etc/ssl/certs/mysql/

ENTRYPOINT ["/bin/bash", "/mysql_entrypoint.sh"]
