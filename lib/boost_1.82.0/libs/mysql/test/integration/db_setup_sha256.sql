--
-- Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
--
-- Distributed under the Boost Software License, Version 1.0. (See accompanying
-- file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
--

-- Setup that requires the presence of SHA256 functionality
DROP USER IF EXISTS 'csha2p_user'@'%';
CREATE USER 'csha2p_user'@'%' IDENTIFIED WITH 'caching_sha2_password';
ALTER USER 'csha2p_user'@'%' IDENTIFIED BY 'csha2p_password';
GRANT ALL PRIVILEGES ON boost_mysql_integtests.* TO 'csha2p_user'@'%';

DROP USER IF EXISTS 'csha2p_empty_password_user'@'%';
CREATE USER 'csha2p_empty_password_user'@'%' IDENTIFIED WITH 'caching_sha2_password';
ALTER USER 'csha2p_empty_password_user'@'%' IDENTIFIED BY '';
GRANT ALL PRIVILEGES ON boost_mysql_integtests.* TO 'csha2p_empty_password_user'@'%';

-- This one is used for unknown authentication plugin tests
DROP USER IF EXISTS 'sha2p_user'@'%';
CREATE USER 'sha2p_user'@'%' IDENTIFIED WITH 'sha256_password';
ALTER USER 'sha2p_user'@'%' IDENTIFIED BY 'sha2p_password';
GRANT ALL PRIVILEGES ON boost_mysql_integtests.* TO 'sha2p_user'@'%';

FLUSH PRIVILEGES;