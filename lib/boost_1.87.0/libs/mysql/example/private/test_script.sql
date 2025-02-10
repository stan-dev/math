--
-- Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
--
-- Distributed under the Boost Software License, Version 1.0. (See accompanying
-- file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
--

USE boost_mysql_examples;

CREATE TEMPORARY TABLE products (
    id VARCHAR(50) PRIMARY KEY,
    description VARCHAR(256)
);

INSERT INTO products VALUES ('PTT', 'Potatoes'), ('CAR', NULL);
SELECT * FROM products;
DROP TABLE products;
