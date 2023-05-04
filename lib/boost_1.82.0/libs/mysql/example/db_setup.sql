--
-- Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
--
-- Distributed under the Boost Software License, Version 1.0. (See accompanying
-- file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
--

-- Connection system variables
SET NAMES utf8;

-- Database
DROP DATABASE IF EXISTS boost_mysql_examples;
CREATE DATABASE boost_mysql_examples;
USE boost_mysql_examples;

-- Tables
CREATE TABLE company(
    id CHAR(10) NOT NULL PRIMARY KEY,
    name VARCHAR(100) NOT NULL
);
CREATE TABLE employee(
    id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
    first_name VARCHAR(100) NOT NULL,
    last_name VARCHAR(100) NOT NULL,
    salary DOUBLE,
    company_id CHAR(10) NOT NULL,
    FOREIGN KEY (company_id) REFERENCES company(id)
);
    
INSERT INTO company (name, id) VALUES
    ("Award Winning Company, Inc.", "AWC"),
    ("Sector Global Leader Plc", "SGL"),
    ("High Growth Startup, Ltd", "HGS")
;
INSERT INTO employee (first_name, last_name, salary, company_id) VALUES
    ("Efficient", "Developer", 30000, "AWC"),
    ("Lazy", "Manager", 80000, "AWC"),
    ("Good", "Team Player", 35000, "HGS"),
    ("Enormous", "Slacker", 45000, "SGL"),
    ("Coffee", "Drinker", 30000, "HGS"),
    ("Underpaid", "Intern", 15000, "AWC")
;

-- User
DROP USER IF EXISTS 'example_user'@'%';
CREATE USER 'example_user'@'%' IDENTIFIED WITH 'mysql_native_password';
ALTER USER 'example_user'@'%' IDENTIFIED BY 'example_password';
GRANT ALL PRIVILEGES ON boost_mysql_examples.* TO 'example_user'@'%';
FLUSH PRIVILEGES;
