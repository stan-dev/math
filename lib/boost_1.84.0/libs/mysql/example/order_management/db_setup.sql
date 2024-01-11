--
-- Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
--
-- Distributed under the Boost Software License, Version 1.0. (See accompanying
-- file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
--

-- Connection system variables
SET NAMES utf8;

-- Database
DROP DATABASE IF EXISTS boost_mysql_order_management;
CREATE DATABASE boost_mysql_order_management;
USE boost_mysql_order_management;

-- User
DROP USER IF EXISTS 'orders_user'@'%';
CREATE USER 'orders_user'@'%' IDENTIFIED WITH 'mysql_native_password';
ALTER USER 'orders_user'@'%' IDENTIFIED BY 'orders_password';
GRANT ALL PRIVILEGES ON boost_mysql_order_management.* TO 'orders_user'@'%';
FLUSH PRIVILEGES;

-- Table definitions
CREATE TABLE products (
    id INT PRIMARY KEY AUTO_INCREMENT,
    short_name VARCHAR(100) NOT NULL,
    descr TEXT,
    price INT NOT NULL,
    FULLTEXT(short_name, descr)
);

CREATE TABLE orders(
    id INT PRIMARY KEY AUTO_INCREMENT,
    `status` ENUM('draft', 'pending_payment', 'complete') NOT NULL DEFAULT 'draft'
);

CREATE TABLE order_items(
    id INT PRIMARY KEY AUTO_INCREMENT,
    order_id INT NOT NULL,
    product_id INT NOT NULL,
    quantity INT NOT NULL,
    FOREIGN KEY (order_id) REFERENCES orders(id),
    FOREIGN KEY (product_id) REFERENCES products(id)
);

-- Procedures
DELIMITER //

CREATE DEFINER = 'orders_user'@'%' PROCEDURE get_products(IN p_search VARCHAR(50))
BEGIN
    DECLARE max_products INT DEFAULT 20;
    IF p_search IS NULL THEN
        SELECT id, short_name, descr, price
        FROM products
        LIMIT max_products;
    ELSE
        SELECT id, short_name, descr, price FROM products
        WHERE MATCH(short_name, descr) AGAINST(p_search)
        LIMIT max_products;
    END IF;
END //

CREATE PROCEDURE create_order()
BEGIN
    START TRANSACTION;

    -- Create the order
    INSERT INTO orders () VALUES ();

    -- Return the order
    SELECT id, `status`
    FROM orders
    WHERE id = LAST_INSERT_ID();

    COMMIT;
END //

CREATE DEFINER = 'orders_user'@'%' PROCEDURE get_order(
    IN p_order_id INT
)
BEGIN
    DECLARE order_status TEXT;
    START TRANSACTION READ ONLY;

    -- Check parameters
    IF p_order_id IS NULL THEN
        SIGNAL SQLSTATE '45000' SET MYSQL_ERRNO = 1048, MESSAGE_TEXT = 'get_order: invalid parameters';
    END IF;

    -- Check that the order exists
    SELECT `status`
    INTO order_status
    FROM orders WHERE id = p_order_id;
    IF order_status IS NULL THEN
        SIGNAL SQLSTATE '45000' SET MYSQL_ERRNO = 1329, MESSAGE_TEXT = 'The given order does not exist';
    END IF;

    -- Return the order. The IFNULL statements make MySQL correctly report the fields as non-NULL
    SELECT
        IFNULL(p_order_id, 0) AS id,
        IFNULL(order_status, 'draft') AS `status`;
    SELECT
        item.id AS id,
        item.quantity AS quantity,
        prod.price AS unit_price
    FROM order_items item
    JOIN products prod ON item.product_id = prod.id
    WHERE item.order_id = p_order_id;

    COMMIT;
END //

CREATE DEFINER = 'orders_user'@'%' PROCEDURE get_orders()
BEGIN
    SELECT id, `status` FROM orders;
END //

CREATE DEFINER = 'orders_user'@'%' PROCEDURE add_line_item(
    IN p_order_id INT,
    IN p_product_id INT,
    IN p_quantity INT,
    OUT pout_line_item_id INT
)
BEGIN
    DECLARE product_price INT;
    DECLARE order_status TEXT;
    START TRANSACTION;

    -- Check parameters
    IF p_order_id IS NULL OR p_product_id IS NULL OR p_quantity IS NULL OR p_quantity <= 0 THEN
        SIGNAL SQLSTATE '45000' SET MYSQL_ERRNO = 1048, MESSAGE_TEXT = 'add_line_item: invalid params';
    END IF;

    -- Ensure that the product is valid
    SELECT price INTO product_price FROM products WHERE id = p_product_id;
    IF product_price IS NULL THEN
        SIGNAL SQLSTATE '45000' SET MYSQL_ERRNO = 1329, MESSAGE_TEXT = 'The given product does not exist';
    END IF;

    -- Get the order
    SELECT `status` INTO order_status FROM orders WHERE id = p_order_id;
    IF order_status IS NULL THEN
        SIGNAL SQLSTATE '45000' SET MYSQL_ERRNO = 1329, MESSAGE_TEXT = 'The given order does not exist';
    END IF;
    IF order_status <> 'draft' THEN
        SIGNAL SQLSTATE '45000' SET MYSQL_ERRNO = 1000, MESSAGE_TEXT = 'The given order is not editable';
    END IF;

    -- Insert the new item
    INSERT INTO order_items (order_id, product_id, quantity) VALUES (p_order_id, p_product_id, p_quantity);

    -- Return value
    SET pout_line_item_id = LAST_INSERT_ID();

    -- Return the edited order
    SELECT id, `status`
    FROM orders WHERE id = p_order_id;
    SELECT
        item.id AS id,
        item.quantity AS quantity,
        prod.price AS unit_price
    FROM order_items item
    JOIN products prod ON item.product_id = prod.id
    WHERE item.order_id = p_order_id;

    COMMIT;
END //

CREATE DEFINER = 'orders_user'@'%' PROCEDURE remove_line_item(
    IN p_line_item_id INT
)
BEGIN
    DECLARE order_id INT;
    DECLARE order_status TEXT;
    START TRANSACTION;

    -- Check parameters
    IF p_line_item_id IS NULL THEN
        SIGNAL SQLSTATE '45000' SET MYSQL_ERRNO = 1048, MESSAGE_TEXT = 'remove_line_item: invalid params';
    END IF;

    -- Get the order
    SELECT orders.id, orders.`status`
    INTO order_id, order_status
    FROM orders
    JOIN order_items items ON (orders.id = items.order_id)
    WHERE items.id = p_line_item_id;

    IF order_status IS NULL THEN
        SIGNAL SQLSTATE '45000' SET MYSQL_ERRNO = 1329, MESSAGE_TEXT = 'The given order item does not exist';
    END IF;
    IF order_status <> 'draft' THEN
        SIGNAL SQLSTATE '45000' SET MYSQL_ERRNO = 1000, MESSAGE_TEXT = 'The given order is not editable';
    END IF;

    -- Delete the line item
    DELETE FROM order_items
    WHERE id = p_line_item_id;

    -- Return the edited order
    SELECT id, `status`
    FROM orders WHERE id = order_id;
    SELECT
        item.id AS id,
        item.quantity AS quantity,
        prod.price AS unit_price
    FROM order_items item
    JOIN products prod ON item.product_id = prod.id
    WHERE item.order_id = order_id;

    COMMIT;
END //

CREATE DEFINER = 'orders_user'@'%' PROCEDURE checkout_order(
    IN p_order_id INT,
    OUT pout_order_total INT
)
BEGIN
    DECLARE order_status TEXT;
    START TRANSACTION;

    -- Check parameters
    IF p_order_id IS NULL THEN
        SIGNAL SQLSTATE '45000' SET MYSQL_ERRNO = 1048, MESSAGE_TEXT = 'checkout_order: invalid params';
    END IF;

    -- Get the order
    SELECT `status`
    INTO order_status
    FROM orders WHERE id = p_order_id;

    IF order_status IS NULL THEN
        SIGNAL SQLSTATE '45000' SET MYSQL_ERRNO = 1329, MESSAGE_TEXT = 'The given order does not exist';
    END IF;
    IF order_status <> 'draft' THEN
        SIGNAL SQLSTATE '45000' SET MYSQL_ERRNO = 1000, MESSAGE_TEXT = 'The given order is not in a state that can be checked out';
    END IF;

    -- Update the order
    UPDATE orders SET `status` = 'pending_payment' WHERE id = p_order_id;

    -- Retrieve the total price
    SELECT SUM(prod.price * item.quantity)
    INTO pout_order_total
    FROM order_items item
    JOIN products prod ON item.product_id = prod.id
    WHERE item.order_id = p_order_id;

    -- Return the edited order
    SELECT id, `status`
    FROM orders WHERE id = p_order_id;
    SELECT
        item.id AS id,
        item.quantity AS quantity,
        prod.price AS unit_price
    FROM order_items item
    JOIN products prod ON item.product_id = prod.id
    WHERE item.order_id = p_order_id;

    COMMIT;
END //


CREATE DEFINER = 'orders_user'@'%' PROCEDURE complete_order(
    IN p_order_id INT
)
BEGIN
    DECLARE order_status TEXT;
    START TRANSACTION;

    -- Check parameters
    IF p_order_id IS NULL THEN
        SIGNAL SQLSTATE '45000' SET MYSQL_ERRNO = 1048, MESSAGE_TEXT = 'complete_order: invalid params';
    END IF;

    -- Get the order
    SELECT `status`
    INTO order_status
    FROM orders WHERE id = p_order_id;

    IF order_status IS NULL THEN
        SIGNAL SQLSTATE '45000' SET MYSQL_ERRNO = 1329, MESSAGE_TEXT = 'The given order does not exist';
    END IF;
    IF order_status <> 'pending_payment' THEN
        SIGNAL SQLSTATE '45000' SET MYSQL_ERRNO = 1000, MESSAGE_TEXT = 'The given order is not in a state that can be completed';
    END IF;

    -- Update the order
    UPDATE orders SET `status` = 'complete' WHERE id = p_order_id;

    -- Return the edited order
    SELECT id, `status`
    FROM orders WHERE id = p_order_id;
    SELECT
        item.id AS id,
        item.quantity AS quantity,
        prod.price AS unit_price
    FROM order_items item
    JOIN products prod ON item.product_id = prod.id
    WHERE item.order_id = p_order_id;

    COMMIT;
END //

DELIMITER ;

-- Contents for the products table
INSERT INTO products (price, short_name, descr) VALUES
    (6400, 'A Feast for Odin', 'A Feast for Odin is a points-driven game, with plethora of pathways to victory, with a range of risk balanced against reward. A significant portion of this is your central hall, which has a whopping -86 points of squares and a major part of your game is attempting to cover these up with various tiles. Likewise, long halls and island colonies can also offer large rewards, but they will have penalties of their own.'),
    (1600, 'Railroad Ink',     'The critically acclaimed roll and write game where you draw routes on your board trying to connect the exits at its edges. The more you connect, the more points you make, but beware: each incomplete route will make you lose points!'),
    (4000, 'Catan',            'Catan is a board game for two to four players in which you compete to gather resources and build the biggest settlements on the fictional island of Catan. It takes approximately one hour to play.'),
    (2500, 'Not Alone',        'It is the 25th century. You are a member of an intergalactic expedition shipwrecked on a mysterious planet named Artemia. While waiting for the rescue ship, you begin to explore the planet but an alien entity picks up your scent and begins to hunt you. You are NOT ALONE! Will you survive the dangers of Artemia?'),
    (4500, 'Dice Hospital',    "In Dice Hospital, a worker placement board game, players are tasked with running a local hospital. Each round you'll be admitting new patients, hiring specialists, building new departments, and treating as many incoming patients as you can.")
;
