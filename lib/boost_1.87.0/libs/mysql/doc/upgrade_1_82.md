# Upgrade instructions to 1.82

This document describes how to upgrade from 0.2.x to 1.82 (please remember that 1.82 is the first stable release in this library, and the versioning follows Boost versioning).

This is a major upgrade, since a lot of changes have been made since the review took place. If you encounter any problems, please file an issue against the repository, and we will be happy to help :)

* The `value` class has been replaced by `field_view`. The generic accessors have been replaced by type-specific functions. Conversions, `get_optional` and `get_std_optional` have been removed, in favor of an interface more similar to `json::value`. The class `field` has been added to be able to take ownership of a `field_view`.
    * **Action**: replace `value` by `field_view`.
* Statements are no longer I/O objects. `connection::prepare_statement` works the same, but execution requires going through the connection. 
    * **Action**: replace `stmt.execute(...)` by `conn.execute_statement(stmt, ...)`,
      `stmt.close()` by `conn.close_statement(stmt)`.
* Resultsets as I/O objects have been removed, and executing queries and statements is now simpler. There are now two ways to run queries or statements:
    * `connection::query` and `connection::execute_statement` now read all rows into memory, into a `results` object (which is similar to `resultset`, but is a plain data object and contains all the rows). These functions are equivalent to the old `connection::query`/`statement::execute` plus `resultset::read_all`.
        * **Action**: if you were using `conn.query(...).read_all()` (or similar), replace it by `conn.query(...)`.
    * `connection::start_query` and `statement::start_execution` behave similarly to the old `connection::query` and `statement::execute`. They use a data structure called `execution_state`, which is similar to the old `resultset`.
        * **Action**: if you were reading rows using `read_one` or `read_many`, consider whether your application requires reading row-by-row or not. If it doesn't, apply the point above. If it does, replace `conn.query(...)` by `conn.start_query(...)`, and `result.read_one(...)` by `conn.read_some_rows(...)`. There is no function to read a single row.
* Statement parameters are now passed as a `std::tuple` to `execute_statement`, rather than a collection.
    * **Action**:
        * If you were using `stmt.execute(make_values(a, b), ...)`, replace it by `conn.execute_statement(std::make_tuple(a, b), ...)`
        * If you were using a collection and can't migrate to a `std::tuple` (because the number of parameters is unknown at compile-time), then have a look at [this solution](https://github.com/boostorg/mysql/issues/110).
* `error_info` has been renamed to `diagnostics` and `error_info::message()` to `diagnostics::server_message()`. The message is no longer included by default in the thrown exceptions's `what()`. This is because the message may not be UTF-8 compatible.
    * **Action**: apply the renames. If you're using exceptions, catch the new `error_with_diagnostics` exception type to retrieve the server message.
* The `date` and `datetime` types are now custom types, instead of aliases for `std::chrono`. This enables them to represent zero and invalid dates.
    * **Action**: to get a `time_point` from a `date` or `datetime`, use `as_time_point` or `get_time_point` + `valid`.
* Binary types (`BLOB`, `BINARY`, `VARBINARY`, `GEOMETRY`) are now represented as a special type `blob_view`/`blob`.
    * **Action**: if you handle these types in your application, use `is_blob`, `get_blob` and `as_blob` functions in `field_view`.
* The `collation` enum has been removed in favor of plain integers.
    * **Action**: if you were using collations explicity, replace the enumerator by the collation ID. You can find them in `<boost/mysql/mysql_collations.hpp>` and `<boost/mysql/mariadb_collations.hpp>`.
* `connection_params` has been reverted to `handshake_params`.
    * **Action**: replace occurrences of `connection_params` by `handshake_params`.
* `row` is no longer streamable. The stream operation on `row` wasn't a universal agreement.
    * **Action**: if you were streaming rows, switch to using a loop and streaming individual `field_view`s, instead.
* Metadata strings is no longer retained by default, to save allocations.
    * **Action**: if your code uses metadata strings (e.g. `metadata::column_name()`, use `conn.set_metadata_mode(metadata_mode::full))` before initiating any query/statement execution.
* The library now uses `boost::core::string_view` (under the alias `boost::mysql::string_view`) rather than `boost::string_view`.
    * **Action**: if you were using `boost::string_view` directly, change it to `boost::mysql::string_view`.
* `errc` has been split into `client_errc`, `server_errc` and server-specific error codes (which don't have an enum).
    * **Action**: if you were referencing `errc` values directly, replace `errc` by `client_errc` or `server_errc`.
* `field_type` has been renamed to `column_type`.
    * **Action**: apply the rename.
* The following members in `metadata` have been renamed:
    * `field_name` => `column_name`
    * `original_field_name` => `original_column_name`
    * `character_set` => `column_collation`
    * **Action**: apply the renames.
* `connection::next_layer()` has been renamed to `connection::stream()`.
    * **Action**: apply the rename.
* `no_statement_params` has been removed.
    * **Action**: replace it by an empty `std::make_tuple()`.
