import clang.cindex
import errno
import glob
import os
import re
import sys
import tempfile
import difflib
import argparse
from argparse import ArgumentParser, RawTextHelpFormatter
from clang.cindex import Config, Index, CursorKind
from typing import List, Tuple, Optional


def processCLIArgs() -> argparse.Namespace:
    """
    Define and process the command line interface to the extract_signatures.py script.

    Returns:
        argparse.Namespace: The parsed command-line arguments.
    """
    cli_description = "Generate forward decls for the stan math library."
    cli_epilog = "See more information at: https://github.com/stan-dev/math"
    parser = ArgumentParser(
        description=cli_description,
        epilog=cli_epilog,
        formatter_class=RawTextHelpFormatter,
    )
    # Now define all the rules of the command line args and opts
    parser.add_argument(
        "--input_base_path", type=str, default="./", help="Base path to stan math"
    )
    parser.add_argument(
        "--math_path", type=str, default="stan/math/", help="Base path to stan math"
    )
    parser.add_argument(
        "--debug_functions",
        type=bool,
        default=False,
        help="Enable debug output for functions",
    )
    parser.add_argument(
        "--debug_classes",
        type=bool,
        default=False,
        help="Enable debug output for classes",
    )
    parser.add_argument(
        "--y",
        type=bool,
        default=False,
        help="Automatically accept all changes",
    )
    parser.add_argument(
        "--input_file",
        type=str,
        default="mix.hpp",
        help="Input file to parse",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        default="forward_decl.hpp",
        help="File to store the functions forward declarations",
    )
    # And parse the command line against those rules
    return parser.parse_args()


manual_decl_funcs: List[str] = []

manual_decl_classes: List[str] = []  # ["fvar", "var_value", "vari_value"]


def silentremove(filename: str) -> None:
    """
    Silently remove a file if it exists; do nothing if it does not.

    Args:
        filename (str): The path to the file to remove.
    """
    try:
        os.remove(filename)
    except OSError as e:  # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occurred


def get_template_parameters(cursor: clang.cindex.Cursor) -> List[str]:
    """
    Extract template parameters from a given cursor.

    Args:
        cursor (clang.cindex.Cursor): The cursor from which to extract template parameters.

    Returns:
        List[str]: A list of template parameter strings.
    """
    template_params = []
    for child in cursor.get_children():
        if child.kind == clang.cindex.CursorKind.TEMPLATE_TYPE_PARAMETER:
            # Template type parameter
            param_name = child.spelling
            # Check for parameter pack
            if any(x.spelling == "..." for x in child.get_tokens()):
                template_params.append(f"typename... {param_name}")
            else:
                template_params.append(f"typename {param_name}")
        elif child.kind == clang.cindex.CursorKind.TEMPLATE_NON_TYPE_PARAMETER:
            # Non-type template parameter
            param_type = child.type.spelling
            param_name = child.spelling
            template_params.append(f"{param_type} {param_name}")
        elif child.kind == clang.cindex.CursorKind.TEMPLATE_TEMPLATE_PARAMETER:
            # Template template parameter
            param_name = child.spelling
            template_params.append(f"template <class...> class {param_name}")
    return template_params


def get_function_declaration(
    func_decl_cursor: clang.cindex.Cursor,
) -> Tuple[bool, str, str, str, str, List[str], str]:
    """
    Extract components of a function declaration from a function cursor.

    Args:
        func_decl_cursor (clang.cindex.Cursor): The cursor representing the function declaration.

    Returns:
        Tuple[bool, str, str, str, str, List[str], str]: A tuple containing:
            - has_defaults (bool): True if the function has default parameters.
            - storage_class (str): The storage class of the function.
            - attributes (str): Attributes like 'virtual', 'inline', etc.
            - return_type (str): The return type of the function.
            - func_name (str): The function name.
            - params (List[str]): A list of parameter strings.
            - func_qualifiers (str): Function qualifiers like 'const', 'noexcept'.
    """
    # Sometimes clang makes auto return types into the actual literal return type.
    # so we need to check if the return type is auto and not use
    # the result type if so
    if func_decl_cursor.result_type.kind == clang.cindex.TypeKind.AUTO:
        return_type = "auto"
    else:
        return_type = func_decl_cursor.result_type.spelling

    # Function name
    func_name = func_decl_cursor.spelling

    # Function parameters
    params = []
    for param_cursor in func_decl_cursor.get_children():
        if param_cursor.kind == clang.cindex.CursorKind.PARM_DECL:
            if any([x.spelling == "=" for x in param_cursor.get_tokens()]):
                return True, "", "", "", "", [], ""
            if func_decl_cursor.spelling == "operator==" and False:
                print(
                    " ".join([x.spelling for x in func_decl_cursor.get_tokens()])
                )
                import pdb

                pdb.set_trace()

            param_type = param_cursor.type.spelling
            param_name = param_cursor.spelling
            params.append(f"{param_type} {param_name}")

    # Function qualifiers (const, noexcept)
    func_qualifiers = ""
    if func_decl_cursor.is_const_method():
        func_qualifiers += " const"
    if (
        func_decl_cursor.exception_specification_kind
        == clang.cindex.ExceptionSpecificationKind.BASIC_NOEXCEPT
    ):
        func_qualifiers += " noexcept"
    # Storage class (static, virtual)
    storage_class = ""
    if func_decl_cursor.storage_class == clang.cindex.StorageClass.STATIC:
        storage_class = "static "
    attributes = ""
    valid_attributes = [
        "virtual",
        "explicit",
        "constexpr",
        "inline",
        "STAN_COLD_PATH",
    ]
    for token in func_decl_cursor.get_tokens():
        (kind, spelling) = token.kind, token.spelling
        if kind == clang.cindex.TokenKind.KEYWORD and spelling in valid_attributes:
            attributes += spelling + " "
        # If we hit the function name, stop
        elif (
            kind == clang.cindex.TokenKind.IDENTIFIER
            and spelling == func_decl_cursor.spelling
        ):
            break
    return False, storage_class, attributes, return_type, func_name, params, func_qualifiers


def extract_function_template_declaration(
    cursor: clang.cindex.Cursor,
) -> str:
    """
    Extracts the function declaration from a Cursor object of kind FUNCTION_TEMPLATE.

    Args:
        cursor (clang.cindex.Cursor): The FUNCTION_TEMPLATE cursor.

    Returns:
        str: The reconstructed function declaration as a string.
    """
    if cursor.spelling in manual_decl_funcs:
        return ""
    if any([x.spelling.strip() == "friend" for x in cursor.get_tokens()]):
        return ""
    # Extract template parameters
    template_params = get_template_parameters(cursor)
    # Extract function declaration components
    (
        has_defaults,
        storage_class,
        attributes,
        return_type,
        func_name,
        params,
        func_qualifiers,
    ) = get_function_declaration(cursor)
    if has_defaults:
        return ""
    # Build the function declaration string
    template_str = (
        f"template <{', '.join(template_params)}>\n" if template_params else ""
    )
    params_str = ", ".join(params)
    func_decl_str = (
        f"{storage_class} {attributes} {return_type} {func_name}({params_str})"
        f"{func_qualifiers};"
    )

    # Combine template and function declaration
    full_decl = template_str + func_decl_str

    return full_decl


def extract_class_template_declaration(cursor: clang.cindex.Cursor) -> str:
    """
    Extracts the class forward declaration from a Cursor object of kind CLASS_TEMPLATE, CLASS_DECL, or STRUCT_DECL.

    Args:
        cursor (clang.cindex.Cursor): The class cursor.

    Returns:
        str: The reconstructed class forward declaration as a string.
    """
    # Skip certain classes if needed
    if cursor.spelling in manual_decl_classes:
        return ""
    if cursor.spelling == "":
        # Anonymous structs/unions can be skipped
        return ""
    if cursor.spelling == "append_return_type" and False:
        import pdb

        pdb.set_trace()
    # Extract template parameters
    template_params = get_template_parameters(cursor)
    # Class name
    class_name = cursor.spelling
    # Class kind (class, struct, union)
    # structs that templates are treated as a clang.cindex.CursorKind.CLASS_TEMPLATE
    # so we need to parse the tokens to see if we can find the struct or class keyword
    class_kind = ""
    for token in cursor.get_tokens():
        if token.kind == clang.cindex.TokenKind.KEYWORD:
            if token.spelling == "struct":
                class_kind = "struct"
                break
            elif token.spelling == "class":
                class_kind = "class"
                break
        elif (
            token.kind == clang.cindex.TokenKind.IDENTIFIER
            and token.spelling == class_name
        ):
            break
    # We need to check if this is just a specialization
    # If it is we can skip it
    raw_tokens = [
        (token.spelling, token.kind) for token in cursor.get_tokens()
    ]
    for i in range(len(raw_tokens)):
        if raw_tokens[i][0] == class_name and len(raw_tokens) > i + 1:
            if raw_tokens[i + 1][0] == "<":
                return ""
            else:
                break

    # Attributes (e.g., final)
    attributes = ""
    valid_attributes = ["final", "sealed", "abstract"]
    tokens = list(cursor.get_tokens())
    for token in tokens:
        if token.spelling in valid_attributes:
            attributes += f" {token.spelling}"

    # Build the class declaration string
    template_str = (
        f"template <{', '.join(template_params)}>\n" if template_params else ""
    )
    class_decl_str = f"{class_kind} {class_name};"

    # Combine template and class declaration
    full_decl = template_str + class_decl_str

    return full_decl


def get_functions_and_classes_in_namespace(
    tu: clang.cindex.TranslationUnit,
    input_path: str,
    namespaces: Tuple[str, ...],
    debug_function_nodes: bool = False,
    debug_class_nodes: bool = False,
) -> Tuple[List[clang.cindex.Cursor], List[clang.cindex.Cursor]]:
    """
    Retrieves functions and classes within specified namespaces from a translation unit.

    Args:
        tu (clang.cindex.TranslationUnit): The translation unit to analyze.
        input_path (str): The path to the input file.
        namespaces (Tuple[str, ...]): A tuple of namespace names to search within.
        debug_function_nodes (bool, optional): If True, print debug information for functions.
        debug_class_nodes (bool, optional): If True, print debug information for classes.

    Returns:
        Tuple[List[clang.cindex.Cursor], List[clang.cindex.Cursor]]: A tuple containing:
            - functions (List[clang.cindex.Cursor]): A list of function cursors.
            - classes (List[clang.cindex.Cursor]): A list of class cursors.
    """
    functions: List[clang.cindex.Cursor] = []
    classes: List[clang.cindex.Cursor] = []

    def visit_node(node: clang.cindex.Cursor, namespaces: Tuple[str, ...]) -> None:
        if (
            node.extent.start.file is None
            or node.extent.start.file.name.find("stan/math") == -1
            or node.extent.start.file.name.find("stan/math/prim") >= 0
        ):
            return
        if node.kind == clang.cindex.CursorKind.NAMESPACE:
            if node.spelling == namespaces[0]:
                if len(namespaces) == 1:
                    # We are in the target namespace
                    for child in node.get_children():
                        visit_node(child, namespaces)
                else:
                    # Need to go deeper into nested namespaces
                    for child in node.get_children():
                        visit_node(child, namespaces[1:])
        elif node.kind in [
            clang.cindex.CursorKind.FUNCTION_DECL,
            clang.cindex.CursorKind.FUNCTION_TEMPLATE,
        ]:
            # Check if the function is a free function (not inside a class)
            if (
                node.semantic_parent.kind == clang.cindex.CursorKind.NAMESPACE
                and node.semantic_parent.spelling == namespaces[-1]
            ):
                if debug_function_nodes:
                    print(
                        node.location.file,
                        node.location.line,
                        node.location.column,
                        node.kind,
                        node.spelling,
                    )
                    import pdb

                    pdb.set_trace()
                if str(node.location.file).find(input_path) == -1:
                    return
                functions.append(node)
        elif node.kind in [
            clang.cindex.CursorKind.CLASS_TEMPLATE,
            clang.cindex.CursorKind.CLASS_DECL,
            clang.cindex.CursorKind.STRUCT_DECL,
        ]:
            # Collect classes and structs
            if (
                node.semantic_parent.kind == clang.cindex.CursorKind.NAMESPACE
                and node.semantic_parent.spelling == namespaces[-1]
            ):
                if debug_class_nodes:
                    print(
                        node.location.file,
                        node.location.line,
                        node.location.column,
                        node.kind,
                        node.spelling,
                    )
                    import pdb

                    pdb.set_trace()
                if str(node.location.file).find(input_path) == -1:
                    return
                classes.append(node)
        else:
            # Continue traversing other nodes
            for child in node.get_children():
                visit_node(child, namespaces)

    visit_node(tu.cursor, namespaces)
    return functions, classes


def main(inputs: argparse.Namespace) -> None:
    """
    Main function to generate forward declarations for the Stan Math library.

    Args:
        inputs (argparse.Namespace): The parsed command-line arguments.
    """
    # Parse the command line arguments
    clang.cindex.Config.set_library_path(
        "/opt/homebrew/Cellar/llvm/19.1.2/lib/"
    )
    base_path = inputs.input_base_path
    filename = inputs.math_path + inputs.input_file
    print("Parsing: ", filename)
    comp_args = [
        "-std=c++17",
        "-D_REENTRANT",
        "-DBOOST_DISABLE_ASSERTS",
        "-DSTAN_THREADS",
        "-DSTAN_MATH_FORWARD_DECL_HPP",
        "-fparse-all-comments",
        "-O0",
        "-g",
        "-x",
        "c++",
        "-fno-delayed-template-parsing",
        "-resource-dir",
        "/opt/homebrew//Cellar/llvm/19.1.2/lib/clang/19",
        "-I" + base_path,
        "-I" + base_path + "lib/tbb_2020.3/include",
        "-I" + base_path + "lib/eigen_3.4.0",
        "-I" + base_path + "lib/boost_1.84.0",
        "-I" + base_path + "lib/sundials_6.1.1/include",
        "-I" + base_path + "lib/sundials_6.1.1/src/sundials",
    ]  # Add any necessary compiler flags here
    print("args: ", comp_args)
    # Generate output in a temporary file
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_output_file:
        temp_output_filename = temp_output_file.name
        # Write to the temporary file
        temp_output_file.write(
            """
#ifndef STAN_MATH_FORWARD_DECL_HPP
#define STAN_MATH_FORWARD_DECL_HPP
#include <stan/math/manual_forward_decl.hpp>
namespace stan {
namespace math {
namespace internal {
\n"""
        )
        options = (
            clang.cindex.TranslationUnit.PARSE_DETAILED_PROCESSING_RECORD
            | clang.cindex.TranslationUnit.PARSE_PRECOMPILED_PREAMBLE
            | clang.cindex.TranslationUnit.PARSE_INCLUDE_BRIEF_COMMENTS_IN_CODE_COMPLETION
            | clang.cindex.TranslationUnit.PARSE_CACHE_COMPLETION_RESULTS
        )
        try:
            translation_unit = clang.cindex.TranslationUnit.from_source(
                filename, comp_args, options=options
            )
        except Exception as e:
            import pdb

            pdb.set_trace()
            print(f"Error parsing {filename}: {e}")
            return
        for diag in translation_unit.diagnostics:
            print("Diagnostic:\n")
            print(diag)
        functions, classes = get_functions_and_classes_in_namespace(
            translation_unit,
            inputs.math_path,
            ("stan", "math", "internal"),
            debug_function_nodes=inputs.debug_functions,
            debug_class_nodes=inputs.debug_classes,
        )
        for class_node in classes:
            forward_declaration = extract_class_template_declaration(class_node)
            temp_output_file.write(forward_declaration + "\n\n")
        for function in functions:
            function_declaration = extract_function_template_declaration(function)
            temp_output_file.write(
                function_declaration + "\n\n"
            )  # Add two newlines for readability
        temp_output_file.write(
            """
} // namespace internal
\n"""
        )
        functions, classes = get_functions_and_classes_in_namespace(
            translation_unit,
            inputs.math_path,
            ("stan", "math"),
            debug_function_nodes=inputs.debug_functions,
            debug_class_nodes=inputs.debug_classes,
        )
        for class_node in classes:
            forward_declaration = extract_class_template_declaration(class_node)
            temp_output_file.write(forward_declaration + "\n\n")
        for function in functions:
            function_declaration = extract_function_template_declaration(function)
            temp_output_file.write(
                function_declaration + "\n\n"
            )  # Add two newlines for readability
        temp_output_file.write(f"// Classes Parsed: {len(classes)}\n")
        temp_output_file.write(f"// Functions Parsed: {len(functions)}\n")
        temp_output_file.write(
            """
} // namespace math
} // namespace stan
#endif // STAN_MATH_FORWARD_DECL_HPP
\n"""
        )
        temp_output_file.flush()
    # Read the old and new outputs
    old_output = []
    if os.path.exists(inputs.output_file):
        with open(inputs.output_file, "r") as old_file:
            old_output = old_file.readlines()
    else:
        old_output = []
    with open(temp_output_filename, "r") as new_file:
        new_output = new_file.readlines()
    # Show the diff to the user
    diff = difflib.unified_diff(
        old_output, new_output, fromfile=inputs.output_file, tofile="New Output"
    )
    diff_text = "".join(diff)
    if diff_text:
        print("Differences between the old and new outputs:\n")
        print(diff_text)
        # Ask the user if they want to overwrite the old file
        if inputs.y:
            user_input = "y"
        else:
            user_input = input(
                f"Do you want to overwrite {inputs.output_file} with the new output? (y/n): "
            )
        if user_input.lower() == "y":
            # Overwrite the old file with the new output
            with open(inputs.output_file, "w") as output_file:
                output_file.writelines(new_output)
            print(f"{inputs.output_file} has been updated.")
        else:
            print("No changes were made.")
    else:
        print("No differences found. The output file remains unchanged.")
    # Clean up the temporary file
    os.remove(temp_output_filename)


if __name__ == "__main__":
    inputs = processCLIArgs()
    main(inputs)
