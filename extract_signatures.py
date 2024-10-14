import clang.cindex
import errno
import glob
import os
import re
import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from clang.cindex import Config, Index, CursorKind


def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e:  # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occurred

import clang.cindex

def extract_function_template_declaration(cursor):
    """
    Extracts the function declaration from a Cursor object of kind FUNCTION_TEMPLATE.

    Args:
        cursor (clang.cindex.Cursor): The FUNCTION_TEMPLATE cursor.

    Returns:
        str: The reconstructed function declaration as a string.
    """

    def get_template_parameters(cursor):
        template_params = []
        for child in cursor.get_children():
            if child.kind == clang.cindex.CursorKind.TEMPLATE_TYPE_PARAMETER:
                # Template type parameter
                param_name = child.spelling
                # Check for default type
                default_type = None
                for gc in child.get_children():
                    if gc.kind == clang.cindex.CursorKind.TYPE_REF:
                        default_type = gc.spelling
                if default_type:
                    template_params.append(f"typename {param_name} = {default_type}")
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

    def get_function_decl(cursor):
        if cursor.kind in [clang.cindex.CursorKind.FUNCTION_DECL,
                           clang.cindex.CursorKind.CXX_METHOD]:
            return cursor
        for child in cursor.get_children():
            result = get_function_decl(child)
            if result is not None:
                return result
        return None
    def get_function_declaration(func_decl_cursor):
        # Return type
        return_type = func_decl_cursor.result_type.spelling

        # Function name
        func_name = func_decl_cursor.spelling

        # Function parameters
        params = []
        for param_cursor in func_decl_cursor.get_arguments():
            param_type = param_cursor.type.spelling
            param_name = param_cursor.spelling
            params.append(f"{param_type} {param_name}")

        # Function qualifiers (const, noexcept)
        func_qualifiers = ''
        if func_decl_cursor.is_const_method():
            func_qualifiers += ' const'
        if func_decl_cursor.exception_specification_kind == clang.cindex.ExceptionSpecificationKind.BasicNoexcept:
            func_qualifiers += ' noexcept'

        # Storage class (static, virtual)
        storage_class = ''
        if func_decl_cursor.is_static_method():
            storage_class = 'static '
        elif func_decl_cursor.is_virtual_method():
            storage_class = 'virtual '

        return storage_class, return_type, func_name, params, func_qualifiers

    # Extract template parameters
    template_params = get_template_parameters(cursor)

    # Find the function declaration cursor
    func_decl_cursor = get_function_decl(cursor)
    if func_decl_cursor is None:
        raise ValueError("Function declaration not found in the FUNCTION_TEMPLATE cursor.")

    # Extract function declaration components
    storage_class, return_type, func_name, params, func_qualifiers = get_function_declaration(func_decl_cursor)

    # Build the function declaration string
    template_str = f"template <{', '.join(template_params)}>\n" if template_params else ""
    params_str = ', '.join(params)
    func_decl_str = f"{storage_class}{return_type} {func_name}({params_str}){func_qualifiers};"

    # Combine template and function declaration
    full_decl = template_str + func_decl_str

    return full_decl


def processCLIArgs():
    """
    Define and process the command line interface to the extract_signatures.py script.
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
        "--math_path", type=str, default="./stan/math/", help="Base path to stan math"
    )
    parser.add_argument(
        "--debug_functions", type=bool, default=False, help="Base path to stan math"
    )
    parser.add_argument(
        "--debug_classes", type=bool, default=False, help="Base path to stan math"
    )
    parser.add_argument(
        "--input_file",
        type=str,
        default="stan/math/mix.hpp",
        help="Input file to parse",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        default="./forward_decls.hpp",
        help="File to store the functions forward declarations",
    )
    # And parse the command line against those rules
    return parser.parse_args()


def get_functions_and_classes_in_namespace(
    tu, input_path, debug_function_nodes=False, debug_class_nodes=False
):
    functions = []
    classes = []

    def visit_node(node, namespaces):
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
                if node.extent.start.file.name.find("forward_decls.hpp") > 0 or \
                    str(node.location.file).find(input_path) == -1:
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
                if node.extent.start.file.name.find("forward_decls.hpp") > 0 or \
                    str(node.location.file).find(input_path) == -1:
                    return
                classes.append(node)
        else:
            # Continue traversing other nodes
            for child in node.get_children():
                visit_node(child, namespaces)

    visit_node(tu.cursor, ("stan", "math"))
    return functions, classes


def generate_forward_declaration(function_node):
    # Skip functions with default values or other manually specified forward declarations
    if function_node.spelling in ["init_threadpool_tbb", 
                                  "finite_diff_grad_hessian",
                                    "grad_reg_inc_gamma", 
                                    "grad_reg_lower_inc_gamma",
                                    "grad_pFq",
                                    "grad_2F1"]:
        return ""
    # Get the start of the function declaration
    start = function_node.extent.start
    # Find the location where the function body starts ('{' or ';')
    tokens = list(function_node.get_tokens())
    body_start = None
    do_debug = False
    # Get the start of the class/struct declaration
    start = function_node.extent.start
    # Find the location where the class/struct body starts ('{' or ';')
    tokens = list(function_node.get_tokens())
    body_start = None
    token_full = ""
    do_debug = False
    if function_node.spelling.find("unnamed") > 0:
        return ""
    
    for i in range(len(tokens)):
        token = tokens[i]
        token_full += token.spelling
        # Look for the first '{', ';', or ' : ' to show end of class/struct declaration
        if (
            token.spelling == "{"
            or token.spelling == ";"
            or (token.spelling == ":" and tokens[i + 1].spelling != ":")
            or (i > 0 and token.spelling == ":" and tokens[i - 1].spelling != ":")
        ):
            body_start = token.location
            break
#    print("Token Full: ", token_full)
    if body_start is None:
        # If no body is found, use the end of the class node
        body_start = function_node.extent.end
    # Ensure that the start and body_start are in the same file
    if start.file.name != body_start.file.name:
        raise Exception("Start and body_start are in different files")
    # Read the source code between start and body_start
    file_name = start.file.name
    with open(file_name, "r") as f:
        code = f.read()
    start_offset = start.offset
    end_offset = body_start.offset
    # Add comment for file and line number
    while code[start_offset] != "\n":
        start_offset -= 1
    while code[end_offset] != "{":
        end_offset -= 1
    forward_declaration = code[start_offset:end_offset].strip()
#    import pdb; pdb.set_trace()
    if not forward_declaration.endswith(";"):
        forward_declaration += ";"
    raw_tokens = [(token.spelling, token.kind) for token in tokens]
    has_default_value = any([x[0].strip() == "=" for x in raw_tokens])
    forward_declaration = forward_declaration.replace("* = nullptr", "*")
    forward_declaration = forward_declaration.replace(" = void", "")
    forward_declaration = forward_declaration.replace(" = double", "")
    forward_declaration = forward_declaration.replace(" = std::tuple<>", "")
    base_pattern = r"((?:[^<>]++|<[^<>]*>)*?)"
    forward_declaration = re.sub(r"= return_type_t<"+ base_pattern + r">,", ",", forward_declaration)
    forward_declaration = re.sub(r"= scalar_type_t<"+ base_pattern + r">,", ",", forward_declaration)
    forward_declaration = re.sub(r"= boost::optional<"+ base_pattern + r">,", ",", forward_declaration)
    forward_declaration = re.sub(r"= promote_scalar_t<"+ base_pattern + r">,", ",", forward_declaration)
    forward_declaration = re.sub(r"= !is_constant<"+ base_pattern + r">::value,", ",", forward_declaration)
    forward_declaration = re.sub(r"= return_type_t<"+ base_pattern + r">,", ",", forward_declaration)
    forward_declaration = f"// {file_name}:{start.line}\n" + forward_declaration
        # Regular expression pattern for matching `return_type_t<...>` with nested brackets
    # Replace all matches of `return_type_t<...>` with an empty string
    # Add a semicolon if necessary
    if do_debug:# or function_node.spelling == "zero_adjoints":
        print("Function Decl: \n", forward_declaration)
        import pdb
        pdb.set_trace()
    return forward_declaration


def generate_class_forward_declaration(class_node):
    # Get the start of the class/struct declaration
    start = class_node.extent.start
    # Find the location where the class/struct body starts ('{' or ';')
    tokens = list(class_node.get_tokens())
    body_start = None
    token_full = ""
    do_debug = False
    if class_node.spelling.find("unnamed") > 0:
        return ""    
    for i in range(len(tokens)):
        token = tokens[i]
        token_full += token.spelling
        # Look for the first '{', ';', or ' : ' to show end of class/struct declaration
        if (
            token.spelling == "{"
            or token.spelling == ";"
            or (token.spelling == ":" and tokens[i + 1].spelling != ":")
            or (i > 0 and token.spelling == ":" and tokens[i - 1].spelling != ":")
        ):
            body_start = token.location
            break
#    print("Token Full: ", token_full)
    if body_start is None:
        # If no body is found, use the end of the class node
        body_start = class_node.extent.end
    # Ensure that the start and body_start are in the same file
    if start.file.name != body_start.file.name:
        raise Exception("Start and body_start are in different files")
    # Read the source code between start and body_start
    file_name = start.file.name
    with open(file_name, "r") as f:
        code = f.read()
    start_offset = start.offset
    end_offset = body_start.offset
    # Add comment for file and line number
    forward_declaration = f"// {file_name}:{start.line}\n"
    forward_declaration += code[start_offset:end_offset].strip()
    if not forward_declaration.endswith(";"):
        forward_declaration += ";"
    raw_tokens = [(token.spelling, token.kind) for token in tokens]
    has_default_value = any([x[0].strip() == "=" for x in raw_tokens])
    forward_declaration = forward_declaration.replace(" = void", "")
    forward_declaration = forward_declaration.replace(" = double", "")
    forward_declaration = forward_declaration.replace(" = std::tuple<>", "")
    forward_declaration = re.sub(r"= scalar_type_t<[^>]+?>", "", forward_declaration)
    forward_declaration = forward_declaration.replace(" final;", ";")
    
    # Add a semicolon if necessary
    if do_debug:
        import pdb
        pdb.set_trace()
    return forward_declaration


def main(inputs):
    # Parse the command line arguments
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
        "-x", "c++",
        "-fno-delayed-template-parsing",
        "-I" + base_path,
        "-I" + base_path + "lib/tbb_2020.3/include",
        "-I" + base_path + "lib/eigen_3.4.0",
        "-I" + base_path + "lib/boost_1.84.0",
        "-I" + base_path + "lib/sundials_6.1.1/include",
        "-I" + base_path + "lib/sundials_6.1.1/src/sundials",
    ]  # Add any necessary compiler flags here
    print("args: ", comp_args)
    silentremove(inputs.output_file)
    #index = clang.cindex.Index.create(excludeDecls=False)
    options = (
    clang.cindex.TranslationUnit.PARSE_INCOMPLETE |
    clang.cindex.TranslationUnit.PARSE_DETAILED_PROCESSING_RECORD |
    clang.cindex.TranslationUnit.PARSE_PRECOMPILED_PREAMBLE |
    clang.cindex.TranslationUnit.PARSE_INCLUDE_BRIEF_COMMENTS_IN_CODE_COMPLETION |
    clang.cindex.TranslationUnit.PARSE_CACHE_COMPLETION_RESULTS
    )
    try:
        translation_unit = clang.cindex.TranslationUnit.from_source(filename, comp_args, options=options)
    except Exception as e:
        import pdb; pdb.set_trace()
        print(f"Error parsing {filename}: {e}")
        return
    for diag in translation_unit.diagnostics:
        print("Diagnostic:\n")
        print(diag)
    functions, classes = get_functions_and_classes_in_namespace(
        translation_unit,
        inputs.math_path,
        debug_function_nodes=inputs.debug_functions,
        debug_class_nodes=inputs.debug_classes,
    )
    # Open the output file and write the forward declarations
    with open(inputs.output_file, "w") as output_file:
        output_file.write(
            """
#ifndef STAN_MATH_FORWARD_DECL_HPP
#define STAN_MATH_FORWARD_DECL_HPP
#include <stan/math/manual_forward_decls.hpp>
namespace stan {
namespace math {
\n"""
        )
        for class_node in classes:
            forward_decl = generate_class_forward_declaration(class_node)
            output_file.write(forward_decl + "\n\n")
        for function in functions:
            try:
                function_declaration = extract_function_template_declaration(function)
                print(function_declaration)
            except Exception as e:
                forward_decl = generate_forward_declaration(function)
                pass
            output_file.write(forward_decl + "\n\n")  # Add two newlines for readability
        output_file.write(f"// Functions Parsed: {len(functions)}\n")
        output_file.write(
            """
} // namespace math
} // namespace stan
#endif // STAN_MATH_FORWARD_DECL_HPP
\n"""
        )
    print(f"Forward declarations have been written to {inputs.output_file}")


if __name__ == "__main__":
    inputs = processCLIArgs()
    main(inputs)

"""
# Your list of tokens
tokens = [
    ("template", clang.cindex.TokenKind.KEYWORD),
    ("<", clang.cindex.TokenKind.PUNCTUATION),
    ("typename", clang.cindex.TokenKind.KEYWORD),
    ("T_y", clang.cindex.TokenKind.IDENTIFIER),
    (",", clang.cindex.TokenKind.PUNCTUATION),
    ("require_any_t", clang.cindex.TokenKind.IDENTIFIER),
    ("<", clang.cindex.TokenKind.PUNCTUATION),
    ("is_matrix", clang.cindex.TokenKind.IDENTIFIER),
    ("<", clang.cindex.TokenKind.PUNCTUATION),
    ("T_y", clang.cindex.TokenKind.IDENTIFIER),
    (">", clang.cindex.TokenKind.PUNCTUATION),
    (",", clang.cindex.TokenKind.PUNCTUATION),
    ("is_prim_or_rev_kernel_expression", clang.cindex.TokenKind.IDENTIFIER),
    ("<", clang.cindex.TokenKind.PUNCTUATION),
    ("T_y", clang.cindex.TokenKind.IDENTIFIER),
    (">>", clang.cindex.TokenKind.PUNCTUATION),
    ("*", clang.cindex.TokenKind.PUNCTUATION),
    ("=", clang.cindex.TokenKind.PUNCTUATION),
    ("nullptr", clang.cindex.TokenKind.KEYWORD),
    (">", clang.cindex.TokenKind.PUNCTUATION),
    ("inline", clang.cindex.TokenKind.KEYWORD),
    ("void", clang.cindex.TokenKind.KEYWORD),
    ("check_square", clang.cindex.TokenKind.IDENTIFIER),
    ("(", clang.cindex.TokenKind.PUNCTUATION),
    ("const", clang.cindex.TokenKind.KEYWORD),
    ("char", clang.cindex.TokenKind.KEYWORD),
    ("*", clang.cindex.TokenKind.PUNCTUATION),
    ("function", clang.cindex.TokenKind.IDENTIFIER),
    (",", clang.cindex.TokenKind.PUNCTUATION),
    ("const", clang.cindex.TokenKind.KEYWORD),
    ("char", clang.cindex.TokenKind.KEYWORD),
    ("*", clang.cindex.TokenKind.PUNCTUATION),
    ("name", clang.cindex.TokenKind.IDENTIFIER),
    (",", clang.cindex.TokenKind.PUNCTUATION),
    ("const", clang.cindex.TokenKind.KEYWORD),
    ("T_y", clang.cindex.TokenKind.IDENTIFIER),
    ("&", clang.cindex.TokenKind.PUNCTUATION),
    ("y", clang.cindex.TokenKind.IDENTIFIER),
    (")", clang.cindex.TokenKind.PUNCTUATION),
    ("{", clang.cindex.TokenKind.PUNCTUATION),
    # Function body tokens...
]

# File information
file_name = "./stan/math/prim/err/check_square.hpp"
line_number = 20

# Generate forward declaration
forward_declaration = generate_forward_declaration_from_tokens(
    tokens, file_name, line_number
)
print(forward_declaration)
"""
