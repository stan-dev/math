import sys
import clang.cindex
import glob
import errno
import os
from argparse import ArgumentParser, RawTextHelpFormatter


def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e:  # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occurred


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
        "--input_base_path", type=str, default="./math/", help="Base path to stan math"
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
    tu, debug_function_nodes=False, debug_class_nodes=False
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
                if node.extent.start.file.name.find("forward_decls.hpp") > 0:
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
                if node.extent.start.file.name.find("forward_decls.hpp") > 0:
                    return
                classes.append(node)
        else:
            # Continue traversing other nodes
            for child in node.get_children():
                visit_node(child, namespaces)

    visit_node(tu.cursor, ("stan", "math"))
    return functions, classes



def generate_forward_declaration_from_tokens(tokens):
    """
    Generate a forward declaration from a list of tokens.

    Args:
        tokens (list of tuples): List of (spelling, kind) tuples representing tokens.
        file_name (str): The name of the file where the function is declared.
        line_number (int): The line number where the function is declared.

    Returns:
        str: The forward declaration as a formatted string.
    """
    code_tokens = []
    nesting_level = 0
    in_template_parameter_list = False
    in_function_parameter_list = False
    i = 0
    # Collect tokens up to the first '{' or ';' (excluding it)
    while i < len(tokens):
        spelling, kind = tokens[i]
        if spelling == '{' or spelling == ';':
            break
        code_tokens.append((spelling, kind))
        i += 1
    # Now, reconstruct the code from code_tokens with appropriate formatting
    code_lines = []
    current_line = ''
    indent = ''
    nesting_level = 0
    in_template_parameter_list = False
    in_default_value = False
    default_nesting_level = 0
    i = 0
    template_str = ""
    signature_str = ""
    debug_function_decl = False
    if debug_function_decl:
        print("BEGIN\n\n\n")
    first_entry = True
    while i < len(code_tokens):
        if debug_function_decl:
            print("-------------------")
            print("I: ", i)
            print("Code Tokens: ", code_tokens[i])
            print("template_str: ", template_str)
            print("signature_str: ", signature_str)
            print("in_template_parameter_list: ", in_template_parameter_list)
            print("nesting_level: ", nesting_level)
            print("in_default_value: ", in_default_value)
            print("default_nesting_level: ", default_nesting_level)
            if (code_tokens[i][0] == ">"):
              import pdb; pdb.set_trace()
        if not in_template_parameter_list and first_entry:
          first_entry = False
          if debug_function_decl:
            print("----\n\nENTERED FUNCTION DECLARATION----\n\n")

        spelling, kind = code_tokens[i]
        # Handle 'template' keyword
        if spelling == 'template':
            template_str += 'template '
            i += 1
            continue
        # Handle entering template parameter list
        if spelling == '<' and (i > 0 and code_tokens[i - 1][0] == 'template'):
            in_template_parameter_list = True
            nesting_level += 1
            template_str += '<'
            i += 1
            continue
        # Handle template parameter list
        if in_template_parameter_list:
            if in_default_value:
              if spelling == '<':
                  default_nesting_level += 1
                  i += 1
                  continue
              elif spelling == '>':
                  # If we never nest (* = nullptr>) then we are done with the default value
                  if default_nesting_level == 0:
                    in_default_value = False
                    nesting_level -= 1
                  else:
                    default_nesting_level -= 1
                  if nesting_level == 0:
                      template_str += '>'
                      in_template_parameter_list = False
                  i += 1
                  continue
              elif spelling == '>>':
                  if default_nesting_level == 0:
                    in_default_value = False
                    nesting_level -= 2
                    if nesting_level == 0:
                        template_str += '>>'
                        in_template_parameter_list = False
                  else:
                    default_nesting_level -= 1
                    # One of them could be in the default value
                    if default_nesting_level == 0:
                      in_default_value = False                 
                      nesting_level -= 1
                      if nesting_level == 0:
                          template_str += '>'
                          in_template_parameter_list = False
                  i += 1
                  continue
              elif spelling == ">" and default_nesting_level == 0:
                  in_default_value = False
                  template_str += '>'
                  nesting_level -= 1
                  if nesting_level == 0:
                      in_template_parameter_list = False
                  i += 1
                  continue
              else:
                  i += 1
                  continue
            if spelling == '<':
                nesting_level += 1
                template_str += '<'
            elif spelling == '>':
                nesting_level -= 1
                template_str += '>'
                if nesting_level == 0:
                    in_template_parameter_list = False
            elif spelling == '>>':
                nesting_level -= 2
                template_str += '>>'
                if nesting_level == 0:
                    in_template_parameter_list = False
            elif spelling == ',':
                # Add comma and line break
                template_str += ',\n    '
            else:
                if kind == clang.cindex.TokenKind.PUNCTUATION and spelling == '=':
                    in_default_value = True
                elif kind == clang.cindex.TokenKind.PUNCTUATION or kind == clang.cindex.TokenKind.IDENTIFIER:
                    template_str += spelling
                else:
                    template_str += spelling + ' '
            i += 1
            continue
        # Handle function declaration
        else:
            if spelling == '(':
                # Start of function parameter list
                in_function_parameter_list = True
                signature_str += '('
            elif spelling == ')':
                # End of function parameter list
                in_function_parameter_list = False
                signature_str += ')'
            elif spelling == ',' and in_function_parameter_list:
                # Function parameter separator
                signature_str += ', '
            else:
                if kind == clang.cindex.TokenKind.KEYWORD or kind == clang.cindex.TokenKind.IDENTIFIER:
                    signature_str += spelling + ' '
                else:
                    signature_str += spelling
            i += 1
            continue
    # Clean up extra spaces
    if debug_function_decl:
        print("Template Str: ", template_str)
        print("Signature Str: ", signature_str)
    current_line = template_str + "\n" + signature_str
    # Ensure it ends with a semicolon
    if not current_line.endswith(';'):
        current_line += ';'
    # Add comment for file and line number
    forward_declaration = current_line
    return forward_declaration


def generate_forward_declaration(function_node):
    # Get the start of the function declaration
    start = function_node.extent.start
    # Find the location where the function body starts ('{' or ';')
    tokens = list(function_node.get_tokens())
    body_start = None
    do_debug = False
    if function_node.spelling == "check_square":
        do_debug = True
        import pdb

        pdb.set_trace()
    sentinal = 0
    raw_tokens = [(token.spelling, token.kind) for token in tokens]
    forward_decl = generate_forward_declaration_from_tokens(raw_tokens)
    file_name = start.file.name
    forward_declaration = f"// {file_name}:{start.line}\n" + forward_decl
    # Add a semicolon if necessary
    if do_debug:
        print("Forward Declaration: ", forward_declaration)
        import pdb

        pdb.set_trace()
    if not forward_declaration.endswith(";"):
        forward_declaration += ";"
    return forward_declaration


def generate_class_forward_declaration(class_node):
    # Get the start of the class/struct declaration
    start = class_node.extent.start
    # Find the location where the class/struct body starts ('{' or ';')
    tokens = list(class_node.get_tokens())
    body_start = None
    token_full = ""
    do_debug = False
    if class_node.spelling == "index_type" and False:
        do_debug = True
        import pdb

        pdb.set_trace()
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
    print("Token Full: ", token_full)
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
    # Add a semicolon if necessary
    if do_debug:
        import pdb

        pdb.set_trace()
    if not forward_declaration.endswith(";"):
        forward_declaration += ";"
    return forward_declaration


def main(inputs):
    # Parse the command line arguments
    base_path = inputs.input_base_path
    filename = base_path + inputs.input_file
    print("Parsing: ", filename)
    args = [
        "-std=c++17",
        "-D_REENTRANT",
        "-DBOOST_DISABLE_ASSERTS",
        "-DSTAN_THREADS",
        "-O0",
        "-g",
        "-I" + base_path,
        "-I" + base_path + "/lib/tbb_2020.3/include",
        "-I" + base_path + "/lib/eigen_3.4.0",
        "-I" + base_path + "/lib/boost_1.84.0",
        "-I" + base_path + "/lib/sundials_6.1.1/include",
        "-I" + base_path + "/lib/sundials_6.1.1/src/sundials",
    ]  # Add any necessary compiler flags here
    index = clang.cindex.Index.create()
    translation_unit = index.parse(filename, args=args)
    functions, classes = get_functions_and_classes_in_namespace(
        translation_unit,
        debug_function_nodes=inputs.debug_functions,
        debug_class_nodes=inputs.debug_classes,
    )
    # Open the output file and write the forward declarations
    silentremove(inputs.output_file)
    with open(inputs.output_file, "w") as output_file:
        output_file.write(
            """
#ifndef STAN_MATH_FORWARD_DECL_HPP
#define STAN_MATH_FORWARD_DECL_HPP
#include <type_traits>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
namespace stan {
namespace math {
namespace internal {
struct nonexisting_adjoint;
template <typename T, typename S, typename Enable>
class empty_broadcast_array;
template <typename T, typename F>
struct callback_vari;
template <typename ReturnType, typename Enable, typename... Ops>
class partials_propagator;
} // namespace internal
\n"""
        )
        for class_node in classes:
            forward_decl = generate_class_forward_declaration(class_node)
            output_file.write(forward_decl + "\n\n")
        for function in functions:
            forward_decl = generate_forward_declaration(function)
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


def generate_forward_declaration_from_tokens(tokens, file_name, line_number):
    """
    Generate a forward declaration from a list of tokens.

    Args:
        tokens (list of tuples): List of (spelling, kind) tuples representing tokens.
        file_name (str): The name of the file where the function is declared.
        line_number (int): The line number where the function is declared.

    Returns:
        str: The forward declaration as a formatted string.
    """
    code_tokens = []
    nesting_level = 0
    in_template_parameter_list = False
    in_function_parameter_list = False
    i = 0
    # Collect tokens up to the first '{' or ';' (excluding it)
    while i < len(tokens):
        spelling, kind = tokens[i]
        if spelling == '{' or spelling == ';':
            break
        code_tokens.append((spelling, kind))
        i += 1
    # Now, reconstruct the code from code_tokens with appropriate formatting
    code_lines = []
    current_line = ''
    indent = ''
    nesting_level = 0
    in_template_parameter_list = False
    i = 0
    template_str = ""
    signature_str = ""
    while i < len(code_tokens):
        spelling, kind = code_tokens[i]
        # Handle 'template' keyword
        if spelling == 'template':
            template_str += 'template '
            i += 1
            continue
        # Handle entering template parameter list
        if spelling == '<' and (i > 0 and code_tokens[i - 1][0] == 'template'):
            in_template_parameter_list = True
            nesting_level += 1
            template_str += '<'
            i += 1
            continue
        # Handle template parameter list
        if in_template_parameter_list:
            if spelling == '<':
                nesting_level += 1
                template_str += '<'
            elif spelling == '>':
                nesting_level -= 1
                template_str += '>'
                if nesting_level == 0:
                    in_template_parameter_list = False
            elif spelling == ',':
                # Add comma and line break
                template_str += ',\n    '
            else:
                if kind == clang.cindex.TokenKind.PUNCTUATION:
                    template_str += spelling
                else:
                    template_str += spelling + ' '
            i += 1
            continue
        # Handle function declaration
        else:
            if spelling == '(':
                # Start of function parameter list
                in_function_parameter_list = True
                signature_str += '('
            elif spelling == ')':
                # End of function parameter list
                in_function_parameter_list = False
                signature_str += ')'
            elif spelling == ',' and in_function_parameter_list:
                # Function parameter separator
                signature_str += ', '
            else:
                if kind == clang.cindex.TokenKind.KEYWORD or kind == clang.cindex.TokenKind.IDENTIFIER:
                    signature_str += spelling + ' '
                else:
                    signature_str += spelling
            i += 1
            continue
    # Clean up extra spaces
    print("Template Str: ", template_str)
    print("Signature Str: ", signature_str)
    current_line = template_str + "\n" + signature_str
    # Ensure it ends with a semicolon
    if not current_line.endswith(';'):
        current_line += ';'
    # Add comment for file and line number
    forward_declaration = f"// {file_name}:{line_number}\n"
    forward_declaration += current_line
    return forward_declaration

# Your list of tokens
tokens = [
    ('template', clang.cindex.TokenKind.KEYWORD),
    ('<', clang.cindex.TokenKind.PUNCTUATION),
    ('typename', clang.cindex.TokenKind.KEYWORD),
    ('T_y', clang.cindex.TokenKind.IDENTIFIER),
    (',', clang.cindex.TokenKind.PUNCTUATION),
    ('require_any_t', clang.cindex.TokenKind.IDENTIFIER),
    ('<', clang.cindex.TokenKind.PUNCTUATION),
    ('is_matrix', clang.cindex.TokenKind.IDENTIFIER),
    ('<', clang.cindex.TokenKind.PUNCTUATION),
    ('T_y', clang.cindex.TokenKind.IDENTIFIER),
    ('>', clang.cindex.TokenKind.PUNCTUATION),
    (',', clang.cindex.TokenKind.PUNCTUATION),
    ('is_prim_or_rev_kernel_expression', clang.cindex.TokenKind.IDENTIFIER),
    ('<', clang.cindex.TokenKind.PUNCTUATION),
    ('T_y', clang.cindex.TokenKind.IDENTIFIER),
    ('>>', clang.cindex.TokenKind.PUNCTUATION),
    ('*', clang.cindex.TokenKind.PUNCTUATION),
    ('=', clang.cindex.TokenKind.PUNCTUATION),
    ('nullptr', clang.cindex.TokenKind.KEYWORD),
    ('>', clang.cindex.TokenKind.PUNCTUATION),
    ('inline', clang.cindex.TokenKind.KEYWORD),
    ('void', clang.cindex.TokenKind.KEYWORD),
    ('check_square', clang.cindex.TokenKind.IDENTIFIER),
    ('(', clang.cindex.TokenKind.PUNCTUATION),
    ('const', clang.cindex.TokenKind.KEYWORD),
    ('char', clang.cindex.TokenKind.KEYWORD),
    ('*', clang.cindex.TokenKind.PUNCTUATION),
    ('function', clang.cindex.TokenKind.IDENTIFIER),
    (',', clang.cindex.TokenKind.PUNCTUATION),
    ('const', clang.cindex.TokenKind.KEYWORD),
    ('char', clang.cindex.TokenKind.KEYWORD),
    ('*', clang.cindex.TokenKind.PUNCTUATION),
    ('name', clang.cindex.TokenKind.IDENTIFIER),
    (',', clang.cindex.TokenKind.PUNCTUATION),
    ('const', clang.cindex.TokenKind.KEYWORD),
    ('T_y', clang.cindex.TokenKind.IDENTIFIER),
    ('&', clang.cindex.TokenKind.PUNCTUATION),
    ('y', clang.cindex.TokenKind.IDENTIFIER),
    (')', clang.cindex.TokenKind.PUNCTUATION),
    ('{', clang.cindex.TokenKind.PUNCTUATION),
    # Function body tokens...
]

# File information
file_name = './stan/math/prim/err/check_square.hpp'
line_number = 20

# Generate forward declaration
forward_declaration = generate_forward_declaration_from_tokens(tokens, file_name, line_number)
print(forward_declaration)
