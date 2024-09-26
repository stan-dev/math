import sys
import clang.cindex
import glob
import errno
import os
from argparse import ArgumentParser, RawTextHelpFormatter



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
    in_default_value = False
    default_nesting_level = 0
    i = 0
    template_str = ""
    signature_str = ""
    print("BEGIN\n\n\n")
    first_entry = True
    while i < len(code_tokens):
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
    print("Template Str: ", template_str)
    print("Signature Str: ", signature_str)
    current_line = template_str + "\n" + signature_str
    # Ensure it ends with a semicolon
    if not current_line.endswith(';'):
        current_line += ';'
    # Add comment for file and line number
    forward_declaration = current_line
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
