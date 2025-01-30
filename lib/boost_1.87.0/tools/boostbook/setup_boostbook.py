#!/usr/bin/env python
#
# Copyright (c) 2002 Douglas Gregor <doug.gregor -at- gmail.com>
# Copyright (c) 2024 Andrey Semashev
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

import os
import re
import sys

# User configuration

DOCBOOK_XSL_VERSION = "1.79.2"
DOCBOOK_XSL_SHA256 = bytes.fromhex("ee8b9eca0b7a8f89075832a2da7534bce8c5478fc8fc2676f512d5d87d832102")
DOCBOOK_XSL_DIR_NAME = "docbook-xsl-nons-%s" % DOCBOOK_XSL_VERSION
DOCBOOK_XSL_TARBALL = "%s.tar.bz2" % DOCBOOK_XSL_DIR_NAME
DOCBOOK_XSL_URL = "https://github.com/docbook/xslt10-stylesheets/releases/download/release%%2F%s/%s" % (DOCBOOK_XSL_VERSION, DOCBOOK_XSL_TARBALL)
#DOCBOOK_XSL_URL = "http://sourceforge.net/projects/docbook/files/docbook-xsl/%s/%s/download" % (DOCBOOK_XSL_VERSION, DOCBOOK_XSL_TARBALL)

DOCBOOK_DTD_VERSION = "4.5"
DOCBOOK_DTD_SHA256 = bytes.fromhex("4e4e037a2b83c98c6c94818390d4bdd3f6e10f6ec62dd79188594e26190dc7b4")
DOCBOOK_DTD_DIR_NAME = "docbook-dtd-%s" % DOCBOOK_DTD_VERSION
DOCBOOK_DTD_ZIP = "docbook-xml-%s.zip" % DOCBOOK_DTD_VERSION
DOCBOOK_DTD_URL = "https://docbook.org/xml/%s/%s" % (DOCBOOK_DTD_VERSION, DOCBOOK_DTD_ZIP)
#DOCBOOK_DTD_URL = "http://www.oasis-open.org/docbook/xml/%s/%s" % (DOCBOOK_DTD_VERSION, DOCBOOK_DTD_ZIP)

FOP_VERSION = "0.94"
FOP_JDK_VERSION = "1.4"
FOP_SHA256 = bytes.fromhex("972b24b0dd2586433881bae421d92ef977c749e9ba064cacaeedc67751d035d4")
FOP_DIR_NAME = "fop-%s" % FOP_VERSION
FOP_TARBALL = "fop-%s-bin-jdk%s.tar.gz" % (FOP_VERSION, FOP_JDK_VERSION)
FOP_URL = "https://archive.apache.org/dist/xmlgraphics/fop/binaries/%s" % FOP_TARBALL
#FOP_URL = "http://mirrors.ibiblio.org/pub/mirrors/apache/xmlgraphics/fop/binaries/%s" % FOP_TARBALL

# No user configuration below this point-------------------------------------

import optparse
import shutil
import hashlib
import urllib.request
import tarfile
import zipfile

def accept_args(args):
    parser = optparse.OptionParser()
    parser.add_option("-t", "--tools", dest="tools", help=("directory downloaded tools will be installed into."
        " Optional. Used by release scripts to put the tools separately from the tree to be archived."))
    parser.usage = "setup_boostbook [options]"

    (options, args) = parser.parse_args(args)
    if options.tools is None:
        options.tools = os.getcwd()

    return options.tools

def to_posix(path):
    if path is None:
        return None
    return path.replace("\\", "/")

def unzip(archive_path, result_dir):
    with zipfile.ZipFile(archive_path, 'r', zipfile.ZIP_DEFLATED) as z:
        for f in z.infolist():
            print(f.filename)
            if not os.path.exists(os.path.join(result_dir, os.path.dirname(f.filename))):
                os.makedirs(os.path.join(result_dir, os.path.dirname(f.filename)))
            with open(os.path.join(result_dir, f.filename), 'wb') as result:
                result.write(z.read(f.filename))

def gunzip(archive_path, result_dir):
    with tarfile.open(archive_path, 'r:*') as tar:
        for tarinfo in tar:
            tar.extract(tarinfo, result_dir)

def http_get(file, url):
    with open(file, "wb") as f:
        f.write(urllib.request.urlopen(url).read())

def verify_checksum(file, sha256):
    with open(file, "rb") as f:
        digest = hashlib.file_digest(f, "sha256").digest()
    if digest != sha256:
        print("    ERROR, file checksum mismatch:\n    Local file:      %s\n    Expected SHA256: %s\n    Actual SHA256:   %s" % (file, sha256.hex(), digest.hex()))
        sys.exit(1)

def find_executable(executable_name, env_variable, error_message):
    print("Looking for %s ..." % executable_name)
    resolved_path = os.getenv(env_variable)
    if resolved_path is not None:
        print("    Trying %s specified in env. variable %s" % (resolved_path, env_variable))
        if not os.path.isfile(specified) or not os.access(resolved_path, os.X_OK):
            print("    Cannot find executable %s specified in env. variable %s" % (resolved_path, env_variable))
            resolved_path = None

    if resolved_path is None:
        resolved_path = shutil.which(executable_name)

    if resolved_path is None:
        print("    Cannot find %s executable" % executable_name)
        print(error_message)
        return None

    resolved_path = os.path.abspath(resolved_path)
    print("    Found: %s" % resolved_path)

    return resolved_path

def adjust_user_config(config_file, docbook_xsl_dir, docbook_dtd_dir, xsltproc, doxygen, fop, java):
    print("Modifying user-config.jam ...")

    os.replace(config_file, config_file + ".bak")
    with open(config_file + ".bak", "r") as in_file, open(config_file, "w") as out_file:
        boostbook_written = 0
        xsltproc_written = 0
        doxygen_written = 0
        fop_written = 0

        def write_boostbook():
            nonlocal boostbook_written, out_file, docbook_xsl_dir, docbook_dtd_dir
            if boostbook_written == 0:
                out_file.write("using boostbook\n  : \"%s\"\n  : \"%s\"\n  ;\n" % (docbook_xsl_dir, docbook_dtd_dir))
                boostbook_written = 1
        def write_xsltproc():
            nonlocal xsltproc_written, out_file, xsltproc
            if xsltproc_written == 0:
                out_file.write("using xsltproc : \"%s\" ;\n" % xsltproc)
                xsltproc_written = 1
        def write_doxygen():
            nonlocal doxygen_written, out_file, doxygen
            if doxygen_written == 0:
                if doxygen is not None:
                    out_file.write("using doxygen : \"%s\" ;\n" % doxygen)
                doxygen_written = 1
        def write_fop():
            nonlocal fop_written, out_file, fop, java
            if fop_written == 0:
                if fop is not None:
                    out_file.write("using fop\n  : \"%s\"\n  :\n  : \"%s\"\n  ;\n" % (fop, java))
                fop_written = 1

        skip_statement = 0
        for line in in_file.readlines():
            stripped_line = line
            match = re.match(r"^(.*?)#", stripped_line)
            if match is not None:
                stripped_line = match.group(1)
            stripped_line = stripped_line.strip()

            if skip_statement == 0:
                if re.match(r"^using\s+boostbook\b", stripped_line):
                    skip_statement = 1
                    write_boostbook()
                elif re.match(r"^using\s+xsltproc\b", stripped_line):
                    skip_statement = 1
                    write_xsltproc()
                elif re.match(r"^using\s+doxygen\b", stripped_line):
                    skip_statement = 1
                    write_doxygen()
                elif re.match(r"^using\s+fop\b", stripped_line):
                    skip_statement = 1
                    write_fop()

            if skip_statement == 0:
                out_file.write(line)
            elif stripped_line.endswith(";"):
                skip_statement = 0

        write_boostbook()
        write_xsltproc()
        write_doxygen()
        write_fop()

    print("    done.")

def setup_docbook_xsl(tools_directory):
    print("DocBook XSLT Stylesheets ...")

    DOCBOOK_XSL_TARBALL_PATH = os.path.join(tools_directory, DOCBOOK_XSL_TARBALL)
    if os.path.exists(DOCBOOK_XSL_TARBALL_PATH):
        print("    Using existing DocBook XSLT Stylesheets (version %s)." % DOCBOOK_XSL_VERSION)
    else:
        print("    Downloading DocBook XSLT Stylesheets version %s...\n        from %s" % (DOCBOOK_XSL_VERSION, DOCBOOK_XSL_URL))
        http_get(DOCBOOK_XSL_TARBALL_PATH, DOCBOOK_XSL_URL)

    verify_checksum(DOCBOOK_XSL_TARBALL_PATH, DOCBOOK_XSL_SHA256)

    DOCBOOK_XSL_DIR = to_posix(os.path.join(tools_directory, DOCBOOK_XSL_DIR_NAME))

    if not os.path.exists(DOCBOOK_XSL_DIR):
        print("    Expanding DocBook XSLT Stylesheets into %s..." % DOCBOOK_XSL_DIR)
        gunzip(DOCBOOK_XSL_TARBALL_PATH, tools_directory)
        print("    done.")

    return DOCBOOK_XSL_DIR

def setup_docbook_dtd(tools_directory):
    print("DocBook DTD ...")
    DOCBOOK_DTD_ZIP_PATH = to_posix(os.path.join(tools_directory, DOCBOOK_DTD_ZIP))
    if os.path.exists(DOCBOOK_DTD_ZIP_PATH):
        print("    Using existing DocBook XML DTD (version %s)." % DOCBOOK_DTD_VERSION)
    else:
        print("    Downloading DocBook XML DTD version %s..." % DOCBOOK_DTD_VERSION)
        http_get(DOCBOOK_DTD_ZIP_PATH, DOCBOOK_DTD_URL)

    verify_checksum(DOCBOOK_DTD_ZIP_PATH, DOCBOOK_DTD_SHA256)

    DOCBOOK_DTD_DIR = to_posix(os.path.join(tools_directory, DOCBOOK_DTD_DIR_NAME))
    if not os.path.exists(DOCBOOK_DTD_DIR):
        print("    Expanding DocBook XML DTD into %s... " % DOCBOOK_DTD_DIR)
        unzip(DOCBOOK_DTD_ZIP_PATH, DOCBOOK_DTD_DIR)
        print("    done.")

    return DOCBOOK_DTD_DIR

def find_xsltproc():
    return to_posix(find_executable("xsltproc", "XSLTPROC",
                                ("If you have already installed xsltproc, please set the environment\n"
                                 "variable XSLTPROC to the xsltproc executable. If you do not have\n"
                                 "xsltproc, you may download it from http://xmlsoft.org/XSLT/.")))

def find_doxygen():
    return to_posix(find_executable("doxygen", "DOXYGEN",
                                ("Warning: unable to find Doxygen executable. You will not be able to\n"
                                 "  use Doxygen to generate BoostBook documentation. If you have Doxygen,\n"
                                 "  please set the DOXYGEN environment variable to the path of the doxygen\n"
                                 "  executable.")))


def find_java():
    return to_posix(find_executable("java", "JAVA",
                                ("Warning: unable to find Java executable. You will not be able to\n"
                                 "  generate PDF documentation. If you have Java, please set the JAVA\n"
                                 "  environment variable to the path of the java executable.")))

def setup_fop(tools_directory):
    print("FOP ...")
    FOP_TARBALL_PATH = os.path.join(tools_directory, FOP_TARBALL)
    if sys.platform == 'win32':
        fop_driver = "fop.bat"
    else:
        fop_driver = "fop"

    FOP_DIR = to_posix(os.path.join(tools_directory, FOP_DIR_NAME))
    FOP = to_posix(os.path.join(FOP_DIR, fop_driver))

    if os.path.exists(FOP_TARBALL_PATH):
        print("    Using existing FOP distribution (version %s)." % FOP_VERSION)
    else:
        print("    Downloading FOP distribution version %s..." % FOP_VERSION)
        http_get(FOP_TARBALL_PATH, FOP_URL)

    verify_checksum(FOP_TARBALL_PATH, FOP_SHA256)

    if not os.path.exists(FOP_DIR):
        print("    Expanding FOP distribution into %s... " % FOP_DIR)
        gunzip(FOP_TARBALL_PATH, tools_directory)
        print("    done.")

    return FOP

def find_user_config():
    print("Looking for user-config.jam ...")
    config_dir = os.getenv("HOME")
    if config_dir is not None:
        jam_config = os.path.join(config_dir, "user-config.jam")
        if os.path.isfile(jam_config):
            print("    Found user-config.jam in HOME directory: %s" % jam_config)
            return jam_config

    config_dir = os.getenv("BOOST_ROOT")
    if config_dir is not None:
        jam_config = os.path.join(config_dir, "tools/build/user-config.jam")
        if os.path.isfile(jam_config):
            print("    Found user-config.jam in BOOST_ROOT directory: %s" % jam_config)
            return jam_config

    print("    Not found")
    return None

def setup_boostbook(tools_directory):
    print(("Setting up boostbook tools...\n"
           "-----------------------------\n"))

    user_config     = find_user_config()
    if user_config is None:
        print(("ERROR: Please set the BOOST_ROOT environment variable to refer to your\n"
               "Boost installation or copy user-config.jam into your home directory."))
        sys.exit(1)

    DOCBOOK_XSL_DIR = setup_docbook_xsl(tools_directory)
    DOCBOOK_DTD_DIR = setup_docbook_dtd(tools_directory)
    XSLTPROC        = find_xsltproc()
    DOXYGEN         = find_doxygen()
    JAVA            = find_java()

    FOP             = None
    if JAVA is not None:
        print("Java is present.")
        FOP = setup_fop(tools_directory)

    adjust_user_config(config_file     = user_config,
                       docbook_xsl_dir = DOCBOOK_XSL_DIR,
                       docbook_dtd_dir = DOCBOOK_DTD_DIR,
                       xsltproc        = XSLTPROC,
                       doxygen         = DOXYGEN,
                       fop             = FOP,
                       java            = JAVA)

    print(("Done! Execute \"b2\" in a documentation directory to generate\n"
           "documentation with BoostBook. If you have not already, you will need\n"
           "to compile Boost.Jam."))

def main():
    (tools_directory) = accept_args(sys.argv[1:])
    setup_boostbook(tools_directory)

if __name__ == "__main__":
    main()
