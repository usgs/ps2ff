#!/usr/bin/env python
"""
Program to generate API docs from doc strings and put the results
into the docs directory.
"""

import argparse
import os.path
import sys
import pathlib

from ps2ff.cmd import get_command_output


def main(args):
    """
    Generate API docs.

    Args:
        args: Output of argparse. Currently only holds the verbose flag.

    Returns:
        Nothing. Function will exit upon success or failure.

    """
    # -------------------------------------------------------------
    # Some useful directories
    # -------------------------------------------------------------
    REPO_DIR = os.path.dirname(os.path.abspath(__file__))
    DOCS_DIR = os.path.join(REPO_DIR, 'docs')
    API_DIR = os.path.join(REPO_DIR, 'doc_source')
    PACKAGE_DIR = os.path.join(REPO_DIR, 'ps2ff')

    # -------------------------------------------------------------
    # what is the package called and who are the authors
    # -------------------------------------------------------------
    PACKAGE = "ps2ff"
    AUTHORS = 'Eric Thompson, Bruce Worden'
    verstr = '1.1'

    # -------------------------------------------------------------
    # run the api doc command; this creates the .rst files
    # -------------------------------------------------------------
    sys.stderr.write('Building ps2ff API documentation (REST)...\n')
    sphinx_cmd = 'sphinx-apidoc -o %s -f -e -l -M -d 12 -H %s -A "%s"'\
                 ' -V %s %s' % (API_DIR, PACKAGE, AUTHORS, verstr,
                                PACKAGE_DIR)
    res, stdout, stderr = get_command_output(sphinx_cmd)

    if not res:
        raise Exception('Could not build ps2ff API documentation'
                        ' - error "%s".' % stderr)

    if args.verbose:
        print(stdout.decode('utf-8'))
        print(stderr.decode('utf-8'))

    # --------------------------------------------
    # try to clean up some of the excess labeling
    # --------------------------------------------
    clean_cmd = "sed -i '' -e 's/ module//g' `find %s/*.rst -type f "\
                "-maxdepth 0 -print`" % API_DIR
    res, stdout, stderr = get_command_output(clean_cmd)
    clean_cmd = "sed -i '' -e 's/ package//g' `find %s/*.rst -type f "\
                "-maxdepth 0 -print`" % API_DIR
    res, stdout, stderr = get_command_output(clean_cmd)
    clean_cmd = "sed -i '' -e '/Subpackages/d' `find %s/*.rst -type f "\
                "-maxdepth 0 -print`" % API_DIR
    res, stdout, stderr = get_command_output(clean_cmd)
    clean_cmd = "sed -i '' -e '/-.*-/d' `find %s/*.rst -type f "\
                "-maxdepth 0 -print`" % API_DIR
    res, stdout, stderr = get_command_output(clean_cmd)

    # -------------------------------------------------------------
    # Build the html
    # -------------------------------------------------------------
    sys.stderr.write('Building ps2ff pages (HTML)...\n')
    res, stdout, stderr = get_command_output(
        'sphinx-build -a -E doc_source docs')
    if not res:
        raise Exception('Could not build HTML for API documentation. - '
                        'error "%s"' % stderr.decode())
    if args.verbose:
        print(stdout.decode('utf-8'))
        print(stderr.decode('utf-8'))

    pathlib.Path(os.path.join(DOCS_DIR, '.nojekyll')).touch(exist_ok=True)


if __name__ == '__main__':
    desc = 'Create API documentation for ps2ff.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='Produce more output to the screen. ')

    pargs = parser.parse_args()
    main(pargs)
