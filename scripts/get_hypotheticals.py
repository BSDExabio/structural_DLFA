#!/usr/bin/env python3
""" This will extract any hypothetical proteins from a given Genbank file.

"""
import argparse
import logging
import string
import sys
from pathlib import Path

from rich import pretty
pretty.install()

from rich.logging import RichHandler
rich_handler = RichHandler(rich_tracebacks=True,
                           markup=True)
logging.basicConfig(level='INFO', format='%(message)s',
                    datefmt="[%Y/%m/%d %H:%M:%S]",
                    handlers=[rich_handler])
logger = logging.getLogger(__name__)

from rich.traceback import install
install()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='find hypothetical proteins in Genbank files')
    parser.add_argument('genbank_file', help='Genbank file in which we want to find hypotheticals')

    args = parser.parse_args()
