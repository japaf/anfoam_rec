#!/usr/bin/env python
__author__ = 'jiri1kolar'
import logging

def init_logging():
    logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
            datefmt='%H:%M:%S')