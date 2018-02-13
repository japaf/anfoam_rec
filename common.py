#!/usr/bin/env python
__author__ = 'jiri1kolar'
import logging
import distutils.spawn as ds
import os
import sys
import fnmatch

def init_logging():
    logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
            datefmt='%H:%M:%S')

def critical_error_exit(msg):
    logging.critical(msg)
    sys.exit()

def check_exec(exec_name):
    return ds.find_executable(exec_name)

def create_parent_directory(path):
    """Create directory if it does not exists.

    :param path:
    :return:
    """
    dir_name = os.path.dirname(path)
    if not os.path.exists(dir_name) and not dir_name == "":
        os.makedirs(dir_name)

def find_files_recursively(directory, pattern):
    matches = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return matches

def find_files(directory,pattern_path):
    matches=[]
    filenames=os.listdir(directory)
    for file in filenames:
        file_path=os.path.join(directory,file)
        if fnmatch.fnmatch(file_path,pattern_path):
            matches.append(file_path)
    return matches

def getparrentdir(path):
    return path[:-len(path.split('/')[-1])]


