import os
import os.path
import sys


def check_path_with_msg(this_path):
    """
    If the specified path exists on the file system, return True.
    If it does not exist print a message to stdout and return False.
    """
    try:
        if not(os.path.exists(this_path)):
            sys.stdout.write('path not found: {}\n'.format(this_path))
            sys.stdout.flush()
            return(False)
        else:
            return(True)
    except TypeError:
        if this_path is None:
            sys.stdout.write('some paths are None\n')
            sys.stdout.flush()
            return(False)
