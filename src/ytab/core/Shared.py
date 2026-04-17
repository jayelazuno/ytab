"""A utility module, containing various shared functions.
"""

import os
import collections
import collections.abc
import functools


def get_script_dir():
    """Get the base folder of this script.

    Returns
    -------
    str
        The base folder of this script.
    """

    return os.path.dirname(os.path.realpath(__file__))


def get_repo_root():
    """Return the YTAB repository root.

    Shared.py lives at:
        <repo>/src/ytab/core/Shared.py

    so going up three levels from this script directory gives <repo>.
    """

    return os.path.dirname(os.path.dirname(os.path.dirname(get_script_dir())))


def get_resources_dir():
    """Return the absolute path to the species resource root.

    Resolution order
    ----------------
    1. YTAB_RESOURCES_DIR environment variable
    2. <repo>/resources/species
    """

    env_path = os.environ.get("YTAB_RESOURCES_DIR")
    if env_path:
        return os.path.realpath(os.path.expanduser(env_path))

    return os.path.join(get_repo_root(), "resources", "species")


def get_species_dir(species_name):
    """Return the absolute path to one species resource directory.

    Parameters
    ----------
    species_name : str
        Species folder name, e.g. 'glabrata'.

    Returns
    -------
    str
        Absolute path to the species resource directory.
    """

    return os.path.join(get_resources_dir(), species_name)


def get_dependency(*path_components):
    """Given path components under the species resource root, return the
    absolute path.

    Examples
    --------
    get_dependency("glabrata", "C_glabrata_BG2_S_cerevisiae_orthologs.txt")
    get_dependency("glabrata", "reference_genome",
                   "GCA_014217725.1_ASM1421772v1_genomic.fna")
    """

    return os.path.join(get_resources_dir(), *path_components)


def flatten(seq_of_seqs):
    """Flatten a sequence of sequences into a single sequence.

    Returns
    -------
    list
        The flattened list.
    """

    return [item for sublist in seq_of_seqs for item in sublist]


def make_dir(dir_path):
    """Makes sure the specified path exists as a folder.

    If the path doesn't exist, it recursively goes back in the folder tree
    and creates each subpath.

    Parameters
    ----------
    dir_path : str
        The folder to create (if necessary).
    """

    if not dir_path:
        return

    if not os.path.isdir(dir_path):
        make_dir(os.path.dirname(dir_path))
        os.mkdir(dir_path)


# This code was borrowed from
# https://wiki.python.org/moin/PythonDecoratorLibrary#Memoize

class memoized(object):
   '''Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
   (not reevaluated).
   '''
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      if not isinstance(args, collections.abc.Hashable):
         # uncacheable. a list, for instance.
         # better to not cache than blow up.
         return self.func(*args)
      if args in self.cache:
         return self.cache[args]
      else:
         value = self.func(*args)
         self.cache[args] = value
         return value
   def __repr__(self):
      '''Return the function's docstring.'''
      return self.func.__doc__
   def __get__(self, obj, objtype):
      '''Support instance methods.'''
      return functools.partial(self.__call__, obj)