#!/usr/bin/env python
'''
    Utilities for zoo files.
'''
from io import BytesIO
import bz2
import gzip
import os
import logging
logger = logging.getLogger(__name__)


# File opening. This is based on the example on SO here:
# http://stackoverflow.com/a/26986344
fmagic = {b'\x1f\x8b\x08': gzip.open,
          b'\x42\x5a\x68': bz2.BZ2File}


def auto_open(fn, *args, **kwargs):
    """function to open regular or compressed files for read / write.

    This function opens files based on their "magic bytes". Supports bz2
    and gzip. If it finds neither of these, presumption is it is a
    standard, uncompressed file.

    Example::

        with auto_open("/path/to/file/maybe/compressed", mode="rb") as fh:
            fh.read()

        with auto_open("/tmp/test.txt.gz", mode="wb") as fh:
            fh.write("my big testfile")

    :param fn: either a string of an existing or new file path, or
        a BytesIO handle
    :param \*\*kwargs: additional arguments that are understood by the
        underlying open handler
    :returns: a file handler
    """
    if isinstance(fn, BytesIO):
        return fn

    if os.path.isfile(fn) and os.stat(fn).st_size > 0:
        with open(fn, 'rb') as fp:
            fs = fp.read(max([len(x) for x in fmagic]))
        for (magic, _open) in fmagic.items():
            if fs.startswith(magic):
                return _open(fn, *args, **kwargs)
    else:
        if fn.endswith('gz'):
            return gzip.open(fn, *args, **kwargs)
        elif fn.endswith('bz2'):
            return bz2.BZ2File(fn, *args, **kwargs)

    return open(fn, *args, **kwargs)


class LazyProperty(object):
    """Decorator to evaluate a property only on access.

    Compute the attribute value and caches it in the instance.
    Python Cookbook (Denis Otkidach) http://stackoverflow.com/users/168352/denis-otkidach
    This decorator allows you to create a property which can be computed once and
    accessed many times.

    Example::

        class Circle:
            def __init__(self, radius):
                self.radius = radius

            @LazyProperty
            def area(self):
                print("computing area")
                return 3.14 * self.radius ** 2

        >>> c = Circle(4)
        >>> c.area
        computing area
        50.24
        >>> c.area
        50.24

    You can see that the property method is only executed once.
    """

    def __init__(self, method, name=None):
        # record the unbound-method and the name
        self.method = method
        self.name = name or method.__name__
        self.__doc__ = method.__doc__

    def __get__(self, inst, cls):
        if inst is None:
            return self
        # compute, cache and return the instance's attribute value
        result = self.method(inst)
        # setattr redefines the instance's attribute so this doesn't get called again
        setattr(inst, self.name, result)
        return result


def unique(seq):
    """Return the elements of a list uniquely while preserving the order

    :param list seq: a list of hashable elements
    :returns: new list with first occurence of elements of seq"""
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]



