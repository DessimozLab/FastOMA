from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future.builtins import dict
from future.builtins import zip
from future.builtins import range
from future import standard_library
standard_library.install_hooks()

try:
    from progressbar import ProgressBar, Percentage, Timer, ETA, Bar
    PROGRESSBAR = True
except ImportError:
    PROGRESSBAR = False

from collections import deque

def setup_progressbar(msg, size):
    if not msg.endswith(': '):
        msg += ': '

    widgets = [msg,
               Percentage(), ' ',
               Bar(), ' ',
               Timer(), ' ',
               ETA()]

    pbar = ProgressBar(widgets=widgets, maxval=size)
    return pbar

def enum(*sequential, **named):
    """creates an Enum type with given values"""
    enums = dict(zip(sequential, range(len(sequential))), **named)
    enums['reverse'] = dict((value, key) for key, value in enums.items())
    return type('Enum', (object, ), enums)


class IterableClassException(Exception):
    pass

def py2_iterable(Class):
    """
    Use as a class decorator to make a class that has a python 3 next method --
    __next__() -- also iterable with python 2, which uses next(). Also checks
    for an __iter__ method -- if this is missing the class won't be iterable anyway.


    e.g.
    @py2_iterable
    class Py2and3Iterator(object):
        def __init__(self):
            self.data = list('somestuff')
            self._pos = 0

        def __iter__(self):
            return self

        def __next__(self):
            if self._pos == len(self.data):
                self._pos = 0
                raise StopIteration
            char = self.data[self._pos]
            self._pos += 1
            return char


    :param Class: the class being decorated
    :return: Class: the decorated class, which is iterable in py2 and py3
    """
    if not hasattr(Class, '__iter__'):
        raise IterableClassException('Class "{}" has no __iter__ method and will not be iterable'
                                     .format(Class.__class__.__name__))

    if hasattr(Class, '__next__'):
        next_method = getattr(Class, '__next__')
        setattr(Class, 'next', next_method)

    return Class


@py2_iterable
class Queue(object):

    def __init__(self):
        self.__queue = deque()

    def __iter__(self):
        return self

    def __len__(self):
        return len(self.__queue)

    def __next__(self):
        if self.isempty():
            raise StopIteration
        return self.dequeue()

    def enqueue(self, item):
        self.__queue.append(item)

    def dequeue(self):
        if self.isempty():
            raise Exception('empty queue')
        return self.__queue.popleft()

    def isempty(self):
        return len(self.__queue) == 0
