import os
import collections
import re



__all__ = ['tail', 'fall_back_tail', 'grep']


def tail(fh, lines=20, block_size=1024):
    """Returns the last n lines from a file

    This function returns the last n lines from an file-like
    object. It does this efficiently without reading the whole
    file, but rather by loading blocks from the end of the file.

    .. note::

        If the file is opened in text mode, i.e. open('/path', 'rt'),
        python3 cannot efficiently move in the file. In this case,
        the function fall back to a slow method that goes through
        the whole file.

    Example:

    >>> with open("/etc/passwd", 'rb') as f:
    ...     last_lines = tail(f, 2)
    ...
    >>> print(last_lines)

    :param fh: file-like object to read from
    :param int lines: number of lines to be returned
    :param int block_size: size of block to be read at once.
        intended for optimisation.
    :returns: The last lines as a list of bytes/str object"""

    if lines <= 0:
        raise ValueError('invalid lines value %r' % lines)

    encoded = getattr(fh, 'encoding', False)
    if encoded:
        return fall_back_tail(fh, lines)
    CR = '\n' if encoded else b'\n'
    data = '' if encoded else b''
    fh.seek(0, os.SEEK_END)
    fsize = fh.tell()
    block = -1
    loaded_enough_data = False
    while not loaded_enough_data:
        step = (block * block_size)
        if abs(step) >= fsize:
            fh.seek(0)
            newdata = fh.read(block_size - (abs(step) - fsize))
            loaded_enough_data = True
        else:
            fh.seek(step, os.SEEK_END)
            newdata = fh.read(block_size)
        data = newdata + data
        if data.count(CR) > lines:
            break
        else:
            block -= 1
    return data.splitlines()[-lines:]


def fall_back_tail(fh, lines):
    fh.seek(0)
    data = collections.deque(fh, maxlen=lines)
    return [e.rstrip('\n') for e in data]


def grep(fh, pat):
    """Yields lines matching a pattern

    This function yields all the lines that match a given pattern.
    The pattern can be either a simple str/bytes, or a compiled
    regex expression. The newline character is not removed.

    Example:
        >>> with open('/etc/hosts', 'rb') as fh:
        ...    for line in grep(fh, b'127.0.0.1'):
        ...        print(line)
        127.0.0.1       localhost

    :param fh: file-like object
    :param pat: search pattern, either str, bytes or compiled regex
    :returns: generator yielding lines matching pattern.

    """
    if isinstance(pat, (str, bytes)):
        encoded = getattr(fh, 'encoding', False)
        if encoded and isinstance(pat, bytes):
            pat = re.compile(pat.decode())
        elif not encoded and isinstance(pat, str):
            pat = re.compile(pat.encode('utf-8'))
        else:
            pat = re.compile(pat)
    fh.seek(0)
    for line in fh:
        if pat.search(line):
            yield line
