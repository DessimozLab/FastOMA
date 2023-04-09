from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future.builtins import next
from future.builtins import str
from future import standard_library
standard_library.install_hooks()
from past.builtins import basestring

import io
import collections
from .tools import enum, py2_iterable


class LexError(Exception):

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


@py2_iterable
class Streamer(object):

    """ Wraps an io.StringIO and iterates a byte at a time,
    instead of a line at a time """

    def __init__(self, stream):
        """ _peek always looks ahead 1 position """
        if isinstance(stream, str) or isinstance(stream, basestring):
            stream = io.StringIO(u'{}'.format(stream))
        self.stream = stream
        self._peek = self.stream.read(1)

    def __iter__(self):
        return self

    def __next__(self):
        char = self._peek

        self._peek = self.stream.read(1)

        if self.stream.closed:
            raise StopIteration

        if char == '':
            self.stream.close()
            raise StopIteration

        return char
    def peek(self):
        return self._peek

    def isclosed(self):
        return self.stream.closed

Token = collections.namedtuple('Token', 'typ val')


@py2_iterable
class NewickLexer(object):

    """ Breaks newick stream into lexing tokens:
    Works as a state machine, like Rob Pike's Go text template parser """

    tokens = enum("EOF", "TREE", "LEAF", "SUBTREE", "LABEL",
                  "LENGTH", "SUPPORT", "ENDSUB", "ENDTREE")

    def __init__(self, streamer):
        self.streamer = streamer
        self.token = None
        self.token_buffer = bytearray()
        self.state = self.lex_tree

    def buffer(self, char):
        """ Adds a streamed character to the token buffer """
        self.token_buffer.append(ord(char))

    def eat_spaces(self):
        while self.streamer.peek().isspace():
            next(self.streamer)

    def emit(self, item):
        """ Emits the token buffer's contents as a token; clears the buffer """
        self.token = item
        self.empty_buffer()

    def empty_buffer(self):
        """ Clears the token buffer (python 2 has no bytearray.clear()
        method ) """
        self.token_buffer = self.token_buffer[0:0]

    def __iter__(self):
        return self

    def __next__(self):
        """ Each iteration returns a token. While a token isn't ready,
        advance the state machine one state. """
        while not self.token:
            self.state = self.state()
        token, self.token = self.token, None
        return token

    def pos(self):
        """ Returns position in input stream """
        return self.streamer.stream.tell()

    def truncated_string(self, s, length=60, ellipsis='...'):
        """ Returns a string `s` truncated to maximum length `length`.
        If `s` is longer than `length` it is truncated and `ellipsis` is
        appended to the end. The ellipsis is included in the length.
        If `s` is shorter than `length` `s` is returned unchanged. """
        l = length - len(ellipsis)
        return s[:l] + (s[l:] and ellipsis)

    def lex_tree(self):
        for x in self.streamer:
            if x == '(':
                break

        if self.streamer.isclosed():
            self.emit(Token(self.tokens.EOF, -1))
            raise StopIteration

        self.emit(Token(self.tokens.TREE, ''))
        return self.lex_subtree_start

    def lex_subtree_start(self):
        self.eat_spaces()
        char = self.streamer.peek()

        if char == '(':
            self.emit(Token(self.tokens.SUBTREE, next(self.streamer)))
            return self.lex_subtree_start

        else:
            self.emit(Token(self.tokens.LEAF, None))
            return self.lex_label

    def lex_label(self):
        self.eat_spaces()
        char = self.streamer.peek()
        if char in ('"', "'"):
            next(self.streamer)  # throw away opening quote
            self._match_delimited(char)
        else:
            despacer = {' ': '_'}
            self._match_run(str.isalnum, accepted_chars='-_|.',
                            denied_chars=':,;', replacements=despacer)
        label = self.token_buffer.decode()
        if label == '':
            label = None
        self.emit(Token(self.tokens.LABEL, label))

        return self.lex_length

    def lex_length(self):
        char = self.streamer.peek()
        if char == ':':
            self.streamer.next()  # throw away colon
            self._match_number()
            if len(self.token_buffer) == 0:
                num = None
            else:
                num = float(self.token_buffer)
        else:
            num = None
        self.emit(Token(self.tokens.LENGTH, num))
        return self.lex_subtree_end

    def lex_subtree_end(self):
        self.eat_spaces()
        char = self.streamer.peek()

        if char == ';':
            next(self.streamer)
            self.emit(Token(self.tokens.ENDTREE, ';'))
            return self.lex_tree

        elif char == ',':
            next(self.streamer)
            return self.lex_subtree_start

        elif char == ')':
            next(self.streamer)
            self.emit(Token(self.tokens.ENDSUB, ')'))
            peek = self.streamer.peek()  # is a label or a support value next?
            if peek.isdigit() or peek == '.':
                return self.lex_support
            return self.lex_label

        else:
            raise LexError('Don\'t know how to lex this: {0} ({1})'.format(
                char, self.streamer.stream.tell()))

    def lex_support(self):
        self._match_number()
        if len(self.token_buffer) == 0:
            num = 0.0
        else:
            num = float(self.token_buffer)
        self.emit(Token(self.tokens.SUPPORT, num))
        return self.lex_length

    def _match_delimited(self, delimiter):
        pos = self.pos() - 2  # stream is 2 chars ahead of the opening delimiter
        for char in self.streamer:
            if char == delimiter:
                return
            self.buffer(char)

        buf = self.truncated_string(self.token_buffer.decode())
        msg = ''.join((
            'Unterminated {0}-delimited string starting at '.format(delimiter),
            'position {0}:\n{1}{2}'.format(pos, delimiter, buf)
            ))
        raise LexError(msg)

    def _match(self, predicate, accepted_chars='', denied_chars='',
               replacements=None):
        """ Checks next character in stream. If predicate returns True, or char
        is in `accepted_chars`, advances the stream and returns 1. Else, or if
        the char is in `denied_chars`, doesn't advance the stream and returns 0.
        Replacements is an optional dictionary that can be used to replace the
        streamed character with an alternative (e.g. replace spaces with
        underscores). """

        replacements = (replacements or {})
        char = self.streamer.peek()
        char = replacements.get(char, char)

        if predicate(char) or char in accepted_chars:
            if len(char) == 1:
                self.buffer(char)
            next(self.streamer)  # advance stream
            return 1
        elif char in denied_chars:
            return 0
        else:
            return 0

    def _match_run(self, predicate, **kwargs):
        """ kwargs are `accepted_chars`, `denied_chars` and `replacements` """
        nchars = 0
        try:
            while True:
                matched = self._match(predicate, **kwargs)
                nchars += matched
                if matched == 0:
                    return nchars
        except StopIteration:
            raise LexError('Unexpected end of stream')

    def _match_number(self):
        digits = 0
        self._match(lambda x: False, '-+')
        digits += self._match_run(str.isdigit)
        if self._match(lambda x: False, '.'):
            digits += self._match_run(str.isdigit)

        if digits > 0:
            if self._match(lambda x: False, 'eE'):
                self._match(lambda x: False, '-+')
                self._match_run(str.isdigit)

        else:
            self.empty_buffer()
