"""
Various useful utilities, mostly taken from the ASPN Python cookbook.
"""

import types

seq_types = type(()), type([])


def flatten(*args):
    for arg in args:
        if type(arg) in seq_types:
            for elem in arg:
                yield from flatten(elem)
        else:
            yield arg


def cross_lists(*sets):
    """Return the cross product of the arguments"""
    wheels = [iter(_) for _ in sets]
    digits = [next(it) for it in wheels]
    while True:
        yield digits[:]
        for i in range(len(digits) - 1, -1, -1):
            try:
                digits[i] = next(wheels[i])
                break
            except StopIteration:
                wheels[i] = iter(sets[i])
                digits[i] = next(wheels[i])
        else:
            break


# Cached / memoized methods


def cachedmethod(function):
    return types.MethodType(Memoize(function), None)


class Memoize:
    def __init__(self, function):
        self._cache = {}
        self._callable = function

    def __call__(self, *args, **kwds):
        cache = self._cache
        key = self._getKey(*args, **kwds)
        try:
            return cache[key]
        except KeyError:
            cachedValue = cache[key] = self._callable(*args, **kwds)
            return cachedValue

    def _getKey(self, *args, **kwds):
        return kwds and (args, ImmutableDict(kwds)) or args


class memoized:
    """Decorator that caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned, and
    not re-evaluated.
    """

    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        try:
            return self.cache[args]
        except KeyError:
            self.cache[args] = value = self.func(*args)
            return value
        except TypeError:
            # uncachable -- for instance, passing a list as an argument.
            # Better to not cache than to blow up entirely.
            return self.func(*args)

    def __repr__(self):
        """Return the function's docstring."""
        return self.func.__doc__


class ImmutableDict(dict):
    """A hashable dict."""

    def __init__(self, *args, **kwds):
        dict.__init__(self, *args, **kwds)

    def __setitem__(self, key, value):
        raise NotImplementedError("dict is immutable")

    def __delitem__(self, key):
        raise NotImplementedError("dict is immutable")

    def clear(self):
        raise NotImplementedError("dict is immutable")

    def setdefault(self, k, default=None):
        raise NotImplementedError("dict is immutable")

    def popitem(self):
        raise NotImplementedError("dict is immutable")

    def update(self, other):
        raise NotImplementedError("dict is immutable")

    def __hash__(self):
        return hash(tuple(self.items()))
