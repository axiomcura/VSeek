# base classes for expcetions
class Error(Exception):
    """Base class error for other exceptions"""
    pass

class DependencyError(Exception):
    """ Base class for any dependency related error"""
    pass


# user-based expcetions
class MissingDependencyError(DependencyError):
    """ Raised when any missing dependencies are"""
    pass

class UnsupportedDependencyError(DependencyError):
    """ Raised when calling an unspported dependency"""
    pass

class ExecutionError(DependencyError):
    """Raised when a dependency executable failed"""