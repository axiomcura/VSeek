# base classes for expcetions
class Error(Exception):
    """Base class error for other exceptions"""

    pass


class InvalidFormats(Exception):
    """Base class for invalid formats exceptions"""

    pass


class DependencyError(Exception):
    """Base class for any dependency related error"""

    pass


# user-based expcetions
class MissingDependencyError(DependencyError):
    """Raised when any missing dependencies are"""

    pass


class UnsupportedDependencyError(DependencyError):
    """Raised when calling an unsupported dependency"""

    pass


class ExecutionError(DependencyError):
    """Raised when a dependency executable failed"""

    pass


class InvalidFileError(InvalidFormats):
    """Raised when an invalid file is used for loading or processing"""

    pass


class InvalidFastaFormatError(InvalidFormats):
    """Raised if FASTA contents possess an invalid format"""

    pass
