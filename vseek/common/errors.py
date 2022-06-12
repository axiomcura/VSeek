# base classes for expcetions
class FileNotFoundError(Exception):
    """Base exceptions for missing files"""

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


class ProfileNotFound(FileNotFoundError):
    """Raised if genome profile is missing"""

    pass


class GenomeDirectoryNotFound(FileNotFoundError):
    """Raised if the Genome directory is not found"""

    pass


class FastaFileNotFound(FileNotFoundError):
    """Raised if the FASTA file is not found"""

    pass


class SequenceFormatNotSupported(InvalidFormats):
    """Raised when attempting to create a specific sequence format
    but it is not supported"""

    pass
