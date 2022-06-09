import shutil

# program imports
from vseek.common.errors import UnsupportedDependencyError, MissingDependencyError

# dowloading the data module
def dependency_check(prog: str) -> bool:
    """Checks if the executable program is available

    Parameters
    ----------
    prog : str
        name of the program or executable

    Returns
    -------
    Bool
        Returns True if the program exists

    Raises:
    -------
    MissingDependencyError:
        Raised If the program is not installed to local machine
    UnsupportedDependencyError:
        Raised if attempt to call unsupported executable. Please look at supported
        executables in the README.md file.
    """
    callables = ["prefetch", "fasterq-dump"]

    # checking if the callable is supported
    if prog not in callables:
        raise UnsupportedDependencyError(
            f"{prog} is not supported. Supported programs {callables}"
        )

    # checking if supported callable is installed
    check = shutil.which(prog)
    if check is None:
        raise MissingDependencyError(f"{prog} is not installed")

    return True

