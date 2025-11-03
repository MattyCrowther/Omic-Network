class ParseError(ValueError):
    """Raised when a parser fails to read or validate an omics data file."""
    def __init__(self, message: str, *, path: str | None = None, cause: Exception | None = None):
        if path:
            message = f"{message} [file={path}]"
        super().__init__(message)
        self.path = path
        self.cause = cause