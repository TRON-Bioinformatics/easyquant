"""
File handler module.
"""

def write_line_to_file(file_handle: object, values: list) -> bool:
    """Writes a specific line to an open file handle."""
    try:
        file_handle.write(
            "{}\n".format(
                "\t".join([str(e) for e in values])
            ).encode()
        )
        return True
    except IOError:
        return False
