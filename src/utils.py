import os


def check_and_create_dir(path: str) -> bool:
    """
    Check if directory exists, create it if not
    :param path: Path to directory
    :return: True if directory was created, False if directory already exists
    """
    path = os.path.dirname(path)
    if not os.path.exists(path):
        os.makedirs(path)
        return True
    return False
