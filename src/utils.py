import os


def check_and_create_dir(path):
    path = os.path.dirname(path)
    if not os.path.exists(path):
        os.makedirs(path)
